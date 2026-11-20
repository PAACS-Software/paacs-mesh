// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/PolyhedralModel.hpp"
#include "paacs/predicates/Predicates.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>

namespace paacs::geo {
namespace {
bool finite(Vec3 p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
double length(Vec3 p) {
    return std::hypot(p.x, p.y, p.z);
}
using EdgeKey = std::pair<std::uint32_t, std::uint32_t>;
EdgeKey edgeKey(std::uint32_t a, std::uint32_t b) {
    return std::minmax(a, b);
}
struct Stop {
    PolyhedralReason reason;
};
struct Work {
    std::uint64_t count = 0, limit = 0;
    void tick(std::uint64_t n = 1) {
        if (n > limit - count)
            throw Stop{PolyhedralReason::BudgetExceeded};
        count += n;
    }
};
int orient(Vec2 a, Vec2 b, Vec2 c) {
    return static_cast<int>(predicates::orient2d(a, b, c));
}
bool onSegment(Vec2 a, Vec2 b, Vec2 p) {
    return orient(a, b, p) == 0 && p.u >= std::min(a.u, b.u) && p.u <= std::max(a.u, b.u) &&
           p.v >= std::min(a.v, b.v) && p.v <= std::max(a.v, b.v);
}
bool intersects(Vec2 a, Vec2 b, Vec2 c, Vec2 d) {
    const int x = orient(a, b, c), y = orient(a, b, d), z = orient(c, d, a), w = orient(c, d, b);
    return (x * y < 0 && z * w < 0) || (x == 0 && onSegment(a, b, c)) ||
           (y == 0 && onSegment(a, b, d)) || (z == 0 && onSegment(c, d, a)) ||
           (w == 0 && onSegment(c, d, b));
}
// -1 outside, 0 boundary, 1 inside. Half-open crossing rule avoids counting
// a ring vertex twice; exact orientation determines the crossing side.
int inRing(const std::vector<Vec2> &ring, Vec2 p) {
    bool inside = false;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const Vec2 a = ring[i], b = ring[(i + 1) % ring.size()];
        if (onSegment(a, b, p))
            return 0;
        if ((a.v > p.v) != (b.v > p.v)) {
            const int s = orient(a, b, p);
            if ((b.v > a.v && s > 0) || (b.v < a.v && s < 0))
                inside = !inside;
        }
    }
    return inside ? 1 : -1;
}
long double signedArea(const std::vector<Vec2> &ring) {
    long double a = 0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const auto p = ring[i], q = ring[(i + 1) % ring.size()];
        a += static_cast<long double>(p.u) * q.v - static_cast<long double>(p.v) * q.u;
    }
    return a / 2;
}
class Line final : public GEdge {
    const GVertex *a_, *b_;

  public:
    Line(int tag, const GVertex *a, const GVertex *b) : GEdge(tag), a_(a), b_(b) {}
    CurveKind kind() const noexcept override { return CurveKind::Line; }
    Vec3 pointAt(double t) const override {
        return a_->position() + (b_->position() - a_->position()) * t;
    }
    Vec3 tangentAt(double) const override {
        const auto d = b_->position() - a_->position();
        return d / geo::length(d);
    }
    double length() const override { return geo::length(b_->position() - a_->position()); }
    int startVertexTag() const override { return a_->tag(); }
    int endVertexTag() const override { return b_->tag(); }
    double closestParameter(const Vec3 &p) const override {
        const auto d = b_->position() - a_->position();
        const auto size = geo::length(d);
        return std::clamp((p - a_->position()).dot(d / size) / size, 0.0, 1.0);
    }
};
class Face final : public GFace {
  public:
    Vec3 origin, u, v, normal;
    Vec2 lo{}, hi{};
    std::vector<std::vector<Vec2>> rings;
    std::vector<std::vector<int>> edges;
    explicit Face(int tag) : GFace(tag) {}
    SurfaceKind kind() const noexcept override { return SurfaceKind::Plane; }
    Vec3 pointAt(const Vec2 &uv) const override { return origin + u * uv.u + v * uv.v; }
    Vec2 paramAt(const Vec3 &p) const override {
        const auto d = p - origin;
        return {d.dot(u), d.dot(v)};
    }
    Vec3 normalAt(const Vec2 &) const override { return normal; }
    void paramBounds(Vec2 &a, Vec2 &b) const override {
        a = lo;
        b = hi;
    }
    std::vector<int> boundaryEdgeTags() const override { return edges.front(); }
    std::vector<std::vector<int>> boundaryLoops() const override { return edges; }
    double curvatureAt(const Vec2 &) const override { return 0; }
    bool containsParam(const Vec2 &p) const override {
        if (inRing(rings.front(), p) < 0)
            return false;
        for (std::size_t i = 1; i < rings.size(); ++i)
            if (inRing(rings[i], p) > 0)
                return false;
        return true;
    }
};
class Solid final : public GSolid {
  public:
    std::vector<const Face *> faces;
    Vec3 lo, hi;
    double tolerance = 0;
    explicit Solid(int tag) : GSolid(tag) {}
    std::vector<int> boundaryFaceTags() const override {
        std::vector<int> out;
        out.reserve(faces.size());
        for (const auto *f : faces)
            out.push_back(f->tag());
        return out;
    }
    void bbox(Vec3 &a, Vec3 &b) const override {
        a = lo;
        b = hi;
    }
    bool isInside(const Vec3 &p) const override {
        if (!finite(p) || p.x < lo.x - tolerance || p.y < lo.y - tolerance ||
            p.z < lo.z - tolerance || p.x > hi.x + tolerance || p.y > hi.y + tolerance ||
            p.z > hi.z + tolerance)
            return false;
        for (const auto *f : faces)
            if (std::abs((p - f->origin).dot(f->normal)) <= tolerance &&
                f->containsParam(f->paramAt(p)))
                return true;
        // Retry directions if a ray hits a topological edge. Oriented crossings
        // naturally handle non-convex bodies; this is a geometric query, not
        // an exact volumetric-intersection certificate.
        const std::array<Vec3, 3> rays{{{1, .3713906763541037, .5291502622129182},
                                        {.1936491673103709, 1, .4170288281141495},
                                        {.6139406135149205, .2718281828459045, 1}}};
        for (const auto ray : rays) {
            int winding = 0;
            bool ambiguous = false;
            for (const auto *f : faces) {
                const double den = ray.dot(f->normal);
                if (std::abs(den) < 1e-14)
                    continue;
                const double t = (f->origin - p).dot(f->normal) / den;
                if (t <= 0)
                    continue;
                const auto uv = f->paramAt(p + ray * t);
                bool boundary = false;
                for (const auto &ring : f->rings)
                    boundary |= inRing(ring, uv) == 0;
                if (boundary) {
                    ambiguous = true;
                    break;
                }
                if (f->containsParam(uv))
                    winding += den > 0 ? 1 : -1;
            }
            if (!ambiguous)
                return winding != 0;
        }
        return false;
    }
};
class OwnedModel final : public Model {
  public:
    std::vector<std::unique_ptr<GVertex>> vs;
    std::vector<std::unique_ptr<Line>> es;
    std::vector<std::unique_ptr<Face>> fs;
    std::unique_ptr<Solid> body;
    template <class T, class U>
    static std::vector<const T *> ptrs(const std::vector<std::unique_ptr<U>> &in) {
        std::vector<const T *> out;
        out.reserve(in.size());
        for (const auto &p : in)
            out.push_back(p.get());
        return out;
    }
    std::vector<const GVertex *> vertices() const override { return ptrs<GVertex>(vs); }
    std::vector<const GEdge *> edges() const override { return ptrs<GEdge>(es); }
    std::vector<const GFace *> faces() const override { return ptrs<GFace>(fs); }
    std::vector<const GSolid *> solids() const override { return {body.get()}; }
    const GVertex *vertex(int t) const override {
        return t > 0 && std::size_t(t) <= vs.size() ? vs[std::size_t(t - 1)].get() : nullptr;
    }
    const GEdge *edge(int t) const override {
        return t > 0 && std::size_t(t) <= es.size() ? es[std::size_t(t - 1)].get() : nullptr;
    }
    const GFace *face(int t) const override {
        return t > 0 && std::size_t(t) <= fs.size() ? fs[std::size_t(t - 1)].get() : nullptr;
    }
    const GSolid *solid(int t) const override {
        return body && body->tag() == t ? body.get() : nullptr;
    }
};
} // namespace

PolyhedralBuild makePolyhedralModel(const PolyhedralInput &input,
                                    const PolyhedralLimits &limits) noexcept {
    PolyhedralBuild out;
    Work work{0, limits.maximumWork};
    try {
        if (input.vertices.size() < 4 || input.faces.size() < 4 || input.solidTag < 0)
            throw Stop{PolyhedralReason::Empty};
        if (input.vertices.size() > limits.maximumVertices ||
            input.faces.size() > limits.maximumFaces ||
            input.vertices.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
            input.faces.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw Stop{PolyhedralReason::BudgetExceeded};
        if (!(limits.relativeTolerance > 0 && limits.relativeTolerance < 1e-3))
            throw Stop{PolyhedralReason::Degenerate};
        Vec3 lo = input.vertices.front(), hi = lo;
        for (const auto p : input.vertices) {
            work.tick();
            if (!finite(p))
                throw Stop{PolyhedralReason::Nonfinite};
            lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
            hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
        }
        const double scale = length(hi - lo);
        if (!(scale > 0) || !std::isfinite(scale))
            throw Stop{PolyhedralReason::Degenerate};
        const double coordinate = std::max({std::abs(lo.x), std::abs(lo.y), std::abs(lo.z),
                                            std::abs(hi.x), std::abs(hi.y), std::abs(hi.z)});
        const double tolerance = scale * limits.relativeTolerance +
                                 coordinate * 32 * std::numeric_limits<double>::epsilon();
        auto model = std::make_unique<OwnedModel>();
        for (std::size_t i = 0; i < input.vertices.size(); ++i)
            model->vs.push_back(
                std::make_unique<GVertex>(static_cast<int>(i + 1), input.vertices[i]));
        struct Use {
            std::size_t face;
            std::uint32_t a, b;
        };
        std::map<EdgeKey, std::vector<Use>> uses;
        std::vector<bool> usedVertices(input.vertices.size(), false);
        long double volume = 0;
        for (std::size_t fi = 0; fi < input.faces.size(); ++fi) {
            const auto &in = input.faces[fi];
            if (in.loops.empty() || in.loops.front().size() < 3)
                throw Stop{PolyhedralReason::InvalidLoop};
            for (const auto &loop : in.loops) {
                if (loop.size() < 3)
                    throw Stop{PolyhedralReason::InvalidLoop};
                std::set<std::uint32_t> ids;
                for (auto id : loop) {
                    work.tick();
                    if (id >= input.vertices.size())
                        throw Stop{PolyhedralReason::InvalidIndex};
                    if (!ids.insert(id).second)
                        throw Stop{PolyhedralReason::InvalidLoop};
                    usedVertices[id] = true;
                }
            }
            auto face = std::make_unique<Face>(static_cast<int>(fi + 1));
            const auto &ring = in.loops.front();
            face->origin = input.vertices[ring.front()];
            auto u = input.vertices[ring[1]] - face->origin;
            if (length(u) <= tolerance)
                throw Stop{PolyhedralReason::Degenerate};
            u = u / length(u);
            Vec3 normal{};
            for (std::size_t k = 2; k < ring.size(); ++k) {
                work.tick();
                const auto n = u.cross((input.vertices[ring[k]] - face->origin) / scale);
                if (length(n) > limits.relativeTolerance) {
                    normal = n / length(n);
                    break;
                }
            }
            if (length(normal) == 0)
                throw Stop{PolyhedralReason::Degenerate};
            face->u = u;
            face->v = normal.cross(u);
            face->normal = normal;
            bool first = true;
            for (const auto &loop : in.loops) {
                std::vector<Vec2> uv;
                uv.reserve(loop.size());
                for (std::size_t j = 0; j < loop.size(); ++j) {
                    work.tick();
                    const auto a = loop[j], b = loop[(j + 1) % loop.size()];
                    const Vec3 p = input.vertices[a];
                    if (length(input.vertices[b] - p) <= tolerance)
                        throw Stop{PolyhedralReason::Degenerate};
                    if (std::abs((p - face->origin).dot(normal)) > tolerance)
                        throw Stop{PolyhedralReason::Nonplanar};
                    const auto q = face->paramAt(p);
                    uv.push_back(q);
                    if (first) {
                        face->lo = face->hi = q;
                        first = false;
                    } else {
                        face->lo = {std::min(face->lo.u, q.u), std::min(face->lo.v, q.v)};
                        face->hi = {std::max(face->hi.u, q.u), std::max(face->hi.v, q.v)};
                    }
                    uses[edgeKey(a, b)].push_back({fi, a, b});
                }
                face->rings.push_back(std::move(uv));
            }
            const auto outerArea = signedArea(face->rings.front());
            if (outerArea == 0)
                throw Stop{PolyhedralReason::Degenerate};
            // Normal follows the outer ring, even for concave starts.
            if (outerArea < 0) {
                face->normal = face->normal * -1;
                face->v = face->v * -1;
                for (auto &r : face->rings)
                    for (auto &p : r)
                        p.v = -p.v;
                const double old = face->lo.v;
                face->lo.v = -face->hi.v;
                face->hi.v = -old;
            }
            for (std::size_t r = 0; r < face->rings.size(); ++r) {
                const auto &a = face->rings[r];
                const auto area = signedArea(a);
                if (r && (!(area < 0) || inRing(face->rings.front(), a.front()) != 1))
                    throw Stop{PolyhedralReason::InvalidLoop};
                if (r)
                    for (std::size_t s = 1; s < r; ++s) {
                        work.tick();
                        if (inRing(face->rings[s], a.front()) >= 0 ||
                            inRing(a, face->rings[s].front()) >= 0)
                            throw Stop{PolyhedralReason::InvalidLoop};
                    }
                for (std::size_t i = 0; i < a.size(); ++i)
                    for (std::size_t s = r; s < face->rings.size(); ++s) {
                        const auto &b = face->rings[s];
                        for (std::size_t j = 0; j < b.size(); ++j) {
                            if (r == s && (j <= i || j == i + 1 || (i == 0 && j + 1 == a.size())))
                                continue;
                            work.tick();
                            if (intersects(a[i], a[(i + 1) % a.size()], b[j],
                                           b[(j + 1) % b.size()]))
                                throw Stop{PolyhedralReason::InvalidLoop};
                        }
                    }
                volume +=
                    area * static_cast<long double>((face->origin - lo).dot(face->normal)) / 3;
            }
            model->fs.push_back(std::move(face));
        }
        if (std::find(usedVertices.begin(), usedVertices.end(), false) != usedVertices.end())
            throw Stop{PolyhedralReason::InvalidIndex};
        std::vector<std::vector<std::size_t>> adjacent(input.faces.size());
        for (const auto &[key, rows] : uses) {
            work.tick();
            if (rows.size() != 2)
                throw Stop{PolyhedralReason::OpenOrNonmanifold};
            if (rows[0].a != rows[1].b || rows[0].b != rows[1].a)
                throw Stop{PolyhedralReason::InconsistentOrientation};
            adjacent[rows[0].face].push_back(rows[1].face);
            adjacent[rows[1].face].push_back(rows[0].face);
            model->es.push_back(std::make_unique<Line>(static_cast<int>(model->es.size() + 1),
                                                       model->vs[key.first].get(),
                                                       model->vs[key.second].get()));
        }
        std::vector<bool> visited(input.faces.size());
        std::vector<std::size_t> todo{0};
        visited[0] = true;
        for (std::size_t i = 0; i < todo.size(); ++i)
            for (auto j : adjacent[todo[i]])
                if (!visited[j]) {
                    visited[j] = true;
                    todo.push_back(j);
                }
        if (todo.size() != input.faces.size())
            throw Stop{PolyhedralReason::DisconnectedShell};
        if (volume == 0 || !std::isfinite(static_cast<double>(volume)))
            throw Stop{PolyhedralReason::Degenerate};
        std::map<EdgeKey, int> tags;
        int next = 1;
        for (const auto &[k, ignored] : uses) {
            (void)ignored;
            tags.emplace(k, next++);
        }
        model->body = std::make_unique<Solid>(input.solidTag);
        model->body->lo = lo;
        model->body->hi = hi;
        model->body->tolerance = tolerance;
        for (std::size_t fi = 0; fi < model->fs.size(); ++fi) {
            auto &face = *model->fs[fi];
            for (const auto &ring : input.faces[fi].loops) {
                std::vector<int> e;
                for (std::size_t i = 0; i < ring.size(); ++i)
                    e.push_back(tags.at(edgeKey(ring[i], ring[(i + 1) % ring.size()])));
                face.edges.push_back(std::move(e));
            }
            if (volume < 0)
                face.normal = face.normal * -1;
            model->body->faces.push_back(&face);
        }
        out.volume = std::abs(static_cast<double>(volume));
        out.model = std::move(model);
        out.reason = PolyhedralReason::None;
    } catch (const Stop &s) {
        out.reason = s.reason;
    } catch (const std::bad_alloc &) {
        out.reason = PolyhedralReason::AllocationFailed;
    } catch (...) {
        out.reason = PolyhedralReason::InternalError;
    }
    out.work = work.count;
    return out;
}

PolyhedralInput extrudedPolygon(std::span<const std::vector<Vec2>> loops, Vec3 sweep, Vec3 origin,
                                Vec3 axisU, Vec3 axisV) {
    PolyhedralInput out;
    std::vector<std::vector<std::uint32_t>> rings;
    for (std::size_t k = 0; k < loops.size(); ++k) {
        auto p = loops[k];
        const auto a = signedArea(p);
        if ((k == 0 && a < 0) || (k > 0 && a > 0))
            std::reverse(p.begin(), p.end());
        std::vector<std::uint32_t> ring;
        for (const auto uv : p) {
            if (out.vertices.size() >= std::numeric_limits<std::uint32_t>::max() / 2U)
                throw std::length_error("extrudedPolygon: vertex count");
            ring.push_back(static_cast<std::uint32_t>(out.vertices.size()));
            out.vertices.push_back(origin + axisU * uv.u + axisV * uv.v);
        }
        rings.push_back(std::move(ring));
    }
    const auto n = static_cast<std::uint32_t>(out.vertices.size());
    for (std::uint32_t i = 0; i < n; ++i)
        out.vertices.push_back(out.vertices[i] + sweep);
    PolyhedralFace bottom, top;
    for (const auto &ring : rings) {
        auto b = ring;
        std::reverse(b.begin(), b.end());
        bottom.loops.push_back(std::move(b));
        auto t = ring;
        for (auto &i : t)
            i += n;
        top.loops.push_back(std::move(t));
    }
    out.faces.push_back(std::move(bottom));
    out.faces.push_back(std::move(top));
    for (const auto &ring : rings)
        for (std::size_t i = 0; i < ring.size(); ++i) {
            const auto a = ring[i], b = ring[(i + 1) % ring.size()];
            out.faces.push_back({{{a, b, b + n, a + n}}});
        }
    return out;
}
const char *polyhedralReasonName(PolyhedralReason r) noexcept {
    switch (r) {
    case PolyhedralReason::None:
        return "none";
    case PolyhedralReason::Empty:
        return "empty";
    case PolyhedralReason::InvalidIndex:
        return "invalid_index";
    case PolyhedralReason::Nonfinite:
        return "nonfinite";
    case PolyhedralReason::Degenerate:
        return "degenerate";
    case PolyhedralReason::Nonplanar:
        return "nonplanar";
    case PolyhedralReason::InvalidLoop:
        return "invalid_loop";
    case PolyhedralReason::OpenOrNonmanifold:
        return "open_or_nonmanifold";
    case PolyhedralReason::InconsistentOrientation:
        return "inconsistent_orientation";
    case PolyhedralReason::DisconnectedShell:
        return "disconnected_shell";
    case PolyhedralReason::BudgetExceeded:
        return "budget_exceeded";
    case PolyhedralReason::AllocationFailed:
        return "allocation_failed";
    case PolyhedralReason::InternalError:
        return "internal_error";
    }
    return "unknown";
}
} // namespace paacs::geo
