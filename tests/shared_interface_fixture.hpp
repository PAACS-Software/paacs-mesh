// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/geometry/Kernel.hpp"
#include <algorithm>
#include <array>
#include <map>
#include <memory>
#include <vector>

namespace geometry_test {
// Two tetrahedral solids use one identical, oppositely owned whole face.
class SharedFaceTetrahedraModel final : public paacs::geo::Model {
    using Vec3 = paacs::Vec3;
    using Vec2 = paacs::Vec2;
    struct Edge final : paacs::geo::GEdge {
        const paacs::geo::GVertex *a, *b;
        Edge(int tag, const paacs::geo::GVertex *first, const paacs::geo::GVertex *last)
            : GEdge(tag), a(first), b(last) {}
        paacs::geo::CurveKind kind() const noexcept override { return paacs::geo::CurveKind::Line; }
        Vec3 pointAt(double t) const override { return a->position() + (b->position() - a->position()) * t; }
        Vec3 tangentAt(double) const override { return (b->position() - a->position()).normalized(); }
        double length() const override { return (b->position() - a->position()).norm(); }
        int startVertexTag() const override { return a->tag(); }
        int endVertexTag() const override { return b->tag(); }
    };
    struct Face final : paacs::geo::GFace {
        Vec3 origin, u, v, first, last;
        Vec2 lower, upper;
        std::vector<int> edges;
        Face(int tag, Vec3 a, Vec3 b, Vec3 c, std::vector<int> loop)
            : GFace(tag), origin(a), first(b - a), last(c - a), edges(std::move(loop)) {
            // Put the UV bounding-box midpoint inside these isosceles faces,
            // allowing ownership classification to sample the actual interior.
            u = (b - c).normalized();
            v = first.cross(last).normalized().cross(u);
            const Vec2 x{first.dot(u), first.dot(v)}, y{last.dot(u), last.dot(v)};
            lower = {std::min({0.0, x.u, y.u}), std::min({0.0, x.v, y.v})};
            upper = {std::max({0.0, x.u, y.u}), std::max({0.0, x.v, y.v})};
        }
        paacs::geo::SurfaceKind kind() const noexcept override { return paacs::geo::SurfaceKind::Plane; }
        Vec3 pointAt(const Vec2 &p) const override { return origin + u * p.u + v * p.v; }
        Vec2 paramAt(const Vec3 &p) const override {
            const auto d = p - origin;
            return {d.dot(u), d.dot(v)};
        }
        Vec3 normalAt(const Vec2 &) const override { return u.cross(v).normalized(); }
        void paramBounds(Vec2 &lo, Vec2 &hi) const override { lo = lower; hi = upper; }
        std::vector<int> boundaryEdgeTags() const override { return edges; }
        double curvatureAt(const Vec2 &) const override { return 0; }
        Vec2 metricScale() const override { return {u.norm(), v.norm()}; }
        bool containsParam(const Vec2 &p) const override {
            const auto d = u * p.u + v * p.v;
            const double aa = first.norm2(), bb = last.norm2(), ab = first.dot(last);
            const double den = aa * bb - ab * ab;
            const double x = (d.dot(first) * bb - d.dot(last) * ab) / den;
            const double y = (d.dot(last) * aa - d.dot(first) * ab) / den;
            return x >= -1e-12 && y >= -1e-12 && x + y <= 1 + 1e-12;
        }
    };
    struct Solid final : paacs::geo::GSolid {
        bool upper;
        Solid(int tag, bool positive) : GSolid(tag), upper(positive) {}
        bool isInside(const Vec3 &p) const override {
            const double z = upper ? p.z : -p.z;
            return p.x >= 0 && p.y >= 0 && z >= 0 && p.x + p.y + z <= 1;
        }
        std::vector<int> boundaryFaceTags() const override {
            return upper ? std::vector<int>{4, 5, 6, 7} : std::vector<int>{1, 2, 3, 7};
        }
        void bbox(Vec3 &lo, Vec3 &hi) const override {
            lo = {0, 0, upper ? 0.0 : -1.0};
            hi = {1, 1, upper ? 1.0 : 0.0};
        }
    };
    std::vector<std::unique_ptr<paacs::geo::GVertex>> vertices_;
    std::vector<std::unique_ptr<Edge>> edges_;
    std::vector<std::unique_ptr<Face>> faces_;
    std::array<Solid, 2> solids_{{Solid{1, false}, Solid{2, true}}};
    std::map<std::array<unsigned, 2>, int> edgeTags_;
    void addFace(std::array<unsigned, 3> triangle) {
        std::vector<int> loop;
        for (unsigned i = 0; i < 3; ++i) {
            std::array<unsigned, 2> key{triangle[i], triangle[(i + 1) % 3]};
            std::sort(key.begin(), key.end());
            auto found = edgeTags_.find(key);
            if (found == edgeTags_.end()) {
                const int tag = static_cast<int>(edges_.size() + 1);
                edges_.push_back(std::make_unique<Edge>(tag, vertices_[key[0]].get(), vertices_[key[1]].get()));
                found = edgeTags_.emplace(key, tag).first;
            }
            loop.push_back(found->second);
        }
        faces_.push_back(std::make_unique<Face>(static_cast<int>(faces_.size() + 1),
            vertices_[triangle[0]]->position(), vertices_[triangle[1]]->position(),
            vertices_[triangle[2]]->position(), std::move(loop)));
    }
    template <class T, class U>
    static std::vector<const T *> refs(const std::vector<std::unique_ptr<U>> &owned) {
        std::vector<const T *> out;
        for (const auto &p : owned) out.push_back(p.get());
        return out;
    }
    template <class T>
    static const T *get(const std::vector<std::unique_ptr<T>> &owned, int tag) {
        for (const auto &p : owned) if (p->tag() == tag) return p.get();
        return nullptr;
    }
public:
    SharedFaceTetrahedraModel() {
        for (Vec3 p : {Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}, Vec3{0, 0, -1}})
            vertices_.push_back(std::make_unique<paacs::geo::GVertex>(static_cast<int>(vertices_.size() + 1), p));
        addFace({2, 1, 4}); addFace({0, 4, 1}); addFace({0, 2, 4});
        addFace({1, 2, 3}); addFace({0, 3, 2}); addFace({0, 1, 3});
        addFace({0, 1, 2});
    }
    std::vector<const paacs::geo::GVertex *> vertices() const override { return refs<paacs::geo::GVertex>(vertices_); }
    std::vector<const paacs::geo::GEdge *> edges() const override { return refs<paacs::geo::GEdge>(edges_); }
    std::vector<const paacs::geo::GFace *> faces() const override { return refs<paacs::geo::GFace>(faces_); }
    std::vector<const paacs::geo::GSolid *> solids() const override { return {&solids_[0], &solids_[1]}; }
    const paacs::geo::GVertex *vertex(int tag) const override { return get(vertices_, tag); }
    const paacs::geo::GEdge *edge(int tag) const override { return get(edges_, tag); }
    const paacs::geo::GFace *face(int tag) const override { return get(faces_, tag); }
    const paacs::geo::GSolid *solid(int tag) const override {
        return tag == 1 ? &solids_[0] : tag == 2 ? &solids_[1] : nullptr;
    }
};
} // namespace geometry_test
