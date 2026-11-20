// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/geometry/Kernel.hpp"
#include <algorithm>

namespace geometry_test {
// Test-only adapter: immutable bodies with independent entity namespaces.
class CombinedModel final : public paacs::geo::Model {
    using Vec3 = paacs::Vec3;
    using Vec2 = paacs::Vec2;
    struct Edge final : paacs::geo::GEdge {
        const GEdge *base;
        int offset;
        Edge(const GEdge *p, int n) : GEdge(p->tag() + n), base(p), offset(n) {}
        paacs::geo::CurveKind kind() const noexcept override { return base->kind(); }
        Vec3 pointAt(double t) const override { return base->pointAt(t); }
        Vec3 tangentAt(double t) const override { return base->tangentAt(t); }
        double length() const override { return base->length(); }
        int startVertexTag() const override { return base->startVertexTag() + offset; }
        int endVertexTag() const override { return base->endVertexTag() + offset; }
    };
    struct Face final : paacs::geo::GFace {
        const GFace *base;
        int offset;
        Face(const GFace *p, int n) : GFace(p->tag() + n), base(p), offset(n) {}
        paacs::geo::SurfaceKind kind() const noexcept override { return base->kind(); }
        Vec3 pointAt(const Vec2 &p) const override { return base->pointAt(p); }
        Vec2 paramAt(const Vec3 &p) const override { return base->paramAt(p); }
        Vec3 normalAt(const Vec2 &p) const override { return base->normalAt(p); }
        void paramBounds(Vec2 &l, Vec2 &h) const override { base->paramBounds(l, h); }
        std::vector<int> boundaryEdgeTags() const override {
            auto e = base->boundaryEdgeTags();
            for (auto &i : e)
                i += offset;
            return e;
        }
        std::vector<std::vector<int>> boundaryLoops() const override {
            auto e = base->boundaryLoops();
            for (auto &r : e)
                for (auto &i : r)
                    i += offset;
            return e;
        }
        double curvatureAt(const Vec2 &p) const override { return base->curvatureAt(p); }
        bool containsParam(const Vec2 &p) const override { return base->containsParam(p); }
    };
    struct Solid final : paacs::geo::GSolid {
        const GSolid *base;
        int offset;
        Solid(const GSolid *p, int n) : GSolid(p->tag() + n), base(p), offset(n) {}
        bool isInside(const Vec3 &p) const override { return base->isInside(p); }
        std::vector<int> boundaryFaceTags() const override {
            auto e = base->boundaryFaceTags();
            for (auto &i : e)
                i += offset;
            return e;
        }
        void bbox(Vec3 &l, Vec3 &h) const override { base->bbox(l, h); }
    };
    std::vector<std::unique_ptr<paacs::geo::GVertex>> vs;
    std::vector<std::unique_ptr<Edge>> es;
    std::vector<std::unique_ptr<Face>> fs;
    std::vector<std::unique_ptr<Solid>> ss;
    template <class T, class U>
    static std::vector<const T *> refs(const std::vector<std::unique_ptr<U>> &v, bool reversed) {
        std::vector<const T *> r;
        for (const auto &p : v)
            r.push_back(p.get());
        if (reversed)
            std::reverse(r.begin(), r.end());
        return r;
    }
    template <class T> static const T *get(const std::vector<std::unique_ptr<T>> &v, int tag) {
        for (const auto &p : v)
            if (p->tag() == tag)
                return p.get();
        return nullptr;
    }

  public:
    bool reversed = false;
    void add(const paacs::geo::Model &model, int offset) {
        for (const auto *p : model.vertices())
            vs.push_back(std::make_unique<paacs::geo::GVertex>(p->tag() + offset, p->position()));
        for (const auto *p : model.edges())
            es.push_back(std::make_unique<Edge>(p, offset));
        for (const auto *p : model.faces())
            fs.push_back(std::make_unique<Face>(p, offset));
        for (const auto *p : model.solids())
            ss.push_back(std::make_unique<Solid>(p, offset));
    }
    std::vector<const paacs::geo::GVertex *> vertices() const override {
        return refs<paacs::geo::GVertex>(vs, reversed);
    }
    std::vector<const paacs::geo::GEdge *> edges() const override {
        return refs<paacs::geo::GEdge>(es, reversed);
    }
    std::vector<const paacs::geo::GFace *> faces() const override {
        return refs<paacs::geo::GFace>(fs, reversed);
    }
    std::vector<const paacs::geo::GSolid *> solids() const override {
        return refs<paacs::geo::GSolid>(ss, reversed);
    }
    const paacs::geo::GVertex *vertex(int t) const override { return get(vs, t); }
    const paacs::geo::GEdge *edge(int t) const override { return get(es, t); }
    const paacs::geo::GFace *face(int t) const override { return get(fs, t); }
    const paacs::geo::GSolid *solid(int t) const override { return get(ss, t); }
};
} // namespace geometry_test
