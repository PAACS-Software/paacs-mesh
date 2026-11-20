// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/Kernel.hpp"

#include <map>
#include <stdexcept>

namespace paacs::geo {
namespace {

// ------------------------------------------------------------ line edge ---
class LineEdge final : public GEdge {
public:
    CurveKind kind() const noexcept override { return CurveKind::Line; }
    LineEdge(int tag, const GVertex* a, const GVertex* b)
        : GEdge(tag), a_(a), b_(b) {}
    Vec3 pointAt(double t) const override {
        return a_->position() + (b_->position() - a_->position()) * t;
    }
    Vec3 tangentAt(double) const override {
        return (b_->position() - a_->position()).normalized();
    }
    double length() const override {
        return (b_->position() - a_->position()).norm();
    }
    int startVertexTag() const override { return a_->tag(); }
    int endVertexTag()   const override { return b_->tag(); }
private:
    const GVertex *a_, *b_;
};

// ---------------------------------------------------------- planar face ---
// Plane through origin_ spanned by orthonormal du_, dv_ with extents.
class PlanarFace final : public GFace {
public:
    SurfaceKind kind() const noexcept override { return SurfaceKind::Plane; }
    PlanarFace(int tag, Vec3 origin, Vec3 du, Vec3 dv, Vec2 lo, Vec2 hi,
               std::vector<int> loop)
        : GFace(tag), o_(origin), du_(du), dv_(dv), lo_(lo), hi_(hi),
          loop_(std::move(loop)) {}
    Vec3 pointAt(const Vec2& uv) const override {
        return o_ + du_ * uv.u + dv_ * uv.v;
    }
    Vec2 paramAt(const Vec3& p) const override {
        Vec3 d = p - o_;
        return {d.dot(du_), d.dot(dv_)};
    }
    Vec3 normalAt(const Vec2&) const override { return du_.cross(dv_).normalized(); }
    void paramBounds(Vec2& lo, Vec2& hi) const override { lo = lo_; hi = hi_; }
    std::vector<int> boundaryEdgeTags() const override { return loop_; }
    double curvatureAt(const Vec2&) const override { return 0.0; }
private:
    Vec3 o_, du_, dv_;
    Vec2 lo_, hi_;
    std::vector<int> loop_;
};

// ------------------------------------------------------------ box solid ---
class BoxSolid final : public GSolid {
public:
    BoxSolid(int tag, Vec3 lo, Vec3 hi, std::vector<int> faces)
        : GSolid(tag), lo_(lo), hi_(hi), faces_(std::move(faces)) {}
    bool isInside(const Vec3& p) const override {
        return p.x >= lo_.x && p.x <= hi_.x &&
               p.y >= lo_.y && p.y <= hi_.y &&
               p.z >= lo_.z && p.z <= hi_.z;
    }
    std::vector<int> boundaryFaceTags() const override { return faces_; }
    void bbox(Vec3& lo, Vec3& hi) const override { lo = lo_; hi = hi_; }
private:
    Vec3 lo_, hi_;
    std::vector<int> faces_;
};

constexpr double kPi = 3.14159265358979323846;

// ------------------------------------------------- circular arc (z-plane) --
class ArcEdge final : public GEdge {
public:
    CurveKind kind() const noexcept override { return CurveKind::Circle; }
    ArcEdge(int tag, Vec3 center, double R, double th0, double th1,
            int vStart, int vEnd)
        : GEdge(tag), c_(center), R_(R), th0_(th0), th1_(th1),
          vs_(vStart), ve_(vEnd) {}
    Vec3 pointAt(double t) const override {
        const double th = th0_ + t * (th1_ - th0_);
        return {c_.x + R_ * std::cos(th), c_.y + R_ * std::sin(th), c_.z};
    }
    Vec3 tangentAt(double t) const override {
        const double th = th0_ + t * (th1_ - th0_);
        const double s = (th1_ > th0_) ? 1.0 : -1.0;
        return Vec3{-std::sin(th) * s, std::cos(th) * s, 0};
    }
    double length()  const override { return R_ * std::abs(th1_ - th0_); }
    double curvatureAt(double) const override { return 1.0 / R_; }
    int startVertexTag() const override { return vs_; }
    int endVertexTag()   const override { return ve_; }
private:
    Vec3 c_; double R_, th0_, th1_; int vs_, ve_;
};

// -------------------------------------------- half-cylinder lateral face ---
// Parameterized (u = angle in [th0, th1], v = height in [0, H]).
class CylFace final : public GFace {
public:
    SurfaceKind kind() const noexcept override { return SurfaceKind::Cylinder; }
    CylFace(int tag, Vec3 base, double R, double H, double th0, double th1,
            std::vector<int> loop)
        : GFace(tag), b_(base), R_(R), H_(H), th0_(th0), th1_(th1),
          loop_(std::move(loop)) {}
    Vec3 pointAt(const Vec2& uv) const override {
        return {b_.x + R_ * std::cos(uv.u), b_.y + R_ * std::sin(uv.u),
                b_.z + uv.v};
    }
    Vec2 paramAt(const Vec3& p) const override {
        double th = std::atan2(p.y - b_.y, p.x - b_.x);
        if (th < th0_ - 1e-9) th += 2 * kPi;   // wrap onto this face's range
        return {th, p.z - b_.z};
    }
    Vec3 normalAt(const Vec2& uv) const override {
        return {std::cos(uv.u), std::sin(uv.u), 0};
    }
    void paramBounds(Vec2& lo, Vec2& hi) const override {
        lo = {th0_, 0}; hi = {th1_, H_};
    }
    std::vector<int> boundaryEdgeTags() const override { return loop_; }
    double curvatureAt(const Vec2&) const override { return 1.0 / R_; }
    Vec2 metricScale() const override { return {R_, 1.0}; }
private:
    Vec3 b_; double R_, H_, th0_, th1_;
    std::vector<int> loop_;
};

// ------------------------------------------------------------- disk face ---
class DiskFace final : public GFace {
public:
    SurfaceKind kind() const noexcept override { return SurfaceKind::Plane; }
    DiskFace(int tag, Vec3 center, double R, double nz, std::vector<int> loop)
        : GFace(tag), c_(center), R_(R), nz_(nz), loop_(std::move(loop)) {}
    Vec3 pointAt(const Vec2& uv) const override {
        return {c_.x + uv.u, c_.y + uv.v, c_.z};
    }
    Vec2 paramAt(const Vec3& p) const override {
        return {p.x - c_.x, p.y - c_.y};
    }
    Vec3 normalAt(const Vec2&) const override { return {0, 0, nz_}; }
    void paramBounds(Vec2& lo, Vec2& hi) const override {
        lo = {-R_, -R_}; hi = {R_, R_};
    }
    std::vector<int> boundaryEdgeTags() const override { return loop_; }
    double curvatureAt(const Vec2&) const override { return 0.0; }
    bool containsParam(const Vec2& uv) const override {
        return uv.u * uv.u + uv.v * uv.v <= R_ * R_ * (1 + 1e-9);
    }
private:
    Vec3 c_; double R_, nz_;
    std::vector<int> loop_;
};

// -------------------------------------------------------- cylinder solid ---
class CylSolid final : public GSolid {
public:
    CylSolid(int tag, Vec3 base, double R, double H, std::vector<int> faces)
        : GSolid(tag), b_(base), R_(R), H_(H), faces_(std::move(faces)) {}
    bool isInside(const Vec3& p) const override {
        const double dx = p.x - b_.x, dy = p.y - b_.y;
        return dx * dx + dy * dy <= R_ * R_ * (1 + 1e-12) &&
               p.z >= b_.z && p.z <= b_.z + H_;
    }
    std::vector<int> boundaryFaceTags() const override { return faces_; }
    void bbox(Vec3& lo, Vec3& hi) const override {
        lo = {b_.x - R_, b_.y - R_, b_.z};
        hi = {b_.x + R_, b_.y + R_, b_.z + H_};
    }
private:
    Vec3 b_; double R_, H_;
    std::vector<int> faces_;
};

// ---------------------------------------------------------------- model ---
class BuiltinModel final : public Model {
public:
    std::vector<const GVertex*> vertices() const override { return ptrs(verts_); }
    std::vector<const GEdge*>   edges()    const override { return ptrs(edges_); }
    std::vector<const GFace*>   faces()    const override { return ptrs(faces_); }
    std::vector<const GSolid*>  solids()   const override { return ptrs(solids_); }

    const GVertex* vertex(int t) const override { return find(verts_, t); }
    const GEdge*   edge(int t)   const override { return find(edges_, t); }
    const GFace*   face(int t)   const override { return find(faces_, t); }
    const GSolid*  solid(int t)  const override { return find(solids_, t); }

    std::vector<std::unique_ptr<GVertex>> verts_;
    std::vector<std::unique_ptr<GEdge>>   edges_;
    std::vector<std::unique_ptr<GFace>>   faces_;
    std::vector<std::unique_ptr<GSolid>>  solids_;

private:
    template <class V>
    static std::vector<const typename V::value_type::element_type*>
    ptrs(const V& v) {
        std::vector<const typename V::value_type::element_type*> out;
        out.reserve(v.size());
        for (const auto& e : v) out.push_back(e.get());
        return out;
    }
    template <class V>
    static const typename V::value_type::element_type*
    find(const V& v, int tag) {
        for (const auto& e : v) if (e->tag() == tag) return e.get();
        return nullptr;
    }
};

// --------------------------------------------------------------- kernel ---
class BuiltinKernel final : public Kernel {
public:
    std::unique_ptr<Model> makeBox(const Vec3& lo, const Vec3& hi) override {
        auto m = std::make_unique<BuiltinModel>();

        // Eight vertices, tagged 1 through 8 in box corner order.
        const Vec3 c[8] = {
            {lo.x, lo.y, lo.z}, {hi.x, lo.y, lo.z}, {hi.x, hi.y, lo.z},
            {lo.x, hi.y, lo.z}, {lo.x, lo.y, hi.z}, {hi.x, lo.y, hi.z},
            {hi.x, hi.y, hi.z}, {lo.x, hi.y, hi.z}};
        for (int i = 0; i < 8; ++i)
            m->verts_.push_back(std::make_unique<GVertex>(i + 1, c[i]));

        // 12 edges, tags 1..12.
        const int ev[12][2] = {{1,2},{2,3},{3,4},{4,1},   // bottom
                               {5,6},{6,7},{7,8},{8,5},   // top
                               {1,5},{2,6},{3,7},{4,8}};  // verticals
        for (int i = 0; i < 12; ++i)
            m->edges_.push_back(std::make_unique<LineEdge>(
                i + 1, m->vertex(ev[i][0]), m->vertex(ev[i][1])));

        // 6 planar faces, tags 1..6: -z +z -y +y +x -x
        auto addFace = [&](int tag, Vec3 origin, Vec3 du, Vec3 dv,
                           Vec2 plo, Vec2 phi, std::vector<int> loop) {
            m->faces_.push_back(std::make_unique<PlanarFace>(
                tag, origin, du, dv, plo, phi, std::move(loop)));
        };
        Vec3 d = hi - lo;
        addFace(1, lo,                {1,0,0}, {0,1,0}, {0,0}, {d.x,d.y}, {1,2,3,4});
        addFace(2, {lo.x,lo.y,hi.z},  {1,0,0}, {0,1,0}, {0,0}, {d.x,d.y}, {5,6,7,8});
        addFace(3, lo,                {1,0,0}, {0,0,1}, {0,0}, {d.x,d.z}, {1,10,5,9});
        addFace(4, {lo.x,hi.y,lo.z},  {1,0,0}, {0,0,1}, {0,0}, {d.x,d.z}, {3,11,7,12});
        addFace(5, {hi.x,lo.y,lo.z},  {0,1,0}, {0,0,1}, {0,0}, {d.y,d.z}, {2,11,6,10});
        addFace(6, lo,                {0,1,0}, {0,0,1}, {0,0}, {d.y,d.z}, {4,12,8,9});

        m->solids_.push_back(std::make_unique<BoxSolid>(
            1, lo, hi, std::vector<int>{1,2,3,4,5,6}));
        return m;
    }

    std::unique_ptr<Model> makeCylinder(const Vec3& base, double R,
                                        double H) override {
        auto m = std::make_unique<BuiltinModel>();
        const Vec3 top = base + Vec3{0, 0, H};

        // 4 vertices at the parametric seams (theta = 0 and pi).
        m->verts_.push_back(std::make_unique<GVertex>(1, base + Vec3{R, 0, 0}));
        m->verts_.push_back(std::make_unique<GVertex>(2, base + Vec3{-R, 0, 0}));
        m->verts_.push_back(std::make_unique<GVertex>(3, top + Vec3{R, 0, 0}));
        m->verts_.push_back(std::make_unique<GVertex>(4, top + Vec3{-R, 0, 0}));

        // 4 half-circle arcs + 2 vertical seam lines.
        m->edges_.push_back(std::make_unique<ArcEdge>(1, base, R, 0, kPi, 1, 2));
        m->edges_.push_back(std::make_unique<ArcEdge>(2, base, R, kPi, 2 * kPi, 2, 1));
        m->edges_.push_back(std::make_unique<ArcEdge>(3, top, R, 0, kPi, 3, 4));
        m->edges_.push_back(std::make_unique<ArcEdge>(4, top, R, kPi, 2 * kPi, 4, 3));
        m->edges_.push_back(std::make_unique<LineEdge>(5, m->vertex(1), m->vertex(3)));
        m->edges_.push_back(std::make_unique<LineEdge>(6, m->vertex(2), m->vertex(4)));

        // 2 half-cylinder lateral faces + 2 disks.
        m->faces_.push_back(std::make_unique<CylFace>(
            1, base, R, H, 0, kPi, std::vector<int>{1, 6, 3, 5}));
        m->faces_.push_back(std::make_unique<CylFace>(
            2, base, R, H, kPi, 2 * kPi, std::vector<int>{2, 5, 4, 6}));
        m->faces_.push_back(std::make_unique<DiskFace>(
            3, base, R, -1.0, std::vector<int>{1, 2}));
        m->faces_.push_back(std::make_unique<DiskFace>(
            4, top, R, 1.0, std::vector<int>{3, 4}));

        m->solids_.push_back(std::make_unique<CylSolid>(
            1, base, R, H, std::vector<int>{1, 2, 3, 4}));
        return m;
    }

    std::unique_ptr<Model> importStep(const std::string&) override {
        return nullptr;  // STEP requires the OCCT backend.
    }
};

} // namespace

std::unique_ptr<Kernel> makeBuiltinKernel() {
    return std::make_unique<BuiltinKernel>();
}

} // namespace paacs::geo
