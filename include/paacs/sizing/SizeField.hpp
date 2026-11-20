// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// Sizing field: a scalar field h(x) giving the target element size anywhere
// in space. All user mesh controls (global size, local refinement balls,
// per-entity sizes, curvature adaptation) compose into one field, and the
// meshers only ever query the composite. This is the layer a UI drives.
//
// Order-awareness: the field carries the target element order; a p=2 mesh is
// generated coarser by `orderFactor` since quadratic elements resolve the
// same geometry/solution with fewer elements.
//
#include <algorithm>
#include <memory>
#include <vector>
#include <vector>
#include "paacs/core/Types.hpp"

namespace paacs::sizing {

class SizeField {
public:
    virtual ~SizeField() = default;
    virtual double eval(const Vec3& p) const = 0;
};

class ConstantField final : public SizeField {
public:
    explicit ConstantField(double h) : h_(h) {}
    double eval(const Vec3&) const override { return h_; }
private:
    double h_;
};

// Refinement ball: size hInner inside radius r, blending linearly back to the
// background over a transition band (classic "sphere of influence" control).
class BallRefinement final : public SizeField {
public:
    BallRefinement(Vec3 c, double r, double hInner, double hOuter,
                   double transition = 1.0)
        : c_(c), r_(r), hi_(hInner), ho_(hOuter), band_(transition * r) {}
    double eval(const Vec3& p) const override {
        double d = (p - c_).norm();
        if (d <= r_) return hi_;
        if (d >= r_ + band_) return ho_;
        double t = (d - r_) / band_;
        return hi_ + t * (ho_ - hi_);
    }
private:
    Vec3 c_; double r_, hi_, ho_, band_;
};

// Composite: minimum of all children (refinement always wins).
class MinField final : public SizeField {
public:
    void add(std::shared_ptr<SizeField> f) { fields_.push_back(std::move(f)); }
    double eval(const Vec3& p) const override {
        double h = 1e300;
        for (const auto& f : fields_) h = std::min(h, f->eval(p));
        return h;
    }
    bool empty() const { return fields_.empty(); }
private:
    std::vector<std::shared_ptr<SizeField>> fields_;
};


// Uniform-grid sampled field with trilinear interpolation.
//
// Two purposes: (1) evaluating a deep composite field once per grid node
// instead of once per candidate point makes point placement O(1) per query;
// (2) a flat float array + trilinear lookup is exactly the representation
// GPU kernels consume - this class IS the accelerator handoff format
// (see accel/Backend.hpp). Baking is parallelized.
class BakedGridField final : public SizeField {
public:
    BakedGridField(const SizeField& src, const Vec3& lo, const Vec3& hi,
                   int nx, int ny, int nz);

    double eval(const Vec3& p) const override;

    // Raw access for accelerator backends.
    const std::vector<float>& data() const { return data_; }
    Vec3 gridLo() const { return lo_; }
    Vec3 gridHi() const { return hi_; }
    int  nx() const { return nx_; }
    int  ny() const { return ny_; }
    int  nz() const { return nz_; }

private:
    Vec3 lo_, hi_;
    int nx_, ny_, nz_;
    std::vector<float> data_;
};

// Top-level mesh specification handed to the meshers.
struct MeshSpec {
    std::shared_ptr<SizeField> field;
    int    order       = 1;    // 1 = linear, 2 = quadratic (tet4 -> tet10)
    double orderFactor = 1.6;  // coarsening factor applied when order == 2
    double gradation   = 1.3;  // max adjacent size ratio (reserved)

    // Curvature adaptation: cap element size so a full circle of curvature
    // kappa gets at least this many elements. 0 disables.
    int elementsPerCircle = 24;

    // Worker threads for all parallel phases. 0 = all hardware cores,
    // 1 = fully serial. Results are bit-identical for any value.
    int threads = 0;

    // Stopgap boundary-sliver peeling in the volume mesher (curved faces
    // only; see PM-ALG-002). Set false to guarantee the represented domain
    // is exactly the Delaunay-filtered volume, at the cost of possible
    // near-flat boundary tets on curved walls. Removed entirely once
    // constrained recovery (M10) lands.
    bool boundaryPeeling = true;

    double h(const Vec3& p) const {
        double base = field ? field->eval(p) : 1.0;
        return order >= 2 ? base * orderFactor : base;
    }
    // Geometry-aware size: field size capped by local curvature.
    double hGeom(const Vec3& p, double kappa) const {
        double hf = h(p);
        if (elementsPerCircle > 0 && kappa > 1e-12) {
            double n = static_cast<double>(elementsPerCircle);
            if (order >= 2) n /= orderFactor;   // curved p2 edges track arcs
            hf = std::min(hf, 2.0 * 3.14159265358979323846 / (n * kappa));
        }
        return hf;
    }
};

} // namespace paacs::sizing
