// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// Geometry kernel abstraction.
//
// Meshing algorithms never talk to a CAD kernel directly; they talk to these
// interfaces. Backends:
//   * BuiltinKernel  - analytic primitives (box, cylinder), zero dependencies.
//   * PolyhedralModel - validated planar loops, including polygonal holes.
// Optional PAACS_WITH_OCCT imports planar STEP/BREP solids into that backend.
//
// Every entity carries a persistent integer tag. Tags are the contract with a
// UI layer: selections, sizing controls, boundary-layer specs, and physical
// groups all reference entities by (dimension, tag).
//
#include <memory>
#include <string>
#include <vector>
#include "paacs/core/Types.hpp"

namespace paacs::geo {

// Analytic type declarations are supplied by the geometry backend. Unknown
// remains the safe default; sampling alone never certifies a straight/planar
// entity for a specialized mesher.
enum class CurveKind : std::uint8_t { Unknown, Line, Circle, Other };
enum class SurfaceKind : std::uint8_t { Unknown, Plane, Cylinder, Cone, Sphere, Torus, Other };

class Entity {
public:
    virtual ~Entity() = default;
    int tag() const { return tag_; }
    const std::string& name() const { return name_; }
    void setName(std::string n) { name_ = std::move(n); }
protected:
    explicit Entity(int tag) : tag_(tag) {}
private:
    int tag_;
    std::string name_;   // e.g. STEP entity name: "inlet", "bolt_hole_3"
};

// -------------------------------------------------------------- 0D vertex --
class GVertex : public Entity {
public:
    GVertex(int tag, Vec3 p) : Entity(tag), pos_(p) {}
    const Vec3& position() const { return pos_; }
private:
    Vec3 pos_;
};

// ---------------------------------------------------------------- 1D edge --
class GEdge : public Entity {
public:
    using Entity::Entity;
    virtual CurveKind kind() const noexcept { return CurveKind::Unknown; }
    virtual Vec3   pointAt(double t)   const = 0;  // t in [0,1]
    virtual Vec3   tangentAt(double t) const = 0;
    virtual double length()            const = 0;
    virtual int    startVertexTag()    const = 0;
    virtual int    endVertexTag()      const = 0;
    // Minimal capability boundary for future CAD adapters. Existing analytic
    // curves use the normalized [0,1] interval; an OCCT adapter may expose a
    // native interval while the mesher remains geometry-kernel independent.
    virtual void parameterBounds(double& lo, double& hi) const {
        lo = 0.0;
        hi = 1.0;
    }
    virtual bool periodic() const { return false; }
    // Curve curvature magnitude at t (0 for straight lines).
    virtual double curvatureAt(double) const { return 0.0; }
    // Closest parameter on the curve. Default: sampled + local refinement,
    // adequate for meshing-tolerance placement; backends may override with a
    // certified projection.
    virtual double closestParameter(const Vec3& p) const {
        double lower = 0.0, upper = 1.0;
        parameterBounds(lower, upper);
        double bestT = lower, bestD = 1e300;
        for (int i = 0; i <= 64; ++i) {
            const double t = lower + (upper - lower) * (i / 64.0);
            const double d = (pointAt(t) - p).norm2();
            if (d < bestD) { bestD = d; bestT = t; }
        }
        const double step = (upper - lower) / 64.0;
        double lo = std::max(lower, bestT - step);
        double hi = std::min(upper, bestT + step);
        for (int it = 0; it < 30; ++it) {
            const double m1 = lo + (hi - lo) / 3.0;
            const double m2 = hi - (hi - lo) / 3.0;
            ((pointAt(m1) - p).norm2() < (pointAt(m2) - p).norm2())
                ? hi = m2
                : lo = m1;
        }
        return (lo + hi) * 0.5;
    }
    virtual Vec3 closestPoint(const Vec3& p) const {
        return pointAt(closestParameter(p));
    }
};

// ---------------------------------------------------------------- 2D face --
class GFace : public Entity {
public:
    using Entity::Entity;
    virtual SurfaceKind kind() const noexcept { return SurfaceKind::Unknown; }
    virtual Vec3 pointAt(const Vec2& uv)  const = 0;   // param -> 3D
    virtual Vec2 paramAt(const Vec3& p)   const = 0;   // 3D -> param (on-face)
    virtual Vec3 normalAt(const Vec2& uv) const = 0;
    virtual void paramBounds(Vec2& lo, Vec2& hi) const = 0;
    // Ordered outer boundary loop. Use boundaryLoops() to include holes.
    virtual std::vector<int> boundaryEdgeTags() const = 0;
    // Ordered rings: outer boundary first, followed by holes. Edge direction
    // is recovered from endpoint identity. Existing backends keep one ring.
    virtual std::vector<std::vector<int>> boundaryLoops() const {
        return {boundaryEdgeTags()};
    }
    // Principal curvature magnitude upper bound at uv (0 for planes).
    virtual double curvatureAt(const Vec2& uv) const = 0;
    // Approximate first-fundamental-form scales (|dP/du|, |dP/dv|) so meshers
    // can work in a near-isometric parametric metric (e.g. (R, 1) for a
    // cylinder parameterized by angle and height).
    virtual Vec2 metricScale() const { return {1.0, 1.0}; }
    virtual bool periodicU() const { return false; }
    virtual bool periodicV() const { return false; }
    virtual bool singularAt(const Vec2&) const { return false; }
    // Explicit projection capability used by surface optimization and future
    // CAD adapters. The default composes the existing inverse/evaluate pair.
    // Failed projections leave outputs unchanged. Allocation failures propagate
    // so callers can preserve their resource-exhaustion diagnostic.
    virtual bool projectPoint(const Vec3& point, Vec2& uv,
                              Vec3& projected) const {
        if (!std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(point.z))
            return false;
        try {
            const auto candidateUv = paramAt(point);
            if (!std::isfinite(candidateUv.u) || !std::isfinite(candidateUv.v) ||
                !containsParam(candidateUv))
                return false;
            const auto candidate = pointAt(candidateUv);
            if (!std::isfinite(candidate.x) || !std::isfinite(candidate.y) ||
                !std::isfinite(candidate.z))
                return false;
            uv = candidateUv;
            projected = candidate;
            return true;
        } catch (const std::bad_alloc&) {
            throw;
        } catch (...) {
            return false;
        }
    }
    // True if uv lies inside the face's trimmed region (e.g. inside a disk).
    virtual bool containsParam(const Vec2&) const { return true; }
};

// --------------------------------------------------------------- 3D solid --
class GSolid : public Entity {
public:
    using Entity::Entity;
    virtual bool isInside(const Vec3& p) const = 0;
    virtual std::vector<int> boundaryFaceTags() const = 0;
    virtual void bbox(Vec3& lo, Vec3& hi) const = 0;
};

// ------------------------------------------------------------------ model --
class Model {
public:
    virtual ~Model() = default;

    virtual std::vector<const GVertex*> vertices() const = 0;
    virtual std::vector<const GEdge*>   edges()    const = 0;
    virtual std::vector<const GFace*>   faces()    const = 0;
    virtual std::vector<const GSolid*>  solids()   const = 0;

    virtual const GVertex* vertex(int tag) const = 0;
    virtual const GEdge*   edge(int tag)   const = 0;
    virtual const GFace*   face(int tag)   const = 0;
    virtual const GSolid*  solid(int tag)  const = 0;
};

// A kernel produces models (from primitives or from CAD files).
class Kernel {
public:
    virtual ~Kernel() = default;
    virtual std::unique_ptr<Model> makeBox(const Vec3& lo, const Vec3& hi) = 0;
    // Z-aligned cylinder: base disk center, radius, height.
    virtual std::unique_ptr<Model> makeCylinder(const Vec3& baseCenter,
                                                double radius, double height) = 0;
    // STEP import: implemented by the OCCT backend; Builtin returns nullptr.
    virtual std::unique_ptr<Model> importStep(const std::string& path) = 0;
};

std::unique_ptr<Kernel> makeBuiltinKernel();
#ifdef PAACS_MESH_HAS_OCCT
std::unique_ptr<Kernel> makeOcctKernel();
#endif

} // namespace paacs::geo
