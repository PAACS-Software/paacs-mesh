// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctModelInternal.hpp"
#include "OcctQueryContext.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepCheck_Analyzer.hxx>
#include "OcctSolidClassifier.hpp"
#include <BRepClass_FaceClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GProp_GProps.hxx>
#include <Geom2d_Curve.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_Sequence.hxx>
#include <Precision.hxx>
#include <STEPControl_Reader.hxx>
#include <Standard_Failure.hxx>
#include <Standard_OutOfMemory.hxx>
#include <TCollection_AsciiString.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>

namespace paacs::geo {
namespace {
using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;
std::atomic<std::uint64_t> nextIdentity{1};
double nan() noexcept {
    return (std::numeric_limits<double>::quiet_NaN)();
}
bool finite(Vec3 p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
bool finite(Vec2 p) noexcept {
    return std::isfinite(p.u) && std::isfinite(p.v);
}
Vec3 vec(const gp_Pnt &p) noexcept {
    return {p.X(), p.Y(), p.Z()};
}
Vec3 vec(const gp_Vec &p) noexcept {
    return {p.X(), p.Y(), p.Z()};
}
Vec3 vec(const gp_Dir &p) noexcept {
    return {p.X(), p.Y(), p.Z()};
}
gp_Pnt point(Vec3 p) {
    return {p.x, p.y, p.z};
}
CadOrientation orientation(TopAbs_Orientation o) noexcept {
    switch (o) {
    case TopAbs_FORWARD:
        return CadOrientation::Forward;
    case TopAbs_REVERSED:
        return CadOrientation::Reversed;
    case TopAbs_INTERNAL:
        return CadOrientation::Internal;
    case TopAbs_EXTERNAL:
        return CadOrientation::External;
    }
    return CadOrientation::External;
}
int sense(TopAbs_Orientation o) noexcept {
    return o == TopAbs_FORWARD ? 1 : o == TopAbs_REVERSED ? -1 : 0;
}
struct Stop {
    CadImportReason reason;
};
void require(bool good, CadImportReason reason = CadImportReason::InvalidShape) {
    if (!good)
        throw Stop{reason};
}
template <class T> const T *lookup(const std::vector<std::unique_ptr<T>> &v, int tag) noexcept {
    return tag > 0 && static_cast<std::size_t>(tag) <= v.size() ? v[tag - 1].get() : nullptr;
}
template <class T> std::vector<const T *> pointers(const std::vector<std::unique_ptr<T>> &v) {
    std::vector<const T *> out;
    out.reserve(v.size());
    for (const auto &p : v)
        out.push_back(p.get());
    return out;
}
template <class T> std::span<const T> row(const std::vector<std::vector<T>> &v, int tag) noexcept {
    return tag > 0 && static_cast<std::size_t>(tag) <= v.size() ? std::span<const T>(v[tag - 1])
                                                                : std::span<const T>{};
}
class NativeEdge final : public GEdge {
  public:
    NativeEdge(int tag, TopoDS_Edge e, int start, int end, Vec3 p)
        : GEdge(tag), edge_(std::move(e)), start_(start), end_(end), pole_(p) {
        degenerated_ = BRep_Tool::Degenerated(edge_);
        if (degenerated_) {
            first_ = 0;
            last_ = 1;
            length_ = 0;
            return;
        }
        BRepAdaptor_Curve curve(edge_);
        first_ = curve.FirstParameter();
        last_ = curve.LastParameter();
        require(std::isfinite(first_) && std::isfinite(last_) && last_ > first_);
        kind_ = curve.GetType() == GeomAbs_Line     ? CurveKind::Line
                : curve.GetType() == GeomAbs_Circle ? CurveKind::Circle
                                                    : CurveKind::Other;
        periodic_ = curve.IsPeriodic();
        length_ = GCPnts_AbscissaPoint::Length(curve, first_, last_);
        require(std::isfinite(length_) && length_ > 0);
    }
    CurveKind kind() const noexcept override { return kind_; }
    Vec3 pointAt(double t) const override {
        if (!std::isfinite(t))
            return {nan(), nan(), nan()};
        return degenerated_ ? pole_ : vec(BRepAdaptor_Curve(edge_).Value(t));
    }
    Vec3 tangentAt(double t) const override {
        if (degenerated_ || !std::isfinite(t))
            return {nan(), nan(), nan()};
        BRepAdaptor_Curve curve(edge_);
        gp_Pnt p;
        gp_Vec d;
        curve.D1(t, p, d);
        const double m = d.Magnitude();
        return m > 0 && std::isfinite(m) ? vec(d / m) : Vec3{nan(), nan(), nan()};
    }
    double length() const override { return length_; }
    int startVertexTag() const override { return start_; }
    int endVertexTag() const override { return end_; }
    void parameterBounds(double &lo, double &hi) const override {
        lo = first_;
        hi = last_;
    }
    bool periodic() const override { return periodic_; }
    double curvatureAt(double t) const override {
        if (degenerated_)
            return 0;
        try {
            BRepAdaptor_Curve curve(edge_);
            gp_Pnt p;
            gp_Vec d, dd;
            curve.D2(t, p, d, dd);
            const double m = d.Magnitude();
            return m > 0 ? d.Crossed(dd).Magnitude() / (m * m * m) : nan();
        } catch (const Standard_OutOfMemory &) {
            throw std::bad_alloc{};
        } catch (const Standard_Failure &) {
            return nan();
        }
    }
    double closestParameter(const Vec3 &input) const override {
        if (!finite(input))
            return nan();
        if (degenerated_)
            return first_;
        try {
            TopLoc_Location location;
            double first, last;
            const auto curve = BRep_Tool::Curve(edge_, location, first, last);
            if (curve.IsNull())
                return nan();
            auto p = point(input);
            p.Transform(location.Transformation().Inverted());
            GeomAPI_ProjectPointOnCurve projection(p, curve, first, last);
            // Bounded extrema omit an endpoint when the unbounded projection
            // lies outside this edge. Endpoints are valid closest candidates.
            double bestParameter = first_, bestDistance = (pointAt(first_) - input).norm2();
            const auto consider = [&](double parameter) {
                const double distance = (pointAt(parameter) - input).norm2();
                if (distance < bestDistance ||
                    (distance == bestDistance && parameter < bestParameter)) {
                    bestDistance = distance;
                    bestParameter = parameter;
                }
            };
            consider(last_);
            for (int i = 1; i <= projection.NbPoints(); ++i)
                consider(projection.Parameter(i));
            return bestParameter;
        } catch (const Standard_OutOfMemory &) {
            throw std::bad_alloc{};
        } catch (const Standard_Failure &) {
            return nan();
        }
    }

  private:
    TopoDS_Edge edge_;
    int start_, end_;
    Vec3 pole_;
    double first_ = 0, last_ = 1, length_ = 0;
    bool degenerated_ = false, periodic_ = false;
    CurveKind kind_ = CurveKind::Other;
};
class NativeFace final : public GFace {
  public:
    NativeFace(int tag, TopoDS_Face face, std::vector<std::vector<int>> loops)
        : GFace(tag), face_(std::move(face)), loops_(std::move(loops)) {
        tolerance_ = std::max(BRep_Tool::Tolerance(face_), Precision::Confusion());
        BRepTools::UVBounds(face_, lo_.u, hi_.u, lo_.v, hi_.v);
        require(finite(lo_) && finite(hi_) && hi_.u > lo_.u && hi_.v > lo_.v);
        BRepAdaptor_Surface s(face_, false);
        switch (s.GetType()) {
        case GeomAbs_Plane:
            kind_ = SurfaceKind::Plane;
            break;
        case GeomAbs_Cylinder:
            kind_ = SurfaceKind::Cylinder;
            break;
        case GeomAbs_Cone:
            kind_ = SurfaceKind::Cone;
            break;
        case GeomAbs_Sphere:
            kind_ = SurfaceKind::Sphere;
            break;
        case GeomAbs_Torus:
            kind_ = SurfaceKind::Torus;
            break;
        default:
            kind_ = SurfaceKind::Other;
        }
    }
    SurfaceKind kind() const noexcept override { return kind_; }
    Vec3 pointAt(const Vec2 &uv) const override {
        if (!finite(uv))
            return {nan(), nan(), nan()};
        return vec(BRepAdaptor_Surface(face_, false).Value(uv.u, uv.v));
    }
    Vec2 paramAt(const Vec3 &p) const override {
        Vec2 uv;
        Vec3 mapped;
        return projectPoint(p, uv, mapped) ? uv : Vec2{nan(), nan()};
    }
    Vec3 normalAt(const Vec2 &uv) const override {
        try {
            BRepAdaptor_Surface s(face_, false);
            gp_Pnt p;
            gp_Vec du, dv;
            s.D1(uv.u, uv.v, p, du, dv);
            auto n = du.Crossed(dv);
            const double m = n.Magnitude();
            if (!(m > 0) || !std::isfinite(m))
                return {nan(), nan(), nan()};
            n /= m;
            if (face_.Orientation() == TopAbs_REVERSED)
                n.Reverse();
            return vec(n);
        } catch (const Standard_OutOfMemory &) {
            throw std::bad_alloc{};
        } catch (const Standard_Failure &) {
            return {nan(), nan(), nan()};
        }
    }
    void paramBounds(Vec2 &lo, Vec2 &hi) const override {
        lo = lo_;
        hi = hi_;
    }
    std::vector<int> boundaryEdgeTags() const override { return loops_.front(); }
    std::vector<std::vector<int>> boundaryLoops() const override { return loops_; }
    bool containsParam(const Vec2 &uv) const override {
        if (!finite(uv))
            return false;
        try {
            // FaceClassifier's 2D tolerance is in parameter coordinates.
            // Applying a physical distance directly would overaccept trims
            // on, for example, a large-radius cylindrical chart. Use the
            // more restrictive source-surface resolution for both axes.
            const BRepAdaptor_Surface surface(face_, false);
            const double uvTolerance =
                std::min(surface.UResolution(tolerance_), surface.VResolution(tolerance_));
            if (!(uvTolerance > 0) || !std::isfinite(uvTolerance))
                return false;
            BRepClass_FaceClassifier c(face_, gp_Pnt2d(uv.u, uv.v), uvTolerance);
            return c.State() == TopAbs_IN || c.State() == TopAbs_ON;
        } catch (const Standard_OutOfMemory &) {
            throw std::bad_alloc{};
        } catch (const Standard_Failure &) {
            return false;
        }
    }
    bool projectPoint(const Vec3 &input, Vec2 &uv, Vec3 &projected) const override {
        if (!finite(input))
            return false;
        try {
            TopLoc_Location location;
            const auto s = BRep_Tool::Surface(face_, location);
            if (s.IsNull())
                return false;
            auto p = point(input);
            p.Transform(location.Transformation().Inverted());
            std::optional<GeomAPI_ProjectPointOnSurf> freshProjection;
            GeomAPI_ProjectPointOnSurf* selectedProjection;
            if(auto* context=detail::OcctQueryContext::active())
                selectedProjection=&context->surfaceProjector(s,lo_.u,hi_.u,lo_.v,hi_.v,p);
            else {
                freshProjection.emplace(p,s,lo_.u,hi_.u,lo_.v,hi_.v);
                selectedProjection=&*freshProjection;
            }
            auto& projection=*selectedProjection;
            double best = (std::numeric_limits<double>::infinity)();
            Vec2 bestUv;
            Vec3 bestPoint;
            bool found = false, ambiguous = false;
            const auto canonical = [](double t, double lo, double hi, bool periodic,
                                      double period) {
                if (!periodic || !(period > 0))
                    return t;
                double result = lo + std::fmod(t - lo, period);
                if (result < lo)
                    result += period;
                return result > hi && t >= lo && t <= hi ? t : result;
            };
            for (int i = 1; i <= projection.NbPoints(); ++i) {
                Vec2 candidate;
                projection.Parameters(i, candidate.u, candidate.v);
                candidate.u = canonical(candidate.u, lo_.u, hi_.u, s->IsUPeriodic(),
                                        s->IsUPeriodic() ? s->UPeriod() : 0);
                candidate.v = canonical(candidate.v, lo_.v, hi_.v, s->IsVPeriodic(),
                                        s->IsVPeriodic() ? s->VPeriod() : 0);
                // Select the nearest underlying stationary branch before
                // applying the trim. Skipping an exterior near branch can
                // otherwise choose a distant interior branch on a folded or
                // periodic surface and pull a query across a hole.
                if (!finite(candidate))
                    continue;
                const auto mapped = pointAt(candidate);
                if (!finite(mapped))
                    continue;
                const double distance = (mapped - input).norm();
                if (!std::isfinite(distance))
                    continue;
                if (!found || distance < best - tolerance_) {
                    best = distance;
                    bestUv = candidate;
                    bestPoint = mapped;
                    found = true;
                    ambiguous = false;
                } else if (std::abs(distance - best) <= tolerance_ &&
                           (mapped - bestPoint).norm() > tolerance_)
                    ambiguous = true;
            }
            if (!found || ambiguous || !containsParam(bestUv))
                return false;
            uv = bestUv;
            projected = bestPoint;
            return true;
        } catch (const Standard_OutOfMemory &) {
            throw std::bad_alloc{};
        } catch (const Standard_Failure &) {
            return false;
        }
    }
    double curvatureAt(const Vec2 &uv) const override {
        try {
            BRepAdaptor_Surface s(face_, false);
            gp_Pnt p;
            gp_Vec du, dv, duu, dvv, duv;
            s.D2(uv.u, uv.v, p, du, dv, duu, dvv, duv);
            const auto n = du.Crossed(dv);
            const double area = n.Magnitude();
            if (!(area > 0))
                return nan();
            const double e = du.Dot(du), f = du.Dot(dv), g = dv.Dot(dv), den = e * g - f * f;
            if (!(den > 0))
                return nan();
            const auto unit = n / area;
            const double l = unit.Dot(duu), m = unit.Dot(duv), r = unit.Dot(dvv);
            const double mean = (e * r - 2 * f * m + g * l) / (2 * den),
                         gauss = (l * r - m * m) / den;
            const double root = std::sqrt(std::max(0., mean * mean - gauss));
            return std::max(std::abs(mean + root), std::abs(mean - root));
        } catch (const Standard_OutOfMemory &) {
            throw std::bad_alloc{};
        } catch (const Standard_Failure &) {
            return nan();
        }
    }
    Vec2 metricScale() const override {
        try {
            BRepAdaptor_Surface s(face_, false);
            gp_Pnt p;
            gp_Vec du, dv;
            s.D1((lo_.u + hi_.u) * .5, (lo_.v + hi_.v) * .5, p, du, dv);
            return {du.Magnitude(), dv.Magnitude()};
        } catch (const Standard_OutOfMemory &) {
            throw std::bad_alloc{};
        } catch (const Standard_Failure &) {
            return {nan(), nan()};
        }
    }
    bool periodicU() const override { return BRepAdaptor_Surface(face_, false).IsUPeriodic(); }
    bool periodicV() const override { return BRepAdaptor_Surface(face_, false).IsVPeriodic(); }
    bool singularAt(const Vec2 &uv) const override {
        try {
            BRepAdaptor_Surface s(face_, false);
            gp_Pnt p;
            gp_Vec du, dv;
            s.D1(uv.u, uv.v, p, du, dv);
            return du.Crossed(dv).Magnitude() <= Precision::Confusion();
        } catch (const Standard_OutOfMemory &) {
            throw std::bad_alloc{};
        } catch (const Standard_Failure &) {
            return true;
        }
    }

  private:
    TopoDS_Face face_;
    std::vector<std::vector<int>> loops_;
    Vec2 lo_, hi_;
    double tolerance_ = 0;
    SurfaceKind kind_ = SurfaceKind::Other;
};
class NativeSolid final : public GSolid {
  public:
    NativeSolid(int tag, TopoDS_Solid solid, std::vector<int> faces, double tolerance)
        : GSolid(tag), solid_(std::move(solid)), faces_(std::move(faces)), tolerance_(tolerance) {
        Bnd_Box box;
        BRepBndLib::AddOptimal(solid_, box, false, false);
        require(!box.IsVoid() && !box.IsOpen());
        box.Get(lo_.x, lo_.y, lo_.z, hi_.x, hi_.y, hi_.z);
        require(finite(lo_) && finite(hi_));
    }
    bool isInside(const Vec3 &p) const override {
        if (!finite(p))
            return false;
        try {
            paacs::geo::detail::OcctSolidClassifier c(solid_);
            c.Perform(point(p), tolerance_);
            return c.State() == TopAbs_IN || c.State() == TopAbs_ON;
        } catch (const Standard_OutOfMemory &) {
            throw std::bad_alloc{};
        } catch (const Standard_Failure &) {
            return false;
        }
    }
    std::vector<int> boundaryFaceTags() const override { return faces_; }
    void bbox(Vec3 &lo, Vec3 &hi) const override {
        lo = lo_;
        hi = hi_;
    }

  private:
    TopoDS_Solid solid_;
    std::vector<int> faces_;
    double tolerance_;
    Vec3 lo_, hi_;
};
OcctSurfaceDescriptor describe(const TopoDS_Face &face, const GFace &geometry) {
    OcctSurfaceDescriptor d;
    d.kind = geometry.kind();
    d.orientation = orientation(face.Orientation());
    geometry.paramBounds(d.parameterLower, d.parameterUpper);
    d.periodicU = geometry.periodicU();
    d.periodicV = geometry.periodicV();
    BRepAdaptor_Surface s(face, false);
    gp_Ax3 position;
    switch (s.GetType()) {
    case GeomAbs_Plane:
        position = s.Plane().Position();
        break;
    case GeomAbs_Cylinder:
        position = s.Cylinder().Position();
        d.radius = s.Cylinder().Radius();
        break;
    case GeomAbs_Cone:
        position = s.Cone().Position();
        d.radius = s.Cone().RefRadius();
        d.semiAngle = s.Cone().SemiAngle();
        break;
    case GeomAbs_Sphere:
        position = s.Sphere().Position();
        d.radius = s.Sphere().Radius();
        break;
    case GeomAbs_Torus:
        position = s.Torus().Position();
        d.radius = s.Torus().MinorRadius();
        d.majorRadius = s.Torus().MajorRadius();
        break;
    default:
        return d;
    }
    d.origin = vec(position.Location());
    d.axis = vec(position.Direction());
    d.xDirection = vec(position.XDirection());
    return d;
}
} // namespace

struct OcctModel::Impl {
    struct Pcurve {
        std::uint64_t id;
        int edgeTag, faceTag;
        double first, last, curveFirst, curveLast, tolerance;
        occ::handle<Geom2d_Curve> curve;
    };
    TopoDS_Shape root;
    ShapeMap vertexMap, edgeMap, faceMap, solidMap;
    std::uint64_t identity = 0;
    std::vector<std::unique_ptr<GVertex>> vertices;
    std::vector<std::unique_ptr<GEdge>> edges;
    std::vector<std::unique_ptr<GFace>> faces;
    std::vector<std::unique_ptr<GSolid>> solids;
    std::vector<double> vertexTolerances, edgeTolerances, faceTolerances, solidTolerances, volumes;
    std::vector<std::vector<OcctWireUse>> wires;
    std::vector<std::vector<OcctFaceUse>> bodyUses;
    std::vector<std::vector<OcctFaceOwner>> owners;
    std::vector<OcctSurfaceDescriptor> descriptors;
    std::vector<Pcurve> pcurves;
    void build(const NativeCadImportOptions &options) {
        TopExp::MapShapes(root, TopAbs_VERTEX, vertexMap);
        TopExp::MapShapes(root, TopAbs_EDGE, edgeMap);
        TopExp::MapShapes(root, TopAbs_FACE, faceMap);
        TopExp::MapShapes(root, TopAbs_SOLID, solidMap);
        identity = nextIdentity.fetch_add(1, std::memory_order_relaxed);
        require(identity != 0);
        for (int i = 1; i <= vertexMap.Extent(); ++i) {
            const auto v = TopoDS::Vertex(vertexMap(i));
            const auto p = vec(BRep_Tool::Pnt(v));
            require(finite(p));
            vertices.push_back(std::make_unique<GVertex>(i, p));
            vertexTolerances.push_back(BRep_Tool::Tolerance(v));
        }
        for (int i = 1; i <= edgeMap.Extent(); ++i) {
            const auto e = TopoDS::Edge(edgeMap(i));
            TopoDS_Vertex a, b;
            TopExp::Vertices(e, a, b, false);
            const int ai = vertexMap.FindIndex(a), bi = vertexMap.FindIndex(b);
            require(ai > 0 && bi > 0);
            edges.push_back(
                std::make_unique<NativeEdge>(i, e, ai, bi, vertices[ai - 1]->position()));
            edgeTolerances.push_back(BRep_Tool::Tolerance(e));
        }
        std::uint64_t useCount = 0;
        owners.resize(faceMap.Extent());
        for (int i = 1; i <= faceMap.Extent(); ++i) {
            const auto face = TopoDS::Face(faceMap(i));
            const auto outer = BRepTools::OuterWire(face);
            require(!outer.IsNull());
            std::vector<TopoDS_Wire> ordered{outer};
            for (TopExp_Explorer w(face, TopAbs_WIRE); w.More(); w.Next())
                if (!w.Current().IsSame(outer))
                    ordered.push_back(TopoDS::Wire(w.Current()));
            std::vector<std::vector<int>> loops;
            std::vector<OcctWireUse> faceWires;
            for (const auto &wire : ordered) {
                OcctWireUse w;
                w.outer = wire.IsSame(outer);
                std::vector<int> tags;
                std::size_t occurrenceCount = 0;
                for (TopExp_Explorer all(wire, TopAbs_EDGE); all.More(); all.Next())
                    ++occurrenceCount;
                for (BRepTools_WireExplorer ex(wire, face); ex.More(); ex.Next()) {
                    require(++useCount <= options.maximumUses, CadImportReason::TopologyBudget);
                    const auto used = ex.Current();
                    const int tag = edgeMap.FindIndex(used);
                    require(tag > 0);
                    TopoDS_Vertex a, b;
                    TopExp::Vertices(used, a, b, true);
                    OcctEdgeUse use;
                    use.useId = useCount;
                    use.edgeTag = tag;
                    use.startVertexTag = vertexMap.FindIndex(a);
                    use.endVertexTag = vertexMap.FindIndex(b);
                    require(use.startVertexTag > 0 && use.endVertexTag > 0);
                    use.orientation = orientation(used.Orientation());
                    use.sense = sense(used.Orientation());
                    use.seam = BRep_Tool::IsClosed(used, face);
                    use.degenerated = BRep_Tool::Degenerated(used);
                    edges[tag - 1]->parameterBounds(use.firstParameter, use.lastParameter);
                    use.tolerance = edgeTolerances[tag - 1];
                    double pf, pl;
                    auto pc = BRep_Tool::CurveOnSurface(used, face, pf, pl);
                    require(!pc.IsNull() && std::isfinite(pf) && std::isfinite(pl) && pl > pf);
                    const double tolerance = use.tolerance + BRep_Tool::Tolerance(face) +
                                             std::max(vertexTolerances[use.startVertexTag - 1],
                                                      vertexTolerances[use.endVertexTag - 1]) +
                                             Precision::Confusion();
                    pcurves.push_back({useCount, tag, i, pf, pl, use.firstParameter,
                                       use.lastParameter, tolerance, std::move(pc)});
                    w.edges.push_back(use);
                    tags.push_back(tag);
                }
                // WireExplorer can omit branching/internal/external uses.
                // Refuse that unsupported boundary contract rather than
                // silently presenting an incomplete topology receipt.
                require(!tags.empty() && tags.size() == occurrenceCount,
                        CadImportReason::UnsupportedCurve);
                loops.push_back(std::move(tags));
                faceWires.push_back(std::move(w));
            }
            faces.push_back(std::make_unique<NativeFace>(i, face, std::move(loops)));
            wires.push_back(std::move(faceWires));
            faceTolerances.push_back(BRep_Tool::Tolerance(face));
            descriptors.push_back(describe(face, *faces.back()));
        }
        for (int i = 1; i <= solidMap.Extent(); ++i) {
            const auto solid = TopoDS::Solid(solidMap(i));
            std::vector<int> tags;
            std::vector<OcctFaceUse> uses;
            double tolerance = Precision::Confusion();
            for (TopExp_Explorer ex(solid, TopAbs_FACE); ex.More(); ex.Next()) {
                require(++useCount <= options.maximumUses, CadImportReason::TopologyBudget);
                const int tag = faceMap.FindIndex(ex.Current());
                require(tag > 0);
                const auto base = faceMap(tag).Orientation(), used = ex.Current().Orientation();
                const int relative =
                    sense(base) != 0 && sense(used) != 0 ? (base == used ? 1 : -1) : 0;
                uses.push_back({tag, orientation(used), relative});
                tags.push_back(tag);
                owners[tag - 1].push_back({i, relative});
                tolerance = std::max(tolerance, faceTolerances[tag - 1]);
            }
            require(!tags.empty());
            // A reversed whole solid can pass BRepCheck while representing
            // the unbounded complement. GSolid describes bounded material;
            // retaining that complement with an absolute volume would give
            // contradictory bounds and inside/outside classification.
            paacs::geo::detail::OcctSolidClassifier classifier(solid);
            classifier.PerformInfinitePoint(tolerance);
            require(classifier.State() == TopAbs_OUT);
            GProp_GProps props;
            BRepGProp::VolumeProperties(solid, props);
            const double volume = std::abs(props.Mass());
            require(std::isfinite(volume) && volume > 0);
            volumes.push_back(volume);
            solidTolerances.push_back(tolerance);
            bodyUses.push_back(std::move(uses));
            solids.push_back(std::make_unique<NativeSolid>(i, solid, std::move(tags), tolerance));
        }
    }
};
OcctModel::OcctModel(std::unique_ptr<Impl> impl) : impl_(std::move(impl)) {}
OcctModel::~OcctModel() = default;
std::vector<const GVertex *> OcctModel::vertices() const {
    return pointers(impl_->vertices);
}
std::vector<const GEdge *> OcctModel::edges() const {
    return pointers(impl_->edges);
}
std::vector<const GFace *> OcctModel::faces() const {
    return pointers(impl_->faces);
}
std::vector<const GSolid *> OcctModel::solids() const {
    return pointers(impl_->solids);
}
const GVertex *OcctModel::vertex(int tag) const {
    return lookup(impl_->vertices, tag);
}
const GEdge *OcctModel::edge(int tag) const {
    return lookup(impl_->edges, tag);
}
const GFace *OcctModel::face(int tag) const {
    return lookup(impl_->faces, tag);
}
const GSolid *OcctModel::solid(int tag) const {
    return lookup(impl_->solids, tag);
}
std::uint64_t OcctModel::identity() const noexcept {
    return impl_->identity;
}
double OcctModel::tolerance(int dimension, int tag) const noexcept {
    const auto *values = dimension == 0   ? &impl_->vertexTolerances
                         : dimension == 1 ? &impl_->edgeTolerances
                         : dimension == 2 ? &impl_->faceTolerances
                         : dimension == 3 ? &impl_->solidTolerances
                                          : nullptr;
    return values && tag > 0 && static_cast<std::size_t>(tag) <= values->size() ? (*values)[tag - 1]
                                                                                : nan();
}
double OcctModel::volume(int tag) const noexcept {
    return tag > 0 && static_cast<std::size_t>(tag) <= impl_->volumes.size()
               ? impl_->volumes[tag - 1]
               : nan();
}
std::span<const OcctFaceUse> OcctModel::faceUses(int tag) const noexcept {
    return row(impl_->bodyUses, tag);
}
std::span<const OcctFaceOwner> OcctModel::faceOwners(int tag) const noexcept {
    return row(impl_->owners, tag);
}
std::span<const OcctWireUse> OcctModel::boundaryUses(int tag) const noexcept {
    return row(impl_->wires, tag);
}
const OcctSurfaceDescriptor *OcctModel::surfaceDescriptor(int tag) const noexcept {
    return tag > 0 && static_cast<std::size_t>(tag) <= impl_->descriptors.size()
               ? &impl_->descriptors[tag - 1]
               : nullptr;
}
bool OcctModel::parameterOnUse(std::uint64_t id, double t, Vec2 &uv) const {
    if (!std::isfinite(t))
        return false;
    const auto it = std::lower_bound(impl_->pcurves.begin(), impl_->pcurves.end(), id,
                                     [](const auto &p, auto value) { return p.id < value; });
    if (it == impl_->pcurves.end() || it->id != id || t < it->curveFirst || t > it->curveLast)
        return false;
    try {
        const double parameter = it->first + (t - it->curveFirst) /
                                                 (it->curveLast - it->curveFirst) *
                                                 (it->last - it->first);
        const auto p = it->curve->Value(parameter);
        const Vec2 candidate{p.X(), p.Y()};
        const auto mapped = face(it->faceTag)->pointAt(candidate),
                   original = edge(it->edgeTag)->pointAt(t);
        if (!finite(candidate) || !finite(mapped) || !finite(original) ||
            (mapped - original).norm() > it->tolerance)
            return false;
        uv = candidate;
        return true;
    } catch (const Standard_OutOfMemory &) {
        throw std::bad_alloc{};
    } catch (const Standard_Failure &) {
        return false;
    }
}
bool detail::OcctModelAccess::fitParameterOnUse(const OcctModel& model, std::uint64_t id,
    double previous, double at, double next, Vec2& uv) {
    const auto& pcurves=model.impl_->pcurves;
    const auto use=std::lower_bound(pcurves.begin(),pcurves.end(),id,
        [](const auto& p, auto key){return p.id<key;});
    if(use==pcurves.end()||use->id!=id||
       !(use->curveFirst<=previous&&previous<at&&at<next&&next<=use->curveLast))return false;
    const auto parameter=[&](double p){
        return use->first+(p-use->curveFirst)/(use->curveLast-use->curveFirst)*(use->last-use->first);
    };
    const auto affine=parameter(at);
        const auto lower=std::midpoint(parameter(previous),affine);
        const auto upper=std::midpoint(affine,parameter(next));
    if(!(lower<affine&&affine<upper))return false;
    try {
        const auto* face=model.face(use->faceTag);
        const auto point=model.edge(use->edgeTag)->pointAt(at);
        if(!finite(point))return false;
        const auto evaluate=[&](double p){
            const auto value=use->curve->Value(p);
            return Vec2{value.X(),value.Y()};
        };
        const auto error=[&](double p){
            const auto candidate=evaluate(p);
            const auto mapped=face->pointAt(candidate);
            if(!finite(candidate)||!finite(mapped))return std::numeric_limits<double>::infinity();
            return (mapped-point).norm();
        };
        auto lo=lower,hi=upper;
        for(unsigned i=0;i<48;++i){
            const auto a=std::lerp(lo,hi,1./3.),b=std::lerp(lo,hi,2./3.);
            if(error(a)<error(b))hi=b;else lo=a;
        }
        const auto fitted=std::midpoint(lo,hi);
        if(!(fitted>lower&&fitted<upper&&error(fitted)<=use->tolerance))return false;
        uv=evaluate(fitted);
        return true;
    } catch(const Standard_OutOfMemory&) {
        throw std::bad_alloc{};
    } catch(const Standard_Failure&) {
        return false;
    }
}
NativeCadImportOutcome importOcctCad(const std::filesystem::path &path,
                                     const NativeCadImportOptions &options) noexcept {
    NativeCadImportOutcome out;
    try {
        require(options.format == CadFileFormat::Step || options.format == CadFileFormat::Brep,
                CadImportReason::InvalidRequest);
        require(options.maximumSolids > 0 && options.maximumFaces > 0 && options.maximumEdges > 0 &&
                    options.maximumVertices > 0 && options.maximumUses > 0,
                CadImportReason::InvalidRequest);
        std::error_code error;
        out.fileBytes = std::filesystem::file_size(path, error);
        require(!error, CadImportReason::FileUnavailable);
        require(out.fileBytes <= options.maximumFileBytes, CadImportReason::FileBudget);
        const auto utf8 = path.u8string();
        const std::string filename(utf8.begin(), utf8.end());
        auto impl = std::make_unique<OcctModel::Impl>();
        if (options.format == CadFileFormat::Step) {
            STEPControl_Reader reader;
            require(reader.ReadFile(filename.c_str()) == IFSelect_RetDone,
                    CadImportReason::ParseFailed);
            reader.SetSystemLengthUnit(1.);
            out.millimetresPerModelUnit = 1.;
            NCollection_Sequence<TCollection_AsciiString> length, angle, solidAngle;
            reader.FileUnits(length, angle, solidAngle);
            for (int i = 1; i <= length.Length(); ++i)
                out.sourceLengthUnits.emplace_back(length(i).ToCString());
            require(reader.TransferRoots() > 0, CadImportReason::ParseFailed);
            impl->root = reader.OneShape();
        } else {
            BRep_Builder builder;
            require(BRepTools::Read(impl->root, filename.c_str(), builder),
                    CadImportReason::ParseFailed);
        }
        require(!impl->root.IsNull(), CadImportReason::ParseFailed);
        ShapeMap solids, faces, edges, vertices;
        TopExp::MapShapes(impl->root, TopAbs_SOLID, solids);
        TopExp::MapShapes(impl->root, TopAbs_FACE, faces);
        TopExp::MapShapes(impl->root, TopAbs_EDGE, edges);
        TopExp::MapShapes(impl->root, TopAbs_VERTEX, vertices);
        out.solids = solids.Extent();
        out.faces = faces.Extent();
        out.edges = edges.Extent();
        out.vertices = vertices.Extent();
        require(out.solids > 0, CadImportReason::UnsupportedSolidCount);
        require(out.solids <= options.maximumSolids && out.faces <= options.maximumFaces &&
                    out.edges <= options.maximumEdges && out.vertices <= options.maximumVertices,
                CadImportReason::TopologyBudget);
        ShapeMap ownedFaces, ownedEdges, ownedVertices;
        for (int i = 1; i <= solids.Extent(); ++i) {
            TopExp::MapShapes(solids(i), TopAbs_FACE, ownedFaces);
            TopExp::MapShapes(solids(i), TopAbs_EDGE, ownedEdges);
            TopExp::MapShapes(solids(i), TopAbs_VERTEX, ownedVertices);
        }
        out.auxiliaryFaces=faces.Extent()-ownedFaces.Extent();
        out.auxiliaryEdges=edges.Extent()-ownedEdges.Extent();
        out.auxiliaryVertices=vertices.Extent()-ownedVertices.Extent();
        const bool auxiliary=out.auxiliaryFaces||out.auxiliaryEdges||out.auxiliaryVertices;
        require(!auxiliary||options.retainAuxiliaryGeometry);
        require(BRepCheck_Analyzer(impl->root).IsValid());
        impl->build(options);
        out.auxiliaryGeometryRetained=auxiliary;
        for (const auto &d : impl->descriptors)
            if (d.kind != SurfaceKind::Plane)
                ++out.curvedFaces;
        for (const auto &e : impl->edges)
            if (e->kind() != CurveKind::Line && e->length() > 0)
                ++out.curvedEdges;
        for (double volume : impl->volumes)
            out.cadVolume += volume;
        require(std::isfinite(out.cadVolume));
        out.model = std::unique_ptr<OcctModel>(new OcctModel(std::move(impl)));
        out.reason = CadImportReason::None;
    } catch (const Stop &stop) {
        out.reason = stop.reason;
    } catch (const std::bad_alloc &) {
        out.reason = CadImportReason::AllocationFailed;
    } catch (const Standard_OutOfMemory &) {
        out.reason = CadImportReason::AllocationFailed;
    } catch (const Standard_Failure &) {
        out.reason = CadImportReason::KernelFailure;
    } catch (...) {
        out.reason = CadImportReason::KernelFailure;
    }
    return out;
}
const TopoDS_Shape &detail::OcctModelAccess::shape(const OcctModel &model) noexcept {
    return model.impl_->root;
}
NativeCadImportOutcome detail::OcctModelAccess::fromShape(const TopoDS_Shape& shape,
                                                        const NativeCadImportOptions& options) noexcept {
    NativeCadImportOutcome out;
    try {
        require(!shape.IsNull()&&options.maximumSolids>0&&options.maximumFaces>0&&options.maximumEdges>0&&
                options.maximumVertices>0&&options.maximumUses>0,CadImportReason::InvalidRequest);
        ShapeMap solids,faces,edges,vertices,ownedFaces,ownedEdges,ownedVertices;
        TopExp::MapShapes(shape,TopAbs_SOLID,solids);TopExp::MapShapes(shape,TopAbs_FACE,faces);
        TopExp::MapShapes(shape,TopAbs_EDGE,edges);TopExp::MapShapes(shape,TopAbs_VERTEX,vertices);
        out.solids=solids.Extent();out.faces=faces.Extent();out.edges=edges.Extent();out.vertices=vertices.Extent();
        require(out.solids>0,CadImportReason::UnsupportedSolidCount);
        require(out.solids<=options.maximumSolids&&out.faces<=options.maximumFaces&&out.edges<=options.maximumEdges&&
                out.vertices<=options.maximumVertices,CadImportReason::TopologyBudget);
        for(int i=1;i<=solids.Extent();++i){TopExp::MapShapes(solids(i),TopAbs_FACE,ownedFaces);
            TopExp::MapShapes(solids(i),TopAbs_EDGE,ownedEdges);TopExp::MapShapes(solids(i),TopAbs_VERTEX,ownedVertices);}
        out.auxiliaryFaces=faces.Extent()-ownedFaces.Extent();out.auxiliaryEdges=edges.Extent()-ownedEdges.Extent();
        out.auxiliaryVertices=vertices.Extent()-ownedVertices.Extent();
        const bool auxiliary=out.auxiliaryFaces||out.auxiliaryEdges||out.auxiliaryVertices;
        require(!auxiliary||options.retainAuxiliaryGeometry);
        require(BRepCheck_Analyzer(shape).IsValid());
        auto impl=std::make_unique<OcctModel::Impl>();impl->root=shape;impl->build(options);
        out.auxiliaryGeometryRetained=auxiliary;
        for(const auto& d:impl->descriptors)if(d.kind!=SurfaceKind::Plane)++out.curvedFaces;
        for(const auto& e:impl->edges)if(e->kind()!=CurveKind::Line&&e->length()>0)++out.curvedEdges;
        for(double volume:impl->volumes)out.cadVolume+=volume;require(std::isfinite(out.cadVolume));
        out.model=std::unique_ptr<OcctModel>(new OcctModel(std::move(impl)));out.reason=CadImportReason::None;
    }catch(const Stop& stop){out.reason=stop.reason;}catch(const Standard_OutOfMemory&){out.reason=CadImportReason::AllocationFailed;}
     catch(const std::bad_alloc&){out.reason=CadImportReason::AllocationFailed;}catch(...){out.reason=CadImportReason::KernelFailure;}
    return out;
}
TopoDS_Solid detail::OcctModelAccess::solid(const OcctModel &model, int tag) {
    if (!model.solid(tag))
        throw std::out_of_range("OCCT solid tag");
    return TopoDS::Solid(model.impl_->solidMap(tag));
}
TopoDS_Face detail::OcctModelAccess::face(const OcctModel &model, int tag) {
    if (!model.face(tag))
        throw std::out_of_range("OCCT face tag");
    return TopoDS::Face(model.impl_->faceMap(tag));
}
TopoDS_Edge detail::OcctModelAccess::edge(const OcctModel &model, int tag) {
    if (!model.edge(tag))
        throw std::out_of_range("OCCT edge tag");
    return TopoDS::Edge(model.impl_->edgeMap(tag));
}
} // namespace paacs::geo
