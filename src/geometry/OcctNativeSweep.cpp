// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctNativeSweep.hpp"
#include "OcctModelInternal.hpp"
#include "paacs/algo/ConvexToTetMesher.hpp"
#include "paacs/algo/LayeredPatchMesher.hpp"
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepCheck_Analyzer.hxx>
#include "OcctSolidClassifier.hpp"
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_Line.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Message_ProgressIndicator.hxx>
#include <NCollection_List.hxx>
#include <Standard_OutOfMemory.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <map>
#include <ostream>
#include <set>
#include <streambuf>

namespace paacs::geo::detail {
namespace {
using Reason = OcctNativeSweepReason;
using Clock = std::chrono::steady_clock;
using Pair = std::pair<Index, Index>;
struct Stop {
    Reason reason;
    const char *diagnostic;
};
void need(bool ok, Reason reason, const char *diagnostic) {
    if (!ok)
        throw Stop{reason, diagnostic};
}
struct Timer {
    double &seconds;
    Clock::time_point start = Clock::now();
    ~Timer() { seconds += std::chrono::duration<double>(Clock::now() - start).count(); }
};
bool finite(Vec3 p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
Pair pairKey(Index a, Index b) {
    return a < b ? Pair{a, b} : Pair{b, a};
}
std::uint64_t add(std::uint64_t a, std::uint64_t b) {
    need(b <= UINT64_MAX - a, Reason::BudgetExceeded, "payload-overflow");
    return a + b;
}
std::uint64_t mul(std::uint64_t a, std::uint64_t b) {
    need(!b || a <= UINT64_MAX / b, Reason::BudgetExceeded, "payload-overflow");
    return a * b;
}
std::uint64_t payload(std::uint64_t n, std::uint64_t c, std::uint64_t cn) {
    return add(add(mul(n, sizeof(Vec3) + sizeof(mesh::Classification)), mul(c, sizeof(mesh::Cell))),
               mul(cn, sizeof(Index)));
}
std::uint64_t payload(const mesh::Mesh &m) {
    return payload(m.nodeCount(), m.cellCount(), m.connectivityCount());
}
struct Context {
    const OcctNativeSweepOptions &options;
    OcctNativeSweepOutcome &out;
    std::uint64_t remaining() const { return options.maximumWork - out.workUsed; }
    void check() {
        need(!options.cancellation || !options.cancellation->load(std::memory_order_relaxed),
             Reason::Cancelled, "cancelled");
    }
    void tick(std::uint64_t n, std::uint64_t &phase) {
        check();
        need(n <= remaining(), Reason::BudgetExceeded, "cumulative-work-budget");
        out.workUsed += n;
        phase += n;
    }
    void collect(std::uint64_t n, std::uint64_t &phase) {
        need(n <= remaining(), Reason::BudgetExceeded, "cumulative-work-budget");
        out.workUsed += n;
        phase += n;
        check();
    }
    void output(std::uint64_t n, std::uint64_t c, std::uint64_t cn) {
        need(n <= options.maximumNodes && c <= options.maximumCells &&
                 payload(n, c, cn) <= options.maximumOutputBytes,
             Reason::BudgetExceeded, "output-budget");
    }
    void working(std::uint64_t bytes) {
        out.maximumWorkingPayloadEstimate = std::max(out.maximumWorkingPayloadEstimate, bytes);
        need(bytes <= options.maximumWorkingBytes, Reason::BudgetExceeded,
             "working-payload-admission");
    }
};
class Buffer final : public std::streambuf {
  public:
    explicit Buffer(std::uint64_t limit) : limit_(limit) {}
    std::string value;

  protected:
    std::streamsize xsputn(const char *p, std::streamsize n) override {
        need(n >= 0 && static_cast<std::uint64_t>(n) <= limit_ - value.size(),
             Reason::BudgetExceeded, "source-snapshot-budget");
        value.append(p, static_cast<std::size_t>(n));
        return n;
    }
    int_type overflow(int_type c) override {
        if (traits_type::eq_int_type(c, traits_type::eof()))
            return traits_type::not_eof(c);
        char ch = traits_type::to_char_type(c);
        xsputn(&ch, 1);
        return c;
    }

  private:
    std::uint64_t limit_;
};
std::string snapshot(const OcctModel &m, std::uint64_t maximum) {
    Buffer buffer(maximum);
    std::ostream stream(&buffer);
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    BRepTools::Write(OcctModelAccess::shape(m), stream);
    return std::move(buffer.value);
}
class Progress final : public Message_ProgressIndicator {
  public:
    explicit Progress(const std::atomic_bool *flag) : flag_(flag) {}

  protected:
    bool UserBreak() override { return flag_ && flag_->load(std::memory_order_relaxed); }
    void Show(const Message_ProgressScope &, bool) override {}

  private:
    const std::atomic_bool *flag_;
};
TopoDS_Shape copy(const TopoDS_Shape &source, Context &context) {
    context.tick(1, context.out.recognitionWork);
    BRepBuilderAPI_Copy operation(source, true, false);
    context.check();
    need(operation.IsDone(), Reason::KernelFailure, "native-copy-failed");
    return operation.Shape();
}
double mass(const TopoDS_Shape &shape, Context &context) {
    double values[2];
    for (int i = 0; i < 2; ++i) {
        context.tick(1, context.out.recognitionWork);
        GProp_GProps properties;
        const auto error = BRepGProp::VolumePropertiesGK(shape, properties, i ? 1e-10 : 1e-8, false,
                                                         true, false, false, false);
        context.check();
        values[i] = properties.Mass();
        need(std::isfinite(error) && error >= 0 && std::isfinite(values[i]) && values[i] > 0,
             Reason::RecognitionFailed, "invalid-native-mass");
    }
    const double stability = std::abs(values[1] - values[0]) / values[1];
    context.out.integrationStability = std::max(context.out.integrationStability, stability);
    need(stability <= context.options.maximumRelativeRecognitionError, Reason::RecognitionFailed,
         "unstable-native-mass");
    return values[1];
}
struct Cap {
    int tag;
    double area;
    gp_Pln plane;
};
struct Recognition {
    int base = -1, top = -1;
    Vec3 direction{};
    double height = 0;
};
Recognition recognize(const OcctModel &model, int body, const std::set<int> &ownedFaces,
                      const std::set<int> &ownedVertices, Context &context) {
    const auto native = OcctModelAccess::solid(model, body);
    std::vector<Cap> caps;
    for (int tag : ownedFaces) {
        context.tick(1, context.out.recognitionWork);
        const auto face = OcctModelAccess::face(model, tag);
        BRepAdaptor_Surface surface(face);
        if (surface.GetType() == GeomAbs_Plane) {
            GProp_GProps p;
            BRepGProp::SurfaceProperties(face, p);
            context.check();
            need(std::isfinite(p.Mass()) && p.Mass() > 0, Reason::RecognitionFailed,
                 "invalid-cap-area");
            caps.push_back({tag, p.Mass(), surface.Plane()});
        }
    }
    std::stable_sort(caps.begin(), caps.end(),
                     [](const auto &a, const auto &b) { return a.area > b.area; });
    need(caps.size() >= 2, Reason::UnsupportedGeometry, "native-solid-has-no-planar-cap-pair");
    context.out.cadVolume = mass(native, context);
    for (std::size_t i = 0; i < caps.size(); ++i)
        for (std::size_t j = i + 1; j < caps.size(); ++j) {
            context.tick(1, context.out.recognitionWork);
            auto base = caps[i], top = caps[j];
            const auto axis = base.plane.Axis().Direction();
            if (std::abs(axis.Dot(top.plane.Axis().Direction())) < 1 - 1e-10 ||
                std::abs(base.area - top.area) / std::max(base.area, top.area) >
                    context.options.maximumRelativeRecognitionError)
                continue;
            double height = gp_Vec(base.plane.Location(), top.plane.Location()).Dot(gp_Vec(axis));
            if (!std::isfinite(height) || std::abs(height) <= std::max(model.tolerance(2, base.tag),
                                                                       model.tolerance(2, top.tag)))
                continue;
            if (height < 0) {
                std::swap(base, top);
                height = -height;
            }
            bool supported = true;
            for (int tag : ownedFaces) {
                if (tag == base.tag || tag == top.tag)
                    continue;
                context.tick(1, context.out.recognitionWork);
                BRepAdaptor_Surface surface(OcctModelAccess::face(model, tag));
                if (surface.GetType() == GeomAbs_Plane)
                    supported = std::abs(surface.Plane().Axis().Direction().Dot(axis)) <= 1e-10;
                else if (surface.GetType() == GeomAbs_Cylinder)
                    supported =
                        std::abs(surface.Cylinder().Axis().Direction().Dot(axis)) >= 1 - 1e-10;
                else if (surface.GetType() == GeomAbs_SurfaceOfExtrusion)
                    supported = std::abs(surface.Direction().Dot(axis)) >= 1 - 1e-10;
                else
                    supported = false;
                if (!supported)
                    break;
            }
            if (!supported)
                continue;
            for (int tag : ownedVertices) {
                context.tick(1, context.out.recognitionWork);
                const auto p = model.vertex(tag)->position();
                const double offset =
                    gp_Vec(base.plane.Location(), gp_Pnt(p.x, p.y, p.z)).Dot(gp_Vec(axis));
                const auto tol = model.tolerance(0, tag);
                if (offset < -tol || offset > height + tol) {
                    supported = false;
                    break;
                }
            }
            if (!supported)
                continue;
            auto cap = TopoDS::Face(copy(OcctModelAccess::face(model, base.tag), context));
            cap.Orientation(TopAbs_FORWARD);
            if (BRepAdaptor_Surface(cap).Plane().Axis().Direction().Dot(axis) < 0)
                cap.Reverse();
            context.tick(1, context.out.recognitionWork);
            BRepPrimAPI_MakePrism prism(cap, gp_Vec(axis) * height, false, true);
            context.check();
            need(prism.IsDone(), Reason::KernelFailure, "native-prism-failed");
            const auto expected = prism.Shape();
            context.tick(1, context.out.recognitionWork);
            need(BRepCheck_Analyzer(expected).IsValid(), Reason::RecognitionFailed,
                 "invalid-native-prism");
            paacs::geo::detail::OcctSolidClassifier classifier(expected);
            classifier.PerformInfinitePoint(model.tolerance(3, body));
            context.check();
            need(classifier.State() == TopAbs_OUT, Reason::RecognitionFailed,
                 "unbounded-native-prism");
            context.out.expectedVolume = mass(expected, context);
            context.out.expectedRelativeError =
                std::abs(context.out.expectedVolume - context.out.cadVolume) /
                context.out.cadVolume;
            if (context.out.expectedRelativeError > context.options.maximumRelativeRecognitionError)
                continue;
            NCollection_List<TopoDS_Shape> arguments, tools;
            arguments.Append(copy(native, context));
            tools.Append(copy(expected, context));
            BRepAlgoAPI_Common operation;
            operation.SetArguments(arguments);
            operation.SetTools(tools);
            operation.SetNonDestructive(true);
            operation.SetRunParallel(false);
            occ::handle<Progress> progress = new Progress(context.options.cancellation);
            context.tick(1, context.out.recognitionWork);
            operation.Build(progress->Start());
            context.check();
            need(operation.IsDone() && !operation.HasErrors(), Reason::KernelFailure,
                 "native-common-failed");
            context.out.commonVolume = mass(operation.Shape(), context);
            context.out.commonRelativeError =
                std::abs(context.out.commonVolume - context.out.cadVolume) / context.out.cadVolume;
            if (context.out.commonRelativeError > context.options.maximumRelativeRecognitionError)
                continue;
            return {base.tag, top.tag, {axis.X(), axis.Y(), axis.Z()}, height};
        }
    throw Stop{Reason::UnsupportedGeometry, "native-solid-is-not-proved-uniform-extrusion"};
}
// Exact native endpoint incidence establishes identity. This separate native
// control-polygon bound verifies straight support without fitting a curve.
void verifyVerticalEdge(const OcctModel &model, int tag, Vec3 origin, Vec3 direction,
                        Context &context) {
    const auto edge = OcctModelAccess::edge(model, tag);
    TopLoc_Location location;
    double lo, hi;
    auto curve = BRep_Tool::Curve(edge, location, lo, hi);
    need(!curve.IsNull(), Reason::UnsupportedGeometry, "missing-native-vertical-curve");
    while (auto trimmed = occ::down_cast<Geom_TrimmedCurve>(curve)) {
        context.tick(1, context.out.recognitionWork);
        curve = trimmed->BasisCurve();
    }
    const double tol = model.tolerance(1, tag);
    auto check = [&](gp_Pnt p) {
        context.tick(1, context.out.recognitionWork);
        p.Transform(location.Transformation());
        const Vec3 d{p.X() - origin.x, p.Y() - origin.y, p.Z() - origin.z};
        need(finite(d) && (d - direction * d.dot(direction)).norm() <= tol,
             Reason::UnsupportedGeometry, "native-vertical-control-polygon-is-not-straight");
    };
    if (auto line = occ::down_cast<Geom_Line>(curve)) {
        check(line->Value(lo));
        check(line->Value(hi));
    } else if (auto spline = occ::down_cast<Geom_BSplineCurve>(curve)) {
        for (int i = 1; i <= spline->NbPoles(); ++i) {
            need(std::isfinite(spline->Weight(i)) && spline->Weight(i) > 0,
                 Reason::UnsupportedGeometry, "invalid-native-vertical-weight");
            check(spline->Pole(i));
        }
    } else if (auto bezier = occ::down_cast<Geom_BezierCurve>(curve)) {
        for (int i = 1; i <= bezier->NbPoles(); ++i) {
            need(std::isfinite(bezier->Weight(i)) && bezier->Weight(i) > 0,
                 Reason::UnsupportedGeometry, "invalid-native-vertical-weight");
            check(bezier->Pole(i));
        }
    } else
        throw Stop{Reason::UnsupportedGeometry, "unsupported-native-vertical-curve"};
}
void childBoundary(OcctBoundaryReason reason) {
    if (reason == OcctBoundaryReason::BudgetExceeded)
        throw Stop{Reason::BudgetExceeded, "surface-budget"};
    if (reason == OcctBoundaryReason::Cancelled)
        throw Stop{Reason::Cancelled, "surface-cancelled"};
    if (reason == OcctBoundaryReason::AllocationFailed)
        throw Stop{Reason::AllocationFailed, "surface-allocation-failed"};
    throw Stop{Reason::SurfaceFailed, occtBoundaryReasonName(reason)};
}
void childSweep(algo::LayeredPatchReason reason) {
    using R = algo::LayeredPatchReason;
    if (reason == R::SourceBudget || reason == R::OutputBudget || reason == R::WorkingBudget)
        throw Stop{Reason::BudgetExceeded, "core-sweep-budget"};
    if (reason == R::Cancelled)
        throw Stop{Reason::Cancelled, "core-sweep-cancelled"};
    if (reason == R::AllocationFailed)
        throw Stop{Reason::AllocationFailed, "core-sweep-allocation-failed"};
    throw Stop{Reason::SweepFailed, "retained-core-layered-sweep-refused"};
}
void childConversion(algo::ConvexToTetReason reason) {
    using R = algo::ConvexToTetReason;
    if (reason == R::OutputBudget || reason == R::WorkingBudget)
        throw Stop{Reason::BudgetExceeded, "core-conversion-budget"};
    if (reason == R::Cancelled)
        throw Stop{Reason::Cancelled, "core-conversion-cancelled"};
    if (reason == R::AllocationFailed)
        throw Stop{Reason::AllocationFailed, "core-conversion-allocation-failed"};
    if (reason == R::QualityRequirement)
        throw Stop{Reason::QualityRequirement, "quality-requirement-not-met"};
    throw Stop{Reason::SweepFailed, algo::convexToTetReasonName(reason)};
}
} // namespace

OcctNativeSweepOutcome buildOcctNativeSweepVolume(const OcctModel &model, int body,
                                                  const OcctNativeSweepOptions &options,
                                                  mesh::Mesh &destination) noexcept {
    const auto start = Clock::now();
    OcctNativeSweepOutcome out;
    Context context{options, out};
    try {
        need(model.solid(body) && options.layers >= 1 && options.layers <= 4096 &&
                 std::isfinite(options.maximumLayerThickness) &&
                 options.maximumLayerThickness >= 0 &&
                 std::isfinite(options.maximumRelativeRecognitionError) &&
                 options.maximumRelativeRecognitionError > 0 &&
                 options.maximumRelativeRecognitionError <= 1e-6 &&
                 std::isfinite(options.maximumRelativeVolumeError) &&
                 options.maximumRelativeVolumeError > 0 && options.maximumRelativeVolumeError < 1 &&
                 std::isfinite(options.minimumQuality) && options.minimumQuality >= 0 &&
                 options.minimumQuality <= 1 && std::isfinite(options.warningQuality) &&
                 options.warningQuality >= 0 && options.warningQuality <= 1 &&
                 std::isfinite(options.boundary.linearDeflection) &&
                 options.boundary.linearDeflection > 0 &&
                 std::isfinite(options.boundary.maximumChordError) &&
                 options.boundary.maximumChordError > 0,
             Reason::InvalidRequest, "invalid-request");
        need(options.maximumWork && options.maximumNodes && options.maximumCells &&
                 options.maximumOutputBytes && options.maximumWorkingBytes &&
                 options.maximumSnapshotBytes,
             Reason::BudgetExceeded, "zero-admission-budget");
        context.check();
        need(options.maximumSnapshotBytes <= options.maximumWorkingBytes / 8,
             Reason::BudgetExceeded, "source-snapshot-working-admission");
        context.tick(1, out.recognitionWork);
        const auto before = snapshot(model, options.maximumSnapshotBytes);
        const auto fixedWorking = mul(before.size(), 4);
        context.working(fixedWorking);
        std::set<int> ownedFaces, ownedEdges, ownedVertices;
        std::map<int, std::set<int>> edgeFaces;
        for (const auto &use : model.faceUses(body)) {
            context.tick(1, out.recognitionWork);
            ownedFaces.insert(use.faceTag);
            for (const auto &wire : model.boundaryUses(use.faceTag))
                for (const auto &e : wire.edges) {
                    context.tick(1, out.recognitionWork);
                    ownedEdges.insert(e.edgeTag);
                    edgeFaces[e.edgeTag].insert(use.faceTag);
                    if (e.startVertexTag > 0)
                        ownedVertices.insert(e.startVertexTag);
                    if (e.endVertexTag > 0)
                        ownedVertices.insert(e.endVertexTag);
                }
        }
        const auto topologyWorking =
            add(fixedWorking,
                mul(add(add(ownedFaces.size(), ownedEdges.size()), ownedVertices.size()), 512));
        context.working(topologyWorking);
        Recognition recognition;
        {
            Timer timer{out.recognitionSeconds};
            recognition = recognize(model, body, ownedFaces, ownedVertices, context);
        }
        out.baseCapTag = recognition.base;
        out.targetCapTag = recognition.top;
        out.direction = recognition.direction;
        out.height = recognition.height;
        unsigned layers = options.layers;
        if (options.maximumLayerThickness > 0) {
            const double required = std::ceil(recognition.height / options.maximumLayerThickness);
            need(std::isfinite(required) && required <= 4096, Reason::BudgetExceeded,
                 "normal-layer-admission");
            layers = std::max(layers, static_cast<unsigned>(std::max(1., required)));
        }
        out.layers = layers;
        std::set<int> baseEdges, topEdges, baseVertices, topVertices;
        auto capTopology = [&](int face, std::set<int> &edges, std::set<int> &vertices) {
            for (const auto &wire : model.boundaryUses(face))
                for (const auto &e : wire.edges) {
                    context.tick(1, out.recognitionWork);
                    need(e.startVertexTag > 0 && e.endVertexTag > 0, Reason::ProvenanceLost,
                         "native-cap-edge-has-no-endpoint-identity");
                    edges.insert(e.edgeTag);
                    vertices.insert(e.startVertexTag);
                    vertices.insert(e.endVertexTag);
                }
        };
        capTopology(recognition.base, baseEdges, baseVertices);
        capTopology(recognition.top, topEdges, topVertices);
        std::map<int, int> sideForEdge, topForEdge, verticalForVertex, topForVertex;
        {
            Timer timer{out.recognitionSeconds};
            for (int tag : baseEdges) {
                std::set<int> sides;
                for (int owner : edgeFaces.at(tag)) {
                    context.tick(1, out.recognitionWork);
                    if (owner != recognition.base)
                        sides.insert(owner);
                }
                need(sides.size() == 1, Reason::ProvenanceLost,
                     "ambiguous-native-base-edge-side-face");
                const int side = *sides.begin();
                sideForEdge[tag] = side;
                std::set<int> targets;
                for (const auto &wire : model.boundaryUses(side))
                    for (const auto &e : wire.edges) {
                        context.tick(1, out.recognitionWork);
                        if (topEdges.contains(e.edgeTag))
                            targets.insert(e.edgeTag);
                    }
                need(targets.size() == 1, Reason::ProvenanceLost,
                     "ambiguous-native-translated-cap-edge");
                topForEdge[tag] = *targets.begin();
            }
            for (int vertex : baseVertices) {
                std::vector<std::pair<int, int>> candidates;
                for (int tag : ownedEdges) {
                    context.tick(1, out.recognitionWork);
                    if (baseEdges.contains(tag) || topEdges.contains(tag))
                        continue;
                    const auto *e = model.edge(tag);
                    const int a = e->startVertexTag(), b = e->endVertexTag();
                    if (a == vertex && topVertices.contains(b))
                        candidates.push_back({tag, b});
                    if (b == vertex && topVertices.contains(a))
                        candidates.push_back({tag, a});
                }
                need(candidates.size() == 1, Reason::ProvenanceLost,
                     "ambiguous-native-vertical-edge-incidence");
                const auto [edge, target] = candidates.front();
                verticalForVertex[vertex] = edge;
                topForVertex[vertex] = target;
                const auto p = model.vertex(vertex)->position(),
                           q = model.vertex(target)->position();
                need((q - p - recognition.direction * recognition.height).norm() <=
                         std::max(model.tolerance(0, vertex), model.tolerance(0, target)),
                     Reason::GeometryError, "native-vertical-endpoint-mismatch");
                verifyVerticalEdge(model, edge, p, recognition.direction, context);
            }
        }
        mesh::Mesh boundary;
        need(context.remaining() > 0, Reason::BudgetExceeded, "surface-work-budget");
        auto bo = options.boundary;
        bo.cancellation = options.cancellation;
        bo.maximumWork = std::min(bo.maximumWork, context.remaining());
        bo.maximumNodes = std::min(bo.maximumNodes, options.maximumNodes);
        bo.maximumTriangles = std::min(bo.maximumTriangles, options.maximumCells);
        bo.maximumOutputBytes = std::min({bo.maximumOutputBytes, options.maximumOutputBytes,
                                          (options.maximumWorkingBytes - topologyWorking) / 3});
        OcctBoundaryOutcome br;
        {
            Timer timer{out.surfaceSeconds};
            br = buildOcctBoundary(model, body, bo, boundary);
        }
        out.boundaryReason = br.reason;
        context.collect(br.work, out.surfaceWork);
        if (!br.published)
            childBoundary(br.reason);
        need(br.sourceModelUnmodified, Reason::SourceMutation, "surface-source-audit-missing");
        context.working(add(topologyWorking, mul(payload(boundary), 3)));
        mesh::Mesh patch;
        std::vector<Index> remap(static_cast<std::size_t>(boundary.nodeCount()), kInvalidIndex);
        {
            Timer timer{out.sweepSeconds};
            for (Index i = 0; i < boundary.cellCount(); ++i) {
                context.tick(1, out.sweepWork);
                const auto &cell = boundary.cell(i);
                if (cell.onEntity != mesh::Classification{2, recognition.base})
                    continue;
                need(cell.type == mesh::CellType::Tri3, Reason::InvalidMesh,
                     "nontriangle-native-cap");
                const auto *ids = boundary.cellNodes(i);
                std::array<Index, 3> p{};
                for (int j = 0; j < 3; ++j) {
                    auto &id = remap[static_cast<std::size_t>(ids[j])];
                    if (id == kInvalidIndex) {
                        context.output(patch.nodeCount() + 1, patch.cellCount(),
                                       patch.connectivityCount());
                        id = patch.addNode(boundary.node(ids[j]), boundary.nodeClass(ids[j]));
                    }
                    p[j] = id;
                }
                if ((patch.node(p[1]) - patch.node(p[0]))
                        .cross(patch.node(p[2]) - patch.node(p[0]))
                        .dot(recognition.direction) < 0)
                    std::swap(p[1], p[2]);
                context.output(patch.nodeCount(), patch.cellCount() + 1,
                               patch.connectivityCount() + 3);
                patch.addCell(mesh::CellType::Tri3, p, {2, recognition.base});
            }
            need(patch.nodeCount() > 0 && patch.cellCount() > 0, Reason::InvalidMesh,
                 "empty-native-cap-triangulation");
            out.capNodes = patch.nodeCount();
            out.capTriangles = patch.cellCount();
        }
        std::map<Pair, unsigned> incidence;
        for (Index i = 0; i < patch.cellCount(); ++i) {
            context.tick(1, out.sweepWork);
            const auto *ids = patch.cellNodes(i);
            for (int j = 0; j < 3; ++j)
                ++incidence[pairKey(ids[j], ids[(j + 1) % 3])];
        }
        std::map<Pair, int> sideByPair;
        for (const auto &[pair, count] : incidence) {
            context.tick(1, out.sweepWork);
            need(count == 1 || count == 2, Reason::InvalidMesh, "nonmanifold-native-cap");
            if (count == 2)
                continue;
            const auto a = patch.nodeClass(pair.first), b = patch.nodeClass(pair.second);
            std::set<int> candidates;
            for (int edge : baseEdges) {
                context.tick(1, out.sweepWork);
                const auto *e = model.edge(edge);
                auto belongs = [&](mesh::Classification c) {
                    return c == mesh::Classification{1, edge} ||
                           (c.dim == 0 &&
                            (c.tag == e->startVertexTag() || c.tag == e->endVertexTag()));
                };
                if (belongs(a) && belongs(b))
                    candidates.insert(edge);
            }
            need(candidates.size() == 1, Reason::ProvenanceLost,
                 "ambiguous-native-cap-segment-identity");
            sideByPair[pair] = sideForEdge.at(*candidates.begin());
        }
        const auto patchWork =
            add(topologyWorking, add(add(payload(boundary), payload(patch)),
                                     add(mul(remap.size(), sizeof(Index)),
                                         mul(add(incidence.size(), sideByPair.size()), 128))));
        context.working(patchWork);
        const auto plannedNodes = mul(patch.nodeCount(), layers + 1);
        const auto plannedTets = mul(mul(patch.cellCount(), layers), 3);
        const auto plannedTriangles =
            add(mul(patch.cellCount(), 2), mul(mul(sideByPair.size(), layers), 2));
        context.output(plannedNodes, add(plannedTets, plannedTriangles),
                       add(mul(plannedTets, 4), mul(plannedTriangles, 3)));
        const auto plannedPayload = payload(plannedNodes, add(plannedTets, plannedTriangles),
                                            add(mul(plannedTets, 4), mul(plannedTriangles, 3)));
        context.working(add(patchWork, mul(plannedPayload, 5)));
        algo::LayeredPatchRequest request;
        request.direction = recognition.direction;
        request.volumeClassification = {3, body};
        request.targetClassification = {2, recognition.top};
        request.sideClassification = {2, 0};
        request.schedule.wallLayers = layers;
        request.schedule.firstLayerThickness = recognition.height / layers;
        request.schedule.growthRatio = 1;
        request.schedule.totalDepth = recognition.height;
        request.schedule.maximumCoreThickness = recognition.height / layers;
        request.schedule.maximumLayers = layers + 1;
        request.maximumSourceCells = options.maximumCells;
        request.maximumSourceNodes = options.maximumNodes;
        request.maximumSourcePairChecks = context.remaining();
        request.maximumNodes = options.maximumNodes;
        request.maximumCells = options.maximumCells;
        request.maximumOutputBytes = options.maximumOutputBytes;
        request.maximumWorkingBytes = options.maximumWorkingBytes - patchWork;
        request.auditLimits = options.auditLimits;
        request.auditLimits.cancellation = options.cancellation;
        request.cancellation = options.cancellation;
        const auto plannedSweepCells =
            add(add(mul(patch.cellCount(), layers), mul(patch.cellCount(), 2)),
                mul(sideByPair.size(), layers));
        // Reserve the retained engine's known emitted-node/cell visits before
        // construction; its separately reported exact source-pair checks then
        // spend only the still-available shared work admission.
        context.tick(add(plannedNodes, plannedSweepCells), out.sweepWork);
        request.maximumSourcePairChecks = context.remaining();
        mesh::Mesh swept;
        algo::LayeredPatchResult sr;
        {
            Timer timer{out.sweepSeconds};
            sr = algo::meshLayeredPatch(patch, request, swept);
        }
        context.collect(sr.sourcePairChecks, out.sweepWork);
        if (!sr.published)
            childSweep(sr.reason);
        need(sr.schedule.offsets.size() == layers + 1 &&
                 swept.nodeCount() == static_cast<Index>(plannedNodes),
             Reason::InvalidMesh, "unexpected-core-sweep-node-layout");
        mesh::Mesh classified;
        {
            Timer timer{out.sweepSeconds};
            for (Index i = 0; i < swept.nodeCount(); ++i) {
                context.tick(1, out.sweepWork);
                const auto index = i % patch.nodeCount(), row = i / patch.nodeCount();
                const auto base = patch.nodeClass(index);
                mesh::Classification cls{3, body};
                Vec3 p = swept.node(i);
                need((p -
                      (patch.node(index) +
                       recognition.direction * sr.schedule.offsets[static_cast<std::size_t>(row)]))
                             .norm() <= 1e-10,
                     Reason::InvalidMesh, "unexpected-core-sweep-coordinate-layout");
                if (row == 0)
                    cls = base;
                else if (row == layers) {
                    if (base.dim == 0) {
                        cls = {0, topForVertex.at(base.tag)};
                        p = model.vertex(cls.tag)->position();
                    } else if (base.dim == 1)
                        cls = {1, topForEdge.at(base.tag)};
                    else {
                        need(base == mesh::Classification{2, recognition.base},
                             Reason::ProvenanceLost, "foreign-base-cap-node");
                        cls = {2, recognition.top};
                    }
                } else if (base.dim == 0)
                    cls = {1, verticalForVertex.at(base.tag)};
                else if (base.dim == 1)
                    cls = {2, sideForEdge.at(base.tag)};
                classified.addNode(p, cls);
            }
            for (Index i = 0; i < swept.cellCount(); ++i) {
                context.tick(1, out.sweepWork);
                const auto &cell = swept.cell(i);
                const auto *ids = swept.cellNodes(i);
                auto cls = cell.onEntity;
                if (cls == request.sideClassification) {
                    std::set<Index> baseIds;
                    for (int j = 0; j < mesh::nodesPerCell(cell.type); ++j)
                        baseIds.insert(ids[j] % patch.nodeCount());
                    need(baseIds.size() == 2, Reason::InvalidMesh,
                         "unexpected-core-swept-side-topology");
                    cls = {2, sideByPair.at(pairKey(*baseIds.begin(), *baseIds.rbegin()))};
                }
                classified.addCell(
                    cell.type, {ids, static_cast<std::size_t>(mesh::nodesPerCell(cell.type))}, cls);
            }
        }
        algo::ConvexToTetOptions conversion;
        conversion.maximumNodes = options.maximumNodes;
        conversion.maximumCells = options.maximumCells;
        conversion.maximumOutputBytes = options.maximumOutputBytes;
        const auto conversionStaging = add(patchWork, add(payload(swept), payload(classified)));
        context.working(conversionStaging);
        conversion.maximumWorkingBytes = options.maximumWorkingBytes - conversionStaging;
        conversion.minimumQuality = options.minimumQuality;
        conversion.warningQuality = options.warningQuality;
        conversion.cancellation = options.cancellation;
        conversion.auditLimits = options.auditLimits;
        conversion.auditLimits.cancellation = options.cancellation;
        context.tick(add(classified.nodeCount(), classified.cellCount()), out.sweepWork);
        mesh::Mesh result;
        algo::ConvexToTetResult tr;
        {
            Timer timer{out.sweepSeconds};
            tr = algo::convertConvexCellsToTet4(classified, conversion, result);
        }
        context.check();
        if (!tr.published)
            childConversion(tr.reason);
        context.output(result.nodeCount(), result.cellCount(), result.connectivityCount());
        mesh::Mesh exterior;
        {
            Timer timer{out.validationSeconds};
            const auto validationStaging = add(conversionStaging, payload(result));
            context.working(add(validationStaging,
                                add(payload(result), mul(result.nodeCount(), sizeof(Index)))));
            std::vector<Index> exteriorIds(static_cast<std::size_t>(result.nodeCount()),
                                           kInvalidIndex);
            for (Index i = 0; i < result.cellCount(); ++i) {
                context.tick(1, out.validationWork);
                if (result.cell(i).type != mesh::CellType::Tri3)
                    continue;
                std::array<Index, 3> triangle{};
                const auto *ids = result.cellNodes(i);
                for (int j = 0; j < 3; ++j) {
                    auto &index = exteriorIds[static_cast<std::size_t>(ids[j])];
                    if (index == kInvalidIndex)
                        index = exterior.addNode(result.node(ids[j]), result.nodeClass(ids[j]));
                    triangle[j] = index;
                }
                exterior.addCell(mesh::CellType::Tri3, triangle, result.cell(i).onEntity);
            }
            const auto validationWorking = add(
                validationStaging, add(payload(exterior), mul(exteriorIds.size(), sizeof(Index))));
            context.working(validationWorking);
            need(context.remaining() > 0, Reason::BudgetExceeded, "exterior-work-budget");
            OcctStarOptions co;
            co.maximumWork = context.remaining();
            co.maximumBoundaryNodes = options.maximumNodes;
            co.maximumBoundaryTriangles = options.maximumCells;
            co.maximumNodes = options.maximumNodes;
            co.maximumCells = options.maximumCells;
            co.maximumOutputBytes = options.maximumOutputBytes;
            co.maximumWorkingBytes = options.maximumWorkingBytes - validationWorking;
            co.cancellation = options.cancellation;
            co.auditLimits = options.auditLimits;
            co.auditLimits.maximumTemporaryBytes =
                std::min(co.auditLimits.maximumTemporaryBytes, co.maximumWorkingBytes);
            co.auditLimits.cancellation = options.cancellation;
            out.exteriorCertificate = certifyOcctBoundarySurface(exterior, co);
            context.collect(out.exteriorCertificate.workUsed, out.validationWork);
            if (!out.exteriorCertificate.certified) {
                const auto reason = out.exteriorCertificate.reason;
                if (reason == OcctStarReason::WorkBudget || reason == OcctStarReason::OutputBudget)
                    throw Stop{Reason::BudgetExceeded, out.exteriorCertificate.diagnostic};
                if (reason == OcctStarReason::Cancelled)
                    throw Stop{Reason::Cancelled, out.exteriorCertificate.diagnostic};
                if (reason == OcctStarReason::AllocationFailed)
                    throw Stop{Reason::AllocationFailed, out.exteriorCertificate.diagnostic};
                throw Stop{Reason::InvalidMesh, out.exteriorCertificate.diagnostic};
            }
            paacs::geo::detail::OcctSolidClassifier classifier(OcctModelAccess::solid(model, body));
            auto inside = [&](Vec3 p) {
                context.tick(1, out.validationWork);
                ++out.nativeInteriorQueries;
                classifier.Perform(gp_Pnt(p.x, p.y, p.z), model.tolerance(3, body));
                context.check();
                return classifier.State() == TopAbs_IN || classifier.State() == TopAbs_ON;
            };
            for (Index i = 0; i < result.nodeCount(); ++i) {
                context.tick(1, out.validationWork);
                const auto p = result.node(i);
                const auto cls = result.nodeClass(i);
                Vec3 q;
                if (cls.dim == 0) {
                    need(ownedVertices.contains(cls.tag), Reason::ProvenanceLost,
                         "foreign-sweep-vertex");
                    q = model.vertex(cls.tag)->position();
                } else if (cls.dim == 1) {
                    need(ownedEdges.contains(cls.tag), Reason::ProvenanceLost,
                         "foreign-sweep-edge");
                    q = model.edge(cls.tag)->closestPoint(p);
                } else if (cls.dim == 2) {
                    need(ownedFaces.contains(cls.tag), Reason::ProvenanceLost,
                         "foreign-sweep-face");
                    need(closestPointOnOcctFace(model, cls.tag, p, q), Reason::GeometryError,
                         "native-sweep-node-projection-failed");
                } else {
                    need(cls == mesh::Classification{3, body} && inside(p), Reason::GeometryError,
                         "native-sweep-interior-node-outside");
                    continue;
                }
                context.check();
                const auto error = (p - q).norm();
                out.maximumBoundaryNodeError = std::max(out.maximumBoundaryNodeError, error);
                need(std::isfinite(error) &&
                         error <= std::max(8 * model.tolerance(cls.dim, cls.tag), 1e-9),
                     Reason::GeometryError, "native-sweep-node-error");
            }
            std::set<int> covered;
            for (Index i = 0; i < result.cellCount(); ++i) {
                context.tick(1, out.validationWork);
                const auto &cell = result.cell(i);
                const auto *ids = result.cellNodes(i);
                if (cell.type == mesh::CellType::Tet4) {
                    need(cell.onEntity == mesh::Classification{3, body}, Reason::ProvenanceLost,
                         "foreign-sweep-volume");
                    need(inside((result.node(ids[0]) + result.node(ids[1]) + result.node(ids[2]) +
                                 result.node(ids[3])) *
                                .25),
                         Reason::GeometryError, "native-sweep-tet-centroid-outside");
                } else {
                    need(cell.type == mesh::CellType::Tri3 && cell.onEntity.dim == 2 &&
                             ownedFaces.contains(cell.onEntity.tag),
                         Reason::ProvenanceLost, "foreign-sweep-boundary");
                    covered.insert(cell.onEntity.tag);
                    const auto a = result.node(ids[0]), b = result.node(ids[1]),
                               c = result.node(ids[2]);
                    for (const auto p :
                         std::array{(a + b + c) / 3., (a + b) * .5, (b + c) * .5, (c + a) * .5}) {
                        context.tick(1, out.validationWork);
                        Vec3 q;
                        need(closestPointOnOcctFace(model, cell.onEntity.tag, p, q),
                             Reason::GeometryError, "native-sweep-facet-projection-failed");
                        context.check();
                        const auto error = (p - q).norm();
                        out.maximumBoundarySampleError =
                            std::max(out.maximumBoundarySampleError, error);
                        need(std::isfinite(error) &&
                                 error <= options.boundary.maximumChordError +
                                              model.tolerance(2, cell.onEntity.tag),
                             Reason::GeometryError, "native-sweep-facet-error");
                    }
                }
            }
            need(covered == ownedFaces, Reason::ProvenanceLost,
                 "native-sweep-omitted-original-face");
            auto limits = options.auditLimits;
            limits.maximumTemporaryBytes = std::min(
                limits.maximumTemporaryBytes, options.maximumWorkingBytes - validationWorking);
            limits.cancellation = options.cancellation;
            out.audit = mesh::auditLinearVolumeBoundary(result, limits);
            context.check();
            if (out.audit.cancelled)
                throw Stop{Reason::Cancelled, "final-audit-cancelled"};
            if (out.audit.budgetExceeded)
                throw Stop{Reason::BudgetExceeded, "final-audit-budget"};
            if (out.audit.allocationFailed)
                throw Stop{Reason::AllocationFailed, "final-audit-allocation-failed"};
            need(out.audit.valid, Reason::InvalidMesh, "native-sweep-final-audit");
            context.tick(result.cellCount(), out.validationWork);
            context.working(add(validationWorking, mul(result.cellCount(), sizeof(double))));
            out.quality = algo::analyzeMixedMeshQuality(result);
            context.check();
            need(!out.quality.allocationFailed, Reason::AllocationFailed,
                 "native-sweep-quality-allocation-failed");
            need(out.quality.complete && out.quality.valid, Reason::InvalidMesh,
                 "native-sweep-quality-invalid");
            need(out.quality.minimumQuality >= options.minimumQuality, Reason::QualityRequirement,
                 "quality-requirement-not-met");
            out.meshVolume = out.quality.totalVolume;
            out.relativeVolumeError = std::abs(out.meshVolume - out.cadVolume) / out.cadVolume;
            need(out.relativeVolumeError <= options.maximumRelativeVolumeError,
                 Reason::VolumeMismatch, "native-sweep-original-volume-error");
            need(std::abs(out.meshVolume - out.exteriorCertificate.signedVolume) <=
                     out.meshVolume * 1e-10,
                 Reason::VolumeMismatch, "native-sweep-exterior-volume-mismatch");
            context.tick(1, out.validationWork);
            need(before == snapshot(model, options.maximumSnapshotBytes), Reason::SourceMutation,
                 "native-sweep-source-mutated");
            out.sourceImmutabilityVerified = true;
            context.check();
        }
        out.nodes = result.nodeCount();
        out.cells = result.cellCount();
        destination.swap(result);
        out.published = true;
        out.destinationUnchanged = false;
        out.reason = Reason::None;
        out.diagnostic = "complete";
    } catch (const Stop &stop) {
        out.reason = stop.reason;
        out.diagnostic = stop.diagnostic;
    } catch (const Standard_OutOfMemory &) {
        out.reason = Reason::AllocationFailed;
        out.diagnostic = "native-allocation-failed";
    } catch (const std::bad_alloc &) {
        out.reason = Reason::AllocationFailed;
        out.diagnostic = "allocation-failed";
    } catch (const Standard_Failure &) {
        out.reason = Reason::KernelFailure;
        out.diagnostic = "native-kernel-failure";
    } catch (...) {
        out.reason = Reason::InternalError;
        out.diagnostic = "internal-error";
    }
    out.totalSeconds = std::chrono::duration<double>(Clock::now() - start).count();
    return out;
}
const char *occtNativeSweepReasonName(OcctNativeSweepReason r) noexcept {
    switch (r) {
    case Reason::None:
        return "none";
    case Reason::InvalidRequest:
        return "invalid_request";
    case Reason::UnsupportedGeometry:
        return "unsupported_geometry";
    case Reason::RecognitionFailed:
        return "recognition_failed";
    case Reason::SurfaceFailed:
        return "surface_failed";
    case Reason::ProvenanceLost:
        return "provenance_lost";
    case Reason::SweepFailed:
        return "sweep_failed";
    case Reason::GeometryError:
        return "geometry_error";
    case Reason::VolumeMismatch:
        return "volume_mismatch";
    case Reason::InvalidMesh:
        return "invalid_mesh";
    case Reason::QualityRequirement:
        return "quality_requirement_not_met";
    case Reason::BudgetExceeded:
        return "budget_exceeded";
    case Reason::Cancelled:
        return "cancelled";
    case Reason::AllocationFailed:
        return "allocation_failed";
    case Reason::SourceMutation:
        return "source_mutation";
    case Reason::KernelFailure:
        return "kernel_failure";
    case Reason::InternalError:
        return "internal_error";
    }
    return "internal_error";
}
} // namespace paacs::geo::detail
