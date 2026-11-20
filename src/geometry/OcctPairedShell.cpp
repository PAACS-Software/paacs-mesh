// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctPairedShell.hpp"
#include "OcctBoundaryInternal.hpp"
#include "OcctModelInternal.hpp"
#include "OcctNativeValidation.hpp"
#include "paacs/predicates/Predicates.hpp"
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <Standard_OutOfMemory.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <streambuf>

namespace paacs::geo::detail {
namespace {
using Reason = OcctPairedShellReason;
using Tri = std::array<Index, 3>;
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
    const OcctPairedShellOptions &options;
    OcctPairedShellOutcome &out;
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
struct Roles {
    int base = -1, top = -1, rim = -1, basePole = -1, topPole = -1, baseRim = -1, topRim = -1;
};
struct Disk {
    int face = -1, pole = -1, rim = -1, seam = -1;
    CadOrientation orientation;
};
Roles recognize(const OcctModel &m, int body, Context &c, std::set<int> &ownedFaces,
                std::set<int> &ownedEdges, std::set<int> &ownedVertices, std::uint64_t &working) {
    need(m.faceUses(body).size() == 3, Reason::UnsupportedGeometry,
         "paired-shell-needs-three-owned-faces");
    std::vector<Disk> disks;
    const auto sourceWorking = working;
    auto incidenceAdmission = [&] {
        c.working(add(sourceWorking,
                      mul(ownedFaces.size() + ownedEdges.size() + ownedVertices.size() + 8, 64)));
    };
    for (const auto &use : m.faceUses(body)) {
        c.tick(1, c.out.recognitionWork);
        incidenceAdmission();
        ownedFaces.insert(use.faceTag);
        need(m.faceOwners(use.faceTag).size() == 1 && m.faceOwners(use.faceTag)[0].solidTag == body,
             Reason::UnsupportedGeometry, "paired-shell-has-shared-original-face");
        BRepAdaptor_Surface surface(OcctModelAccess::face(m, use.faceTag));
        Disk d;
        d.face = use.faceTag;
        d.orientation = use.orientation;
        bool candidate = true;
        for (const auto &w : m.boundaryUses(use.faceTag))
            for (const auto &e : w.edges) {
                c.tick(1, c.out.recognitionWork);
                incidenceAdmission();
                ownedEdges.insert(e.edgeTag);
                if (e.startVertexTag > 0)
                    ownedVertices.insert(e.startVertexTag);
                if (e.endVertexTag > 0)
                    ownedVertices.insert(e.endVertexTag);
                if (e.degenerated) {
                    candidate = candidate && e.startVertexTag > 0 &&
                                e.startVertexTag == e.endVertexTag &&
                                (d.pole < 0 || d.pole == e.startVertexTag);
                    d.pole = e.startVertexTag;
                } else if (e.seam) {
                    candidate = candidate && (d.seam < 0 || d.seam == e.edgeTag);
                    d.seam = e.edgeTag;
                } else {
                    candidate = candidate && (d.rim < 0 || d.rim == e.edgeTag);
                    d.rim = e.edgeTag;
                }
            }
        candidate = candidate && d.pole > 0 && d.rim > 0 && d.seam > 0 &&
                    m.boundaryUses(d.face).size() == 1 &&
                    (surface.GetType() == GeomAbs_BSplineSurface ||
                     surface.GetType() == GeomAbs_Sphere || surface.GetType() == GeomAbs_Cone);
        if (candidate) {
            const auto *edge = m.edge(d.rim);
            candidate = edge && edge->startVertexTag() > 0 &&
                        edge->startVertexTag() == edge->endVertexTag();
        }
        if (candidate)
            disks.push_back(d);
    }
    need(disks.size() == 2, Reason::UnsupportedGeometry,
         "paired-shell-native-disk-topology-not-supported");
    if (disks[0].orientation == CadOrientation::Forward)
        std::swap(disks[0], disks[1]);
    need(disks[0].orientation == CadOrientation::Reversed &&
             disks[1].orientation == CadOrientation::Forward && disks[0].pole != disks[1].pole &&
             disks[0].rim != disks[1].rim,
         Reason::UnsupportedGeometry, "paired-shell-opposed-sheet-ownership-not-supported");
    Roles r{disks[0].face, disks[1].face, -1,          disks[0].pole,
            disks[1].pole, disks[0].rim,  disks[1].rim};
    for (int face : ownedFaces)
        if (face != r.base && face != r.top)
            r.rim = face;
    BRepAdaptor_Surface lip(OcctModelAccess::face(m, r.rim));
    need(lip.GetType() == GeomAbs_Plane || lip.GetType() == GeomAbs_SurfaceOfExtrusion,
         Reason::UnsupportedGeometry, "paired-shell-connector-support-not-supported");
    std::set<int> rimEdges;
    for (const auto &w : m.boundaryUses(r.rim))
        for (const auto &e : w.edges) {
            c.tick(1, c.out.recognitionWork);
            if (!e.degenerated)
                rimEdges.insert(e.edgeTag);
        }
    need(rimEdges == std::set<int>{r.baseRim, r.topRim}, Reason::UnsupportedGeometry,
         "paired-shell-connector-incidence-not-supported");
    working = add(sourceWorking,
                  mul(ownedFaces.size() + ownedEdges.size() + ownedVertices.size() + 8, 64));
    c.working(working);
    return r;
}
Vec3 curveClosest(const OcctModel &m, int tag, Vec3 p, Context &c, std::uint64_t &phase) {
    c.tick(1, phase);
    BRepExtrema_DistShapeShape d(BRepBuilderAPI_MakeVertex(gp_Pnt(p.x, p.y, p.z)).Vertex(),
                                 OcctModelAccess::edge(m, tag));
    c.check();
    need(d.IsDone() && d.NbSolution() > 0, Reason::GeometryError,
         "paired-shell-native-curve-projection-failed");
    const auto q = d.PointOnShape2(1);
    return {q.X(), q.Y(), q.Z()};
}
std::map<Pair, unsigned> edgeCounts(const std::vector<Tri> &triangles, Context &c) {
    std::map<Pair, unsigned> counts;
    for (const auto &t : triangles)
        for (unsigned j = 0; j < 3; ++j) {
            c.tick(1, c.out.sweepWork);
            auto &n = counts[pairKey(t[j], t[(j + 1) % 3])];
            need(++n <= 2, Reason::InvalidMesh, "paired-shell-cap-edge-incidence");
        }
    return counts;
}
void childFailure(OcctBoundaryReason reason, const char *diagnostic) {
    if (reason == OcctBoundaryReason::Cancelled)
        throw Stop{Reason::Cancelled, diagnostic};
    if (reason == OcctBoundaryReason::AllocationFailed)
        throw Stop{Reason::AllocationFailed, diagnostic};
    if (reason == OcctBoundaryReason::BudgetExceeded)
        throw Stop{Reason::BudgetExceeded, diagnostic};
    if (reason == OcctBoundaryReason::SourceChanged)
        throw Stop{Reason::SourceMutation, diagnostic};
    throw Stop{Reason::SurfaceFailed, diagnostic};
}
} // namespace
OcctPairedShellOutcome buildOcctPairedShellVolume(const OcctModel &m, int body,
                                                  const OcctPairedShellOptions &o,
                                                  mesh::Mesh &destination) noexcept {
    const auto start = Clock::now();
    OcctPairedShellOutcome out;
    Context context{o, out};
    std::string before;
    bool snapshotReady = false;
    try {
        need(m.solid(body) && o.layers >= 1 && o.layers <= 4096 &&
                 o.maximumNativeValidationThreads >= 1 && o.maximumNativeValidationThreads <= 8 &&
                 std::isfinite(o.maximumLayerThickness) && o.maximumLayerThickness >= 0 &&
                 std::isfinite(o.maximumRelativeRecognitionError) &&
                 o.maximumRelativeRecognitionError > 0 &&
                 o.maximumRelativeRecognitionError <= 1e-6 &&
                 std::isfinite(o.maximumRelativeVolumeError) && o.maximumRelativeVolumeError > 0 &&
                 o.maximumRelativeVolumeError < 1 && std::isfinite(o.minimumQuality) &&
                 o.minimumQuality >= 0 && o.minimumQuality <= 1 &&
                 std::isfinite(o.warningQuality) && o.warningQuality >= 0 &&
                 o.warningQuality <= 1 && std::isfinite(o.boundary.linearDeflection) &&
                 o.boundary.linearDeflection > 0 && std::isfinite(o.boundary.maximumChordError) &&
                 o.boundary.maximumChordError > 0,
             Reason::InvalidRequest, "invalid-paired-shell-request");
        need(o.maximumWork && o.maximumNodes && o.maximumCells && o.maximumOutputBytes &&
                 o.maximumWorkingBytes && o.maximumSnapshotBytes,
             Reason::BudgetExceeded, "zero-paired-shell-budget");
        context.check();
        need(o.maximumSnapshotBytes <= o.maximumWorkingBytes / 8, Reason::BudgetExceeded,
             "paired-shell-snapshot-working-admission");
        context.tick(1, out.recognitionWork);
        before = snapshot(m, o.maximumSnapshotBytes);
        snapshotReady = true;
        auto fixedWorking = mul(before.size(), 4);
        context.working(fixedWorking);
        std::set<int> ownedFaces, ownedEdges, ownedVertices;
        Roles roles;
        {
            Timer timer{out.recognitionSeconds};
            roles =
                recognize(m, body, context, ownedFaces, ownedEdges, ownedVertices, fixedWorking);
            out.cadVolume = mass(OcctModelAccess::solid(m, body), context);
        }
        out.baseCapTag = roles.base;
        out.targetCapTag = roles.top;
        auto bo = o.boundary;
        bo.maximumWork = std::min(bo.maximumWork, context.remaining());
        bo.maximumNodes = std::min(bo.maximumNodes, o.maximumNodes);
        bo.maximumTriangles = std::min(bo.maximumTriangles, o.maximumCells);
        bo.maximumOutputBytes = std::min({bo.maximumOutputBytes, o.maximumOutputBytes,
                                          (o.maximumWorkingBytes - fixedWorking) / 3});
        bo.cancellation = o.cancellation;
        need(bo.maximumWork && bo.maximumNodes && bo.maximumTriangles && bo.maximumOutputBytes,
             Reason::BudgetExceeded, "paired-shell-surface-admission-budget");
        mesh::Mesh closed;
        {
            Timer timer{out.surfaceSeconds};
            const auto prepared = buildOcctBoundaryCandidate(m, body, bo, closed);
            context.collect(prepared.work, out.surfaceWork);
            out.boundaryReason = prepared.reason;
            if (!prepared.published || !prepared.sourceModelUnmodified)
                childFailure(prepared.reason, "paired-shell-source-surface-failed");
        }
        mesh::Mesh patch;
        std::vector<Tri> triangles;
        std::map<Index, Index> local;
        std::vector<Vec3> targets;
        std::vector<mesh::Classification> targetClasses;
        std::map<Pair, unsigned> counts;
        Index inserted = -1;
        const int targetLip = m.edge(roles.topRim)->startVertexTag();
        std::uint64_t patchWorking = 0;
        {
            Timer timer{out.sweepSeconds};
            context.working(
                add(fixedWorking,
                    add(payload(closed), mul(closed.nodeCount() + closed.cellCount(), 160))));
            for (Index cell = 0; cell < closed.cellCount(); ++cell) {
                context.tick(1, out.sweepWork);
                if (closed.cell(cell).onEntity == mesh::Classification{2, roles.base}) {
                    Tri t;
                    const auto *p = closed.cellNodes(cell);
                    for (unsigned j = 0; j < 3; ++j) {
                        auto [it, fresh] = local.try_emplace(p[j], patch.nodeCount());
                        if (fresh)
                            patch.addNode(closed.node(p[j]), closed.nodeClass(p[j]));
                        t[j] = it->second;
                    }
                    triangles.push_back(t);
                }
            }
            need(!triangles.empty(), Reason::SurfaceFailed, "paired-shell-empty-base-cap");
            {
                mesh::Mesh empty;
                closed.swap(empty);
            }
            local.clear();
            patchWorking = add(fixedWorking,
                               add(payload(patch), mul(triangles.size() + patch.nodeCount(), 256)));
            context.working(patchWorking);
            counts = edgeCounts(triangles, context);
            const auto desiredTarget = m.vertex(targetLip)->position();
            const auto desiredSource =
                curveClosest(m, roles.baseRim, desiredTarget, context, out.sweepWork);
            const auto *curve = m.edge(roles.baseRim);
            double lo, hi;
            curve->parameterBounds(lo, hi);
            const double period = hi - lo;
            need(std::isfinite(period) && period > 0, Reason::GeometryError,
                 "paired-shell-rim-parameter-range");
            context.tick(1, out.sweepWork);
            const double tp = curve->closestParameter(desiredSource);
            need(std::isfinite(tp), Reason::GeometryError, "paired-shell-rim-preimage-parameter");
            // Native curve parameters select correspondence to an existing
            // source node; original source and target identities stay distinct.
            for (Index i = 0; i < patch.nodeCount() && inserted < 0; ++i) {
                context.tick(1, out.sweepWork);
                auto cls = patch.nodeClass(i);
                if ((cls.dim == 1 && cls.tag == roles.baseRim) ||
                    (cls.dim == 0 && cls.tag == curve->startVertexTag())) {
                    const double t = curve->closestParameter(patch.node(i));
                    double delta = std::abs(t - tp);
                    delta = std::min(delta, std::abs(period - delta));
                    if (delta <= period * 1e-12 &&
                        (patch.node(i) - desiredSource).norm() <= m.tolerance(cls.dim, cls.tag))
                        inserted = i;
                }
            }
            if (inserted < 0) {
                bool split = false;
                for (std::size_t i = 0; i < triangles.size() && !split; ++i) {
                    const auto t = triangles[i];
                    for (unsigned j = 0; j < 3 && !split; ++j) {
                        context.tick(1, out.sweepWork);
                        const auto a = t[j], b = t[(j + 1) % 3], v = t[(j + 2) % 3];
                        if (counts[pairKey(a, b)] != 1)
                            continue;
                        context.tick(2, out.sweepWork);
                        double ta = curve->closestParameter(patch.node(a)),
                               tb = curve->closestParameter(patch.node(b)), q = tp;
                        need(std::isfinite(ta) && std::isfinite(tb), Reason::GeometryError,
                             "paired-shell-boundary-parameter");
                        if (std::abs(ta - tb) > period * .5) {
                            if (ta < tb)
                                ta += period;
                            else
                                tb += period;
                            if (q < std::min(ta, tb))
                                q += period;
                        }
                        if (q <= std::min(ta, tb) + period * 1e-12 ||
                            q >= std::max(ta, tb) - period * 1e-12)
                            continue;
                        context.output(mul(patch.nodeCount() + 1, 2), mul(triangles.size() + 1, 5),
                                       mul(triangles.size() + 1, 24));
                        inserted = patch.addNode(desiredSource, {1, roles.baseRim});
                        triangles[i] = {a, inserted, v};
                        triangles.push_back({inserted, b, v});
                        split = true;
                    }
                }
                need(split, Reason::GeometryError,
                     "paired-shell-target-vertex-preimage-not-resolved");
                counts = edgeCounts(triangles, context);
            }
            targets.reserve(patch.nodeCount());
            targetClasses.reserve(patch.nodeCount());
            double maximumGap = 0;
            for (Index i = 0; i < patch.nodeCount(); ++i) {
                context.tick(1, out.sweepWork);
                const auto cls = patch.nodeClass(i);
                Vec3 p;
                mesh::Classification target;
                if (i == inserted) {
                    p = desiredTarget;
                    target = {0, targetLip};
                } else if (cls == mesh::Classification{0, roles.basePole}) {
                    p = m.vertex(roles.topPole)->position();
                    target = {0, roles.topPole};
                } else if ((cls.dim == 1 && cls.tag == roles.baseRim) ||
                           (cls.dim == 0 && cls.tag != roles.basePole)) {
                    p = curveClosest(m, roles.topRim, patch.node(i), context, out.sweepWork);
                    target = {1, roles.topRim};
                } else {
                    need(closestPointOnOcctFace(m, roles.top, patch.node(i), p),
                         Reason::GeometryError, "paired-shell-target-projection-failed");
                    target = {2, roles.top};
                    context.check();
                }
                need(finite(p), Reason::GeometryError, "paired-shell-nonfinite-target");
                maximumGap = std::max(maximumGap, (p - patch.node(i)).norm());
                targets.push_back(p);
                targetClasses.push_back(target);
            }
            out.height = maximumGap;
            const auto poleDirection =
                m.vertex(roles.topPole)->position() - m.vertex(roles.basePole)->position();
            need(poleDirection.norm() > 0, Reason::UnsupportedGeometry,
                 "paired-shell-coincident-poles");
            out.direction = poleDirection / poleDirection.norm();
            out.layers = o.layers;
            if (o.maximumLayerThickness > 0) {
                const auto required = std::ceil(maximumGap / o.maximumLayerThickness);
                need(std::isfinite(required) && required <= 4096, Reason::BudgetExceeded,
                     "paired-shell-layer-admission");
                out.layers = std::max(out.layers, static_cast<unsigned>(std::max(1., required)));
            }
            out.capNodes = patch.nodeCount();
            out.capTriangles = triangles.size();
        }
        mesh::Mesh result;
        const auto n = patch.nodeCount();
        const auto levels = out.layers + 1;
        const std::uint64_t tetCount = mul(mul(triangles.size(), 3), out.layers);
        std::uint64_t boundaryEdges = 0;
        for (const auto &[e, count] : counts)
            if (count == 1)
                ++boundaryEdges;
        const auto boundaryTriangles =
            add(mul(triangles.size(), 2), mul(mul(boundaryEdges, 2), out.layers));
        const auto totalCells = add(tetCount, boundaryTriangles), totalNodes = mul(n, levels),
                   connectivity = add(mul(tetCount, 4), mul(boundaryTriangles, 3));
        context.output(totalNodes, totalCells, connectivity);
        need(totalNodes <= std::uint64_t((std::numeric_limits<Index>::max)()) &&
                 totalCells <= std::uint64_t((std::numeric_limits<Index>::max)()),
             Reason::BudgetExceeded, "paired-shell-index-admission");
        patchWorking =
            add(fixedWorking,
                add(payload(patch), add(mul(triangles.size(), sizeof(Tri)),
                                        add(mul(counts.size(), 128),
                                            mul(n, sizeof(Vec3) + sizeof(mesh::Classification))))));
        const auto resultWorking = add(patchWorking, payload(totalNodes, totalCells, connectivity));
        context.working(add(resultWorking, mul(mul(tetCount, 4), 160)));
        {
            Timer timer{out.sweepSeconds};
            result.reserveAppend(totalNodes, totalCells, connectivity);
            std::set<Index> rimNodes;
            for (const auto &[e, count] : counts)
                if (count == 1) {
                    rimNodes.insert(e.first);
                    rimNodes.insert(e.second);
                }
            for (unsigned level = 0; level < levels; ++level)
                for (Index i = 0; i < n; ++i) {
                    context.tick(1, out.sweepWork);
                    Vec3 p;
                    mesh::Classification cls;
                    if (level == 0) {
                        p = patch.node(i);
                        cls = patch.nodeClass(i);
                    } else if (level == out.layers) {
                        p = targets[i];
                        cls = targetClasses[i];
                    } else {
                        p = patch.node(i) +
                            (targets[i] - patch.node(i)) * (double(level) / out.layers);
                        cls = {3, body};
                        if (rimNodes.contains(i)) {
                            Vec3 q;
                            need(closestPointOnOcctFace(m, roles.rim, p, q), Reason::GeometryError,
                                 "paired-shell-layer-lip-projection-failed");
                            p = q;
                            cls = {2, roles.rim};
                            context.check();
                        }
                    }
                    result.addNode(p, cls);
                }
            auto tet = [&](std::array<Index, 4> ids) {
                context.tick(1, out.sweepWork);
                const auto sign = predicates::orient3d(result.node(ids[0]), result.node(ids[1]),
                                                       result.node(ids[2]), result.node(ids[3]));
                need(sign != predicates::Sign::Zero, Reason::InvalidMesh,
                     "paired-shell-collapsed-tet");
                if (sign != predicates::Sign::Positive)
                    std::swap(ids[1], ids[2]);
                result.addCell(mesh::CellType::Tet4, ids, {3, body});
            };
            // Consistent pulling diagonals depend only on source node indices.
            // The resulting oriented simplex complex is subsequently audited;
            // no claim of convex/planar Prism6 cells is made for curved sheets.
            for (unsigned level = 0; level < out.layers; ++level)
                for (auto t : triangles) {
                    std::sort(t.begin(), t.end());
                    const auto [a, b, v] = t;
                    const Index low = level * n, high = (level + 1) * n;
                    tet({a + low, b + low, v + low, v + high});
                    tet({a + low, b + low, v + high, b + high});
                    tet({a + low, a + high, b + high, v + high});
                }
            std::map<Tri, std::pair<Tri, unsigned>> faces;
            const std::array<Tri, 4> order{Tri{0, 2, 1}, Tri{0, 1, 3}, Tri{0, 3, 2}, Tri{1, 2, 3}};
            for (Index cell = 0; cell < result.cellCount(); ++cell) {
                const auto *p = result.cellNodes(cell);
                for (auto o : order) {
                    context.tick(1, out.sweepWork);
                    Tri t{p[o[0]], p[o[1]], p[o[2]]}, key = t;
                    std::sort(key.begin(), key.end());
                    auto [it, fresh] = faces.try_emplace(key, std::pair<Tri, unsigned>{t, 0});
                    need(++it->second.second <= 2, Reason::InvalidMesh,
                         "paired-shell-interior-incidence");
                }
            }
            for (const auto &[key, use] : faces)
                if (use.second == 1) {
                    context.tick(1, out.sweepWork);
                    int face;
                    if (key[2] < n)
                        face = roles.base;
                    else if (key[0] >= out.layers * n)
                        face = roles.top;
                    else {
                        std::set<Index> base;
                        for (auto id : key)
                            base.insert(id % n);
                        need(base.size() == 2, Reason::InvalidMesh,
                             "paired-shell-mixed-exterior-facet");
                        auto it = base.begin();
                        const auto a = *it++, b = *it;
                        need(counts[pairKey(a, b)] == 1, Reason::InvalidMesh,
                             "paired-shell-side-facet-not-on-original-rim");
                        face = roles.rim;
                    }
                    result.addCell(mesh::CellType::Tri3, use.first, {2, face});
                }
            need(result.cellCount() == totalCells, Reason::InvalidMesh,
                 "paired-shell-exterior-count-mismatch");
        }
        {
            Timer timer{out.validationSeconds};
            mesh::Mesh exterior;
            std::map<Index, Index> surfaceNodes;
            context.working(add(resultWorking, add(payload(totalNodes, boundaryTriangles,
                                                           mul(boundaryTriangles, 3)),
                                                   mul(totalNodes, 128))));
            for (Index cell = 0; cell < result.cellCount(); ++cell) {
                context.tick(1, out.validationWork);
                if (result.cell(cell).type == mesh::CellType::Tri3) {
                    Tri ids;
                    const auto *p = result.cellNodes(cell);
                    for (unsigned j = 0; j < 3; ++j) {
                        auto [it, fresh] = surfaceNodes.try_emplace(p[j], exterior.nodeCount());
                        if (fresh)
                            exterior.addNode(result.node(p[j]), result.nodeClass(p[j]));
                        ids[j] = it->second;
                    }
                    exterior.addCell(mesh::CellType::Tri3, ids, result.cell(cell).onEntity);
                }
            }
            const auto validationWorking =
                add(resultWorking, add(payload(exterior), mul(surfaceNodes.size(), 128)));
            context.working(validationWorking);
            auto limits = o.auditLimits;
            limits.maximumTemporaryBytes =
                std::min(limits.maximumTemporaryBytes, o.maximumWorkingBytes - validationWorking);
            limits.cancellation = o.cancellation;
            context.tick(add(result.nodeCount(), mul(result.cellCount(), 5)), out.validationWork);
            out.audit = mesh::auditLinearVolumeBoundary(result, limits);
            context.check();
            if (out.audit.cancelled)
                throw Stop{Reason::Cancelled, "paired-shell-audit-cancelled"};
            if (out.audit.budgetExceeded)
                throw Stop{Reason::BudgetExceeded, "paired-shell-audit-budget"};
            if (out.audit.allocationFailed)
                throw Stop{Reason::AllocationFailed, "paired-shell-audit-allocation"};
            need(out.audit.complete && out.audit.valid, Reason::InvalidMesh,
                 "paired-shell-full-oriented-volume-audit");
            OcctStarOptions cert;
            cert.maximumWork = context.remaining();
            cert.maximumBoundaryNodes = o.maximumNodes;
            cert.maximumBoundaryTriangles = o.maximumCells;
            cert.maximumNodes = o.maximumNodes;
            cert.maximumCells = o.maximumCells;
            cert.maximumOutputBytes = o.maximumOutputBytes;
            cert.maximumWorkingBytes = o.maximumWorkingBytes - validationWorking;
            cert.cancellation = o.cancellation;
            cert.auditLimits = limits;
            need(cert.maximumWork && cert.maximumWorkingBytes, Reason::BudgetExceeded,
                 "paired-shell-exterior-certificate-admission");
            out.exteriorCertificate = certifyOcctBoundarySurface(exterior, cert);
            context.collect(out.exteriorCertificate.workUsed, out.validationWork);
            if (!out.exteriorCertificate.certified) {
                const auto reason = out.exteriorCertificate.reason;
                if (reason == OcctStarReason::Cancelled)
                    throw Stop{Reason::Cancelled, out.exteriorCertificate.diagnostic};
                if (reason == OcctStarReason::AllocationFailed)
                    throw Stop{Reason::AllocationFailed, out.exteriorCertificate.diagnostic};
                if (reason == OcctStarReason::WorkBudget || reason == OcctStarReason::OutputBudget)
                    throw Stop{Reason::BudgetExceeded, out.exteriorCertificate.diagnostic};
                throw Stop{Reason::InvalidMesh, out.exteriorCertificate.diagnostic};
            }
            bo.maximumWork = std::min(o.boundary.maximumWork, context.remaining());
            need(bo.maximumWork, Reason::BudgetExceeded, "paired-shell-facet-work-admission");
            const auto facet =
                auditOcctBoundaryFacets(m, body, exterior, bo, o.maximumSnapshotBytes);
            context.collect(facet.work, out.validationWork);
            if (!facet.completed || !facet.sourceModelUnmodified)
                childFailure(facet.reason, "paired-shell-original-facet-audit-failed");
            out.maximumBoundarySampleError = facet.maximumFacetSampleError;
            need(out.maximumBoundarySampleError <= o.boundary.maximumChordError,
                 Reason::GeometryError, "paired-shell-original-facet-error");
            std::set<int> retainedVertices, coveredFaces;
            std::vector<Vec3> interior;
            context.working(add(validationWorking, mul(totalNodes + tetCount, sizeof(Vec3))));
            interior.reserve(totalNodes + tetCount);
            for (Index i = 0; i < result.nodeCount(); ++i) {
                context.tick(1, out.validationWork);
                const auto cls = result.nodeClass(i);
                const auto p = result.node(i);
                interior.push_back(p);
                Vec3 q;
                if (cls.dim == 0) {
                    need(ownedVertices.contains(cls.tag), Reason::ProvenanceLost,
                         "paired-shell-foreign-vertex");
                    retainedVertices.insert(cls.tag);
                    q = m.vertex(cls.tag)->position();
                } else if (cls.dim == 1) {
                    need(ownedEdges.contains(cls.tag), Reason::ProvenanceLost,
                         "paired-shell-foreign-edge");
                    q = curveClosest(m, cls.tag, p, context, out.validationWork);
                } else if (cls.dim == 2) {
                    need(ownedFaces.contains(cls.tag), Reason::ProvenanceLost,
                         "paired-shell-foreign-face");
                    need(closestPointOnOcctFace(m, cls.tag, p, q), Reason::GeometryError,
                         "paired-shell-original-node-projection-failed");
                    context.check();
                } else {
                    need(cls == mesh::Classification{3, body}, Reason::ProvenanceLost,
                         "paired-shell-foreign-body-node");
                    continue;
                }
                const auto error = (p - q).norm();
                out.maximumBoundaryNodeError = std::max(out.maximumBoundaryNodeError, error);
                need(std::isfinite(error) && error <= m.tolerance(cls.dim, cls.tag),
                     Reason::GeometryError, "paired-shell-original-node-tolerance");
            }
            need(retainedVertices == ownedVertices, Reason::ProvenanceLost,
                 "paired-shell-original-vertex-coverage");
            for (Index cell = 0; cell < result.cellCount(); ++cell) {
                context.tick(1, out.validationWork);
                const auto &c = result.cell(cell);
                const auto *p = result.cellNodes(cell);
                if (c.type == mesh::CellType::Tet4) {
                    need(c.onEntity == mesh::Classification{3, body}, Reason::ProvenanceLost,
                         "paired-shell-foreign-volume-cell");
                    interior.push_back((result.node(p[0]) + result.node(p[1]) + result.node(p[2]) +
                                        result.node(p[3])) *
                                       .25);
                } else {
                    need(c.type == mesh::CellType::Tri3 && c.onEntity.dim == 2 &&
                             ownedFaces.contains(c.onEntity.tag),
                         Reason::ProvenanceLost, "paired-shell-foreign-boundary-cell");
                    coveredFaces.insert(c.onEntity.tag);
                }
            }
            need(coveredFaces == ownedFaces, Reason::ProvenanceLost,
                 "paired-shell-original-face-coverage");
            OcctNativeValidationOptions native;
            native.maximumThreads = o.maximumNativeValidationThreads;
            native.maximumWork = context.remaining();
            native.maximumPoints = interior.size();
            native.maximumWorkingBytes = o.maximumWorkingBytes - validationWorking -
                                         mul(totalNodes + tetCount, sizeof(Vec3));
            native.maximumSnapshotBytes =
                std::min(o.maximumSnapshotBytes, native.maximumWorkingBytes / 4);
            native.cancellation = o.cancellation;
            need(native.maximumWork && native.maximumWorkingBytes && native.maximumSnapshotBytes,
                 Reason::BudgetExceeded, "paired-shell-native-validation-admission");
            const auto checked = validateOcctNativeInteriorPoints(m, body, interior, native);
            context.collect(checked.workUsed, out.validationWork);
            out.nativeValidationThreads = checked.threadsUsed;
            out.nativeInteriorQueries = checked.nativeQueries + checked.originalReferenceQueries;
            if (!checked.completed || !checked.sourceModelUnmodified) {
                if (checked.reason == OcctNativeValidationReason::Cancelled)
                    throw Stop{Reason::Cancelled, "paired-shell-native-validation-cancelled"};
                if (checked.reason == OcctNativeValidationReason::AllocationFailed)
                    throw Stop{Reason::AllocationFailed,
                               "paired-shell-native-validation-allocation"};
                if (checked.reason == OcctNativeValidationReason::BudgetExceeded)
                    throw Stop{Reason::BudgetExceeded, "paired-shell-native-validation-budget"};
                if (checked.reason == OcctNativeValidationReason::SourceChanged)
                    throw Stop{Reason::SourceMutation,
                               "paired-shell-native-validation-source-changed"};
                if (checked.reason != OcctNativeValidationReason::OutsideMaterial)
                    throw Stop{Reason::KernelFailure, "paired-shell-native-validation-failed"};
                throw Stop{Reason::GeometryError,
                           "paired-shell-original-node-or-tet-centroid-outside"};
            }
            std::vector<Vec3>().swap(interior);
            context.tick(result.cellCount(), out.validationWork);
            context.working(add(validationWorking, mul(result.cellCount(), sizeof(double))));
            out.quality = algo::analyzeMixedMeshQuality(result);
            context.check();
            need(!out.quality.allocationFailed, Reason::AllocationFailed,
                 "paired-shell-quality-allocation-failed");
            need(out.quality.complete && out.quality.valid, Reason::InvalidMesh,
                 "paired-shell-quality-invalid");
            need(out.quality.minimumQuality >= o.minimumQuality, Reason::QualityRequirement,
                 "paired-shell-quality-requirement-not-met");
            out.meshVolume = out.quality.totalVolume;
            out.relativeVolumeError = std::abs(out.meshVolume - out.cadVolume) / out.cadVolume;
            need(out.relativeVolumeError <= o.maximumRelativeVolumeError, Reason::VolumeMismatch,
                 "paired-shell-original-gk-volume-error");
            need(std::abs(out.meshVolume - out.exteriorCertificate.signedVolume) <=
                     out.meshVolume * 1e-10,
                 Reason::VolumeMismatch, "paired-shell-exterior-volume-mismatch");
            context.tick(1, out.validationWork);
            need(before == snapshot(m, o.maximumSnapshotBytes), Reason::SourceMutation,
                 "paired-shell-source-mutated");
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
    } catch (const Stop &s) {
        out.reason = s.reason;
        out.diagnostic = s.diagnostic;
    } catch (const Standard_OutOfMemory &) {
        out.reason = Reason::AllocationFailed;
        out.diagnostic = "paired-shell-native-allocation-failed";
    } catch (const std::bad_alloc &) {
        out.reason = Reason::AllocationFailed;
        out.diagnostic = "paired-shell-allocation-failed";
    } catch (const Standard_Failure &) {
        out.reason = Reason::KernelFailure;
        out.diagnostic = "paired-shell-native-kernel-failure";
    } catch (...) {
        out.reason = Reason::InternalError;
        out.diagnostic = "paired-shell-internal-error";
    }
    if (!out.published && snapshotReady) {
        try {
            out.sourceImmutabilityVerified = before == snapshot(m, o.maximumSnapshotBytes);
            if (!out.sourceImmutabilityVerified) {
                out.reason = Reason::SourceMutation;
                out.diagnostic = "paired-shell-source-mutated";
            }
        } catch (...) {
        }
    }
    out.totalSeconds = std::chrono::duration<double>(Clock::now() - start).count();
    return out;
}
} // namespace paacs::geo::detail
