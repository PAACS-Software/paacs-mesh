// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctTrimTolerance.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <Standard_OutOfMemory.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <cmath>
#include <map>
#include <ostream>
#include <set>
#include <streambuf>

namespace paacs::geo::detail {
namespace {
using Reason = OcctTrimToleranceReason;
struct Stop {
    Reason reason;
};
void need(bool ok, Reason reason) {
    if (!ok)
        throw Stop{reason};
}
class Buffer final : public std::streambuf {
    std::uint64_t limit_;

  public:
    std::string value;
    explicit Buffer(std::uint64_t n) : limit_(n) {}

  protected:
    std::streamsize xsputn(const char *p, std::streamsize n) override {
        need(n >= 0 && std::uint64_t(n) <= limit_ - value.size(), Reason::BudgetExceeded);
        value.append(p, std::size_t(n));
        return n;
    }
    int_type overflow(int_type c) override {
        if (traits_type::eq_int_type(c, traits_type::eof()))
            return traits_type::not_eof(c);
        char ch = traits_type::to_char_type(c);
        xsputn(&ch, 1);
        return c;
    }
};
std::string snapshot(const TopoDS_Shape &s, std::uint64_t limit) {
    Buffer b(limit);
    std::ostream out(&b);
    out.exceptions(std::ios::badbit | std::ios::failbit);
    BRepTools::Write(s, out);
    return std::move(b.value);
}
} // namespace
const char *occtTrimToleranceReasonName(OcctTrimToleranceReason r) noexcept {
    switch (r) {
    case Reason::None:
        return "none";
    case Reason::InvalidRequest:
        return "invalid_request";
    case Reason::SupportMismatch:
        return "copied_support_changed";
    case Reason::SourceCoverageLost:
        return "original_topology_not_retained";
    case Reason::OriginalToleranceExceeded:
        return "original_trim_tolerance_exceeded";
    case Reason::DriftEnvelopeExceeded:
        return "source_tolerance_volume_envelope_exceeded";
    case Reason::NativeQueryFailed:
        return "native_query_failed";
    case Reason::BudgetExceeded:
        return "budget_exceeded";
    case Reason::Cancelled:
        return "cancelled";
    case Reason::AllocationFailed:
        return "allocation_failed";
    case Reason::SourceMutation:
        return "source_mutation";
    }
    return "unknown";
}
OcctTrimToleranceReceipt auditOcctTrimTolerance(const OcctTrimToleranceInput &i,
                                                const OcctTrimToleranceOptions &o) noexcept {
    OcctTrimToleranceReceipt r;
    r.attempted = true;
    r.maximumPhysicalDrift = o.maximumPhysicalDrift;
    try {
        need(std::isfinite(o.maximumPhysicalDrift) && o.maximumPhysicalDrift > 0 &&
                 std::isfinite(i.sourceVolume) && i.sourceVolume > 0 &&
                 std::isfinite(i.derivedVolume) && i.derivedVolume > 0 &&
                 !i.selectedSourceBodies.empty() &&
                 i.parentByFace.size() == std::size_t(i.derivedFaces.Extent() + 1) &&
                 i.parentByEdge.size() == std::size_t(i.derivedEdges.Extent() + 1) &&
                 i.parentByVertex.size() == std::size_t(i.derivedVertices.Extent() + 1) &&
                 i.faceOwners.size() == std::size_t(i.derivedFaces.Extent() + 1),
             Reason::InvalidRequest);
        r.nativeVolumeDrift = std::abs(i.derivedVolume - i.sourceVolume);
        auto tick = [&](std::uint64_t n = 1) {
            need(!o.cancellation || !o.cancellation->load(std::memory_order_relaxed),
                 Reason::Cancelled);
            need(n <= o.maximumWork - r.work, Reason::BudgetExceeded);
            r.work += n;
        };
        tick(0);
        need(o.maximumWork > 0, Reason::BudgetExceeded);
        const std::uint64_t entries = std::uint64_t(i.sourceFaces.Extent()) +
                                      i.sourceEdges.Extent() + i.sourceVertices.Extent() +
                                      i.derivedFaces.Extent() + i.derivedEdges.Extent() +
                                      i.derivedVertices.Extent();
        need(entries <= o.maximumLogicalBytes / 64, Reason::BudgetExceeded);
        const auto snapshotLimit =
            std::min(o.maximumSnapshotBytes, (o.maximumLogicalBytes - entries * 64) / 2);
        const auto before = snapshot(i.originalShape, snapshotLimit);
        std::set<int> expectedFaces, expectedEdges, expectedVertices, observedFaces, observedEdges,
            observedVertices;
        std::set<int> selectedBodies;
        std::map<int, std::set<int>> sourceFacesByBody;
        std::set<int> boundaryEdges, boundaryVertices;
        for (int body : i.selectedSourceBodies) {
            tick();
            need(body > 0 && body <= i.sourceSolids.Extent(), Reason::InvalidRequest);
            need(selectedBodies.insert(body).second, Reason::InvalidRequest);
            OcctWitnessShapeMap children;
            TopExp::MapShapes(i.sourceSolids(body), TopAbs_FACE, children);
            for (int n = 1; n <= children.Extent(); ++n) {
                tick();
                const int tag = i.sourceFaces.FindIndex(children(n));
                need(tag > 0, Reason::SourceCoverageLost);
                expectedFaces.insert(tag);
                sourceFacesByBody[body].insert(tag);
            }
            children.Clear();
            TopExp::MapShapes(i.sourceSolids(body), TopAbs_EDGE, children);
            for (int n = 1; n <= children.Extent(); ++n) {
                tick();
                const int tag = i.sourceEdges.FindIndex(children(n));
                need(tag > 0, Reason::SourceCoverageLost);
                if (BRep_Tool::Degenerated(TopoDS::Edge(children(n))))
                    ++r.degeneratedSourceEdges;
                else
                    expectedEdges.insert(tag);
            }
            children.Clear();
            TopExp::MapShapes(i.sourceSolids(body), TopAbs_VERTEX, children);
            for (int n = 1; n <= children.Extent(); ++n) {
                tick();
                const int tag = i.sourceVertices.FindIndex(children(n));
                need(tag > 0, Reason::SourceCoverageLost);
                expectedVertices.insert(tag);
            }
        }
        r.sourceFaces = expectedFaces.size();
        r.sourceEdges = expectedEdges.size();
        r.sourceVertices = expectedVertices.size();
        for (int tag = 1; tag <= i.derivedFaces.Extent(); ++tag) {
            tick();
            const auto &owners = i.faceOwners[tag];
            need(owners.size() == 1 || owners.size() == 2, Reason::SourceCoverageLost);
            for (const auto owner : owners) {
                tick();
                need(owner.first > 0 && std::size_t(owner.first) < i.parentBodyByPiece.size() &&
                         selectedBodies.contains(i.parentBodyByPiece[owner.first]) &&
                         (owner.second == 1 || owner.second == -1),
                     Reason::SourceCoverageLost);
            }
            if (owners.size() == 2) {
                const auto a = owners[0], b = owners[1];
                need(a.first > 0 && b.first > 0 && a.first != b.first &&
                         std::size_t(a.first) < i.parentBodyByPiece.size() &&
                         std::size_t(b.first) < i.parentBodyByPiece.size() && a.second != 0 &&
                         a.second == -b.second &&
                         i.parentBodyByPiece[a.first] == i.parentBodyByPiece[b.first],
                     Reason::SourceCoverageLost);
                ++r.sharedCutFaces;
                continue;
            }
            ++r.exteriorFaces;
            const auto parent = i.parentByFace[tag];
            const int sourceBody = i.parentBodyByPiece[owners[0].first];
            need(parent.dim == 2 && sourceFacesByBody.at(sourceBody).contains(parent.tag),
                 Reason::SourceCoverageLost);
            observedFaces.insert(parent.tag);
            const auto derived = TopoDS::Face(i.derivedFaces(tag));
            const auto copiedShape = i.originalCopy.ModifiedShape(i.sourceFaces(parent.tag));
            need(!copiedShape.IsNull() && copiedShape.ShapeType() == TopAbs_FACE,
                 Reason::SupportMismatch);
            TopLoc_Location a, b;
            need(BRep_Tool::Surface(derived, a) ==
                         BRep_Tool::Surface(TopoDS::Face(copiedShape), b) &&
                     a == b,
                 Reason::SupportMismatch);
            OcctWitnessShapeMap children;
            TopExp::MapShapes(derived, TopAbs_VERTEX, children);
            for (int n = 1; n <= children.Extent(); ++n) {
                tick();
                const auto id = i.derivedVertices.FindIndex(children(n));
                need(id > 0, Reason::SourceCoverageLost);
                boundaryVertices.insert(id);
            }
            children.Clear();
            TopExp::MapShapes(derived, TopAbs_EDGE, children);
            for (int n = 1; n <= children.Extent(); ++n) {
                tick();
                const auto id = i.derivedEdges.FindIndex(children(n));
                need(id > 0, Reason::SourceCoverageLost);
                boundaryEdges.insert(id);
            }
        }
        r.supportsUnchanged = true;
        r.originalFacesRetained = observedFaces == expectedFaces;
        for (int tag : boundaryEdges) {
            const auto p = i.parentByEdge[tag];
            if (p.dim == 1)
                observedEdges.insert(p.tag);
        }
        for (int tag : boundaryVertices) {
            const auto p = i.parentByVertex[tag];
            if (p.dim == 0)
                observedVertices.insert(p.tag);
        }
        r.originalEdgesRetained = std::includes(observedEdges.begin(), observedEdges.end(),
                                                expectedEdges.begin(), expectedEdges.end());
        r.originalVerticesRetained =
            std::includes(observedVertices.begin(), observedVertices.end(),
                          expectedVertices.begin(), expectedVertices.end());
        need(r.originalFacesRetained && r.originalEdgesRetained && r.originalVerticesRetained,
             Reason::SourceCoverageLost);
        for (int tag : expectedFaces) {
            tick();
            const auto face = TopoDS::Face(i.sourceFaces(tag));
            double tol = BRep_Tool::Tolerance(face);
            OcctWitnessShapeMap children;
            TopExp::MapShapes(face, TopAbs_EDGE, children);
            for (int n = 1; n <= children.Extent(); ++n) {
                tick();
                tol = std::max(tol, BRep_Tool::Tolerance(TopoDS::Edge(children(n))));
            }
            children.Clear();
            TopExp::MapShapes(face, TopAbs_VERTEX, children);
            for (int n = 1; n <= children.Extent(); ++n) {
                tick();
                tol = std::max(tol, BRep_Tool::Tolerance(TopoDS::Vertex(children(n))));
            }
            GProp_GProps area;
            BRepGProp::SurfaceProperties(face, area, 1e-8, false);
            tick(0);
            need(std::isfinite(area.Mass()) && area.Mass() > 0 && std::isfinite(tol) && tol > 0,
                 Reason::NativeQueryFailed);
            r.sourceToleranceVolumeEnvelope += area.Mass() * std::min(tol, o.maximumPhysicalDrift);
        }
        auto point = [&](const gp_Pnt &p, mesh::Classification parent) {
            tick();
            ++r.samples;
            double error = 0, tol = 0;
            if (parent.dim == 0) {
                need(expectedVertices.contains(parent.tag), Reason::SourceCoverageLost);
                const auto vertex = TopoDS::Vertex(i.sourceVertices(parent.tag));
                error = p.Distance(BRep_Tool::Pnt(vertex));
                tol = BRep_Tool::Tolerance(vertex);
                r.maximumOriginalVertexShift = std::max(r.maximumOriginalVertexShift, error);
            } else {
                need((parent.dim == 1 && expectedEdges.contains(parent.tag)) ||
                         (parent.dim == 2 && expectedFaces.contains(parent.tag)),
                     Reason::SourceCoverageLost);
                const auto &target =
                    parent.dim == 1 ? i.sourceEdges(parent.tag) : i.sourceFaces(parent.tag);
                BRepExtrema_DistShapeShape distance(BRepBuilderAPI_MakeVertex(p).Vertex(), target);
                tick(0);
                need(distance.IsDone() && distance.NbSolution() > 0, Reason::NativeQueryFailed);
                error = distance.Value();
                tol = parent.dim == 1 ? BRep_Tool::Tolerance(TopoDS::Edge(target))
                                      : BRep_Tool::Tolerance(TopoDS::Face(target));
                // Returned supports are authoritative native nearest solutions,
                // including trimming edges/vertices of the owning source face.
                for (int n = 1; n <= distance.NbSolution(); ++n) {
                    tick();
                    const auto support = distance.SupportOnShape2(n);
                    if (support.ShapeType() == TopAbs_EDGE)
                        tol = std::max(tol, BRep_Tool::Tolerance(TopoDS::Edge(support)));
                    if (support.ShapeType() == TopAbs_VERTEX)
                        tol = std::max(tol, BRep_Tool::Tolerance(TopoDS::Vertex(support)));
                }
            }
            need(std::isfinite(error) && error >= 0 && std::isfinite(tol) && tol > 0,
                 Reason::NativeQueryFailed);
            tol = std::min(tol, o.maximumPhysicalDrift);
            r.maximumResidual = std::max(r.maximumResidual, error);
            r.maximumToleranceRatio = std::max(r.maximumToleranceRatio, error / tol);
            if (error > tol) {
                ++r.outsideOriginalTolerance;
                if (r.firstFailedParentTag < 0) {
                    r.firstFailedParentDimension = parent.dim;
                    r.firstFailedParentTag = parent.tag;
                }
            }
        };
        for (int tag : boundaryVertices)
            point(BRep_Tool::Pnt(TopoDS::Vertex(i.derivedVertices(tag))), i.parentByVertex[tag]);
        for (int tag : boundaryEdges) {
            tick();
            const auto edge = TopoDS::Edge(i.derivedEdges(tag));
            if (BRep_Tool::Degenerated(edge))
                continue;
            BRepAdaptor_Curve curve(edge);
            const auto lo = curve.FirstParameter(), hi = curve.LastParameter();
            need(std::isfinite(lo) && std::isfinite(hi) && hi >= lo, Reason::NativeQueryFailed);
            for (int n = 0; n <= 4; ++n)
                point(curve.Value(lo + (hi - lo) * n / 4), i.parentByEdge[tag]);
        }
        r.nativeVolumeDrift = std::abs(i.derivedVolume - i.sourceVolume);
        need(std::isfinite(r.sourceToleranceVolumeEnvelope) && r.sourceToleranceVolumeEnvelope > 0,
             Reason::NativeQueryFailed);
        r.sourceUnmodified = before == snapshot(i.originalShape, snapshotLimit);
        need(r.sourceUnmodified, Reason::SourceMutation);
        r.completed = true;
        need(r.outsideOriginalTolerance == 0, Reason::OriginalToleranceExceeded);
        need(r.nativeVolumeDrift <= r.sourceToleranceVolumeEnvelope, Reason::DriftEnvelopeExceeded);
        tick(0);
        r.admitted = true;
        r.reason = Reason::None;
    } catch (const Stop &s) {
        r.reason = s.reason;
    } catch (const Standard_OutOfMemory &) {
        r.reason = Reason::AllocationFailed;
    } catch (const std::bad_alloc &) {
        r.reason = Reason::AllocationFailed;
    } catch (...) {
        r.reason = Reason::NativeQueryFailed;
    }
    return r;
}
} // namespace paacs::geo::detail
