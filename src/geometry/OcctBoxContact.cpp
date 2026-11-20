// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctBoxContact.hpp"
#include "OcctModelInternal.hpp"
#include "paacs/predicates/Predicates.hpp"
#include <BRepAlgoAPI_Common.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <Message_ProgressIndicator.hxx>
#include <NCollection_List.hxx>
#include <Standard_OutOfMemory.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Compound.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <ostream>
#include <set>
#include <streambuf>

namespace paacs::geo::detail {
namespace {
using Reason = OcctBoxContactReason;
struct Stop { Reason reason; };
void need(bool value, Reason reason) { if (!value) throw Stop{reason}; }
bool finite(Vec3 p) { return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z); }
double coord(Vec3 p, int axis) { return axis == 0 ? p.x : axis == 1 ? p.y : p.z; }
struct Context {
    const OcctBoxContactOptions& options;
    std::uint64_t& work;
    void tick(std::uint64_t count = 1) {
        need(!options.cancellation || !options.cancellation->load(std::memory_order_relaxed), Reason::Cancelled);
        need(count <= options.maximumWork - work, Reason::BudgetExceeded);
        work += count;
    }
};
class Progress final : public Message_ProgressIndicator {
public:
    explicit Progress(const std::atomic<bool>* flag) : flag_(flag) {}
protected:
    bool UserBreak() override { return flag_ && flag_->load(std::memory_order_relaxed); }
    void Show(const Message_ProgressScope&, bool) override {}
private:
    const std::atomic<bool>* flag_;
};
class Buffer final : public std::streambuf {
public:
    explicit Buffer(std::uint64_t maximum) : maximum_(maximum) {}
    std::string value;
protected:
    std::streamsize xsputn(const char* data, std::streamsize count) override {
        need(count >= 0 && static_cast<std::uint64_t>(count) <= maximum_ - value.size(), Reason::BudgetExceeded);
        value.append(data, static_cast<std::size_t>(count));
        return count;
    }
    int_type overflow(int_type c) override {
        if (traits_type::eq_int_type(c, traits_type::eof())) return traits_type::not_eof(c);
        const char b = traits_type::to_char_type(c);
        xsputn(&b, 1);
        return c;
    }
private:
    std::uint64_t maximum_;
};
std::string snapshot(const OcctModel& model, std::uint64_t maximum) {
    Buffer buffer(maximum);
    std::ostream stream(&buffer);
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    BRepTools::Write(OcctModelAccess::shape(model), stream);
    return std::move(buffer.value);
}
void validate(const OcctBoxContactOptions& options, Context& context) {
    context.tick();
    need(options.maximumSnapshotBytes > 0, Reason::InvalidRequest);
    // Two snapshots plus string growth and bounded topology/staging storage.
    need(options.maximumSnapshotBytes <= options.maximumWorkingBytes / 8, Reason::BudgetExceeded);
}
TopoDS_Shape copy(const TopoDS_Shape& source, Context& context) {
    context.tick();
    BRepBuilderAPI_Copy operation;
    operation.Perform(source, true, false);
    context.tick();
    need(operation.IsDone(), Reason::NativeBooleanFailed);
    return operation.Shape();
}
TopoDS_Shape common(const TopoDS_Shape& first, const TopoDS_Shape& second, Context& context) {
    context.tick();
    BRepAlgoAPI_Common operation;
    NCollection_List<TopoDS_Shape> arguments, tools;
    arguments.Append(copy(first, context));
    tools.Append(copy(second, context));
    operation.SetArguments(arguments);
    operation.SetTools(tools);
    operation.SetNonDestructive(true);
    operation.SetRunParallel(false);
    occ::handle<Progress> progress = new Progress(context.options.cancellation);
    operation.Build(progress->Start());
    context.tick();
    need(operation.IsDone() && !operation.HasErrors(), Reason::NativeBooleanFailed);
    return operation.Shape();
}
double volume(const TopoDS_Shape& shape) {
    GProp_GProps properties;
    BRepGProp::VolumeProperties(shape, properties);
    return std::abs(properties.Mass());
}
double area(const TopoDS_Shape& shape) {
    GProp_GProps properties;
    BRepGProp::SurfaceProperties(shape, properties);
    return std::abs(properties.Mass());
}
bool has(const TopoDS_Shape& shape, TopAbs_ShapeEnum kind) { return TopExp_Explorer(shape, kind).More(); }
struct Box {
    int tag = -1;
    Vec3 lower{1e300, 1e300, 1e300}, upper{-1e300, -1e300, -1e300};
    double cadVolume = 0;
    std::array<const GVertex*, 8> corners{};
    std::set<int> faceTags;
    TopoDS_Shape native;
};
Box certify(const OcctModel& model, int tag, Context& context) {
    const auto* solid = model.solid(tag);
    need(solid != nullptr, Reason::InvalidRequest);
    Box result;
    result.tag = tag;
    result.native = OcctModelAccess::solid(model, tag);
    result.cadVolume = model.volume(tag);
    need(std::isfinite(result.cadVolume) && result.cadVolume > 0, Reason::UnsupportedBody);
    std::set<int> edges, vertices;
    const auto faces = solid->boundaryFaceTags();
    need(faces.size() == 6, Reason::UnsupportedBody);
    for (int faceTag : faces) {
        context.tick();
        need(result.faceTags.insert(faceTag).second, Reason::UnsupportedBody);
        const auto* face = model.face(faceTag);
        need(face && face->kind() == SurfaceKind::Plane, Reason::UnsupportedBody);
        const auto loops = model.boundaryUses(faceTag);
        need(loops.size() == 1 && loops.front().outer && loops.front().edges.size() == 4, Reason::UnsupportedBody);
        for (const auto& use : loops.front().edges) {
            context.tick();
            const auto* edge = model.edge(use.edgeTag);
            need(edge && edge->kind() == CurveKind::Line && !use.seam && !use.degenerated, Reason::UnsupportedBody);
            edges.insert(use.edgeTag);
            vertices.insert(use.startVertexTag);
            vertices.insert(use.endVertexTag);
        }
    }
    need(edges.size() == 12 && vertices.size() == 8, Reason::UnsupportedBody);
    for (int vertexTag : vertices) {
        context.tick();
        const auto* vertex = model.vertex(vertexTag);
        need(vertex && finite(vertex->position()), Reason::UnsupportedBody);
        const auto p = vertex->position();
        result.lower = {std::min(result.lower.x, p.x), std::min(result.lower.y, p.y), std::min(result.lower.z, p.z)};
        result.upper = {std::max(result.upper.x, p.x), std::max(result.upper.y, p.y), std::max(result.upper.z, p.z)};
    }
    const auto span = result.upper - result.lower;
    need(finite(span) && span.x > 0 && span.y > 0 && span.z > 0, Reason::UnsupportedBody);
    for (int vertexTag : vertices) {
        context.tick();
        const auto* vertex = model.vertex(vertexTag);
        const auto p = vertex->position();
        unsigned mask = 0;
        for (int axis = 0; axis < 3; ++axis) {
            // No coordinate tolerance: actual native vertices must be exactly
            // the eight corners of the common world-aligned box frame.
            need(coord(p, axis) == coord(result.lower, axis) || coord(p, axis) == coord(result.upper, axis), Reason::UnsupportedBody);
            if (coord(p, axis) == coord(result.upper, axis)) mask |= 1U << axis;
        }
        need(result.corners[mask] == nullptr, Reason::UnsupportedBody);
        result.corners[mask] = vertex;
    }
    std::set<int> sides;
    for (int faceTag : result.faceTags) {
        context.tick();
        const auto* descriptor = model.surfaceDescriptor(faceTag);
        need(descriptor && finite(descriptor->axis) && finite(descriptor->origin), Reason::UnsupportedBody);
        int normalAxis = -1;
        // Native 90-degree directions contain sin/cos roundoff (e.g. 6e-17).
        // This recognition guard never changes corner/contact coordinates;
        // complete expected-box Common and native volume prove the body below.
        for (int axis = 0; axis < 3; ++axis)
            if (std::abs(std::abs(coord(descriptor->axis, axis)) - 1) <= 1e-14) normalAxis = axis;
        need(normalAxis >= 0, Reason::UnsupportedBody);
        for (int axis = 0; axis < 3; ++axis)
            if (axis != normalAxis) need(std::abs(coord(descriptor->axis, axis)) <= 1e-14, Reason::UnsupportedBody);
        const auto at = coord(descriptor->origin, normalAxis);
        const bool high = std::abs(at - coord(result.upper, normalAxis)) < std::abs(at - coord(result.lower, normalAxis));
        need(std::abs(at - coord(high ? result.upper : result.lower, normalAxis)) <= model.tolerance(2, faceTag), Reason::UnsupportedBody);
        need(sides.insert(2 * normalAxis + high).second, Reason::UnsupportedBody);
    }
    const double expected = span.x * span.y * span.z;
    need(std::isfinite(expected) && expected > 0 && std::abs(expected - result.cadVolume) <= expected * 1e-10, Reason::UnsupportedBody);
    context.tick();
    const auto ideal = BRepPrimAPI_MakeBox(gp_Pnt(result.lower.x, result.lower.y, result.lower.z), span.x, span.y, span.z).Shape();
    const double commonVolume = volume(common(result.native, ideal, context));
    need(std::isfinite(commonVolume) && std::abs(commonVolume - expected) <= expected * 1e-10, Reason::UnsupportedBody);
    return result;
}
TopoDS_Shape skin(const TopoDS_Shape& source, Context& context) {
    const auto owned = copy(source, context);
    BRep_Builder builder;
    TopoDS_Compound result;
    builder.MakeCompound(result);
    for (TopExp_Explorer face(owned, TopAbs_FACE); face.More(); face.Next()) {
        context.tick();
        builder.Add(result, face.Current());
    }
    return result;
}
double fullFaceArea(const Box& first, const Box& second) {
    int normalAxis = -1;
    for (int axis = 0; axis < 3; ++axis) {
        if (coord(first.upper, axis) == coord(second.lower, axis) || coord(second.upper, axis) == coord(first.lower, axis)) {
            if (normalAxis >= 0) return 0;
            normalAxis = axis;
        }
    }
    if (normalAxis < 0) return 0;
    double result = 1;
    for (int axis = 0; axis < 3; ++axis) if (axis != normalAxis) {
        if (coord(first.lower, axis) != coord(second.lower, axis) || coord(first.upper, axis) != coord(second.upper, axis)) return 0;
        result *= coord(first.upper, axis) - coord(first.lower, axis);
    }
    return std::isfinite(result) && result > 0 ? result : 0;
}
bool unchanged(const OcctModel& model, const std::string& before, const OcctBoxContactOptions& options) noexcept {
    if (before.empty()) return false;
    try { return snapshot(model, options.maximumSnapshotBytes) == before; }
    catch (...) { return false; }
}
} // namespace

struct OcctBoxCertificateSession::Impl {
    struct Entry {
        Box box;
        std::uint64_t work = 0;
    };
    static constexpr std::size_t maximumEntries = 256;
    static constexpr std::uint64_t reservedEntryBytes = maximumEntries * 1024;
    Impl(const OcctModel& model, const std::string& bytes)
        : model(&model), identity(model.identity()), source(bytes) {}

    const OcctModel* model;
    std::uint64_t identity;
    std::string source;
    std::map<int,Entry> entries;
};

OcctBoxCertificateSession::OcctBoxCertificateSession() noexcept = default;
OcctBoxCertificateSession::~OcctBoxCertificateSession() = default;

struct OcctBoxCertificateAccess {
    using Cache = OcctBoxCertificateSession::Impl;

    static Cache* bind(const OcctModel& model, const std::string& before,
                       const OcctBoxContactOptions& options) {
        if (!options.certificateSession) return nullptr;
        auto& cache = options.certificateSession->impl_;
        // Reserve a bounded share of the existing snapshot/staging allowance.
        const auto allowance = options.maximumWorkingBytes / 8;
        if (allowance < Cache::reservedEntryBytes ||
            before.capacity() > allowance - Cache::reservedEntryBytes) {
            cache.reset();
            return nullptr;
        }
        if (cache && (cache->model != &model || cache->identity != model.identity() ||
                      cache->source != before))
            cache.reset();
        if (!cache) {
            try { cache = std::make_unique<Cache>(model, before); }
            catch (const std::bad_alloc&) { return nullptr; }
        }
        return cache.get();
    }

    static Box certifiedBox(Cache* cache, const OcctModel& model, int tag, Context& context) {
        if (!cache) return certify(model, tag, context);
        const auto found = cache->entries.find(tag);
        if (found != cache->entries.end() &&
            found->second.work <= context.options.maximumWork - context.work) {
            context.tick(found->second.work);
            return found->second.box;
        }
        const auto start = context.work;
        auto box = certify(model, tag, context);
        if (cache->entries.size() < Cache::maximumEntries) {
            try { cache->entries.try_emplace(tag, Cache::Entry{box, context.work - start}); }
            catch (const std::bad_alloc&) { /* Retain the certificate for this call only. */ }
        }
        return box;
    }

    static void discard(const OcctBoxContactOptions& options) noexcept {
        if (options.certificateSession) options.certificateSession->impl_.reset();
    }
};

OcctBoxFaceContactOutcome admitOcctBoxFaceContact(const OcctModel& model, int firstBodyTag,
    int secondBodyTag, const OcctBoxContactOptions& options) noexcept {
    OcctBoxFaceContactOutcome result;
    result.firstBodyTag = firstBodyTag;
    result.secondBodyTag = secondBodyTag;
    Context context{options, result.workUsed};
    std::string before;
    try {
        validate(options, context);
        need(firstBodyTag != secondBodyTag, Reason::InvalidRequest);
        context.tick(); before = snapshot(model, options.maximumSnapshotBytes);
        auto* certificates = OcctBoxCertificateAccess::bind(model, before, options);
        const auto first = OcctBoxCertificateAccess::certifiedBox(certificates, model, firstBodyTag, context);
        const auto second = OcctBoxCertificateAccess::certifiedBox(certificates, model, secondBodyTag, context);
        const auto material = common(first.native, second.native, context);
        need(!has(material, TopAbs_SOLID), Reason::MaterialOverlap);
        result.expectedFullFaceArea = fullFaceArea(first, second);
        need(result.expectedFullFaceArea > 0, Reason::PartialFaceContact);
        const auto nativeSkin = common(skin(first.native, context), skin(second.native, context), context);
        need(has(nativeSkin, TopAbs_FACE), Reason::PartialFaceContact);
        context.tick(); result.nativeSkinArea = area(nativeSkin);
        need(std::isfinite(result.nativeSkinArea) && result.nativeSkinArea > 0 &&
            std::abs(result.nativeSkinArea - result.expectedFullFaceArea) <= 1e-10 * std::max(1.0, result.expectedFullFaceArea), Reason::GeometryError);
        context.tick(); result.sourceImmutabilityVerified = unchanged(model, before, options);
        need(result.sourceImmutabilityVerified, Reason::SourceChanged);
        context.tick(); result.admitted = true; result.reason = Reason::None;
    } catch (const Stop& stopped) { result.reason = stopped.reason; }
    catch (const Standard_OutOfMemory&) { result.reason = Reason::AllocationFailed; }
    catch (const std::bad_alloc&) { result.reason = Reason::AllocationFailed; }
    catch (...) { result.reason = Reason::InternalError; }
    if (!result.admitted && !before.empty()) {
        result.sourceImmutabilityVerified = unchanged(model, before, options);
        if (result.reason == Reason::SourceChanged || !result.sourceImmutabilityVerified)
            OcctBoxCertificateAccess::discard(options);
    }
    return result;
}

OcctBoxCornerEvaluationOutcome evaluateOcctBoxContactCorners(const OcctModel& model,
    int bodyTag, const OcctBoxContactOptions& options, mesh::Mesh& destination) noexcept {
    OcctBoxCornerEvaluationOutcome result;
    Context context{options, result.workUsed};
    std::string before;
    try {
        validate(options, context);
        need(destination.nodeCount() == 8 && destination.cellCount() == 18 &&
            destination.countCells(mesh::CellType::Tet4) == 6 && destination.countCells(mesh::CellType::Tri3) == 12, Reason::UnsupportedResolution);
        context.tick(); before = snapshot(model, options.maximumSnapshotBytes);
        auto* certificates = OcctBoxCertificateAccess::bind(model, before, options);
        const auto box = OcctBoxCertificateAccess::certifiedBox(certificates, model, bodyTag, context);
        result.cadVolume = box.cadVolume;
        auto local = destination;
        std::set<int> evaluated;
        for (Index i = 0; i < local.nodeCount(); ++i) {
            context.tick();
            const auto cls = local.nodeClass(i);
            need(cls.dim == 0 && finite(local.node(i)), Reason::UnsupportedResolution);
            const GVertex* vertex = nullptr;
            for (const auto* candidate : box.corners) if (candidate->tag() == cls.tag) vertex = candidate;
            need(vertex && evaluated.insert(cls.tag).second, Reason::InvalidMesh);
            const double correction = (local.node(i) - vertex->position()).norm();
            need(std::isfinite(correction) && correction <= std::max(8 * model.tolerance(0, cls.tag), 1e-9), Reason::GeometryError);
            result.maximumEndpointCorrection = std::max(result.maximumEndpointCorrection, correction);
            local.setNode(i, vertex->position());
            ++result.evaluatedVertices;
            context.tick(); need(model.solid(bodyTag)->isInside(local.node(i)), Reason::GeometryError);
            ++result.nativeInclusionSamples;
        }
        std::set<int> coveredFaces;
        for (Index i = 0; i < local.cellCount(); ++i) {
            context.tick();
            const auto& cell = local.cell(i);
            const auto* ids = local.cellNodes(i);
            const int count = mesh::nodesPerCell(cell.type);
            need(count == 3 || count == 4, Reason::InvalidMesh);
            Vec3 center{};
            for (int j = 0; j < count; ++j) {
                need(ids[j] >= 0 && ids[j] < local.nodeCount(), Reason::InvalidMesh);
                center = center + local.node(ids[j]);
            }
            center = center / static_cast<double>(count);
            if (cell.type == mesh::CellType::Tet4) {
                need(cell.onEntity == mesh::Classification{3, bodyTag}, Reason::InvalidMesh);
                context.tick();
                need(predicates::orient3d(local.node(ids[0]), local.node(ids[1]), local.node(ids[2]), local.node(ids[3])) == predicates::Sign::Positive, Reason::InvalidMesh);
                context.tick(); need(model.solid(bodyTag)->isInside(center), Reason::GeometryError);
                ++result.nativeInclusionSamples;
            } else {
                need(cell.type == mesh::CellType::Tri3 && cell.onEntity.dim == 2 && box.faceTags.contains(cell.onEntity.tag), Reason::InvalidMesh);
                coveredFaces.insert(cell.onEntity.tag);
                const auto* face = model.face(cell.onEntity.tag);
                const double bound = std::max(8 * model.tolerance(2, cell.onEntity.tag), 1e-9);
                Vec2 uv; Vec3 point;
                context.tick(); need(face->projectPoint(center, uv, point), Reason::GeometryError);
                const double error = (point - center).norm();
                need(std::isfinite(error) && error <= bound, Reason::GeometryError);
                result.maximumBoundarySampleError = std::max(result.maximumBoundarySampleError, error);
                for (int j = 0; j < 3; ++j) {
                    context.tick(); need(face->projectPoint(local.node(ids[j]), uv, point), Reason::GeometryError);
                    const double nodeError = (point - local.node(ids[j])).norm();
                    need(std::isfinite(nodeError) && nodeError <= bound, Reason::GeometryError);
                    result.maximumBoundaryNodeError = std::max(result.maximumBoundaryNodeError, nodeError);
                }
            }
        }
        need(coveredFaces == box.faceTags, Reason::InvalidMesh);
        mesh::LinearVolumeAuditLimits limits;
        limits.maximumFaceUses = 36;
        limits.maximumTemporaryBytes = options.maximumWorkingBytes / 8;
        limits.cancellation = options.cancellation;
        context.tick(36); result.audit = mesh::auditLinearVolumeBoundary(local, limits);
        need(!result.audit.cancelled, Reason::Cancelled);
        need(!result.audit.budgetExceeded, Reason::BudgetExceeded);
        need(!result.audit.allocationFailed, Reason::AllocationFailed);
        need(result.audit.complete && result.audit.valid, Reason::InvalidMesh);
        context.tick(18); result.quality = algo::analyzeMixedMeshQuality(local);
        need(!result.quality.allocationFailed, Reason::AllocationFailed);
        need(result.quality.complete && result.quality.valid, Reason::InvalidMesh);
        result.meshVolume = result.quality.totalVolume;
        need(std::isfinite(result.meshVolume) && std::abs(result.meshVolume - result.cadVolume) <= result.cadVolume * 1e-10, Reason::GeometryError);
        context.tick(); result.sourceImmutabilityVerified = unchanged(model, before, options);
        need(result.sourceImmutabilityVerified, Reason::SourceChanged);
        context.tick(); destination.swap(local);
        result.published = true; result.destinationUnchanged = false; result.reason = Reason::None;
    } catch (const Stop& stopped) { result.reason = stopped.reason; }
    catch (const Standard_OutOfMemory&) { result.reason = Reason::AllocationFailed; }
    catch (const std::bad_alloc&) { result.reason = Reason::AllocationFailed; }
    catch (...) { result.reason = Reason::InternalError; }
    if (!result.published && !before.empty()) {
        result.sourceImmutabilityVerified = unchanged(model, before, options);
        if (result.reason == Reason::SourceChanged || !result.sourceImmutabilityVerified)
            OcctBoxCertificateAccess::discard(options);
    }
    return result;
}

const char* occtBoxContactReasonName(OcctBoxContactReason reason) noexcept {
    switch (reason) {
    case Reason::None: return "none";
    case Reason::InvalidRequest: return "invalid_request";
    case Reason::UnsupportedBody: return "unsupported_body";
    case Reason::PartialFaceContact: return "partial_face_contact";
    case Reason::MaterialOverlap: return "material_overlap";
    case Reason::NativeBooleanFailed: return "native_boolean_failed";
    case Reason::GeometryError: return "geometry_error";
    case Reason::UnsupportedResolution: return "unsupported_resolution";
    case Reason::InvalidMesh: return "invalid_mesh";
    case Reason::BudgetExceeded: return "budget_exceeded";
    case Reason::Cancelled: return "cancelled";
    case Reason::AllocationFailed: return "allocation_failed";
    case Reason::SourceChanged: return "source_changed";
    case Reason::InternalError: return "internal_error";
    }
    return "internal_error";
}
} // namespace paacs::geo::detail
