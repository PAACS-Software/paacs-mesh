// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/OcctBodyRelations.hpp"
#include "OcctContactAssembly.hpp"
#include "OcctCancellationProgress.hpp"
#include "OcctModelInternal.hpp"
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <NCollection_List.hxx>
#include <Precision.hxx>
#include <Standard_OutOfMemory.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Compound.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <string>
#include <tuple>

namespace paacs::geo {
namespace {
using detail::OcctModelAccess;
struct Stop { NativeCadMeshReason reason; const char* diagnostic; };
void need(bool value, NativeCadMeshReason reason, const char* diagnostic) {
    if (!value) throw Stop{reason, diagnostic};
}
std::uint64_t sum(std::uint64_t a, std::uint64_t b) {
    need(b <= std::numeric_limits<std::uint64_t>::max() - a,
         NativeCadMeshReason::BudgetExceeded, "body_relations_storage_overflow");
    return a + b;
}
std::uint64_t product(std::uint64_t a, std::uint64_t b) {
    need(a == 0 || b <= std::numeric_limits<std::uint64_t>::max() / a,
         NativeCadMeshReason::BudgetExceeded, "body_relations_storage_overflow");
    return a * b;
}
struct Budget {
    const OcctBodyRelationOptions& options;
    OcctBodyRelationReport& report;
    void tick(std::uint64_t work = 1) {
        need(!options.cancellation || !options.cancellation->load(std::memory_order_relaxed),
             NativeCadMeshReason::Cancelled, "body_relations_cancelled");
        need(work <= options.maximumWork - report.work,
             NativeCadMeshReason::BudgetExceeded, "body_relations_work_budget");
        report.work += work;
    }
    void storage(std::uint64_t bytes) {
        bytes = sum(bytes, options.cancellation ? detail::OcctCancellationProgress::retainedPayloadBytes() : 0);
        report.maximumLogicalWorkingBytes = std::max(report.maximumLogicalWorkingBytes, bytes);
        need(bytes <= options.maximumRetainedBytes,
             NativeCadMeshReason::BudgetExceeded, "body_relations_storage_budget");
        tick(0);
    }
};
bool has(const TopoDS_Shape& shape, TopAbs_ShapeEnum kind) {
    return TopExp_Explorer(shape, kind).More();
}
TopoDS_Shape common(const TopoDS_Shape& first, const TopoDS_Shape& second,
                    OcctBodyPairRelation& pair, Budget& budget) {
    budget.tick(); ++pair.nativeBooleanCalls;
    BRepAlgoAPI_Common query;
    NCollection_List<TopoDS_Shape> a, b;
    a.Append(first); b.Append(second);
    query.SetArguments(a); query.SetTools(b);
    query.SetNonDestructive(true); query.SetRunParallel(false);
    const auto progress = detail::makeOcctCancellationProgress(budget.options.cancellation);
    query.Build(Message_ProgressIndicator::Start(progress));
    budget.tick(0);
    if (!query.IsDone() || query.HasErrors()) return {};
    return query.Shape();
}
TopoDS_Shape section(const TopoDS_Shape& first, const TopoDS_Shape& second,
                     OcctBodyPairRelation& pair, Budget& budget) {
    budget.tick(); ++pair.nativeSectionCalls;
    BRepAlgoAPI_Section query;
    NCollection_List<TopoDS_Shape> a, b;
    a.Append(first); b.Append(second);
    query.SetArguments(a); query.SetTools(b);
    query.SetNonDestructive(true); query.SetRunParallel(false);
    query.Approximation(false);
    const auto progress = detail::makeOcctCancellationProgress(budget.options.cancellation);
    query.Build(Message_ProgressIndicator::Start(progress)); budget.tick(0);
    if (!query.IsDone() || query.HasErrors()) return {};
    return query.Shape();
}
TopoDS_Shape skin(const TopoDS_Shape& shape, Budget& budget) {
    TopoDS_Compound result;
    BRep_Builder builder; builder.MakeCompound(result);
    for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next()) {
        budget.tick(); builder.Add(result, it.Current());
    }
    return result;
}
bool separated(Vec3 a, Vec3 b, Vec3 c, Vec3 d, double tolerance) {
    // Expanding the combined separation margin outwards prevents a rounding
    // equality from becoming a broad-phase proof of separation.
    const auto beyond = [tolerance](double low, double high) {
        return low > std::nextafter(high + tolerance, std::numeric_limits<double>::infinity());
    };
    return beyond(a.x, d.x) || beyond(c.x, b.x) || beyond(a.y, d.y) ||
           beyond(c.y, b.y) || beyond(a.z, d.z) || beyond(c.z, b.z);
}
void inspectPair(const OcctModel& model, OcctBodyPairRelation& pair, Budget& budget) {
    budget.tick();
    const auto* a = model.solid(pair.firstSolidTag);
    const auto* b = model.solid(pair.secondSolidTag);
    const auto ta = model.tolerance(3, pair.firstSolidTag);
    const auto tb = model.tolerance(3, pair.secondSolidTag);
    need(a && b && std::isfinite(ta) && ta >= 0 && std::isfinite(tb) && tb >= 0,
         NativeCadMeshReason::InvalidRequest, "body_relations_invalid_body");
    pair.classificationTolerance = std::max({Precision::Confusion(), ta + tb, budget.options.distanceTolerance});
    need(std::isfinite(pair.classificationTolerance), NativeCadMeshReason::InvalidRequest,
         "body_relations_invalid_tolerance");
    Vec3 al, ah, bl, bh; a->bbox(al, ah); b->bbox(bl, bh);
    if (separated(al, ah, bl, bh, pair.classificationTolerance)) {
        pair.relation = OcctBodyRelation::Disconnected;
        pair.evidence = OcctBodyRelationEvidence::SeparatedBounds;
        pair.issue = OcctBodyRelationIssue::None;
        return;
    }
    pair.issue = OcctBodyRelationIssue::CopyFailed;
    budget.tick(); ++pair.nativeCopies;
    BRepBuilderAPI_Copy ca(OcctModelAccess::solid(model, pair.firstSolidTag), true, false);
    budget.tick(0);
    budget.tick(); ++pair.nativeCopies;
    BRepBuilderAPI_Copy cb(OcctModelAccess::solid(model, pair.secondSolidTag), true, false);
    budget.tick(0);
    if (!ca.IsDone() || !cb.IsDone() || ca.Shape().IsNull() || cb.Shape().IsNull()) return;
    pair.issue = OcctBodyRelationIssue::MaterialQueryFailed;
    const auto material = common(ca.Shape(), cb.Shape(), pair, budget);
    if (material.IsNull()) return;
    pair.issue = OcctBodyRelationIssue::InvalidMaterialResult;
    budget.tick(); ++pair.nativeShapeValidations;
    const bool valid = BRepCheck_Analyzer(material).IsValid(); budget.tick(0);
    if (!valid) return;
    pair.materialQueryComplete = true;
    if (has(material, TopAbs_SOLID)) {
        pair.issue = OcctBodyRelationIssue::InvalidOverlapVolume;
        GProp_GProps coarse, fine;
        budget.tick(2); ++pair.nativeVolumeIntegrations;
        const auto e1 = BRepGProp::VolumePropertiesGK(material, coarse, 1e-8, false, true, false, false, false);
        budget.tick(0);
        budget.tick(2); ++pair.nativeVolumeIntegrations;
        const auto e2 = BRepGProp::VolumePropertiesGK(material, fine, 1e-10, false, true, false, false, false);
        budget.tick(0);
        if (!(std::isfinite(e1) && e1 >= 0 && std::isfinite(e2) && e2 >= 0 &&
              std::isfinite(coarse.Mass()) && coarse.Mass() > 0 &&
              std::isfinite(fine.Mass()) && fine.Mass() > 0 &&
              std::abs(coarse.Mass() - fine.Mass()) <= fine.Mass() * 1e-6)) return;
        pair.nativeOverlapVolume = fine.Mass();
        pair.relation = OcctBodyRelation::MaterialOverlap;
        pair.evidence = OcctBodyRelationEvidence::NativeMaterialCommon;
        pair.issue = OcctBodyRelationIssue::None;
        return;
    }
    // A successful solid Common with no solids does not distinguish a gap
    // from touching. Distance and boundary topology are separate evidence.
    pair.issue = OcctBodyRelationIssue::DistanceQueryFailed;
    budget.tick(); ++pair.nativeDistanceCalls;
    BRepExtrema_DistShapeShape distance(ca.Shape(), cb.Shape()); budget.tick(0);
    if (!distance.IsDone() || !std::isfinite(distance.Value()) || distance.Value() < 0) return;
    pair.minimumDistance = distance.Value();
    if (*pair.minimumDistance > pair.classificationTolerance) {
        pair.relation = OcctBodyRelation::Disconnected;
        pair.evidence = OcctBodyRelationEvidence::NativeDistance;
        pair.issue = OcctBodyRelationIssue::None;
        return;
    }
    pair.issue = OcctBodyRelationIssue::BoundaryQueryFailed;
    const auto firstSkin = skin(ca.Shape(), budget), secondSkin = skin(cb.Shape(), budget);
    auto boundary = common(firstSkin, secondSkin, pair, budget);
    if (boundary.IsNull()) return;
    pair.issue = OcctBodyRelationIssue::InvalidBoundaryResult;
    budget.tick(); ++pair.nativeShapeValidations;
    const bool boundaryValid = BRepCheck_Analyzer(boundary).IsValid(); budget.tick(0);
    if (!boundaryValid) return;
    pair.boundaryQueryComplete = true;
    pair.evidence = OcctBodyRelationEvidence::NativeBoundaryCommon;
    // Face Common deliberately removes lower-dimensional intersections of
    // faces. A near distance alone cannot replace that missing evidence.
    // Section supplies native intersection edges/vertices (including common
    // blocks) without joining solids or changing either copied source shape.
    if (!has(boundary, TopAbs_FACE) && !has(boundary, TopAbs_EDGE) && !has(boundary, TopAbs_VERTEX)) {
        pair.issue = OcctBodyRelationIssue::SectionQueryFailed;
        boundary = section(firstSkin, secondSkin, pair, budget);
        if (boundary.IsNull()) return;
        pair.issue = OcctBodyRelationIssue::InvalidSectionResult;
        budget.tick(); ++pair.nativeShapeValidations;
        const bool sectionValid = BRepCheck_Analyzer(boundary).IsValid(); budget.tick(0);
        if (!sectionValid || has(boundary, TopAbs_FACE) || has(boundary, TopAbs_SOLID)) return;
        pair.evidence = OcctBodyRelationEvidence::NativeBoundarySection;
    }
    pair.issue = OcctBodyRelationIssue::InvalidContactMeasure;
    if (has(boundary, TopAbs_FACE)) {
        GProp_GProps properties; budget.tick(); ++pair.nativeBoundaryMeasures;
        BRepGProp::SurfaceProperties(boundary, properties, true, false); budget.tick(0);
        if (!std::isfinite(properties.Mass()) || properties.Mass() <= 0) return;
        pair.nativeContactArea = properties.Mass(); pair.relation = OcctBodyRelation::FaceContact;
    } else if (has(boundary, TopAbs_EDGE)) {
        GProp_GProps properties; budget.tick(); ++pair.nativeBoundaryMeasures;
        BRepGProp::LinearProperties(boundary, properties, true, false); budget.tick(0);
        if (!std::isfinite(properties.Mass()) || properties.Mass() <= 0) return;
        pair.nativeContactLength = properties.Mass(); pair.relation = OcctBodyRelation::EdgeContact;
    } else if (has(boundary, TopAbs_VERTEX)) pair.relation = OcctBodyRelation::PointContact;
    else { pair.issue = OcctBodyRelationIssue::NearButUnproved; return; }
    pair.issue = OcctBodyRelationIssue::None;
}
} // namespace

OcctBodyRelationReport inspectOcctBodyPairs(const OcctModel& model,
    std::span<const OcctBodyPairRequest> selected, const OcctBodyRelationOptions& options) noexcept {
    OcctBodyRelationReport report;
    Budget budget{options, report};
    std::string before;
    bool finalSourceAttempted = false;
    std::uint64_t retained = 0;
    const auto sourceLimit = std::min(options.maximumSourceBytes, options.maximumRetainedBytes / 4);
    try {
        need(options.maximumPairs > 0 && options.maximumWork > 0 &&
             options.maximumRetainedBytes > 0 && options.maximumSourceBytes > 0 &&
             std::isfinite(options.distanceTolerance) && options.distanceTolerance >= 0,
             NativeCadMeshReason::InvalidRequest, "body_relations_invalid_limits");
        need(sourceLimit > 0, NativeCadMeshReason::BudgetExceeded, "body_relations_storage_budget");
        need(selected.size() <= options.maximumPairs, NativeCadMeshReason::BudgetExceeded,
             "body_relations_pair_budget");
        retained = product(selected.size(), sizeof(OcctBodyPairRequest) + sizeof(OcctBodyPairRelation));
        budget.storage(sum(retained, product(sourceLimit, 2)));
        budget.tick();
        std::vector<OcctBodyPairRequest> pairs;
        pairs.reserve(selected.size()); report.pairs.reserve(selected.size());
        for (auto pair : selected) {
            budget.tick();
            if (pair.firstSolidTag > pair.secondSolidTag) std::swap(pair.firstSolidTag, pair.secondSolidTag);
            need(pair.firstSolidTag > 0 && pair.firstSolidTag != pair.secondSolidTag &&
                 model.solid(pair.firstSolidTag) && model.solid(pair.secondSolidTag),
                 NativeCadMeshReason::InvalidRequest, "body_relations_invalid_pair");
            pairs.push_back(pair);
        }
        const auto less = [&](const auto& a, const auto& b) {
            budget.tick(); return std::tie(a.firstSolidTag, a.secondSolidTag) < std::tie(b.firstSolidTag, b.secondSolidTag);
        };
        std::sort(pairs.begin(), pairs.end(), less);
        for (std::size_t i = 1; i < pairs.size(); ++i) {
            budget.tick();
            need(pairs[i - 1].firstSolidTag != pairs[i].firstSolidTag ||
                 pairs[i - 1].secondSolidTag != pairs[i].secondSolidTag,
                 NativeCadMeshReason::InvalidRequest, "body_relations_duplicate_pair");
        }
        // Reserve both source audits before any geometry query. The second is
        // still performed during cancellation/budget cleanup and stays charged.
        budget.tick(2); ++report.sourceSnapshots;
        const auto start = detail::snapshotOcctContactSource(model, sourceLimit, before); budget.tick(0);
        need(start == NativeCadMeshReason::None, start, "body_relations_source_audit");
        for (const auto request : pairs) {
            report.pairs.emplace_back(); auto& pair = report.pairs.back();
            pair.firstSolidTag = request.firstSolidTag; pair.secondSolidTag = request.secondSolidTag;
            try { inspectPair(model, pair, budget); }
            catch (const Standard_OutOfMemory&) { throw; }
            catch (const Standard_Failure&) { /* Preserve stage-specific unresolved evidence. */ }
        }
        budget.tick(0);
        std::string after;
        finalSourceAttempted = true; ++report.sourceSnapshots;
        const auto end = detail::snapshotOcctContactSource(model, sourceLimit, after);
        report.sourceAuditAvailable = end == NativeCadMeshReason::None;
        report.sourceUnmodified = report.sourceAuditAvailable && before == after;
        need(end == NativeCadMeshReason::None, end, "body_relations_source_audit");
        need(report.sourceUnmodified, NativeCadMeshReason::InternalError, "body_relations_source_changed");
        budget.tick(0);
        report.complete = true; report.reason = NativeCadMeshReason::None;
        report.diagnostic = "body_relations_measured";
    } catch (const Stop& stop) { report.reason = stop.reason; report.diagnostic = stop.diagnostic; }
    catch (const Standard_OutOfMemory&) { report.reason = NativeCadMeshReason::AllocationFailed; report.diagnostic = "body_relations_allocation_failed"; }
    catch (const std::bad_alloc&) { report.reason = NativeCadMeshReason::AllocationFailed; report.diagnostic = "body_relations_allocation_failed"; }
    catch (...) { report.reason = NativeCadMeshReason::InternalError; report.diagnostic = "body_relations_internal_error"; }
    if (!finalSourceAttempted && !before.empty()) {
        // Diagnostic cancellation is not permission to omit the final source
        // guard. This read-only cleanup does not admit another geometry query.
        std::string after;
        ++report.sourceSnapshots;
        const auto end = detail::snapshotOcctContactSource(model, sourceLimit, after);
        report.sourceAuditAvailable = end == NativeCadMeshReason::None;
        report.sourceUnmodified = report.sourceAuditAvailable && before == after;
        if (report.sourceAuditAvailable && !report.sourceUnmodified) {
            report.reason = NativeCadMeshReason::InternalError;
            report.diagnostic = "body_relations_source_changed";
        }
    }
    return report;
}

const char* occtBodyRelationName(OcctBodyRelation value) noexcept {
    switch (value) {
    case OcctBodyRelation::Disconnected: return "disconnected";
    case OcctBodyRelation::FaceContact: return "face_contact";
    case OcctBodyRelation::EdgeContact: return "edge_contact";
    case OcctBodyRelation::PointContact: return "point_contact";
    case OcctBodyRelation::MaterialOverlap: return "material_overlap";
    default: return "unresolved";
    }
}
const char* occtBodyRelationEvidenceName(OcctBodyRelationEvidence value) noexcept {
    switch (value) {
    case OcctBodyRelationEvidence::SeparatedBounds: return "separated_bounds";
    case OcctBodyRelationEvidence::NativeDistance: return "native_distance";
    case OcctBodyRelationEvidence::NativeMaterialCommon: return "native_material_common";
    case OcctBodyRelationEvidence::NativeBoundaryCommon: return "native_boundary_common";
    case OcctBodyRelationEvidence::NativeBoundarySection: return "native_boundary_section";
    default: return "none";
    }
}
const char* occtBodyRelationIssueName(OcctBodyRelationIssue value) noexcept {
    switch (value) {
    case OcctBodyRelationIssue::None: return "none";
    case OcctBodyRelationIssue::NotStarted: return "not_started";
    case OcctBodyRelationIssue::CopyFailed: return "copy_failed";
    case OcctBodyRelationIssue::MaterialQueryFailed: return "material_query_failed";
    case OcctBodyRelationIssue::InvalidMaterialResult: return "invalid_material_result";
    case OcctBodyRelationIssue::InvalidOverlapVolume: return "invalid_overlap_volume";
    case OcctBodyRelationIssue::DistanceQueryFailed: return "distance_query_failed";
    case OcctBodyRelationIssue::BoundaryQueryFailed: return "boundary_query_failed";
    case OcctBodyRelationIssue::InvalidBoundaryResult: return "invalid_boundary_result";
    case OcctBodyRelationIssue::InvalidContactMeasure: return "invalid_contact_measure";
    case OcctBodyRelationIssue::NearButUnproved: return "near_but_unproved";
    case OcctBodyRelationIssue::SectionQueryFailed: return "section_query_failed";
    case OcctBodyRelationIssue::InvalidSectionResult: return "invalid_section_result";
    default: return "unknown";
    }
}
} // namespace paacs::geo
