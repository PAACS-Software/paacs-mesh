// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/OcctBodyRelations.hpp"
#include "OcctModelInternal.hpp"
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <NCollection_List.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <new>
#include <sstream>
#include <stdexcept>

namespace { std::atomic_bool failNext{false}, failEvery{false}; }
void* operator new(std::size_t n) {
    if (failEvery.load() || failNext.exchange(false)) throw std::bad_alloc{};
    if (void* p = std::malloc(n ? n : 1)) return p;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }

using namespace paacs;
using namespace paacs::geo;
using paacs::geo::detail::OcctModelAccess;
namespace {
void check(bool value, const char* message) { if (!value) throw std::runtime_error(message); }
std::string snapshot(const OcctModel& model) {
    std::ostringstream output; BRepTools::Write(OcctModelAccess::shape(model), output); return output.str();
}
TopoDS_Compound compound(std::initializer_list<TopoDS_Shape> shapes) {
    TopoDS_Compound result; BRep_Builder builder; builder.MakeCompound(result);
    for (const auto& shape : shapes) builder.Add(result, shape);
    return result;
}
NativeCadImportOutcome importOwned(TopoDS_Shape shape) {
    auto imported = OcctModelAccess::fromShape(shape);
    check(imported.imported(), "owned body-relation source imports");
    return imported;
}
OcctBodyRelationReport pairReport(const TopoDS_Shape& first, const TopoDS_Shape& second,
                                 OcctBodyRelationOptions options = {}) {
    auto imported = importOwned(compound({first, second}));
    check(imported.model->solids().size() == 2, "two placed bodies retained");
    const auto original = snapshot(*imported.model);
    const OcctBodyPairRequest selected[]{ {2, 1} };
    auto report = inspectOcctBodyPairs(*imported.model, selected, options);
    check(report.complete && report.reason == NativeCadMeshReason::None &&
          report.sourceAuditAvailable && report.sourceUnmodified && report.sourceSnapshots == 2,
          "completed selected-pair query guards original source");
    check(report.pairs.size() == 1 && report.pairs[0].firstSolidTag == 1 &&
          report.pairs[0].secondSolidTag == 2 && snapshot(*imported.model) == original,
          "canonical body tags and original BREP identity");
    check(report.work <= options.maximumWork && report.maximumLogicalWorkingBytes <= options.maximumRetainedBytes,
          "query receipts stay within admitted diagnostic bounds");
    return report;
}
void analyticRelations() {
    const auto box = BRepPrimAPI_MakeBox(1., 1., 1.).Shape();
    const auto shifted = [](double x, double y, double z) {
        return BRepPrimAPI_MakeBox(gp_Pnt(x, y, z), 1., 1., 1.).Shape();
    };
    auto report = pairReport(box, shifted(3, 0, 0));
    const auto& disconnected = report.pairs.front();
    check(disconnected.relation == OcctBodyRelation::Disconnected &&
          disconnected.evidence == OcctBodyRelationEvidence::SeparatedBounds &&
          disconnected.nativeCopies == 0 && disconnected.nativeBooleanCalls == 0 &&
          !disconnected.materialQueryComplete && !disconnected.minimumDistance,
          "strict separated bounds prove only separation without native queries");
    report = pairReport(box, shifted(1, 0, 0));
    const auto face = report.pairs.front();
    check(face.relation == OcctBodyRelation::FaceContact && face.materialQueryComplete &&
          face.boundaryQueryComplete && face.nativeContactArea &&
          std::abs(*face.nativeContactArea - 1.) < 1e-10 && !face.nativeOverlapVolume,
          "unit face contact is not disconnected or material overlap");
    report = pairReport(box, shifted(1, 1, 0));
    if (report.pairs[0].relation != OcctBodyRelation::EdgeContact ||
        !report.pairs[0].nativeContactLength ||
        std::abs(*report.pairs[0].nativeContactLength - 1.) >= 1e-10) {
        const auto& pair = report.pairs[0];
        std::cerr << "edge fixture: relation=" << occtBodyRelationName(pair.relation)
                  << " evidence=" << occtBodyRelationEvidenceName(pair.evidence)
                  << " issue=" << occtBodyRelationIssueName(pair.issue)
                  << " distance=" << pair.minimumDistance.value_or(-1)
                  << " area=" << pair.nativeContactArea.value_or(-1)
                  << " length=" << pair.nativeContactLength.value_or(-1) << '\n';
    }
    check(report.pairs[0].relation == OcctBodyRelation::EdgeContact &&
          report.pairs[0].evidence == OcctBodyRelationEvidence::NativeBoundarySection &&
          report.pairs[0].nativeSectionCalls == 1 &&
          report.pairs[0].nativeContactLength &&
          std::abs(*report.pairs[0].nativeContactLength - 1.) < 1e-10,
          "unit edge contact has native length evidence");
    report = pairReport(box, shifted(1, 1, 1));
    check(report.pairs[0].relation == OcctBodyRelation::PointContact &&
          (report.pairs[0].evidence == OcctBodyRelationEvidence::NativeBoundaryCommon ||
           report.pairs[0].evidence == OcctBodyRelationEvidence::NativeBoundarySection) &&
          report.pairs[0].minimumDistance && *report.pairs[0].minimumDistance == 0,
          "point contact requires native boundary topology");
    report = pairReport(box, shifted(.5, 0, 0));
    check(report.pairs[0].relation == OcctBodyRelation::MaterialOverlap &&
          report.pairs[0].nativeOverlapVolume &&
          std::abs(*report.pairs[0].nativeOverlapVolume - .5) < 1e-10 &&
          report.pairs[0].nativeVolumeIntegrations == 2,
          "positive native common volume reports interference without fusion");
    report = pairReport(box, BRepPrimAPI_MakeBox(gp_Pnt(.25, .25, .25), .5, .5, .5).Shape());
    check(report.pairs[0].relation == OcctBodyRelation::MaterialOverlap &&
          report.pairs[0].nativeOverlapVolume &&
          std::abs(*report.pairs[0].nativeOverlapVolume - .125) < 1e-10,
          "contained material remains overlap despite separated boundary surfaces");
    OcctBodyRelationOptions ambiguity; ambiguity.distanceTolerance = .01;
    report = pairReport(box, shifted(1.001, 0, 0), ambiguity);
    check(report.pairs[0].relation == OcctBodyRelation::Unresolved &&
          report.pairs[0].issue == OcctBodyRelationIssue::NearButUnproved &&
          report.pairs[0].minimumDistance &&
          std::abs(*report.pairs[0].minimumDistance - .001) < 1e-10 &&
          !report.pairs[0].nativeContactArea,
          "near bounding boxes and a small distance do not fabricate contact");
    const auto cylinder = BRepPrimAPI_MakeCylinder(1., 1.).Shape();
    const auto diagonal = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1.5, 1.5, 0), gp_Dir(0, 0, 1)), 1., 1.).Shape();
    report = pairReport(cylinder, diagonal);
    check(report.pairs[0].relation == OcctBodyRelation::Disconnected &&
          report.pairs[0].evidence == OcctBodyRelationEvidence::NativeDistance &&
          report.pairs[0].minimumDistance &&
          std::abs(*report.pairs[0].minimumDistance - (std::sqrt(4.5) - 2.)) < 1e-7,
          "overlapping bounds do not imply contact: native distance resolves the gap");
}
void placedAndSharedTopology() {
    const auto box = BRepPrimAPI_MakeBox(1., 1., 1.).Shape();
    gp_Trsf translate; translate.SetTranslation(gp_Vec(1, 0, 0));
    auto report = pairReport(box, box.Moved(TopLoc_Location(translate)));
    check(report.pairs[0].relation == OcctBodyRelation::FaceContact,
          "repeated underlying shape retains distinct placed body identities");
    gp_Trsf placement;
    placement.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), .37);
    placement.SetTranslationPart(gp_Vec(5, -2, 3));
    auto placed = importOwned(compound({box, box.Moved(TopLoc_Location(translate))}).Moved(TopLoc_Location(placement)));
    const auto original = snapshot(*placed.model);
    const OcctBodyPairRequest selected[]{ {1, 2} };
    report = inspectOcctBodyPairs(*placed.model, selected);
    check(report.complete && report.pairs[0].relation == OcctBodyRelation::FaceContact &&
          report.pairs[0].nativeContactArea && std::abs(*report.pairs[0].nativeContactArea - 1.) < 1e-9 &&
          report.sourceUnmodified && snapshot(*placed.model) == original,
          "rotated translated assembly retains contact area and native source");
    NCollection_List<TopoDS_Shape> arguments, tools;
    arguments.Append(BRepPrimAPI_MakeBox(2., 1., 1.).Shape());
    tools.Append(BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(1, 0, 0), gp_Dir(1, 0, 0)), -3., 3., -3., 3.).Shape());
    BRepAlgoAPI_Splitter split; split.SetArguments(arguments); split.SetTools(tools);
    split.SetNonDestructive(true); split.SetRunParallel(false); split.Build();
    check(split.IsDone() && !split.HasErrors(), "shared native face fixture splits");
    auto shared = importOwned(split.Shape());
    check(shared.model->solids().size() == 2, "shared source retains two solids");
    report = inspectOcctBodyPairs(*shared.model, selected);
    check(report.complete && report.pairs[0].relation == OcctBodyRelation::FaceContact &&
          report.pairs[0].nativeContactArea && std::abs(*report.pairs[0].nativeContactArea - 1.) < 1e-9,
          "shared native topology is contact evidence, not automatic node bonding");
}
void batchAndBounds() {
    auto imported = importOwned(compound({BRepPrimAPI_MakeBox(1., 1., 1.).Shape(),
        BRepPrimAPI_MakeBox(gp_Pnt(1, 0, 0), 1., 1., 1.).Shape(),
        BRepPrimAPI_MakeBox(gp_Pnt(4, 0, 0), 1., 1., 1.).Shape(),
        BRepPrimAPI_MakeBox(gp_Pnt(.5, 0, 0), 1., 1., 1.).Shape()}));
    const auto& model = *imported.model; const auto original = snapshot(model);
    const OcctBodyPairRequest selected[]{ {4, 1}, {2, 1}, {1, 3} };
    const auto report = inspectOcctBodyPairs(model, selected);
    check(report.complete && report.pairs.size() == 3 && report.sourceSnapshots == 2 &&
          report.pairs[0].secondSolidTag == 2 && report.pairs[0].relation == OcctBodyRelation::FaceContact &&
          report.pairs[1].secondSolidTag == 3 && report.pairs[1].relation == OcctBodyRelation::Disconnected &&
          report.pairs[2].secondSolidTag == 4 && report.pairs[2].relation == OcctBodyRelation::MaterialOverlap,
          "selected batch has canonical order, no unrequested all-pairs graph, and one source guard");
    auto options = OcctBodyRelationOptions{}; options.maximumPairs = 2;
    auto refused = inspectOcctBodyPairs(model, selected, options);
    check(!refused.complete && refused.reason == NativeCadMeshReason::BudgetExceeded &&
          refused.pairs.empty() && refused.sourceSnapshots == 0, "pair count refuses before geometry work");
    const OcctBodyPairRequest duplicate[]{ {1, 2}, {2, 1} }, missing[]{ {1, 99} }, self[]{ {1, 1} };
    for (const auto input : {std::span<const OcctBodyPairRequest>(duplicate),
                             std::span<const OcctBodyPairRequest>(missing),
                             std::span<const OcctBodyPairRequest>(self)}) {
        refused = inspectOcctBodyPairs(model, input);
        check(!refused.complete && refused.reason == NativeCadMeshReason::InvalidRequest &&
              refused.sourceSnapshots == 0, "invalid or duplicate pair refuses before native queries");
    }
    options = {}; options.maximumRetainedBytes = 128;
    refused = inspectOcctBodyPairs(model, selected, options);
    check(!refused.complete && refused.reason == NativeCadMeshReason::BudgetExceeded &&
          refused.sourceSnapshots == 0, "bounded diagnostic storage refuses before native queries");
    options = {}; options.maximumSourceBytes = 1;
    refused = inspectOcctBodyPairs(model, selected, options);
    check(!refused.complete && refused.reason == NativeCadMeshReason::BudgetExceeded,
          "bounded BREP guard refuses an oversized source");
    options = {}; options.maximumWork = report.work - 1;
    refused = inspectOcctBodyPairs(model, selected, options);
    check(!refused.complete && refused.reason == NativeCadMeshReason::BudgetExceeded &&
          refused.work <= options.maximumWork && refused.sourceAuditAvailable &&
          refused.sourceUnmodified && refused.sourceSnapshots == 2,
          "one-less work refuses without losing reserved final source audit");
    std::atomic_bool cancellation{true}; options = {}; options.cancellation = &cancellation;
    refused = inspectOcctBodyPairs(model, selected, options);
    check(!refused.complete && refused.reason == NativeCadMeshReason::Cancelled && refused.sourceSnapshots == 0,
          "pre-cancelled batch does no native geometry work");
    failNext.store(true); refused = inspectOcctBodyPairs(model, selected); failNext.store(false);
    check(!refused.complete && refused.reason == NativeCadMeshReason::AllocationFailed,
          "single allocation failure becomes an explicit bounded query refusal");
    failEvery.store(true); refused = inspectOcctBodyPairs(model, selected); failEvery.store(false);
    check(!refused.complete && refused.reason == NativeCadMeshReason::AllocationFailed,
          "persistent allocation failure does not escape the diagnostic API");
    check(snapshot(model) == original, "all successful and failed relation queries preserve exact source bytes");
    const auto empty = inspectOcctBodyPairs(model, {});
    check(empty.complete && empty.pairs.empty() && empty.sourceSnapshots == 2,
          "explicit empty selection does not expand to an implicit full graph");
}
} // namespace
int main() {
    try {
        analyticRelations(); placedAndSharedTopology(); batchAndBounds();
        std::cout << "Selected native body relations: independent measures, contact dimensions, placed identities, source guards and bounded refusal passed\n";
        return 0;
    } catch (const std::exception& error) {
        failNext.store(false); failEvery.store(false); std::cerr << error.what() << '\n'; return 1;
    } catch (...) {
        failNext.store(false); failEvery.store(false); std::cerr << "unexpected body-relation failure\n"; return 2;
    }
}
