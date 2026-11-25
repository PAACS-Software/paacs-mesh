// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctModelInternal.hpp"
#include "paacs/geometry/OcctMesher.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <gp_Pln.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <TopoDS_Compound.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Trsf.hxx>
#include <gp_Ax1.hxx>
#include <sstream>
#include <iostream>
#include <set>
#include <atomic>
#include <bit>
#include <cstdint>
#include <string_view>
using namespace paacs;using namespace paacs::geo;using namespace paacs::geo::detail;
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
std::string snapshot(const OcctModel& m){std::ostringstream s;BRepTools::Write(OcctModelAccess::shape(m),s);return s.str();}
namespace {
enum class CancelChildAt { Never, SecondBegin, SecondEnd, InsideSecond };
struct PartTrace {
    std::atomic_bool* cancellation = nullptr;
    CancelChildAt cancelAt = CancelChildAt::Never;
    int activeOriginalPart = -1;
    unsigned begun = 0, ended = 0;
    bool fired = false;
};
thread_local PartTrace* partTrace = nullptr;

void observePart(const NativeCadStageEvent& event) noexcept {
    if(!partTrace || !event.stage) return;
    auto& state = *partTrace;
    const std::string_view stage(event.stage);
    bool cancel = false;
    if(stage == "contact_part_meshing") {
        if(event.beginning) {
            state.activeOriginalPart = event.entityTag;
            ++state.begun;
            cancel = event.entityTag == 2 && state.cancelAt == CancelChildAt::SecondBegin;
        } else {
            ++state.ended;
            cancel = event.entityTag == 2 && state.cancelAt == CancelChildAt::SecondEnd;
        }
    } else if(stage == "cad_validation" && event.beginning &&
              state.activeOriginalPart == 2 && state.cancelAt == CancelChildAt::InsideSecond) {
        cancel = true;
    }
    if(cancel && !state.fired) {
        state.fired = true;
        state.cancellation->store(true, std::memory_order_relaxed);
    }
}

class PartTraceScope {
public:
    explicit PartTraceScope(PartTrace& state) noexcept : previous_(partTrace) { partTrace = &state; }
    ~PartTraceScope() { partTrace = previous_; }
    PartTraceScope(const PartTraceScope&) = delete;
    PartTraceScope& operator=(const PartTraceScope&) = delete;
private:
    PartTrace* previous_;
};

// Compare the complete semantic payload bit for bit, without reading struct
// padding or substituting a hash for the coordinate/connectivity comparison.
bool sameMeshPayload(const mesh::Mesh& a, const mesh::Mesh& b) {
    if(a.nodeCount() != b.nodeCount() || a.cellCount() != b.cellCount() ||
       a.connectivityCount() != b.connectivityCount()) return false;
    for(Index i = 0; i < a.nodeCount(); ++i) {
        const auto p = a.node(i), q = b.node(i);
        if(std::bit_cast<std::uint64_t>(p.x) != std::bit_cast<std::uint64_t>(q.x) ||
           std::bit_cast<std::uint64_t>(p.y) != std::bit_cast<std::uint64_t>(q.y) ||
           std::bit_cast<std::uint64_t>(p.z) != std::bit_cast<std::uint64_t>(q.z) ||
           a.nodeClass(i) != b.nodeClass(i)) return false;
    }
    for(Index i = 0; i < a.cellCount(); ++i) {
        const auto& p = a.cell(i);
        const auto& q = b.cell(i);
        if(p.type != q.type || p.first != q.first || p.onEntity != q.onEntity) return false;
        for(int j = 0; j < mesh::nodesPerCell(p.type); ++j)
            if(a.cellNodes(i)[j] != b.cellNodes(i)[j]) return false;
    }
    return true;
}

void checkCancelledPartReceipts() {
    BRep_Builder builder;
    TopoDS_Compound shape;
    builder.MakeCompound(shape);
    for(int i = 0; i < 3; ++i)
        builder.Add(shape, BRepPrimAPI_MakeBox(gp_Pnt(3. * i, 0., 0.), 1., 1., 1.).Shape());
    auto imported = OcctModelAccess::fromShape(shape);
    check(imported.imported() && imported.model->solids().size() == 3, "cancelled child fixture import");
    const auto& model = *imported.model;
    const auto original = snapshot(model);
    NativeCadMeshOptions options;
    options.assemblyPolicy = NativeCadAssemblyPolicy::IndependentParts;
    options.useAdvancingInteriorFront = false; // Isolate wrapper telemetry from AUTO's aggregate reason policy.
    options.geometry.targetSize = .5;
    options.geometry.threads = 1;
    options.maximumNativeValidationThreads = 1;
    options.maximumValidationWork = 10'000'000;
    options.geometry.maximumMeshingWork = options.maximumValidationWork;

    mesh::Mesh reference;
    const auto referenceResult = meshNativeCad(model, options, reference);
    check(referenceResult.published && referenceResult.validatedParts == 3 && referenceResult.bodies.size() == 3 &&
          referenceResult.audit.valid && referenceResult.nativeContactSourceGuard && snapshot(model) == original,
          "three-part reference publication and source identity");

    std::atomic_bool cancellation{false};
    auto recording = options;
    recording.geometry.cancellation = &cancellation;
    recording.stageObserver = observePart;
    PartTrace noCancellation{&cancellation};
    mesh::Mesh observed;
    NativeCadMeshOutcome observedResult;
    {
        PartTraceScope tracing(noCancellation);
        observedResult = meshNativeCad(model, recording, observed);
    }
    check(!noCancellation.fired && noCancellation.begun == 3 && noCancellation.ended == 3 &&
          observedResult.published && observedResult.validationWork == referenceResult.validationWork &&
          sameMeshPayload(reference, observed) && snapshot(model) == original,
          "recording observer preserves exact mesh payload, work and original source bytes");
    check(observedResult.interferences.size() == referenceResult.interferences.size(), "observer interference count");
    for(std::size_t i = 0; i < observedResult.interferences.size(); ++i) {
        const auto& p = observedResult.interferences[i];
        const auto& q = referenceResult.interferences[i];
        check(p.firstSolidTag == q.firstSolidTag && p.secondSolidTag == q.secondSolidTag &&
              p.state == q.state && p.nativeOverlapVolume == q.nativeOverlapVolume,
              "observer preserves interference evidence");
    }
    for(std::size_t i = 0; i < observedResult.bodies.size(); ++i) {
        const auto& p = observedResult.bodies[i];
        const auto& q = referenceResult.bodies[i];
        check(p.solidTag == q.solidTag && p.reason == q.reason &&
              p.nativeCentroidsChecked == q.nativeCentroidsChecked &&
              p.nativeCentroidSourceGuard == q.nativeCentroidSourceGuard,
              "observer preserves per-part Native certificate counts");
    }

    // A sentinel destination must survive every refusal exactly.
    mesh::Mesh destination;
    destination.addNode({9., 8., 7.}, {3, 17});
    const mesh::Mesh sentinel = destination;
    auto cancelledAt = [&](CancelChildAt point) {
        cancellation.store(false, std::memory_order_relaxed);
        PartTrace state{&cancellation, point};
        NativeCadMeshOutcome result;
        {
            PartTraceScope tracing(state);
            result = meshNativeCad(model, recording, destination);
        }
        check(state.fired && state.begun == 2 && state.ended == 2,
              "synchronous cancellation entered only two owned children");
        check(!result.published && result.destinationUnchanged && result.reason == NativeCadMeshReason::Cancelled &&
              result.validatedParts == 1 && sameMeshPayload(destination, sentinel) && snapshot(model) == original &&
              result.nativeContactSourceGuard && result.validationWork <= options.maximumValidationWork,
              "child cancellation preserves aggregate resource, source and publication gates");
        check(result.bodies.size() == 2 && result.bodies[0].solidTag == 1 && result.bodies[1].solidTag == 2 &&
              result.bodies[0].nativeCentroidSourceGuard && result.bodies[0].nativeCentroidsChecked > 0,
              "cancelled second child receipt survives beside the validated first part");
        return result;
    };

    const auto beforeSecond = cancelledAt(CancelChildAt::SecondBegin);
    check(beforeSecond.bodies[1].reason == NativeCadMeshReason::Cancelled,
          "child cancelled before construction keeps its refusal");
    const auto afterSecond = cancelledAt(CancelChildAt::SecondEnd);
    check(afterSecond.validationWork > beforeSecond.validationWork &&
          afterSecond.validationWork < referenceResult.validationWork &&
          afterSecond.generationSeconds > 0. &&
          afterSecond.bodies[1].reason == NativeCadMeshReason::None,
          "completed second child work and generation survive cancellation before wrapper validation");
    const auto insideSecond = cancelledAt(CancelChildAt::InsideSecond);
    check(insideSecond.bodies[1].reason == NativeCadMeshReason::Cancelled &&
          insideSecond.validationWork > beforeSecond.validationWork,
          "inside-child cancellation retains the refused child and its completed work");

    cancellation.store(false, std::memory_order_relaxed);
    auto tiny = options;
    tiny.maximumValidationWork = 1;
    tiny.geometry.maximumMeshingWork = 1;
    const auto exhausted = meshNativeCad(model, tiny, destination);
    check(!exhausted.published && exhausted.reason == NativeCadMeshReason::BudgetExceeded &&
          exhausted.validationWork <= 1 && sameMeshPayload(destination, sentinel) && snapshot(model) == original,
          "independent-part aggregate work boundary stays transactional");
}
} // namespace

int main(){try{
 checkCancelledPartReceipts();
 BRep_Builder builder;TopoDS_Compound shape;builder.MakeCompound(shape);builder.Add(shape,BRepPrimAPI_MakeBox(1.,1.,1.).Shape());builder.Add(shape,BRepPrimAPI_MakeBox(gp_Pnt(.5,0,0),1.,1.,1.).Shape());auto imported=OcctModelAccess::fromShape(shape);check(imported.imported()&&imported.model->solids().size()==2,"overlapping assembly import");const auto& model=*imported.model;const auto source=snapshot(model);
 NativeCadMeshOptions o;o.geometry.targetSize=.6;o.maximumValidationWork=10000000;o.geometry.maximumMeshingWork=10000000;o.maximumNativeValidationThreads=2;
 mesh::Mesh output;output.addNode({9,8,7});const auto untouched=mesh::meshHash(output);auto rejected=meshNativeCad(model,o,output);check(!rejected.published&&rejected.reason==NativeCadMeshReason::UnresolvedBodyInterface&&mesh::meshHash(output)==untouched,"default conforming overlap refusal");
 o.assemblyPolicy=NativeCadAssemblyPolicy::IndependentParts;auto result=meshNativeCad(model,o,output);if(!result.published)std::cerr<<nativeCadMeshReasonName(result.reason)<<' '<<result.nativeContactDiagnostic<<'\n';check(result.published&&result.assemblyPolicy==NativeCadAssemblyPolicy::IndependentParts&&result.validatedParts==2&&result.bodies.size()==2,"independent assembly publication");check(result.interferences.size()==1&&result.interferences[0].state==NativeCadInterferenceState::MaterialOverlap&&std::abs(result.interferences[0].nativeOverlapVolume-.5)<1e-10,"native interference report");check(result.audit.valid&&std::abs(result.quality.totalVolume-2.)<1e-10,"part sum is not union volume");check(result.nativeContactSourceGuard&&snapshot(model)==source,"assembly source preserved");
 std::set<Index> nodes[2];for(Index i=0;i<output.cellCount();++i){const auto& c=output.cell(i);if(c.type!=mesh::CellType::Tet4)continue;check(c.onEntity.tag==1||c.onEntity.tag==2,"original body tag");for(unsigned j=0;j<4;++j)nodes[c.onEntity.tag-1].insert(output.cellNodes(i)[j]);}for(auto i:nodes[0])check(!nodes[1].contains(i),"parts were welded");for(const auto& b:result.bodies)check(b.nativeCentroidSourceGuard&&b.nativeCentroidsChecked>0,"original body centroid certificate");
 auto adaptive=o;
 adaptive.useAdaptiveTetRecovery=true;
 adaptive.reuseNativeQueries=true;
 mesh::Mesh adaptiveMesh;
 const auto adaptiveResult=meshNativeCad(model,adaptive,adaptiveMesh);
 check(adaptiveResult.published&&adaptiveResult.validatedParts==2&&mesh::meshHash(adaptiveMesh)==mesh::meshHash(output),
       "adaptive assembly preserves independently meshed body identities");
 check(adaptiveResult.interferences.size()==1&&adaptiveResult.interferences[0].state==NativeCadInterferenceState::MaterialOverlap&&
       adaptiveResult.interferences[0].nativeOverlapVolume==result.interferences[0].nativeOverlapVolume&&snapshot(model)==source,
       "adaptive assembly preserves native interference reporting");
 const auto saved=mesh::meshHash(output);auto low=o;low.maximumValidationWork=1;auto budget=meshNativeCad(model,low,output);check(!budget.published&&budget.reason==NativeCadMeshReason::BudgetExceeded&&mesh::meshHash(output)==saved,"assembly budget transaction");
 std::atomic_bool cancellation{true};low=o;low.geometry.cancellation=&cancellation;auto cancelled=meshNativeCad(model,low,output);check(!cancelled.published&&cancelled.reason==NativeCadMeshReason::Cancelled&&mesh::meshHash(output)==saved,"assembly cancellation transaction");
 low=o;low.geometry.minimumElementQuality=1.;auto strict=meshNativeCad(model,low,output);check(!strict.published&&strict.reason==NativeCadMeshReason::QualityRequirementNotMet&&mesh::meshHash(output)==saved,"assembly strict quality gate");
 low=o;low.faceSizes.push_back({999999,.1,1.3});auto bad=meshNativeCad(model,low,output);check(!bad.published&&bad.reason==NativeCadMeshReason::InvalidRequest&&mesh::meshHash(output)==saved&&snapshot(model)==source,"unknown selection not silently discarded");
 gp_Trsf transform;transform.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,0,1)),.4);transform.SetTranslationPart(gp_Vec(5,2,3));TopoDS_Compound located;builder.MakeCompound(located);builder.Add(located,BRepPrimAPI_MakeBox(1.,1.,1.).Shape());builder.Add(located,BRepPrimAPI_MakeBox(1.,1.,1.).Shape().Moved(TopLoc_Location(transform)));auto placed=OcctModelAccess::fromShape(located);check(placed.imported(),"located assembly import");mesh::Mesh placedMesh;const auto placementBytes=snapshot(*placed.model);auto placedResult=meshNativeCad(*placed.model,o,placedMesh);if(!placedResult.published)std::cerr<<nativeCadMeshReasonName(placedResult.reason)<<' '<<placedResult.nativeContactDiagnostic<<'\n';check(placedResult.published&&placedResult.validatedParts==2&&placedResult.interferences.size()==1&&placedResult.interferences[0].state==NativeCadInterferenceState::NoMaterialOverlap&&snapshot(*placed.model)==placementBytes,"located part geometry/history");
 TopoDS_Compound cylinders;builder.MakeCompound(cylinders);builder.Add(cylinders,BRepPrimAPI_MakeCylinder(1.,2.).Shape());builder.Add(cylinders,BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(5,0,0),gp_Dir(0,0,1)),1.,2.).Shape());auto sized=OcctModelAccess::fromShape(cylinders);check(sized.imported(),"sized assembly import");int selectedCap=-1;for(const auto* face:sized.model->faces())if(face->kind()==SurfaceKind::Plane){for(auto owner:sized.model->faceOwners(face->tag()))if(owner.solidTag==1)selectedCap=face->tag();if(selectedCap>0)break;}check(selectedCap>0,"original selected cap identity");auto sizing=o;sizing.minimumSize=.2;sizing.faceSizes.push_back({selectedCap,.25,1.3});mesh::Mesh sizedMesh;const auto sizedSource=snapshot(*sized.model);const auto sizedResult=meshNativeCad(*sized.model,sizing,sizedMesh);if(!sizedResult.published)std::cerr<<nativeCadMeshReasonName(sizedResult.reason)<<' '<<sizedResult.nativeContactDiagnostic<<'\n';check(sizedResult.published&&sizedResult.validatedParts==2&&snapshot(*sized.model)==sizedSource,"selected sizing must not reject unselected parts");unsigned cellsByPart[2]{};for(Index i=0;i<sizedMesh.cellCount();++i)if(sizedMesh.cell(i).type==mesh::CellType::Tet4){const int part=sizedMesh.cell(i).onEntity.tag;check(part==1||part==2,"sized part identity");++cellsByPart[part-1];}const auto& selectedProfile=sizedResult.bodies.at(0).primitive;check(cellsByPart[0]>0&&cellsByPart[1]>0&&selectedProfile.profiledVolume&&!sizedResult.bodies.at(1).primitive.profiledVolume&&selectedProfile.maximumSectionNodes>selectedProfile.minimumSectionNodes,"selected-face sizing applies to its original part only");const auto sizedHash=mesh::meshHash(sizedMesh);sizing.faceSizes.clear();const auto missingField=meshNativeCad(*sized.model,sizing,sizedMesh);check(!missingField.published&&missingField.reason==NativeCadMeshReason::UnsupportedControl&&mesh::meshHash(sizedMesh)==sizedHash,"minimum-size request without a face control remains explicit");
 TopoDS_Compound auxiliary;builder.MakeCompound(auxiliary);builder.Add(auxiliary,shape);builder.Add(auxiliary,BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(10,10,10),gp_Dir(0,0,1)),-1.,1.,-1.,1.).Face());NativeCadImportOptions auxiliaryOptions;auxiliaryOptions.retainAuxiliaryGeometry=true;auto withAuxiliary=OcctModelAccess::fromShape(auxiliary,auxiliaryOptions);check(withAuxiliary.imported(),"assembly with retained auxiliary face");int auxiliaryFace=-1;for(const auto* face:withAuxiliary.model->faces())if(withAuxiliary.model->faceOwners(face->tag()).empty())auxiliaryFace=face->tag();check(auxiliaryFace>0,"retained auxiliary face identity");auto auxiliaryControl=o;auxiliaryControl.faceSizes.push_back({auxiliaryFace,.25,1.3});const auto auxiliaryResult=meshNativeCad(*withAuxiliary.model,auxiliaryControl,sizedMesh);check(!auxiliaryResult.published&&auxiliaryResult.reason==NativeCadMeshReason::InvalidRequest&&mesh::meshHash(sizedMesh)==sizedHash,"auxiliary sizing selection must not be silently discarded");
 std::cout<<"Independent part identity, interference, sizing, native validation and transactional gates passed\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
