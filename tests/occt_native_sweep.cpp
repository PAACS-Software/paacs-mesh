// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctModelInternal.hpp"
#include "OcctNativeSweep.hpp"
#include "paacs/geometry/OcctMesher.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <Geom_BSplineCurve.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <cmath>
#include <array>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <gp_Ax1.hxx>
#include <gp_Trsf.hxx>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>
using namespace paacs;
namespace {
using namespace geo::detail;
void require(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
std::string snapshot(const geo::OcctModel &m) {
    std::ostringstream stream;
    BRepTools::Write(OcctModelAccess::shape(m), stream);
    return stream.str();
}
mesh::Mesh sentinel() {
    mesh::Mesh m;
    m.addNode({17, 19, 23}, {0, 777});
    return m;
}
void unchanged(const mesh::Mesh &m) {
    require(m.nodeCount() == 1 && m.cellCount() == 0 && m.node(0).x == 17 && m.node(0).y == 19 &&
                m.node(0).z == 23 && m.nodeClass(0) == mesh::Classification{0, 777},
            "destination rollback");
}
void validate(const geo::OcctModel &model, const OcctNativeSweepOutcome &out,
              const mesh::Mesh &mesh) {
    require(out.published && !out.destinationUnchanged && out.reason == OcctNativeSweepReason::None,
            out.diagnostic);
    require(out.sourceImmutabilityVerified && out.audit.valid && out.quality.valid &&
                out.exteriorCertificate.certified,
            "all final gates");
    require(out.expectedRelativeError <= 1e-6 && out.commonRelativeError <= 1e-6 &&
                out.relativeVolumeError <= .01,
            "native recognition and mesh volume");
    require(out.nativeInteriorQueries >= mesh.countCells(mesh::CellType::Tet4),
            "every tet centroid native checked");
    require(out.maximumBoundaryNodeError <= 1e-6 && out.maximumBoundarySampleError <= .05,
            "original CAD node and facet samples");
    require(out.workUsed ==
                out.recognitionWork + out.surfaceWork + out.sweepWork + out.validationWork,
            "shared work receipt");
    for (Index i = 0; i < mesh.nodeCount(); ++i) {
        const auto c = mesh.nodeClass(i);
        require(c.dim == 0   ? model.vertex(c.tag) != nullptr
                : c.dim == 1 ? model.edge(c.tag) != nullptr
                : c.dim == 2 ? model.face(c.tag) != nullptr
                             : c.dim == 3 && model.solid(c.tag) != nullptr,
                "original entity ownership");
    }
}

void sizedRotatedHollowCaps() {
    constexpr double wall=1.,length=8.,size=2.;
    for(const auto dimensions:std::array{std::array{12.,12.},std::array{18.,8.}}){
        const auto [width,height]=dimensions;
        const auto outer=BRepPrimAPI_MakeBox(gp_Pnt(-width/2,-height/2,0),width,height,length).Shape();
        const auto inner=BRepPrimAPI_MakeBox(gp_Pnt(wall-width/2,wall-height/2,-1),
            width-2*wall,height-2*wall,length+2).Shape();
        BRepAlgoAPI_Cut cut(outer,inner);
        require(cut.IsDone(),"owned hollow cap cut failed");
        gp_Trsf placement;
        placement.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(1,2,3)),.61);
        placement.SetTranslationPart(gp_Vec(11,-7,19));
        BRepBuilderAPI_Transform moved(cut.Shape(),placement,true);
        auto imported=OcctModelAccess::fromShape(moved.Shape());
        require(imported.imported(),"owned rotated hollow cap import failed");
        const auto before=snapshot(*imported.model);
        OcctNativeSweepOptions request;
        request.maximumWork=50'000'000;
        request.maximumLayerThickness=size;
        request.boundary.linearDeflection=.01;
        request.boundary.maximumChordError=.05;
        mesh::Mesh coarse,sized;
        const auto original=buildOcctNativeSweepVolume(*imported.model,1,request,coarse);
        validate(*imported.model,original,coarse);
        require(!original.capSizingApplied&&original.requestedCapSize==0&&original.capSizingWork==0&&
            std::string_view(original.capSourceProvider)=="native_tessellation",
            "untouched cap sweep no longer retains the historical source");
        request.maximumCapSize=size;
        request.boundary.construction.resourcePolicy=geo::CadResourcePolicy::Normal;
        const auto refined=buildOcctNativeSweepVolume(*imported.model,1,request,sized);
        validate(*imported.model,refined,sized);
        require(refined.requestedCapSize==size&&refined.capSizingApplied&&
            std::string_view(refined.capSourceProvider)=="regular_chart"&&
            std::string_view(refined.capSizingDiagnostic)=="regular_cap_source_applied"&&
            refined.capSizingWork==refined.surfaceWork&&refined.capSizingSeconds==refined.surfaceSeconds&&
            refined.capSizingWork>0,
            "regular cap source lacks a separate, nonduplicated sizing/work receipt");
        require(refined.capTriangles>original.capTriangles&&refined.capNodes>original.capNodes&&
            refined.maximumCapEdgeLength>0&&refined.maximumCapEdgeLength<=2*size,
            "physically sized cap retained the original oversized tangential triangles");
        require(refined.quality.minimumQuality>=original.quality.minimumQuality&&
            refined.quality.medianQuality>original.quality.medianQuality+.1,
            "owned cap sizing did not improve the median while preserving its low tail");
        const auto expected=(width*height-(width-2*wall)*(height-2*wall))*length;
        require(std::abs(refined.meshVolume-expected)<=expected*1e-10&&
            refined.exteriorCertificate.eulerCharacteristic==0,
            "physically sized cap filled its bore or changed the analytic material volume");
        const auto inverse=placement.Inverted();
        for(Index i=0;i<sized.cellCount();++i)if(sized.cell(i).type==mesh::CellType::Tet4){
            Vec3 center{};const auto* ids=sized.cellNodes(i);
            for(unsigned j=0;j<4;++j)center=center+sized.node(ids[j])*.25;
            gp_Pnt local(center.x,center.y,center.z);local.Transform(inverse);
            require(!(std::abs(local.X())<width/2-wall-1e-7&&
                      std::abs(local.Y())<height/2-wall-1e-7),
                    "a sized-cap tetrahedron centroid occupies the analytic through-hole");
        }
        require(snapshot(*imported.model)==before,"regular cap source mutated the original CAD");
        auto limited=request;limited.boundary.maximumWork=1;
        auto destination=sentinel();
        const auto refused=buildOcctNativeSweepVolume(*imported.model,1,limited,destination);
        require(!refused.published&&refused.destinationUnchanged&&!refused.capSizingApplied&&
            refused.reason==OcctNativeSweepReason::BudgetExceeded&&
            std::string_view(refused.capSourceProvider)=="regular_chart"&&
            refused.capSizingWork==refused.surfaceWork&&
            refused.workUsed==refused.recognitionWork+refused.surfaceWork+refused.sweepWork+refused.validationWork,
            "regular cap source refusal lost its charged cost or silently used a coarse provider");
        unchanged(destination);
        require(snapshot(*imported.model)==before,"refused regular cap source changed CAD");
        std::cout<<"sized rotated hollow "<<width<<'x'<<height<<": cap "<<original.capTriangles
            <<" -> "<<refined.capTriangles<<", median "<<original.quality.medianQuality<<" -> "
            <<refined.quality.medianQuality<<", min "<<original.quality.minimumQuality<<" -> "
            <<refined.quality.minimumQuality<<", work="<<refined.workUsed
            <<", seconds="<<refined.totalSeconds<<'\n';
    }
}

// The observer fires after the private cap call returns, before the caller
// accounts for that call and polls cancellation. This tests the public native
// dispatch receipt rather than cancellation inside the private source builder.
thread_local std::atomic_bool* cancelAtCapCompletion=nullptr;
thread_local unsigned completedSizedCaps=0;
void cancelCompletedSizedCap(const geo::NativeCadStageEvent& event) noexcept {
    if(!event.beginning&&std::string_view(event.stage)=="sized_cap_sweep"&&
       std::string_view(event.diagnostic)=="complete"){
        ++completedSizedCaps;
        if(cancelAtCapCompletion)cancelAtCapCompletion->store(true,std::memory_order_relaxed);
    }
}

void wholeNativeCapReceipts() {
    constexpr double width=60,height=30,wall=1,length=40,size=4;
    const auto outer=BRepPrimAPI_MakeBox(gp_Pnt(-width/2,-height/2,0),width,height,length).Shape();
    const auto inner=BRepPrimAPI_MakeBox(gp_Pnt(wall-width/2,wall-height/2,-1),
        width-2*wall,height-2*wall,length+2).Shape();
    BRepAlgoAPI_Cut cut(outer,inner);
    require(cut.IsDone(),"whole native cap owned annulus cut failed");
    gp_Trsf placement;
    placement.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(1,2,3)),.73);
    placement.SetTranslationPart(gp_Vec(11,-7,5));
    BRepBuilderAPI_Transform moved(cut.Shape(),placement,true);
    // Keep the owned representation that exercises the imported native cap
    // portfolio. A directly constructed BREP can take the separate planar
    // route, so its success would not test the cap accounting boundary.
    struct OwnedStepFile {
        std::filesystem::path path;
        bool owned=false;
        ~OwnedStepFile(){if(owned){std::error_code ignored;std::filesystem::remove(path,ignored);}}
    } step{std::filesystem::temp_directory_path()/
        ("paacs-native-cap-receipt-"+std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count())+".step")};
    require(!std::filesystem::exists(step.path),"owned native cap STEP temporary name already exists");
    step.owned=true;
    STEPControl_Writer writer;
    require(writer.Transfer(moved.Shape(),STEPControl_AsIs)==IFSelect_RetDone&&
        writer.Write(step.path.string().c_str())==IFSelect_RetDone,
        "owned native cap STEP round-trip write failed");
    STEPControl_Reader reader;
    require(reader.ReadFile(step.path.string().c_str())==IFSelect_RetDone&&reader.TransferRoots()>0,
        "owned native cap STEP round-trip read failed");
    auto imported=OcctModelAccess::fromShape(reader.OneShape());
    require(imported.imported(),"whole native cap owned annulus import failed");
    const auto before=snapshot(*imported.model);
    auto request=geo::normalNativeCadMeshOptions();
    request.geometry.targetSize=size;
    request.maximumChordError=.1;
    request.maximumPrimitiveWorkingBytes=512ULL*1024*1024;
    request.geometry.maximumMeshingWork=request.maximumValidationWork=50'000'000;
    // Isolate the cap portfolio: a hollow annulus has no strict star kernel,
    // while the legacy planar route cannot represent this oblique source.
    request.useAdvancingInteriorFront=false;
    request.useNativePartitionFallback=false;
    request.useNativeKernelAdmission=true;
    std::atomic_bool cancellation{false};
    request.geometry.cancellation=&cancellation;
    request.stageObserver=cancelCompletedSizedCap;
    auto destination=sentinel();
    completedSizedCaps=0;
    cancelAtCapCompletion=&cancellation;
    const auto cancelled=geo::meshNativeCad(*imported.model,request,destination);
    cancelAtCapCompletion=nullptr;
    if(completedSizedCaps!=1){
        std::cerr<<"whole cap diagnostic: reason="<<int(cancelled.reason)<<", published="<<cancelled.published;
        for(const auto& body:cancelled.bodies)std::cerr<<", route="<<body.route
            <<", cap="<<body.nativeCapSizingDiagnostic<<", attempts="<<body.nativeCapSizingAttempts;
        std::cerr<<'\n';
    }
    require(completedSizedCaps==1&&cancellation.load(std::memory_order_relaxed),
        "whole native annulus did not reach the completed sized cap stage");
    require(!cancelled.published&&cancelled.destinationUnchanged&&
        cancelled.reason==geo::NativeCadMeshReason::Cancelled&&cancelled.bodies.size()==1,
        "cap completion cancellation did not roll back the whole request");
    const auto& stopped=cancelled.bodies.front();
    require(stopped.nativeCapRequestedSize==size&&stopped.nativeCapSizingAttempts==1&&
        stopped.nativeCapSizingWork>0&&stopped.nativeCapSizingSeconds>0&&
        std::string_view(stopped.nativeCapSizingDiagnostic)=="complete",
        "cancellation erased the completed cap attempt receipt");
    require(stopped.nativeSurfaceWork>0&&stopped.nativeVolumeWork>0&&
        stopped.nativeSurfaceSeconds>0&&stopped.nativeVolumeSeconds>0&&
        stopped.nativeSurfaceWork+stopped.nativeVolumeWork>=stopped.nativeCapSizingWork&&
        cancelled.validationWork>=stopped.nativeCapSizingWork&&
        cancelled.generationSeconds>=stopped.nativeCapSizingSeconds,
        "cancelled cap work or surface/volume stage costs were lost");
    unchanged(destination);
    require(snapshot(*imported.model)==before,"cancelled cap completion changed original CAD");

    // A whole-source node ceiling can reject the regular cap provider even
    // though the established coarse cap and its ten axial layers fit. Both
    // attempts must stay in the same hard work and output admissions.
    request.geometry.cancellation=nullptr;
    request.stageObserver=nullptr;
    request.geometry.maximumNodes=100;
    destination=sentinel();
    const auto fallback=geo::meshNativeCad(*imported.model,request,destination);
    require(fallback.published&&!fallback.destinationUnchanged&&
        fallback.reason==geo::NativeCadMeshReason::None&&fallback.bodies.size()==1&&
        fallback.audit.valid&&fallback.quality.valid,
        "bounded regular-source refusal did not preserve the valid legacy cap");
    const auto& chosen=fallback.bodies.front();
    require(std::string_view(chosen.route)=="native_cap_sweep"&&
        chosen.nativeCapRequestedSize==size&&chosen.nativeCapSizingAttempts==2&&
        !chosen.nativeCapSizingApplied&&chosen.nativeCapNodes>0&&chosen.nativeCapTriangles>0&&
        std::string_view(chosen.nativeCapSizingDiagnostic)!="complete"&&
        std::string_view(chosen.nativeCapSizingDiagnostic)!="not_attempted",
        "coarse fallback hid the sized source refusal or claimed regular sizing");
    require(chosen.nativeCapSizingWork>0&&chosen.nativeCapSizingSeconds>0&&
        chosen.nativeSurfaceWork>0&&chosen.nativeVolumeWork>0&&
        chosen.nativeSurfaceWork+chosen.nativeVolumeWork>=chosen.nativeCapSizingWork&&
        fallback.validationWork>=chosen.nativeCapSizingWork&&
        fallback.validationWork<=request.maximumValidationWork&&
        destination.nodeCount()<=request.geometry.maximumNodes&&
        fallback.generationSeconds>=chosen.nativeCapSizingSeconds,
        "legacy cap fallback lost retry costs or exceeded shared hard admissions");
    require(snapshot(*imported.model)==before,"regular-source fallback changed original CAD");
}
} // namespace
int main(int argc, char **argv) {
    try {
        using namespace geo::detail;
        sizedRotatedHollowCaps();
        wholeNativeCapReceipts();
        OcctNativeSweepOptions options;
        options.boundary.linearDeflection = .01;
        options.layers = 3;
        auto cube = OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(2, 3, 4).Shape());
        require(cube.imported(), "cube snapshot");
        const auto cubeBefore = snapshot(*cube.model);
        mesh::Mesh cubeMesh;
        auto result = buildOcctNativeSweepVolume(*cube.model, 1, options, cubeMesh);
        validate(*cube.model, result, cubeMesh);
        require(std::abs(result.meshVolume - 24) < 1e-10, "cube exact volume");
        require(result.layers == 3, "explicit multilayers");
        require(snapshot(*cube.model) == cubeBefore, "cube immutable");
        auto automatic = options;
        automatic.layers = 1;
        automatic.maximumLayerThickness = .6;
        mesh::Mesh automaticMesh;
        auto automaticResult = buildOcctNativeSweepVolume(*cube.model, 1, automatic, automaticMesh);
        validate(*cube.model, automaticResult, automaticMesh);
        require(automaticResult.layers == std::ceil(automaticResult.height / .6) &&
                    automaticResult.height / automaticResult.layers <= .6,
                "derived normal layers honor admitted thickness");
        // Native STEP exporters can encode straight edges as splines. Keep
        // the original planar faces and pcurves, replacing only equal-value
        // 3D line representations on an owned fixture copy.
        BRepBuilderAPI_Copy splineCopy(OcctModelAccess::shape(*cube.model), true, false);
        auto splineBox = splineCopy.Shape();
        BRep_Builder builder;
        for (TopExp_Explorer explorer(splineBox, TopAbs_EDGE); explorer.More(); explorer.Next()) {
            auto edge = TopoDS::Edge(explorer.Current());
            BRepAdaptor_Curve curve(edge);
            if (curve.GetType() != GeomAbs_Line)
                continue;
            NCollection_Array1<gp_Pnt> poles(1, 2);
            NCollection_Array1<double> knots(1, 2);
            NCollection_Array1<int> multiplicities(1, 2);
            poles(1) = curve.Value(curve.FirstParameter());
            poles(2) = curve.Value(curve.LastParameter());
            knots(1) = curve.FirstParameter();
            knots(2) = curve.LastParameter();
            multiplicities(1) = multiplicities(2) = 2;
            occ::handle<Geom_BSplineCurve> spline =
                new Geom_BSplineCurve(poles, knots, multiplicities, 1);
            builder.UpdateEdge(edge, spline, 1e-7);
            builder.Range(edge, knots(1), knots(2));
        }
        auto splineModel = OcctModelAccess::fromShape(splineBox);
        require(splineModel.imported(), "native straight spline representation");
        mesh::Mesh splineMesh;
        result = buildOcctNativeSweepVolume(*splineModel.model, 1, options, splineMesh);
        validate(*splineModel.model, result, splineMesh);
        require(std::abs(result.meshVolume - 24) < 1e-10,
                "spline line support preserves exact box");
        const auto box = BRepPrimAPI_MakeBox(gp_Pnt(-3, -3, 0), 6, 6, 2).Shape();
        const auto cylinder =
            BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -1), gp_Dir(0, 0, 1)), 1, 4).Shape();
        BRepAlgoAPI_Cut cut(box, cylinder);
        require(cut.IsDone(), "through hole fixture");
        gp_Trsf transform;
        transform.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), .73);
        transform.SetTranslationPart(gp_Vec(11, -7, 19));
        BRepBuilderAPI_Transform moved(cut.Shape(), transform, true);
        auto hole = OcctModelAccess::fromShape(moved.Shape());
        require(hole.imported(), "rotated through bore");
        const auto holeBefore = snapshot(*hole.model);
        mesh::Mesh holeMesh;
        result = buildOcctNativeSweepVolume(*hole.model, 1, options, holeMesh);
        validate(*hole.model, result, holeMesh);
        require(result.exteriorCertificate.eulerCharacteristic == 0, "through-hole genus retained");
        require(snapshot(*hole.model) == holeBefore, "rotated bore immutable");
        auto cone = OcctModelAccess::fromShape(BRepPrimAPI_MakeCone(2, 1, 3).Shape());
        require(cone.imported(), "cone snapshot");
        auto old = sentinel();
        auto failed = buildOcctNativeSweepVolume(*cone.model, 1, options, old);
        require(!failed.published && failed.reason == OcctNativeSweepReason::UnsupportedGeometry,
                "nonuniform cone refused");
        unchanged(old);
        auto sphere = OcctModelAccess::fromShape(BRepPrimAPI_MakeSphere(2).Shape());
        require(sphere.imported(), "sphere snapshot");
        old = sentinel();
        failed = buildOcctNativeSweepVolume(*sphere.model, 1, options, old);
        require(!failed.published && failed.reason == OcctNativeSweepReason::UnsupportedGeometry,
                "missing planar caps refused");
        unchanged(old);
        auto refusal = [&](OcctNativeSweepOptions request, OcctNativeSweepReason expected) {
            auto destination = sentinel();
            const auto receipt = buildOcctNativeSweepVolume(*cube.model, 1, request, destination);
            require(!receipt.published && receipt.destinationUnchanged &&
                        receipt.reason == expected,
                    receipt.diagnostic);
            unchanged(destination);
            require(snapshot(*cube.model) == cubeBefore, "refusal source immutable");
        };
        auto bounded = options;
        bounded.maximumWork = 0;
        refusal(bounded, OcctNativeSweepReason::BudgetExceeded);
        bounded = options;
        bounded.maximumWork = 1;
        refusal(bounded, OcctNativeSweepReason::BudgetExceeded);
        bounded = options;
        bounded.maximumNodes = 1;
        refusal(bounded, OcctNativeSweepReason::BudgetExceeded);
        bounded = options;
        bounded.maximumCells = 1;
        refusal(bounded, OcctNativeSweepReason::BudgetExceeded);
        bounded = options;
        bounded.maximumWorkingBytes = 1;
        refusal(bounded, OcctNativeSweepReason::BudgetExceeded);
        bounded = options;
        bounded.maximumSnapshotBytes = 1;
        refusal(bounded, OcctNativeSweepReason::BudgetExceeded);
        bounded = options;
        bounded.auditLimits.maximumFaceUses = 1;
        refusal(bounded, OcctNativeSweepReason::BudgetExceeded);
        bounded = options;
        bounded.minimumQuality = .999999;
        // Warning and required quality are independent thresholds.
        bounded.warningQuality = .25;
        refusal(bounded, OcctNativeSweepReason::QualityRequirement);
        bounded = options;
        bounded.layers = 0;
        refusal(bounded, OcctNativeSweepReason::InvalidRequest);
        bounded = options;
        bounded.layers = 4097;
        refusal(bounded, OcctNativeSweepReason::InvalidRequest);
        bounded = options;
        bounded.maximumLayerThickness = 1e-12;
        refusal(bounded, OcctNativeSweepReason::BudgetExceeded);
        bounded = options;
        bounded.maximumLayerThickness = -1;
        refusal(bounded, OcctNativeSweepReason::InvalidRequest);
        bounded = options;
        bounded.maximumCapSize = -1;
        refusal(bounded, OcctNativeSweepReason::InvalidRequest);
        bounded.maximumCapSize = std::numeric_limits<double>::infinity();
        refusal(bounded, OcctNativeSweepReason::InvalidRequest);
        std::atomic_bool cancelled{true};
        bounded = options;
        bounded.cancellation = &cancelled;
        refusal(bounded, OcctNativeSweepReason::Cancelled);
        if (argc == 2) {
            auto imported = geo::importOcctCad(argv[1]);
            require(imported.imported(), "actual STEP import");
            const auto before = snapshot(*imported.model);
            options.layers = 1;
            options.boundary.linearDeflection = .02;
            mesh::Mesh actual;
            result = buildOcctNativeSweepVolume(*imported.model, 1, options, actual);
            validate(*imported.model, result, actual);
            require(snapshot(*imported.model) == before, "actual STEP immutable");
            std::cout << "actual sweep nodes=" << result.nodes << " cells=" << result.cells
                      << " work=" << result.workUsed << " native-volume=" << result.cadVolume
                      << " mesh-volume-error=" << result.relativeVolumeError
                      << " seconds=" << result.totalSeconds << '\n';
        }
        std::cout << "native sweep ownership, holes, rotated geometry, nonuniform refusals, "
                     "budgets, cancellation and rollback PASS\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
