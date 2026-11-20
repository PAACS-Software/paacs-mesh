// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctModelInternal.hpp"
#include "OcctNativeSweep.hpp"
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
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <cmath>
#include <gp_Ax1.hxx>
#include <gp_Trsf.hxx>
#include <iostream>
#include <sstream>
#include <stdexcept>
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
} // namespace
int main(int argc, char **argv) {
    try {
        using namespace geo::detail;
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
