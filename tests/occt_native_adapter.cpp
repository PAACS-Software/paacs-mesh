// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/OcctModel.hpp"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <Interface_Static.hxx>
#include <STEPControl_Writer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS_Compound.hxx>
#include <chrono>
#include <fstream>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <stdexcept>
using namespace paacs;
namespace {
void require(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
bool close(double a, double b, double tolerance = 1e-7) {
    return std::abs(a - b) <= tolerance * std::max(1., std::abs(b));
}
struct UnitGuard {
    std::string old;
    UnitGuard() : old(Interface_Static::CVal("xstep.cascade.unit")) {
        require(Interface_Static::SetCVal("xstep.cascade.unit", "M"), "hostile global STEP units");
    }
    ~UnitGuard() { Interface_Static::SetCVal("xstep.cascade.unit", old.c_str()); }
};
geo::NativeCadImportOutcome importBrep(const TopoDS_Shape &shape,
                                       const std::filesystem::path &path) {
    require(BRepTools::Write(shape, path.string().c_str()), "write BREP fixture");
    geo::NativeCadImportOptions options;
    options.format = geo::CadFileFormat::Brep;
    auto result = geo::importOcctCad(path, options);
    require(result.imported(), geo::cadImportReasonName(result.reason));
    require(result.millimetresPerModelUnit == 0 && result.sourceLengthUnits.empty(),
            "BREP unspecified units");
    return result;
}
void inspectTopology(const geo::OcctModel &model) {
    require(model.identity() != 0, "import identity");
    for (const auto *edge : model.edges()) {
        double lo, hi;
        edge->parameterBounds(lo, hi);
        require(std::isfinite(lo) && std::isfinite(hi) && hi > lo, "native edge interval");
        require(model.vertex(edge->startVertexTag()) && model.vertex(edge->endVertexTag()),
                "edge endpoints resolve");
        require((edge->pointAt(lo) - model.vertex(edge->startVertexTag())->position()).norm() <
                    1e-6,
                "first endpoint matches native curve");
        require((edge->pointAt(hi) - model.vertex(edge->endVertexTag())->position()).norm() < 1e-6,
                "last endpoint matches native curve");
        require(model.edge(edge->tag()) == edge, "stable pointer identity");
        require(std::isfinite(model.tolerance(1, edge->tag())), "edge tolerance");
    }
    for (const auto *face : model.faces()) {
        require(model.face(face->tag()) == face && model.surfaceDescriptor(face->tag()),
                "stable face and analytic receipt");
        require(!model.boundaryUses(face->tag()).empty(), "face wires preserved");
        require(model.boundaryUses(face->tag()).front().outer, "outer wire first");
        for (const auto &wire : model.boundaryUses(face->tag()))
            for (const auto &use : wire.edges) {
                require(model.edge(use.edgeTag) && model.vertex(use.startVertexTag) &&
                            model.vertex(use.endVertexTag),
                        "oriented use identities resolve");
                Vec2 uv{91, 92};
                require(model.parameterOnUse(use.useId,
                                             (use.firstParameter + use.lastParameter) * .5, uv),
                        "pcurve agrees with native curve");
                require((face->pointAt(uv) -
                         model.edge(use.edgeTag)
                             ->pointAt((use.firstParameter + use.lastParameter) * .5))
                                .norm() < 1e-5,
                        "pcurve geometric agreement");
            }
    }
    for (const auto *body : model.solids()) {
        require(model.volume(body->tag()) > 0 && std::isfinite(model.volume(body->tag())),
                "exact positive body volume");
        require(!model.faceUses(body->tag()).empty(), "body face owners");
        for (const auto &use : model.faceUses(body->tag())) {
            require(model.face(use.faceTag) && std::abs(use.sense) == 1, "body face sense");
            bool owner = false;
            for (const auto &entry : model.faceOwners(use.faceTag))
                owner = owner || (entry.solidTag == body->tag() && entry.sense == use.sense);
            require(owner, "reciprocal face owner metadata");
        }
    }
}
} // namespace
int main(int argc, char **argv) {
    const auto directory =
        argc == 2
            ? std::filesystem::path(argv[1])
            : std::filesystem::temp_directory_path() /
                  ("paacs-native-occt-" +
                   std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    bool ownsDirectory = false;
    try {
        if (argc == 2)
            std::filesystem::create_directories(directory);
        else {
            ownsDirectory = std::filesystem::create_directory(directory);
            require(ownsDirectory, "exclusive fixtures");
        }
        const auto box = BRepPrimAPI_MakeBox(2., 3., 4.).Shape();
        auto basic = importBrep(box, directory / "box.brep");
        inspectTopology(*basic.model);
        require(close(basic.model->volume(1), 24), "box volume");
        const auto *segment = basic.model->edges().front();
        double segmentLo, segmentHi;
        segment->parameterBounds(segmentLo, segmentHi);
        const auto afterEnd = segment->pointAt(segmentHi) + segment->tangentAt(segmentHi) * 10;
        require(close(segment->closestParameter(afterEnd), segmentHi),
                "curve projection compares segment endpoints");
        std::size_t reversed = 0;
        const Vec3 center{1, 1.5, 2};
        for (const auto *face : basic.model->faces()) {
            Vec2 lo, hi;
            face->paramBounds(lo, hi);
            const Vec2 uv{(lo.u + hi.u) * .5, (lo.v + hi.v) * .5};
            const auto p = face->pointAt(uv), n = face->normalAt(uv);
            require((p - center).dot(n) > 0, "reversed BRep faces retain outward normals");
            reversed += basic.model->surfaceDescriptor(face->tag())->orientation ==
                        geo::CadOrientation::Reversed;
            Vec2 projectedUv{71, 72};
            Vec3 projected{81, 82, 83};
            require(face->projectPoint(p + n, projectedUv, projected) &&
                        (projected - p).norm() < 1e-7,
                    "exact trimmed planar projection");
            const auto oldUv = projectedUv;
            const auto oldPoint = projected;
            require(!face->projectPoint({(std::numeric_limits<double>::quiet_NaN)(), 0, 0},
                                        projectedUv, projected),
                    "nonfinite projection refusal");
            require(projectedUv.u == oldUv.u && projectedUv.v == oldUv.v &&
                        (projected - oldPoint).norm2() == 0,
                    "failed projection output rollback");
        }
        require(reversed > 0, "fixture exercises reversed faces");
        auto complement = box;
        complement.Reverse();
        require(BRepTools::Write(complement, (directory / "reversed-solid.brep").string().c_str()),
                "write reversed whole solid");
        geo::NativeCadImportOptions complementOptions;
        complementOptions.format = geo::CadFileFormat::Brep;
        auto unbounded = geo::importOcctCad(directory / "reversed-solid.brep", complementOptions);
        require(!unbounded.imported() && unbounded.reason == geo::CadImportReason::InvalidShape,
                "unbounded complement refuses without a partial model");
        const auto cylinder = BRepPrimAPI_MakeCylinder(2., 5.).Shape();
        auto curved = importBrep(cylinder, directory / "cylinder.brep");
        inspectTopology(*curved.model);
        require(close(curved.model->volume(1), 20 * std::acos(-1.)), "cylinder exact volume");
        const geo::GFace *lateral = nullptr;
        for (const auto *face : curved.model->faces())
            if (face->kind() == geo::SurfaceKind::Cylinder)
                lateral = face;
        require(lateral && lateral->periodicU(), "native cylinder periodic chart");
        const auto *descriptor = curved.model->surfaceDescriptor(lateral->tag());
        require(close(descriptor->radius, 2) && close(descriptor->axis.z, 1),
                "world cylinder descriptor");
        std::map<int, std::vector<geo::OcctEdgeUse>> seams;
        for (const auto &wire : curved.model->boundaryUses(lateral->tag()))
            for (const auto &use : wire.edges)
                if (use.seam)
                    seams[use.edgeTag].push_back(use);
        require(seams.size() == 1 && seams.begin()->second.size() == 2,
                "seam retains two uses of one edge identity");
        const auto &seam = seams.begin()->second;
        require(seam[0].sense == -seam[1].sense, "opposite seam use orientations");
        Vec2 first, second;
        require(
            curved.model->parameterOnUse(
                seam[0].useId, (seam[0].firstParameter + seam[0].lastParameter) * .5, first) &&
                curved.model->parameterOnUse(
                    seam[1].useId, (seam[1].firstParameter + seam[1].lastParameter) * .5, second),
            "both seam pcurves resolve");
        require(close(std::abs(first.u - second.u), 2 * std::acos(-1.)),
                "seam uses retain distinct parameter branches");
        Vec2 uv;
        Vec3 projection;
        require(lateral->projectPoint({3, 0, 2.5}, uv, projection) &&
                    (projection - Vec3{2, 0, 2.5}).norm() < 1e-6,
                "periodic seam projection");
        require(curved.model->solid(1)->isInside({0, 0, 2.5}) &&
                    !curved.model->solid(1)->isInside({3, 0, 2.5}),
                "cylinder inside/outside");
        auto halfCylinder = importBrep(BRepPrimAPI_MakeCylinder(1e6, 5., std::acos(-1.)).Shape(),
                                       directory / "half-cylinder-large.brep");
        for (const auto *face : halfCylinder.model->faces())
            if (face->kind() == geo::SurfaceKind::Cylinder) {
                Vec2 lo, hi;
                face->paramBounds(lo, hi);
                const Vec2 outsideTrim{lo.u - 5e-8, (lo.v + hi.v) * .5};
                require(!face->containsParam(outsideTrim),
                        "large-radius trim tolerance is measured in physical units");
                const auto outsidePoint = face->pointAt(outsideTrim);
                Vec2 projectedUv;
                Vec3 projectedPoint;
                require(!face->projectPoint(outsidePoint, projectedUv, projectedPoint) ||
                            (outsidePoint - projectedPoint).norm() > .049,
                        "projection cannot accept a point physically outside the trim");
            }
        auto smallCylinder = importBrep(BRepPrimAPI_MakeCylinder(1e-4, 3e-4).Shape(),
                                        directory / "cylinder-small.brep");
        inspectTopology(*smallCylinder.model);
        for (const auto *face : smallCylinder.model->faces())
            if (face->kind() == geo::SurfaceKind::Cylinder) {
                require(face->containsParam({0, 1.5e-4}), "small chart seam remains on face");
                Vec2 projectedUv;
                Vec3 projectedPoint;
                require(face->projectPoint({1.1e-4, 0, 1.5e-4}, projectedUv, projectedPoint) &&
                            (projectedPoint - Vec3{1e-4, 0, 1.5e-4}).norm() < 1e-10,
                        "small chart seam projection");
            }
        const auto hole = BRepAlgoAPI_Cut(BRepPrimAPI_MakeCylinder(2., 5.).Shape(),
                                          BRepPrimAPI_MakeCylinder(1., 5.).Shape())
                              .Shape();
        auto tube = importBrep(hole, directory / "tube.brep");
        inspectTopology(*tube.model);
        require(close(tube.model->volume(1), 15 * std::acos(-1.)), "tube exact volume");
        require(!tube.model->solid(1)->isInside({0, 0, 2.5}) &&
                    tube.model->solid(1)->isInside({1.5, 0, 2.5}),
                "hole inside classification");
        std::size_t trimmed = 0;
        for (const auto *face : tube.model->faces())
            if (face->kind() == geo::SurfaceKind::Plane) {
                require(face->boundaryLoops().size() == 2, "cap hole loops");
                Vec2 rejectedUv{41, 42};
                Vec3 rejected{51, 52, 53};
                require(!face->projectPoint({0, 0, face->pointAt({0, 0}).z}, rejectedUv, rejected),
                        "hole projection outside trim");
                require(rejectedUv.u == 41 && rejectedUv.v == 42 && rejected.x == 51 &&
                            rejected.y == 52 && rejected.z == 53,
                        "trim refusal preserves outputs");
                ++trimmed;
            }
        require(trimmed == 2, "two trimmed tube caps");
        // A blind radial hole removes the nearest stationary branch of the
        // outer cylinder while its antipodal branch remains inside the trim.
        // NURBS conversion reproduces the BSpline projection found in actual
        // CAD without depending on external corpus files or tag ordering.
        const auto blindHole = BRepAlgoAPI_Cut(
            BRepPrimAPI_MakeCylinder(1., 2.).Shape(),
            BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1.6, 0, 1), gp_Dir(-1, 0, 0)),
                                    .25, 1.).Shape()).Shape();
        auto trimmedNurbs = importBrep(BRepBuilderAPI_NurbsConvert(blindHole, true).Shape(),
                                       directory / "trimmed-cylinder-nurbs.brep");
        const geo::GFace* nurbsOuter = nullptr;
        for (const auto* face : trimmedNurbs.model->faces()) {
            Vec2 candidateUv;
            Vec3 candidate;
            if (face->projectPoint({-1.01, 0, 1}, candidateUv, candidate) &&
                (candidate - Vec3{-1, 0, 1}).norm() < 1e-7) {
                require(nurbsOuter == nullptr, "unique NURBS outer face");
                nurbsOuter = face;
            }
        }
        require(nurbsOuter != nullptr, "NURBS nearest valid branch remains projectable");
        Vec2 rejectedNurbsUv{41, 42};
        Vec3 rejectedNurbsPoint{51, 52, 53};
        require(!nurbsOuter->projectPoint({1.01, 0, 1}, rejectedNurbsUv, rejectedNurbsPoint),
                "nearest exterior stationary branch cannot select distant interior branch");
        require(rejectedNurbsUv.u == 41 && rejectedNurbsUv.v == 42 &&
                    rejectedNurbsPoint.x == 51 && rejectedNurbsPoint.y == 52 && rejectedNurbsPoint.z == 53,
                "NURBS trim refusal preserves outputs");
        auto sphere = importBrep(BRepPrimAPI_MakeSphere(2.).Shape(), directory / "sphere.brep");
        inspectTopology(*sphere.model);
        bool degenerate = false;
        for (const auto *face : sphere.model->faces())
            for (const auto &wire : sphere.model->boundaryUses(face->tag()))
                for (const auto &use : wire.edges)
                    degenerate = degenerate || use.degenerated;
        require(degenerate, "sphere preserves degenerate pole edges");
        TopoDS_Compound compound;
        BRep_Builder builder;
        builder.MakeCompound(compound);
        builder.Add(compound, box);
        gp_Trsf move;
        move.SetTranslation(gp_Vec(10, 20, 30));
        builder.Add(compound, box.Moved(TopLoc_Location(move)));
        auto disconnected = importBrep(compound, directory / "disconnected.brep");
        inspectTopology(*disconnected.model);
        require(disconnected.solids == 2 && disconnected.vertices == 16,
                "located disconnected body identities");
        require(disconnected.model->solid(1)->isInside(center) &&
                    disconnected.model->solid(2)->isInside(center + Vec3{10, 20, 30}),
                "located body classification");
        require(!disconnected.model->solid(1)->isInside(center + Vec3{10, 20, 30}),
                "no body identity alias");
        TopTools_ListOfShape arguments, tools;
        arguments.Append(box);
        tools.Append(
            BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(1, 0, 0), gp_Dir(1, 0, 0)), -10., 10., -10., 10.)
                .Shape());
        BRepAlgoAPI_Splitter split;
        split.SetArguments(arguments);
        split.SetTools(tools);
        split.Build();
        require(split.IsDone(), "split box shared-face fixture");
        auto shared = importBrep(split.Shape(), directory / "shared-face.brep");
        inspectTopology(*shared.model);
        require(shared.solids == 2 && close(shared.model->volume(1), 12) &&
                    close(shared.model->volume(2), 12),
                "shared-face body volumes");
        std::size_t interfaces = 0;
        for (const auto *face : shared.model->faces()) {
            const auto owners = shared.model->faceOwners(face->tag());
            if (owners.size() == 2) {
                require(owners[0].solidTag != owners[1].solidTag &&
                            owners[0].sense == -owners[1].sense,
                        "shared face has opposite owner senses");
                ++interfaces;
            }
        }
        require(interfaces == 1, "one exact shared face identity preserved");
        STEPControl_Writer writer;
        require(writer.Transfer(box, STEPControl_AsIs) == IFSelect_RetDone, "STEP transfer");
        const auto mmPath = directory / "box-mm.step";
        require(writer.Write(mmPath.string().c_str()) == IFSelect_RetDone, "STEP write");
        std::ifstream input(mmPath);
        std::string text{std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
        input.close();
        const std::string declaration = "SI_UNIT(.MILLI.,.METRE.)";
        const auto at = text.find(declaration);
        require(at != std::string::npos, "fixture declares mm");
        text.replace(at, declaration.size(), "SI_UNIT($,.METRE.)");
        const auto mPath = directory / "box-m.step";
        std::ofstream(mPath) << text;
        UnitGuard units;
        const auto mm = geo::importOcctCad(mmPath), metres = geo::importOcctCad(mPath);
        require(mm.imported() && metres.imported() && mm.millimetresPerModelUnit == 1 &&
                    metres.millimetresPerModelUnit == 1,
                "explicit mm transfers");
        require(close(mm.cadVolume, 24) && close(metres.cadVolume, 24e9),
                "metre declared dimensions scale to mm");
        require(mm.model->identity() != metres.model->identity() && !mm.sourceLengthUnits.empty() &&
                    !metres.sourceLengthUnits.empty(),
                "instance identity and unit receipt");
        geo::NativeCadImportOptions limit;
        limit.maximumFileBytes = 1;
        require(geo::importOcctCad(mmPath, limit).reason == geo::CadImportReason::FileBudget,
                "input file budget");
        limit.maximumFileBytes = 64000000;
        limit.maximumFaces = 1;
        require(geo::importOcctCad(mmPath, limit).reason == geo::CadImportReason::TopologyBudget,
                "topology budget");
        const auto malformed = directory / "malformed.step";
        std::ofstream(malformed) << "not a STEP document\n";
        require(geo::importOcctCad(malformed).reason == geo::CadImportReason::ParseFailed,
                "malformed STEP refusal");
        geo::NativeCadImportOptions useLimit;
        useLimit.format = geo::CadFileFormat::Brep;
        useLimit.maximumUses = 1;
        require(geo::importOcctCad(directory / "box.brep", useLimit).reason ==
                    geo::CadImportReason::TopologyBudget,
                "oriented-use budget");
        std::cout << "Native OCCT BRep adapter passed geometry, seams, pcurves, holes, reversed "
                     "faces, poles, disconnected locations, units and limits\n";
        if (ownsDirectory)
            std::filesystem::remove_all(directory);
    } catch (const std::exception &error) {
        if (ownsDirectory) {
            std::error_code ignored;
            std::filesystem::remove_all(directory, ignored);
        }
        std::cerr << error.what() << '\n';
        return 1;
    }
}
