// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "geometry_checks.hpp"
#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/geometry/OcctImport.hpp"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <STEPControl_Writer.hxx>
#include <TopoDS_Compound.hxx>
#include <chrono>
#include <fstream>
#include <gp_Ax1.hxx>
#include <gp_Trsf.hxx>
#include <iostream>

using namespace geometry_test;
namespace {
TopoDS_Shape prism(const std::vector<std::vector<Vec2>> &loops, gp_Vec direction) {
    std::unique_ptr<BRepBuilderAPI_MakeFace> face;
    for (std::size_t i = 0; i < loops.size(); ++i) {
        BRepBuilderAPI_MakePolygon polygon;
        for (const auto p : loops[i])
            polygon.Add(gp_Pnt(p.u, p.v, 0));
        polygon.Close();
        require(polygon.IsDone(), "OCCT polygon fixture");
        if (i == 0)
            face = std::make_unique<BRepBuilderAPI_MakeFace>(polygon.Wire(), true);
        else
            face->Add(polygon.Wire());
    }
    require(face && face->IsDone(), "OCCT face fixture");
    return BRepPrimAPI_MakePrism(face->Face(), direction).Shape();
}
struct Fixture {
    std::string name;
    TopoDS_Shape shape;
    double volume;
    bool supported = true;
};
} // namespace
int main(int argc, char **argv) {
    const bool retain = argc == 3 && std::string(argv[1]) == "--fixtures";
    const auto directory =
        retain ? std::filesystem::path(argv[2])
               : std::filesystem::temp_directory_path() /
                     ("paacs-cad-" +
                      std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        std::filesystem::create_directories(directory);
        std::vector<Fixture> fixtures;
        fixtures.push_back({"box", BRepPrimAPI_MakeBox(2., 1., 1.).Shape(), 2});
        fixtures.push_back({"triangle", prism({{{0, 0}, {2, 0}, {.3, 1}}}, {0, 0, 1}), 1});
        const std::vector<std::vector<Vec2>> l{{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}};
        fixtures.push_back({"concave_l", prism(l, {0, 0, 1}), 3});
        fixtures.push_back({"oblique_l", prism(l, {.2, -.1, 1}), 3});
        fixtures.push_back(
            {"through_hole",
             prism({{{0, 0}, {3, 0}, {3, 3}, {0, 3}}, {{1, 1}, {1, 2}, {2, 2}, {2, 1}}}, {0, 0, 1}),
             8});
        fixtures.push_back({"two_holes",
                            prism({{{0, 0}, {4, 0}, {4, 3}, {0, 3}},
                                   {{.5, .5}, {.5, 1.5}, {1.5, 1.5}, {1.5, .5}},
                                   {{2.5, 1.5}, {2.5, 2.5}, {3.5, 2.5}, {3.5, 1.5}}},
                                  {0, 0, 1}),
                            10});
        const auto outer = BRepPrimAPI_MakeBox(3., 3., 1.).Shape();
        const auto cut = BRepPrimAPI_MakeBox(gp_Pnt(1, 1, -1), 1., 1., 3.).Shape();
        fixtures.push_back({"boolean_cut", BRepAlgoAPI_Cut(outer, cut).Shape(), 8});
        gp_Trsf rotation;
        rotation.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), .73);
        fixtures.push_back({"rotated_box",
                            BRepBuilderAPI_Transform(fixtures[0].shape, rotation, true).Shape(),
                            2});
        fixtures.push_back({"thin_plate", BRepPrimAPI_MakeBox(2., 1., .002).Shape(), .004});
        fixtures.push_back({"beam_i",
                            prism({{{0, 0},
                                    {3, 0},
                                    {3, .3},
                                    {1.7, .3},
                                    {1.7, 2.7},
                                    {3, 2.7},
                                    {3, 3},
                                    {0, 3},
                                    {0, 2.7},
                                    {1.3, 2.7},
                                    {1.3, .3},
                                    {0, .3}}},
                                  {0, 0, 2}),
                            5.52});
        fixtures.push_back(
            {"slotted_plate",
             prism({{{0, 0}, {6, 0}, {6, 2}, {0, 2}}, {{.5, .5}, {.5, 1.5}, {5.5, 1.5}, {5.5, .5}}},
                   {0, 0, .05}),
             .35});
        fixtures.push_back(
            {"thin_duct",
             prism({{{0, 0}, {3, 0}, {3, 2}, {0, 2}}, {{.1, .1}, {.1, 1.9}, {2.9, 1.9}, {2.9, .1}}},
                   {.3, 0, 2}),
             1.92});
        std::vector<std::vector<Vec2>> grid{{{0, 0}, {4, 0}, {4, 4}, {0, 4}}};
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i) {
                const double x = .5 + i, y = .5 + j;
                grid.push_back({{x, y}, {x, y + .5}, {x + .5, y + .5}, {x + .5, y}});
            }
        fixtures.push_back({"perforated_plate", prism(grid, {0, 0, .15}), 2.0625});
        fixtures.push_back(
            {"cylinder", BRepPrimAPI_MakeCylinder(1., 2.).Shape(), 2 * std::acos(-1.), false});
        fixtures.push_back(
            {"sphere", BRepPrimAPI_MakeSphere(1.).Shape(), 4 * std::acos(-1.) / 3, false});
        TopoDS_Compound compound;
        BRep_Builder builder;
        builder.MakeCompound(compound);
        builder.Add(compound, fixtures[0].shape);
        builder.Add(compound, BRepPrimAPI_MakeBox(gp_Pnt(4, 0, 0), 1., 1., 1.).Shape());
        fixtures.push_back({"assembly", compound, 3, false});
        int count = 0, graded = 0;
        for (const auto &fixture : fixtures)
            for (const auto format : {geo::CadFileFormat::Step, geo::CadFileFormat::Brep}) {
                const auto path =
                    directory /
                    (fixture.name + (format == geo::CadFileFormat::Step ? ".step" : ".brep"));
                if (format == geo::CadFileFormat::Step) {
                    STEPControl_Writer writer;
                    require(writer.Transfer(fixture.shape, STEPControl_AsIs) == IFSelect_RetDone,
                            "STEP fixture transfer");
                    require(writer.Write(path.string().c_str()) == IFSelect_RetDone,
                            "STEP fixture write");
                } else
                    require(BRepTools::Write(fixture.shape, path.string().c_str()),
                            "BREP fixture write");
                geo::CadImportOptions io;
                io.format = format;
                auto imported = geo::importPlanarCad(path, io);
                if (fixture.supported) {
                    require(imported.imported(),
                            fixture.name + ": " + geo::cadImportReasonName(imported.reason) +
                                " / " + geo::polyhedralReasonName(imported.polyhedralReason));
                    require(std::abs(imported.cadVolume - fixture.volume) < 1e-8 * fixture.volume,
                            "CAD analytic volume");
                    algo::GeometryMeshingOptions options;
                    options.targetSize = .3;
                    mesh::Mesh m;
                    const auto result = algo::meshGeometryModel(*imported.model, options, m);
                    require(result.published,
                            fixture.name + ": " + algo::geometryReasonName(result.reason));
                    require(result.audit.valid, "CAD mesh boundary audit");
                    require(std::abs(positiveVolume(m) - fixture.volume) < 1e-8 * fixture.volume,
                            "CAD mesh analytic volume");
                    checkClassification(*imported.model, m, 1e-7);
                    std::cout << "{\"case\":\"" << path.filename().string()
                              << "\",\"published\":true,\"nodes\":" << m.nodeCount()
                              << ",\"cells\":" << m.cellCount()
                              << ",\"volume\":" << imported.cadVolume << "}\n";
                    if (fixture.name == "thin_plate" || fixture.name == "oblique_l" ||
                        fixture.name == "perforated_plate") {
                        algo::GeometryBodyControl c;
                        c.solidTag = imported.model->solids().at(0)->tag();
                        c.sourceFaceTag = result.analysis.bodies.at(0).targetFaceTag;
                        c.layers = 12;
                        c.endToStartLayerWidthRatio = .25;
                        options.bodyControls = {c};
                        mesh::Mesh layered;
                        const auto gr = algo::meshGeometryModel(*imported.model, options, layered);
                        require(gr.published && gr.audit.valid &&
                                    gr.analysis.bodies.at(0).layers == 12,
                                "CAD graded reverse sweep");
                        require(std::abs(positiveVolume(layered) - fixture.volume) <
                                    1e-8 * fixture.volume,
                                "CAD graded volume");
                        checkClassification(*imported.model, layered, 1e-7);
                        require(algo::analyzeMixedMeshQuality(layered).valid,
                                "CAD graded validity");
                        std::cout << "{\"case\":\"" << path.filename().string()
                                  << ":graded\",\"published\":true,\"layers\":12}\n";
                        ++graded;
                    }
                } else {
                    require(!imported.imported(), "unsupported CAD was accepted");
                    require(imported.reason == (fixture.name == "assembly"
                                                    ? geo::CadImportReason::UnsupportedSolidCount
                                                    : geo::CadImportReason::UnsupportedSurface),
                            "wrong CAD refusal");
                    std::cout << "{\"case\":\"" << path.filename().string()
                              << "\",\"published\":false,\"correct_refusal\":\""
                              << geo::cadImportReasonName(imported.reason) << "\"}\n";
                }
                io.maximumFileBytes = 1;
                require(geo::importPlanarCad(path, io).reason == geo::CadImportReason::FileBudget,
                        "CAD file budget");
                ++count;
            }
        const auto bad = directory / "malformed.brep";
        std::ofstream(bad) << "not a BREP";
        geo::CadImportOptions io;
        io.format = geo::CadFileFormat::Brep;
        require(!geo::importPlanarCad(bad, io).imported(), "malformed BREP accepted");
        require(geo::importPlanarCad(directory / "missing.step").reason ==
                    geo::CadImportReason::FileUnavailable,
                "missing CAD file");
        io.format = static_cast<geo::CadFileFormat>(255);
        require(geo::importPlanarCad(bad, io).reason == geo::CadImportReason::InvalidRequest,
                "invalid CAD format");
        std::cout << "CAD corpus: " << count << " file round trips + " << graded
                  << " graded cases and file/refusal checks passed\n";
        if (!retain)
            std::filesystem::remove_all(directory);
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
