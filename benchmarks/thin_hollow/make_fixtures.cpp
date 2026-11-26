// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepTools.hxx>
#include <Bnd_Box.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <Standard_Version.hxx>
#include <TColStd_SequenceOfAsciiString.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double tolerance = 1e-7;

void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

enum class Profile { Round, Square, Rectangle, RoundedRectangle };

struct Specification {
    Profile profile;
    double width, height, wall, length;
    bool rotated = false;
    double outerCornerRadius = 0;
};

std::string profileName(Profile profile) {
    switch (profile) {
    case Profile::Round: return "round";
    case Profile::Square: return "square";
    case Profile::Rectangle: return "rectangle";
    case Profile::RoundedRectangle: return "rounded-rectangle";
    }
    throw std::runtime_error("Unknown fixture profile");
}

std::string identifier(const Specification& spec) {
    return profileName(spec.profile) + "-t" + std::to_string(int(spec.wall)) +
           "-l" + std::to_string(int(spec.length)) + (spec.rotated ? "-rotated" : "");
}

double expectedVolume(const Specification& spec) {
    if (spec.profile == Profile::Round) {
        const auto radius = spec.width / 2;
        return std::numbers::pi * (radius * radius -
            (radius - spec.wall) * (radius - spec.wall)) * spec.length;
    }
    const auto innerWidth = spec.width - 2 * spec.wall;
    const auto innerHeight = spec.height - 2 * spec.wall;
    double area = spec.width * spec.height - innerWidth * innerHeight;
    if (spec.outerCornerRadius > 0) {
        const auto innerRadius = spec.outerCornerRadius - spec.wall;
        area -= (4 - std::numbers::pi) *
            (spec.outerCornerRadius * spec.outerCornerRadius - innerRadius * innerRadius);
    }
    return area * spec.length;
}

TopoDS_Wire roundedWire(double width, double height, double radius, double z) {
    require(radius > 0 && 2 * radius < std::min(width, height), "Invalid rounded profile");
    const double x = width / 2, y = height / 2, k = radius / std::sqrt(2.);
    const std::array<gp_Pnt, 8> points{{
        {-x + radius, -y, z}, {x - radius, -y, z},
        {x, -y + radius, z}, {x, y - radius, z},
        {x - radius, y, z}, {-x + radius, y, z},
        {-x, y - radius, z}, {-x, -y + radius, z}}};
    const std::array<gp_Pnt, 4> middle{{
        {x - radius + k, -y + radius - k, z},
        {x - radius + k, y - radius + k, z},
        {-x + radius - k, y - radius + k, z},
        {-x + radius - k, -y + radius - k, z}}};
    BRepBuilderAPI_MakeWire wire;
    for (int i = 0; i < 4; ++i) {
        BRepBuilderAPI_MakeEdge straight(points[2 * i], points[2 * i + 1]);
        GC_MakeArcOfCircle curve(points[2 * i + 1], middle[i], points[(2 * i + 2) % 8]);
        require(straight.IsDone() && curve.IsDone(), "Could not construct rounded profile edge");
        BRepBuilderAPI_MakeEdge arc(curve.Value());
        require(arc.IsDone(), "Could not construct rounded profile arc");
        wire.Add(straight.Edge());
        wire.Add(arc.Edge());
    }
    require(wire.IsDone(), "Could not construct rounded profile wire");
    return wire.Wire();
}

TopoDS_Shape roundedPrism(double width, double height, double radius, double z, double length) {
    BRepBuilderAPI_MakeFace face(roundedWire(width, height, radius, z), true);
    require(face.IsDone(), "Could not construct rounded profile face");
    BRepPrimAPI_MakePrism prism(face.Face(), gp_Vec(0, 0, length));
    require(prism.IsDone(), "Could not construct rounded profile extrusion");
    return prism.Shape();
}

gp_Trsf placement(const Specification& spec) {
    gp_Trsf transform;
    if (spec.rotated) {
        gp_Trsf rotate, translate;
        rotate.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), .73);
        translate.SetTranslation(gp_Vec(11, -7, 5));
        transform = translate * rotate;
    }
    return transform;
}

TopoDS_Shape construct(const Specification& spec) {
    require(spec.wall > 0 && spec.length > 0 &&
        2 * spec.wall < std::min(spec.width, spec.height), "Invalid hollow fixture dimensions");
    // The cutter crosses both end planes. A closed cap across the bore would
    // fail the classifier probes below and the independent analytic volume.
    constexpr double extension = 1.;
    TopoDS_Shape outer, inner;
    if (spec.profile == Profile::Round) {
        outer = BRepPrimAPI_MakeCylinder(spec.width / 2, spec.length).Shape();
        inner = BRepPrimAPI_MakeCylinder(
            gp_Ax2(gp_Pnt(0, 0, -extension), gp_Dir(0, 0, 1)),
            spec.width / 2 - spec.wall, spec.length + 2 * extension).Shape();
    } else if (spec.outerCornerRadius > 0) {
        outer = roundedPrism(spec.width, spec.height, spec.outerCornerRadius, 0, spec.length);
        inner = roundedPrism(spec.width - 2 * spec.wall, spec.height - 2 * spec.wall,
            spec.outerCornerRadius - spec.wall, -extension, spec.length + 2 * extension);
    } else {
        outer = BRepPrimAPI_MakeBox(gp_Pnt(-spec.width / 2, -spec.height / 2, 0),
            spec.width, spec.height, spec.length).Shape();
        inner = BRepPrimAPI_MakeBox(
            gp_Pnt(-spec.width / 2 + spec.wall, -spec.height / 2 + spec.wall, -extension),
            spec.width - 2 * spec.wall, spec.height - 2 * spec.wall,
            spec.length + 2 * extension).Shape();
    }
    BRepAlgoAPI_Cut cut(outer, inner);
    cut.SetRunParallel(false);
    cut.Build();
    require(cut.IsDone() && !cut.HasErrors(), "Could not cut the hollow fixture");
    if (!spec.rotated) return cut.Shape();
    BRepBuilderAPI_Transform moved(cut.Shape(), placement(spec), true);
    require(moved.IsDone(), "Could not transform the hollow fixture");
    return moved.Shape();
}

struct Measurement {
    double volume;
    std::array<double, 6> bounds;
    int solids, faces, edges, vertices;
};

Measurement measure(const TopoDS_Shape& shape, const Specification& spec) {
    require(!shape.IsNull() && BRepCheck_Analyzer(shape).IsValid(), "Fixture BRep is not valid");
    TopTools_IndexedMapOfShape solids, faces, edges, vertices;
    TopExp::MapShapes(shape, TopAbs_SOLID, solids);
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    TopExp::MapShapes(shape, TopAbs_EDGE, edges);
    TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
    require(solids.Extent() == 1, "Fixture must contain exactly one solid");
    GProp_GProps properties;
    BRepGProp::VolumeProperties(shape, properties);
    const auto volume = properties.Mass(), expected = expectedVolume(spec);
    require(std::isfinite(volume) && volume > 0 &&
        std::abs(volume - expected) <= expected * 1e-9,
        "Fixture volume disagrees with the analytic hollow profile");
    const auto transform = placement(spec);
    const auto solid = TopoDS::Solid(solids(1));
    const auto state = [&](gp_Pnt point) {
        point.Transform(transform);
        return BRepClass3d_SolidClassifier(solid, point, tolerance).State();
    };
    for (const auto fraction : {0., .125, .5, .875, 1.})
        require(state(gp_Pnt(0, 0, fraction * spec.length)) == TopAbs_OUT,
            "The through-bore or an end opening is not empty");
    require(state(gp_Pnt(spec.width / 2 - spec.wall / 2, 0, spec.length / 2)) == TopAbs_IN,
        "The wall material probe is not inside the solid");
    require(state(gp_Pnt(spec.width / 2 + 1, 0, spec.length / 2)) == TopAbs_OUT,
        "The exterior probe is not outside the solid");
    Bnd_Box box;
    BRepBndLib::AddOptimal(shape, box, false, false);
    require(!box.IsVoid() && !box.IsOpen(), "Fixture bounds are unavailable");
    std::array<double, 6> bounds;
    box.Get(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
    require(std::ranges::all_of(bounds, [](double x) { return std::isfinite(x); }),
        "Fixture bounds must be finite");
    return {volume, bounds, solids.Extent(), faces.Extent(), edges.Extent(), vertices.Extent()};
}

void writeMeasurement(std::ostream& out, const Measurement& value) {
    out << "{\"volume_mm3\":" << value.volume << ",\"solid_count\":" << value.solids
        << ",\"face_count\":" << value.faces << ",\"edge_count\":" << value.edges
        << ",\"vertex_count\":" << value.vertices << ",\"bounds_mm\":[[";
    for (int i = 0; i < 3; ++i) out << (i ? "," : "") << value.bounds[i];
    out << "],[";
    for (int i = 3; i < 6; ++i) out << (i > 3 ? "," : "") << value.bounds[i];
    out << "]]}";
}

std::vector<Specification> specifications(const std::string& group) {
    std::vector<Specification> result;
    for (double length : {40., 160.}) {
        if (group == "basic" && length != 40.) continue;
        for (double wall : {1., 4.}) {
            if (group == "basic" && wall != 1.) continue;
            for (const auto profile : {Profile::Round, Profile::Square, Profile::Rectangle}) {
                Specification spec{profile, profile == Profile::Rectangle ? 60. : 40.,
                    profile == Profile::Rectangle ? 30. : 40., wall, length};
                result.push_back(spec);
                if (group == "full" && wall == 1.) {
                    spec.rotated = true;
                    result.push_back(spec);
                }
            }
        }
        if (group == "full") {
            Specification rounded{Profile::RoundedRectangle, 60., 30., 1., length, false, 4.};
            result.push_back(rounded);
            rounded.rotated = true;
            result.push_back(rounded);
        }
    }
    return result;
}

void generate(const std::filesystem::path& directory, const std::string& group) {
    if (!directory.parent_path().empty()) std::filesystem::create_directories(directory.parent_path());
    require(std::filesystem::create_directory(directory), "Output directory must not already exist");
    const auto temporary = directory / "manifest.json.tmp";
    std::ofstream manifest(temporary, std::ios::binary);
    manifest.exceptions(std::ios::badbit | std::ios::failbit);
    manifest << std::setprecision(17)
        << "{\n  \"schema\":\"paacs-owned-thin-hollow-fixtures-v1\",\n"
        << "  \"complete\":true,\n  \"units\":\"mm\",\n  \"occt_version\":\""
        << OCC_VERSION_COMPLETE << "\",\n  \"group\":\"" << group << "\",\n"
        << "  \"provenance\":\"Constructed by the PAACS Software Ltd generator in benchmarks/thin_hollow; no external CAD assets.\",\n"
        << "  \"cad_meshing_performed\":false,\n"
        << "  \"repair_performed\":false,\n  \"defeaturing_performed\":false,\n"
        << "  \"fixtures\":[\n";
    const auto cases = specifications(group);
    for (std::size_t index = 0; index < cases.size(); ++index) {
        const auto& spec = cases[index];
        const auto id = identifier(spec);
        const auto shape = construct(spec);
        const auto original = measure(shape, spec);
        const auto step = directory / (id + ".step"), brep = directory / (id + ".brep");
        STEPControl_Writer writer;
        require(Interface_Static::SetCVal("write.step.unit", "MM"), "Could not select millimetre STEP units");
        require(writer.Transfer(shape, STEPControl_AsIs) == IFSelect_RetDone &&
            writer.Write(step.string().c_str()) == IFSelect_RetDone, "Could not write STEP fixture");
        require(BRepTools::Write(shape, brep.string().c_str()), "Could not write BRep fixture");
        STEPControl_Reader reader;
        require(reader.ReadFile(step.string().c_str()) == IFSelect_RetDone, "Could not reopen STEP fixture");
        TColStd_SequenceOfAsciiString lengths, angles, solidAngles;
        reader.FileUnits(lengths, angles, solidAngles);
        require(lengths.Length() == 1 && lengths(1).IsEqual("millimetre"),
            "Reopened STEP does not declare the expected millimetre units");
        reader.SetSystemLengthUnit(1.);
        require(reader.TransferRoots() > 0, "Could not transfer reopened STEP fixture");
        const auto reopened = measure(reader.OneShape(), spec);
        require(original.solids == reopened.solids && original.faces == reopened.faces &&
            original.edges == reopened.edges && original.vertices == reopened.vertices,
            "STEP roundtrip changed fixture topology counts");
        for (std::size_t i = 0; i < original.bounds.size(); ++i)
            require(std::abs(original.bounds[i] - reopened.bounds[i]) <= 1e-6,
                "STEP roundtrip changed fixture world bounds");
        if (index) manifest << ",\n";
        manifest << "    {\"id\":\"" << id << "\",\"profile\":\"" << profileName(spec.profile)
            << "\",\"step\":\"" << step.filename().string() << "\",\"step_file\":\"" << step.filename().string() << "\",\"brep_file\":\""
            << brep.filename().string() << "\",\"step_bytes\":" << std::filesystem::file_size(step)
            << ",\"brep_bytes\":" << std::filesystem::file_size(brep)
            << ",\"outer_width_mm\":" << spec.width << ",\"outer_height_mm\":" << spec.height
            << ",\"wall_mm\":" << spec.wall << ",\"length_mm\":" << spec.length
            << ",\"outer_corner_radius_mm\":" << spec.outerCornerRadius
            << ",\"inner_corner_radius_mm\":" << (spec.outerCornerRadius > 0 ? spec.outerCornerRadius - spec.wall : 0)
            << ",\"rotated\":" << (spec.rotated ? "true" : "false")
            << ",\"analytic_volume_mm3\":" << expectedVolume(spec)
            << ",\"through_bore_probes\":5,\"material_probe_inside\":true,\"exterior_probe_outside\":true"
            << ",\"placement_matrix\":[";
        const auto transform = placement(spec);
        for (int row = 1; row <= 3; ++row) {
            if (row > 1) manifest << ',';
            manifest << '[';
            for (int col = 1; col <= 4; ++col) manifest << (col > 1 ? "," : "") << transform.Value(row, col);
            manifest << ']';
        }
        manifest << ",[0,0,0,1]],\"constructed\":";
        writeMeasurement(manifest, original);
        manifest << ",\"step_roundtrip\":";
        writeMeasurement(manifest, reopened);
        manifest << '}';
    }
    manifest << "\n  ]\n}\n";
    manifest.close();
    std::filesystem::rename(temporary, directory / "manifest.json");
    std::cout << "{\"complete\":true,\"fixture_count\":" << cases.size()
        << ",\"cad_meshing_performed\":false}\n";
}
} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2 || (argc == 4 && std::string(argv[2]) == "--group"),
            "Usage: paacs-make-thin-hollow FRESH_DIRECTORY [--group basic|core|full]");
        const std::string group = argc == 4 ? argv[3] : "full";
        require(group == "basic" || group == "core" || group == "full", "Unknown fixture group");
        generate(std::filesystem::absolute(argv[1]), group);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Fixture generation failed: " << error.what() << '\n';
        return 1;
    }
}
