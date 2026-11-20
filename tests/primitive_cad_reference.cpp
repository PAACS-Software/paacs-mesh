// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "primitive_checks.hpp"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Compound.hxx>
#include <chrono>
#include <filesystem>
#include <gp_Ax2.hxx>
#include <iostream>

using namespace primitive_test;
namespace {
TopoDS_Shape shape(const algo::PrimitiveShape &s) {
    using K = algo::PrimitiveKind;
    const gp_Ax2 axes(gp_Pnt(s.origin.x, s.origin.y, s.origin.z),
                      gp_Dir(s.axis.x, s.axis.y, s.axis.z),
                      gp_Dir(s.xDirection.x, s.xDirection.y, s.xDirection.z));
    TopoDS_Shape outer, inner;
    if (s.kind == K::Cylinder || s.kind == K::Tube)
        outer = BRepPrimAPI_MakeCylinder(axes, s.radius, s.height).Shape();
    else if (s.kind == K::Cone || s.kind == K::Frustum)
        outer = BRepPrimAPI_MakeCone(axes, s.radius, s.kind == K::Cone ? 0 : s.endRadius, s.height)
                    .Shape();
    else if (s.kind == K::Torus)
        outer = BRepPrimAPI_MakeTorus(axes, s.majorRadius, s.radius).Shape();
    else
        outer = BRepPrimAPI_MakeTorus(axes, s.majorRadius, s.radius, s.sweepAngle).Shape();
    if (s.innerRadius > 0) {
        if (s.kind == K::Tube)
            inner = BRepPrimAPI_MakeCylinder(axes, s.innerRadius, s.height).Shape();
        else if (s.kind == K::Torus)
            inner = BRepPrimAPI_MakeTorus(axes, s.majorRadius, s.innerRadius).Shape();
        else if (s.kind == K::Elbow)
            inner = BRepPrimAPI_MakeTorus(axes, s.majorRadius, s.innerRadius, s.sweepAngle).Shape();
        if (!inner.IsNull())
            return BRepAlgoAPI_Cut(outer, inner).Shape();
    }
    return outer;
}
gp_Pnt point(Vec3 p) {
    return {p.x, p.y, p.z};
}
} // namespace
int main(int argc, char **argv) {
    const bool retain = argc == 3 && std::string(argv[1]) == "--fixtures";
    const auto directory =
        retain ? std::filesystem::path(argv[2])
               : std::filesystem::temp_directory_path() /
                     ("paacs-curved-cad-" +
                      std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        std::filesystem::create_directories(directory);
        int passed = 0;
        for (const auto &c : cases()) {
            if (c.name != algo::primitiveKindName(c.shape.kind) && c.name != "cylinder_rotated" &&
                c.name != "cone_rotated" && c.name != "torus_rotated" && c.name != "torus_hollow" &&
                c.name != "elbow_hollow" && c.name != "tube_thin_wall" && c.name != "tube_film" &&
                c.name != "torus_thin_wall" && c.name != "elbow_thin_wall")
                continue;
            const auto original = shape(c.shape);
            require(BRepCheck_Analyzer(original).IsValid(), "reference CAD invalid");
            for (auto family : {algo::PrimitiveElements::HexDominant, algo::PrimitiveElements::Tet4,
                                algo::PrimitiveElements::Tet10}) {
                mesh::Mesh m;
                algo::PrimitiveMeshingOptions o;
                o.targetSize = c.size;
                o.elements = family;
                const auto r = algo::meshPrimitive(c.shape, o, m);
                require(r.published, "explicit primitive meshing failed");
                for (bool step : {false, true}) {
                    const auto path = directory / (c.name + (step ? ".step" : ".brep"));
                    TopoDS_Shape cad;
                    if (step) {
                        STEPControl_Writer writer;
                        require(writer.Transfer(original, STEPControl_AsIs) == IFSelect_RetDone,
                                "STEP transfer");
                        require(writer.Write(path.string().c_str()) == IFSelect_RetDone,
                                "STEP write");
                        STEPControl_Reader reader;
                        require(reader.ReadFile(path.string().c_str()) == IFSelect_RetDone &&
                                    reader.TransferRoots() > 0,
                                "STEP read");
                        cad = reader.OneShape();
                    } else {
                        require(BRepTools::Write(original, path.string().c_str()), "BREP write");
                        BRep_Builder builder;
                        require(BRepTools::Read(cad, path.string().c_str(), builder), "BREP read");
                    }
                    require(BRepCheck_Analyzer(cad).IsValid(), "roundtrip CAD invalid");
                    GProp_GProps props;
                    BRepGProp::VolumeProperties(cad, props);
                    require(std::abs(std::abs(props.Mass()) - r.analyticVolume) <=
                                2e-7 * r.analyticVolume,
                            "independent CAD volume");
                    TopoDS_Compound boundary;
                    BRep_Builder builder;
                    builder.MakeCompound(boundary);
                    for (TopExp_Explorer ex(cad, TopAbs_FACE); ex.More(); ex.Next())
                        builder.Add(boundary, ex.Current());
                    std::vector<Index> boundaryNodes, boundaryCells;
                    for (Index i = 0; i < m.nodeCount(); ++i)
                        if (m.nodeClass(i).dim < 3)
                            boundaryNodes.push_back(i);
                    for (Index i = 0; i < m.cellCount(); ++i)
                        if (m.cell(i).onEntity.dim == 2)
                            boundaryCells.push_back(i);
                    double maxNodeDistance = 0, maxCenterDistance = 0;
                    for (std::size_t j = 0; j < 128; ++j) {
                        const auto node = boundaryNodes[j * boundaryNodes.size() / 128];
                        BRepExtrema_DistShapeShape distance(
                            BRepBuilderAPI_MakeVertex(point(m.node(node))).Shape(), boundary);
                        require(distance.IsDone() &&
                                    distance.Value() <= (r.tet10s ? r.chordErrorBound : 0) + 1e-6,
                                "node misses CAD boundary");
                        maxNodeDistance = std::max(maxNodeDistance, distance.Value());
                        const auto id = boundaryCells[j * boundaryCells.size() / 128];
                        Vec3 center{};
                        const auto n = mesh::nodesPerCell(m.cell(id).type);
                        for (int k = 0; k < n; ++k)
                            center = center + m.node(m.cellNodes(id)[k]) / double(n);
                        BRepExtrema_DistShapeShape centerDistance(
                            BRepBuilderAPI_MakeVertex(point(center)).Shape(), boundary);
                        require(centerDistance.IsDone() &&
                                    centerDistance.Value() <= r.chordErrorBound + 1e-6,
                                "facet exceeds CAD chord-error bound");
                        maxCenterDistance = std::max(maxCenterDistance, centerDistance.Value());
                    }
                    BRepClass3d_SolidClassifier classifier(cad);
                    for (std::size_t j = 0; j < 128; ++j) {
                        const Index i = static_cast<Index>(
                            j * (r.hexes + r.pyramids + r.tet4s + r.tet10s) / 128);
                        Vec3 center{};
                        const auto n = mesh::nodesPerCell(m.cell(i).type);
                        for (int k = 0; k < n; ++k)
                            center = center + m.node(m.cellNodes(i)[k]) / double(n);
                        classifier.Perform(point(center), 1e-7);
                        require(classifier.State() == TopAbs_IN || classifier.State() == TopAbs_ON,
                                "cell center outside CAD solid");
                    }
                    ++passed;
                    std::cout << "{\"reference\":\"" << path.filename().string()
                              << "\",\"element_family\":\"" << algo::primitiveElementsName(family)
                              << "\",\"passed\":true,\"node_distance\":" << maxNodeDistance
                              << ",\"facet_distance\":" << maxCenterDistance
                              << ",\"cad_volume\":" << std::abs(props.Mass())
                              << ",\"explicit_primitive_parameters\":true}\n";
                }
            }
        }
        if (!retain)
            std::filesystem::remove_all(directory);
        std::cout << passed
                  << " CAD reference roundtrips passed (parameters supplied explicitly; not a "
                     "curved-CAD importer)\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
