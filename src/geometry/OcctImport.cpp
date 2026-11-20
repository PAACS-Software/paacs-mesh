// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/OcctImport.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_Sequence.hxx>
#include <STEPControl_Reader.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TCollection_AsciiString.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Wire.hxx>
#include <algorithm>
#include <cmath>
#include <limits>

namespace paacs::geo {
namespace {
using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;
struct Stop {
    CadImportReason reason;
};
void require(bool good, CadImportReason reason) {
    if (!good)
        throw Stop{reason};
}
} // namespace
CadImportOutcome importPlanarCad(const std::filesystem::path &path,
                                 const CadImportOptions &options) noexcept {
    CadImportOutcome out;
    try {
        require(options.format == CadFileFormat::Step || options.format == CadFileFormat::Brep,
                CadImportReason::InvalidRequest);
        require(options.polyhedral.relativeTolerance > 0 &&
                    options.polyhedral.relativeTolerance < 1e-3 &&
                    std::isfinite(options.polyhedral.relativeTolerance),
                CadImportReason::InvalidRequest);
        std::error_code error;
        const auto size = std::filesystem::file_size(path, error);
        require(!error, CadImportReason::FileUnavailable);
        out.fileBytes = size;
        require(size <= options.maximumFileBytes, CadImportReason::FileBudget);
        const auto utf8 = path.u8string();
        const std::string filename(utf8.begin(), utf8.end());
        TopoDS_Shape shape;
        if (options.format == CadFileFormat::Step) {
            STEPControl_Reader reader;
            require(reader.ReadFile(filename.c_str()) == IFSelect_RetDone,
                    CadImportReason::ParseFailed);
            // Set the unit on this reader's loaded model, before transfer.
            // Depending on process-global xstep.cascade.unit makes dimensions
            // vary with unrelated OCCT users in a GUI host.
            reader.SetSystemLengthUnit(1.0);
            out.millimetresPerModelUnit = 1.0;
            NCollection_Sequence<TCollection_AsciiString> lengths, angles, solidAngles;
            reader.FileUnits(lengths, angles, solidAngles);
            for (int i = 1; i <= lengths.Length(); ++i)
                out.sourceLengthUnits.emplace_back(lengths.Value(i).ToCString());
            require(reader.TransferRoots() > 0, CadImportReason::ParseFailed);
            shape = reader.OneShape();
        } else {
            BRep_Builder builder;
            require(BRepTools::Read(shape, filename.c_str(), builder),
                    CadImportReason::ParseFailed);
        }
        require(!shape.IsNull(), CadImportReason::ParseFailed);
        ShapeMap solids, faces, edges, vertices;
        TopExp::MapShapes(shape, TopAbs_SOLID, solids);
        TopExp::MapShapes(shape, TopAbs_FACE, faces);
        TopExp::MapShapes(shape, TopAbs_EDGE, edges);
        TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
        out.solids = solids.Extent();
        out.faces = faces.Extent();
        out.edges = edges.Extent();
        out.vertices = vertices.Extent();
        require(out.solids == 1, CadImportReason::UnsupportedSolidCount);
        require(out.vertices <= options.polyhedral.maximumVertices &&
                    out.faces <= options.polyhedral.maximumFaces &&
                    out.edges <= options.polyhedral.maximumWork,
                CadImportReason::TopologyBudget);
        const auto solid = TopoDS::Solid(solids(1));
        ShapeMap solidFaces, solidEdges, solidVertices;
        TopExp::MapShapes(solid, TopAbs_FACE, solidFaces);
        TopExp::MapShapes(solid, TopAbs_EDGE, solidEdges);
        TopExp::MapShapes(solid, TopAbs_VERTEX, solidVertices);
        require(solidFaces.Extent() == faces.Extent() && solidEdges.Extent() == edges.Extent() &&
                    solidVertices.Extent() == vertices.Extent(),
                CadImportReason::InvalidShape);
        require(BRepCheck_Analyzer(shape).IsValid(), CadImportReason::InvalidShape);
        for (int i = 1; i <= faces.Extent(); ++i)
            if (BRepAdaptor_Surface(TopoDS::Face(faces(i))).GetType() != GeomAbs_Plane)
                ++out.curvedFaces;
        for (int i = 1; i <= edges.Extent(); ++i)
            if (BRep_Tool::Degenerated(TopoDS::Edge(edges(i))) ||
                BRepAdaptor_Curve(TopoDS::Edge(edges(i))).GetType() != GeomAbs_Line)
                ++out.curvedEdges;
        require(out.curvedFaces == 0, CadImportReason::UnsupportedSurface);
        require(out.curvedEdges == 0, CadImportReason::UnsupportedCurve);
        PolyhedralInput input;
        input.vertices.reserve(static_cast<std::size_t>(vertices.Extent()));
        for (int i = 1; i <= vertices.Extent(); ++i) {
            const auto p = BRep_Tool::Pnt(TopoDS::Vertex(vertices(i)));
            input.vertices.push_back({p.X(), p.Y(), p.Z()});
        }
        std::uint64_t uses = 0;
        for (int i = 1; i <= solidFaces.Extent(); ++i) {
            const auto face = TopoDS::Face(solidFaces(i));
            require(face.Orientation() == TopAbs_FORWARD || face.Orientation() == TopAbs_REVERSED,
                    CadImportReason::InvalidShape);
            auto forward = face;
            forward.Orientation(TopAbs_FORWARD);
            const auto outer = BRepTools::OuterWire(forward);
            require(!outer.IsNull(), CadImportReason::InvalidShape);
            std::vector<TopoDS_Wire> wires{outer};
            for (TopExp_Explorer ex(forward, TopAbs_WIRE); ex.More(); ex.Next())
                if (!ex.Current().IsSame(outer))
                    wires.push_back(TopoDS::Wire(ex.Current()));
            PolyhedralFace polygon;
            for (const auto &wire : wires) {
                std::vector<std::uint32_t> loop;
                ShapeMap wireEdges;
                TopExp::MapShapes(wire, TopAbs_EDGE, wireEdges);
                for (BRepTools_WireExplorer ex(wire, forward); ex.More(); ex.Next()) {
                    require(++uses <= options.polyhedral.maximumWork,
                            CadImportReason::TopologyBudget);
                    const auto id = vertices.FindIndex(ex.CurrentVertex());
                    require(id > 0, CadImportReason::InvalidShape);
                    loop.push_back(static_cast<std::uint32_t>(id - 1));
                }
                require(loop.size() == static_cast<std::size_t>(wireEdges.Extent()),
                        CadImportReason::InvalidShape);
                if (face.Orientation() == TopAbs_REVERSED)
                    std::reverse(loop.begin(), loop.end());
                polygon.loops.push_back(std::move(loop));
            }
            input.faces.push_back(std::move(polygon));
        }
        auto built = makePolyhedralModel(input, options.polyhedral);
        out.polyhedralReason = built.reason;
        out.polyhedralVolume = built.volume;
        require(built.built(), CadImportReason::PolyhedralValidationFailed);
        GProp_GProps properties;
        BRepGProp::VolumeProperties(solid, properties);
        out.cadVolume = std::abs(properties.Mass());
        require(std::isfinite(out.cadVolume) && out.cadVolume > 0 &&
                    std::abs(out.cadVolume - out.polyhedralVolume) <=
                        100 * options.polyhedral.relativeTolerance *
                            std::max(out.cadVolume, out.polyhedralVolume),
                CadImportReason::VolumeMismatch);
        out.model = std::move(built.model);
        out.reason = CadImportReason::None;
    } catch (const Stop &s) {
        out.reason = s.reason;
    } catch (const std::bad_alloc &) {
        out.reason = CadImportReason::AllocationFailed;
    } catch (const Standard_Failure &) {
        out.reason = CadImportReason::KernelFailure;
    } catch (...) {
        out.reason = CadImportReason::KernelFailure;
    }
    return out;
}
const char *cadImportReasonName(CadImportReason reason) noexcept {
    switch (reason) {
    case CadImportReason::None:
        return "none";
    case CadImportReason::InvalidRequest:
        return "invalid-request";
    case CadImportReason::FileUnavailable:
        return "file-unavailable";
    case CadImportReason::FileBudget:
        return "file-budget";
    case CadImportReason::ParseFailed:
        return "parse-failed";
    case CadImportReason::InvalidShape:
        return "invalid-shape";
    case CadImportReason::UnsupportedSolidCount:
        return "unsupported-solid-count";
    case CadImportReason::UnsupportedSurface:
        return "unsupported-surface";
    case CadImportReason::UnsupportedCurve:
        return "unsupported-curve";
    case CadImportReason::TopologyBudget:
        return "topology-budget";
    case CadImportReason::PolyhedralValidationFailed:
        return "polyhedral-validation-failed";
    case CadImportReason::VolumeMismatch:
        return "volume-mismatch";
    case CadImportReason::AllocationFailed:
        return "allocation-failed";
    case CadImportReason::KernelFailure:
        return "kernel-failure";
    }
    return "unknown";
}
} // namespace paacs::geo
