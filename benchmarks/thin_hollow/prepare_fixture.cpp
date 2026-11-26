// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Defeaturing.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRepTools_History.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <Interface_Static.hxx>
#include <NCollection_Sequence.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <ShapeFix_Shape.hxx>
#include <Standard_Failure.hxx>
#include <TCollection_AsciiString.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Pnt.hxx>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
constexpr std::uintmax_t maxBytes = 8 * 1024 * 1024;
void need(bool value, const char* message) { if (!value) throw std::runtime_error(message); }
std::string jsonString(const std::string& text) {
    std::ostringstream out; out << '"';
    for (unsigned char c : text) {
        if (c == '"' || c == '\\') out << '\\' << char(c);
        else if (c < 32) out << "\\u00" << std::hex << std::setw(2) << std::setfill('0') << unsigned(c);
        else out << char(c);
    }
    out << '"'; return out.str();
}
std::string read(const std::filesystem::path& path) {
    const auto bytes = std::filesystem::file_size(path);
    need(bytes <= maxBytes, "input_byte_limit");
    std::ifstream file(path, std::ios::binary); std::ostringstream out; out << file.rdbuf();
    need(file.good() || file.eof(), "input_read_failed"); need(out.str().size() == bytes, "input_changed");
    return out.str();
}
std::string brep(const TopoDS_Shape& shape) {
    std::ostringstream out; BRepTools::Write(shape, out);
    need(bool(out) && out.str().size() <= maxBytes, "brep_snapshot_limit"); return out.str();
}
void write(const std::filesystem::path& path, const std::string& text) {
    need(text.size() <= maxBytes, "output_byte_limit");
    std::ofstream out(path, std::ios::binary); out << text; out.close(); need(bool(out), "output_write_failed");
}
std::set<int> ids(const std::string& text) {
    std::set<int> result; std::istringstream input(text); std::string token;
    while (std::getline(input, token, ',')) {
        std::size_t end = 0; const auto id = std::stoi(token, &end);
        need(end == token.size() && id > 0 && result.insert(id).second, "invalid_face_selection");
        need(result.size() <= 128, "face_selection_limit");
    }
    need(!result.empty(), "empty_face_selection"); return result;
}
std::string list(const std::set<int>& values) {
    std::ostringstream out; out << '['; bool first = true;
    for (auto value : values) { out << (first ? "" : ",") << value; first = false; }
    out << ']'; return out.str();
}
double real(const std::string& text) {
    std::size_t end = 0; const auto value = std::stod(text, &end);
    need(end == text.size() && std::isfinite(value), "invalid_numeric_argument"); return value;
}
void millimetres(STEPControl_Reader& reader) {
    NCollection_Sequence<TCollection_AsciiString> lengths, angles, solidAngles;
    reader.FileUnits(lengths, angles, solidAngles);
    need(lengths.Length() == 1 && lengths.First().IsEqual("millimetre"), "millimetre_units_required");
    reader.SetSystemLengthUnit(1.);
}
struct Measure {
    double volume, maximumTolerance = 0;
    int edges = 0, vertices = 0;
    std::array<double, 6> bounds;
    TopTools_IndexedMapOfShape faces;
    void tolerance(double value) { need(std::isfinite(value) && value >= 0, "nonfinite_shape_tolerance"); maximumTolerance = std::max(maximumTolerance, value); }
    std::string json() const {
        std::ostringstream out; out << std::setprecision(17) << "{\"volume_mm3\":" << volume
            << ",\"solid_count\":1,\"face_count\":" << faces.Extent() << ",\"edge_count\":" << edges << ",\"vertex_count\":" << vertices << ",\"maximum_tolerance_mm\":" << maximumTolerance
            << ",\"bounds_mm\":[";
        for (int i = 0; i < 6; ++i) out << (i ? "," : "") << bounds[i];
        out << "]}"; return out.str();
    }
};
Measure measure(const TopoDS_Shape& shape) {
    need(!shape.IsNull() && BRepCheck_Analyzer(shape).IsValid(), "invalid_prepared_brep");
    TopTools_IndexedMapOfShape solids, shells, edges, vertices; Measure out;
    TopExp::MapShapes(shape, TopAbs_SOLID, solids); need(solids.Extent() == 1, "one_solid_required");
    TopExp::MapShapes(shape, TopAbs_SHELL, shells); need(shells.Extent() > 0, "closed_shell_required");
    for (int i = 1; i <= shells.Extent(); ++i) need(BRep_Tool::IsClosed(shells(i)), "open_shell_refused");
    TopExp::MapShapes(shape, TopAbs_FACE, out.faces); need(out.faces.Extent() <= 128, "face_limit");
    TopExp::MapShapes(shape, TopAbs_EDGE, edges); need(edges.Extent() <= 1024, "edge_limit");
    TopExp::MapShapes(shape, TopAbs_VERTEX, vertices); need(vertices.Extent() <= 1024, "vertex_limit");
    out.edges = edges.Extent(); out.vertices = vertices.Extent();
    for (int i = 1; i <= out.faces.Extent(); ++i) out.tolerance(BRep_Tool::Tolerance(TopoDS::Face(out.faces(i))));
    for (int i = 1; i <= edges.Extent(); ++i) out.tolerance(BRep_Tool::Tolerance(TopoDS::Edge(edges(i))));
    for (int i = 1; i <= vertices.Extent(); ++i) out.tolerance(BRep_Tool::Tolerance(TopoDS::Vertex(vertices(i))));
    GProp_GProps props; BRepGProp::VolumeProperties(shape, props); out.volume = props.Mass();
    need(std::isfinite(out.volume) && out.volume > 0 && std::isfinite(out.maximumTolerance), "invalid_measurement");
    Bnd_Box box; BRepBndLib::AddOptimal(shape, box, false, false); need(!box.IsVoid() && !box.IsOpen(), "invalid_bounds");
    box.Get(out.bounds[0], out.bounds[1], out.bounds[2], out.bounds[3], out.bounds[4], out.bounds[5]);
    need(std::all_of(out.bounds.begin(), out.bounds.end(), [](double value) { return std::isfinite(value); }), "nonfinite_bounds");
    const auto solid = TopoDS::Solid(solids(1));
    for (double z : {0., 5., 20., 35., 40.})
        need(BRepClass3d_SolidClassifier(solid, gp_Pnt(0, 0, z), 1e-7).State() == TopAbs_OUT, "protected_bore_not_empty");
    need(BRepClass3d_SolidClassifier(solid, gp_Pnt(29.5, 0, 20), 1e-7).State() == TopAbs_IN, "wall_probe_failed");
    return out;
}
std::set<int> mapped(const NCollection_List<TopoDS_Shape>& shapes, const Measure& prepared) {
    std::set<int> result;
    for (NCollection_List<TopoDS_Shape>::Iterator it(shapes); it.More(); it.Next()) {
        TopTools_IndexedMapOfShape faces; TopExp::MapShapes(it.Value(), TopAbs_FACE, faces);
        for (int i = 1; i <= faces.Extent(); ++i) {
            const auto id = prepared.faces.FindIndex(faces(i)); need(id > 0, "history_result_face_missing"); result.insert(id);
        }
    }
    return result;
}
} // namespace

int main(int argc, char** argv) {
    std::filesystem::path input, output; std::string operation = "none", fixture = "rounded-rectangle-t1-l40";
    std::set<int> remove, protect; double maxTolerance = 0, maxSeconds = 30;
    TopoDS_Shape original; std::string sourceBytes, originalSnapshot, stage = "arguments"; bool loaded = false, created = false;
    const auto started = Clock::now();
    const auto elapsed = [&] { return std::chrono::duration<double>(Clock::now() - started).count(); };
    const auto check = [&] { need(elapsed() <= maxSeconds, "cooperative_time_limit"); };
    const auto guard = [&] { return loaded && read(input) == sourceBytes && brep(original) == originalSnapshot; };
    const auto refuse = [&](const std::string& reason) {
        bool unchanged = false; try { unchanged = guard(); } catch (...) {}
        try { if (created) write(output / "report.json", "{\"schema\":\"paacs-preparation-prototype-v1\",\"status\":\"refused\",\"published\":false,\"operation\":" +
            jsonString(operation) + ",\"stage\":" + jsonString(stage) + ",\"reason\":" + jsonString(reason) + ",\"original_guard_available\":" + (loaded ? "true" : "false") +
            ",\"original_unchanged\":" + (unchanged ? "true" : "false") + "}\n"); } catch (...) { std::cerr << "Refusal report could not be saved.\n"; }
        std::cerr << "Preparation refused: " << reason << '\n'; return 2;
    };
    try {
        for (int i = 1; i < argc; i += 2) {
            need(i + 1 < argc, "option_value_required"); const std::string key = argv[i], value = argv[i + 1];
            if (key == "--input") input = value; else if (key == "--output") output = value;
            else if (key == "--operation") operation = value; else if (key == "--fixture-id") fixture = value;
            else if (key == "--remove-faces") remove = ids(value); else if (key == "--protect-faces") protect = ids(value);
            else if (key == "--max-tolerance") maxTolerance = real(value);
            else if (key == "--max-seconds") maxSeconds = real(value); else throw std::runtime_error("unknown_option");
        }
        need(!input.empty() && !output.empty() && fixture == "rounded-rectangle-t1-l40", "owned_rounded_rectangle_fixture_required");
        need(operation == "none" || operation == "defeature" || operation == "repair", "unknown_operation");
        need(std::isfinite(maxSeconds) && maxSeconds > 0 && maxSeconds <= 300, "invalid_time_limit");
        need(maxTolerance >= 0 && maxTolerance <= .01, "invalid_tolerance");
        if (!output.parent_path().empty()) std::filesystem::create_directories(output.parent_path());
        need(std::filesystem::create_directory(output), "fresh_output_required"); created = true;
        stage = "source_import"; sourceBytes = read(input); STEPControl_Reader reader;
        need(reader.ReadFile(input.string().c_str()) == IFSelect_RetDone, "step_import_failed");
        millimetres(reader); need(reader.TransferRoots() > 0, "step_transfer_failed"); original = reader.OneShape();
        stage = "source_validation"; const auto before = measure(original); originalSnapshot = brep(original); loaded = true; check();
        const double expected = (176 - (4 - std::numbers::pi) * 7) * 40;
        need(std::abs(before.volume - expected) <= expected * 1e-9 && before.faces.Extent() == 18, "owned_fixture_geometry_mismatch");
        const std::array<double, 6> expectedBounds{-30, -15, 0, 30, 15, 40};
        for (int i = 0; i < 6; ++i) need(std::isfinite(before.bounds[i]) && std::abs(before.bounds[i] - expectedBounds[i]) <= 1e-6, "owned_fixture_bounds_mismatch");
        for (int id : remove) need(id <= before.faces.Extent(), "remove_face_missing");
        for (int id : protect) need(id <= before.faces.Extent(), "protected_face_missing");
        for (int id : remove) need(!protect.contains(id), "protected_face_selected_for_removal");
        need(operation == "defeature" || remove.empty(), "remove_selection_requires_defeaturing");
        need(operation != "defeature" || !remove.empty(), "defeaturing_requires_explicit_faces");
        need(operation != "repair" || maxTolerance > 0, "repair_requires_max_tolerance");
        need(operation == "repair" || maxTolerance == 0, "tolerance_control_requires_repair");
        need(operation != "repair" || protect.empty(), "repair_history_unavailable_for_protected_faces");
        stage = "owned_copy"; BRepBuilderAPI_Copy copy(original, true, false); need(copy.IsDone(), "owned_copy_failed");
        TopTools_IndexedMapOfShape copiedFaces; TopExp::MapShapes(copy.Shape(), TopAbs_FACE, copiedFaces);
        need(copiedFaces.Extent() == before.faces.Extent(), "copy_face_count_changed");
        std::vector<TopoDS_Shape> lineage(before.faces.Extent() + 1); std::set<int> copyIds;
        for (int id = 1; id <= before.faces.Extent(); ++id) {
            lineage[id] = copy.ModifiedShape(before.faces(id)); const auto copiedId = copiedFaces.FindIndex(lineage[id]);
            need(!lineage[id].IsNull() && copiedId > 0 && copyIds.insert(copiedId).second, "copy_face_lineage_unproved");
        }
        TopoDS_Shape prepared = copy.Shape(); const auto copiedSnapshot = brep(prepared); BRepAlgoAPI_Defeaturing removal;
        bool repairPerformed = false; stage = "explicit_operation"; check();
        if (operation == "defeature") {
            removal.SetShape(prepared); removal.SetRunParallel(false); removal.SetToFillHistory(true);
            for (int id : remove) {
                need(BRepAdaptor_Surface(TopoDS::Face(lineage[id])).GetType() == GeomAbs_Cylinder, "only_selected_cylindrical_fillets_supported");
                removal.AddFaceToRemove(lineage[id]);
            }
            removal.Build(); check(); need(!removal.HasErrors() && !removal.Shape().IsNull() && !removal.History().IsNull(), "defeaturing_failed");
            prepared = removal.Shape();
        } else if (operation == "repair") {
            need(before.maximumTolerance <= maxTolerance, "input_exceeds_tolerance_cap");
            ShapeFix_Shape repair(prepared); repair.SetPrecision(std::min(1e-7, maxTolerance));
            repair.SetMinTolerance(std::min(1e-7, maxTolerance)); repair.SetMaxTolerance(maxTolerance);
            repairPerformed = repair.Perform(); check(); prepared = repair.Shape();
        }
        stage = "prepared_validation"; const auto after = measure(prepared); check();
        need(operation != "repair" || after.maximumTolerance <= maxTolerance, "repair_exceeds_tolerance_cap");
        if (operation == "none") need(after.faces.Extent() == before.faces.Extent() && after.edges == before.edges && after.vertices == before.vertices && std::abs(after.volume - before.volume) <= before.volume * 1e-12, "unchanged_mode_geometry_mismatch");
        if (operation == "none") for (int i = 0; i < 6; ++i) need(std::abs(after.bounds[i] - before.bounds[i]) <= 1e-9, "unchanged_mode_bounds_mismatch");
        std::ostringstream history, catalog; std::set<int> attributed;
        history << '['; catalog << '[';
        for (int id = 1; id <= before.faces.Extent(); ++id) {
            const auto& copiedFace = lineage[id]; const auto retained = after.faces.FindIndex(copiedFace);
            std::set<int> modified, generated, targets; bool deleted = false;
            if (operation == "defeature") { modified = mapped(removal.Modified(copiedFace), after); generated = mapped(removal.Generated(copiedFace), after); deleted = removal.IsDeleted(copiedFace); }
            if (operation != "repair") { if (retained) targets.insert(retained); targets.insert(modified.begin(), modified.end()); targets.insert(generated.begin(), generated.end()); attributed.insert(targets.begin(), targets.end()); }
            if (remove.contains(id)) need(retained == 0 && (deleted || !modified.empty() || !generated.empty()), "requested_feature_not_removed");
            if (protect.contains(id)) need(retained > 0 && after.faces(retained).IsEqual(copiedFace) && !deleted && modified.empty() && generated.empty() && targets.size() == 1, "protected_face_changed_or_ambiguous");
            if (id > 1) { history << ','; catalog << ','; }
            history << "{\"original_face\":" << id << ",\"copy_face\":" << copiedFaces.FindIndex(copiedFace)
                << ",\"original_orientation\":" << int(before.faces(id).Orientation()) << ",\"copy_orientation\":" << int(copiedFace.Orientation())
                << ",\"history_available\":" << (operation == "repair" ? "false" : "true")
                << ",\"retained_prepared_face\":" << (operation != "repair" ? retained : 0)
                << ",\"modified_prepared_faces\":" << list(modified) << ",\"generated_prepared_faces\":" << list(generated)
                << ",\"deleted\":" << (operation == "repair" ? "null" : deleted ? "true" : "false")
                << ",\"prepared_faces\":" << list(targets) << '}';
            BRepAdaptor_Surface surface(TopoDS::Face(before.faces(id))); const bool cylinder = surface.GetType() == GeomAbs_Cylinder;
            catalog << "{\"original_face\":" << id << ",\"surface\":\"" << (cylinder ? "cylinder" : "other")
                << "\",\"radius_mm\":"; if (cylinder) { const auto radius = surface.Cylinder().Radius(); need(std::isfinite(radius) && radius > 0, "nonfinite_catalog_radius"); catalog << std::setprecision(17) << radius; } else catalog << "null"; catalog << '}';
        }
        history << ']'; catalog << ']'; std::set<int> unattributed;
        for (int id = 1; id <= after.faces.Extent(); ++id) if (!attributed.contains(id)) unattributed.insert(id);
        stage = "prepared_step_roundtrip"; const auto staging = output / ".prepared-staging"; need(std::filesystem::create_directory(staging), "staging_directory_failed");
        const auto staged = staging / "prepared.step"; STEPControl_Writer writer;
        need(Interface_Static::SetCVal("write.step.unit", "MM") && writer.Transfer(prepared, STEPControl_AsIs) == IFSelect_RetDone &&
            writer.Write(staged.string().c_str()) == IFSelect_RetDone, "prepared_step_write_failed");
        need(std::filesystem::file_size(staged) <= maxBytes, "prepared_step_byte_limit");
        STEPControl_Reader reopened; need(reopened.ReadFile(staged.string().c_str()) == IFSelect_RetDone, "prepared_step_reopen_failed");
        millimetres(reopened); need(reopened.TransferRoots() > 0, "prepared_step_transfer_failed");
        const auto roundtrip = measure(reopened.OneShape()); check();
        need(roundtrip.faces.Extent() == after.faces.Extent() && std::abs(roundtrip.volume - after.volume) <= after.volume * 1e-9, "prepared_step_roundtrip_mismatch");
        for (int i = 0; i < 6; ++i) need(std::abs(roundtrip.bounds[i] - after.bounds[i]) <= 1e-6, "prepared_step_bounds_mismatch");
        need(guard(), "original_geometry_changed");
        std::ostringstream report; report << std::setprecision(17)
            << "{\"schema\":\"paacs-preparation-prototype-v1\",\"status\":\"succeeded\",\"published\":true,\"prototype_only\":true,\"units\":\"mm\",\"fixture_id\":" << jsonString(fixture) << ",\"operation\":" << jsonString(operation)
            << ",\"input_valid_required\":true,\"original_unchanged\":true,\"body_identity\":\"one owned solid; no assembly or public CAD-tag rebinding\",\"face_ids_are_model_local\":true"
            << ",\"prepared_step_file\":\"prepared.step\",\"atomic_result_directory\":true"
            << ",\"history_scope\":\"copy lineage then provider history; no nearest-face matching\",\"history_available\":" << (operation == "repair" ? "false" : "true")
            << ",\"history_target\":\"prepared BRep before STEP export\",\"step_reimport_face_mapping_available\":false"
            << ",\"remove_faces\":" << list(remove) << ",\"protected_faces\":" << list(protect)
            << ",\"max_tolerance_mm\":" << maxTolerance << ",\"repair_perform_return\":" << (repairPerformed ? "true" : "false")
            << ",\"prepared_differs_from_owned_copy\":" << (brep(prepared) != copiedSnapshot ? "true" : "false")
            << ",\"volume_change_mm3\":" << after.volume - before.volume << ",\"before\":" << before.json()
            << ",\"after\":" << after.json() << ",\"roundtrip\":" << roundtrip.json()
            << ",\"original_face_catalog\":" << catalog.str() << ",\"face_history\":" << history.str()
            << ",\"unattributed_prepared_faces\":" << list(unattributed)
            << ",\"max_seconds\":" << maxSeconds << ",\"seconds\":" << elapsed()
            << ",\"cooperative_deadline_only\":true,\"geometric_hausdorff_deviation_measured\":false}\n";
        stage = "publish"; write(staging / "report.json", report.str()); check(); need(guard(), "original_geometry_changed"); check();
        std::filesystem::rename(staging, output / "result");
        return 0;
    } catch (const Standard_Failure& error) { return refuse(error.GetMessageString() ? error.GetMessageString() : "occt_failure"); }
    catch (const std::exception& error) { return refuse(error.what()); }
}
