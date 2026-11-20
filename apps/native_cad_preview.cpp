// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Optional OCCT display exporter. No PAACS surface/volume mesh is generated.
#include "OcctModelInternal.hpp"
#include "NumberArgument.hpp"
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <Standard_Failure.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace {
constexpr std::uint64_t kDefaultPoints = 1'000'000, kDefaultTriangles = 2'000'000;
constexpr std::uint64_t kDefaultJsonBytes = 256ULL * 1024 * 1024;
template<class T> T number(std::string_view value) {
    T result{};
    if (!paacs::cli::parseNumber(value, result))
        throw std::invalid_argument("invalid numeric argument");
    return result;
}
void require(bool value, const char *message) {
    if (!value) throw std::runtime_error(message);
}
std::string quoted(std::string_view text) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : text) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << unsigned(c) << std::dec;
        else out << c;
    }
    out << '"';
    return out.str();
}
struct OwnedTemporary {
    std::filesystem::path path;
    bool owned = false;
    ~OwnedTemporary() {
        if (owned) {
            std::error_code ignored;
            std::filesystem::remove(path,ignored);
        }
    }
};
void publish(const std::filesystem::path &destination, const std::string &json) {
    OwnedTemporary temporary;
    std::ofstream stream;
    const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
    for (unsigned attempt=0; attempt<32; ++attempt) {
        temporary.path = destination;
        temporary.path += ".preview-" + std::to_string(nonce) + "-" + std::to_string(attempt) + ".tmp";
        stream.open(temporary.path,std::ios::out|std::ios::binary|std::ios::noreplace);
        if (stream) { temporary.owned=true; break; }
        stream.clear();
    }
    require(temporary.owned,"cannot exclusively create display temporary file");
    stream.write(json.data(),static_cast<std::streamsize>(json.size()));
    stream.close();
    require(bool(stream),"display JSON write/close failed");
#ifdef _WIN32
    require(MoveFileExW(temporary.path.c_str(),destination.c_str(),MOVEFILE_REPLACE_EXISTING) != 0,
            "display JSON publication failed");
#else
    std::error_code error;
    std::filesystem::rename(temporary.path,destination,error);
    require(!error,"display JSON publication failed");
#endif
    temporary.owned=false;
}
struct FaceCopy {
    int tag = 0;
    TopoDS_Face original, copied;
    occ::handle<Poly_Triangulation> initialTriangulation;
    int initialNodes = 0, initialTriangles = 0;
};
} // namespace

int main(int argc, char **argv) {
    using namespace paacs;
    if (argc < 4) {
        std::cerr << "usage: paacs-native-cad-preview INPUT.step|brep OUTPUT.json LINEAR_DEFLECTION\n"
                     "  [--max-points N] [--max-triangles N] [--max-json-bytes N] [--max-file-bytes N]\n"
                     "Display geometry only; deflection uses imported units (STEP mm, BREP unspecified).\n";
        return 2;
    }
    try {
        const std::filesystem::path input=argv[1], output=argv[2];
        if (std::filesystem::exists(output))
            require(!std::filesystem::equivalent(input,output),"input and output must differ");
        const double deflection=number<double>(argv[3]);
        require(std::isfinite(deflection) && deflection > 0,"linear deflection must be finite and positive");
        std::uint64_t maxPoints=kDefaultPoints, maxTriangles=kDefaultTriangles, maxJsonBytes=kDefaultJsonBytes;
        geo::NativeCadImportOptions options;
        std::string extension=input.extension().string();
        for (auto &c : extension) if (c >= 'A' && c <= 'Z') c += 'a'-'A';
        options.format = extension == ".brep" ? geo::CadFileFormat::Brep : geo::CadFileFormat::Step;
        for (int i=4; i<argc; i+=2) {
            require(i+1 < argc,"missing option value");
            const std::string_view option=argv[i];
            const auto value=number<std::uint64_t>(argv[i+1]);
            require(value > 0,"display budgets must be positive");
            if (option == "--max-points") maxPoints=value;
            else if (option == "--max-triangles") maxTriangles=value;
            else if (option == "--max-json-bytes") maxJsonBytes=value;
            else if (option == "--max-file-bytes") options.maximumFileBytes=value;
            else throw std::invalid_argument("unknown option");
        }
        const auto cad=geo::importOcctCad(input,options);
        require(cad.imported(),geo::cadImportReasonName(cad.reason));
        const auto identity=cad.model->identity();
        BRepBuilderAPI_Copy copy(geo::detail::OcctModelAccess::shape(*cad.model),true,false);
        require(copy.IsDone(),"native deep copy failed");
        auto copiedRoot=copy.Shape();
        std::vector<FaceCopy> faces;
        for (const auto *entity : cad.model->faces()) {
            FaceCopy face;
            face.tag=entity->tag();
            face.original=geo::detail::OcctModelAccess::face(*cad.model,face.tag);
            face.copied=TopoDS::Face(copy.ModifiedShape(face.original));
            // ModifiedShape retains the copied TShape but can normalize its
            // orientation. Restore the source face's oriented use explicitly.
            face.copied.Orientation(face.original.Orientation());
            require(!face.copied.IsNull() && !face.copied.IsSame(face.original),"copy lost native face identity");
            TopLoc_Location location;
            face.initialTriangulation=BRep_Tool::Triangulation(face.original,location);
            if (!face.initialTriangulation.IsNull()) {
                face.initialNodes=face.initialTriangulation->NbNodes();
                face.initialTriangles=face.initialTriangulation->NbTriangles();
            }
            faces.push_back(std::move(face));
        }
        BRepMesh_IncrementalMesh mesher(copiedRoot,deflection,false,.4,false);
        require(mesher.IsDone(),"OCCT display tessellation failed");
        std::vector<std::array<double,3>> points;
        std::vector<std::array<std::uint64_t,3>> triangles;
        std::vector<int> tags;
        for (const auto &face : faces) {
            TopLoc_Location location;
            const auto triangulation=BRep_Tool::Triangulation(face.copied,location);
            require(!triangulation.IsNull() && triangulation->NbNodes() > 0 && triangulation->NbTriangles() > 0,
                    "copied face lacks display tessellation");
            const auto nodes=static_cast<std::uint64_t>(triangulation->NbNodes());
            const auto count=static_cast<std::uint64_t>(triangulation->NbTriangles());
            require(nodes <= maxPoints-points.size() && count <= maxTriangles-triangles.size(),
                    "display point/triangle output budget exceeded");
            const auto first=static_cast<std::uint64_t>(points.size());
            for (int n=1; n<=triangulation->NbNodes(); ++n) {
                const auto p=triangulation->Node(n).Transformed(location.Transformation());
                require(std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()),
                        "nonfinite display point");
                points.push_back({p.X(),p.Y(),p.Z()});
            }
            for (int t=1; t<=triangulation->NbTriangles(); ++t) {
                int a,b,c;
                triangulation->Triangle(t).Get(a,b,c);
                require(a >= 1 && b >= 1 && c >= 1 && a <= triangulation->NbNodes() &&
                        b <= triangulation->NbNodes() && c <= triangulation->NbNodes() &&
                        a != b && b != c && a != c,"invalid native display triangle");
                if (face.copied.Orientation() == TopAbs_REVERSED) std::swap(b,c);
                triangles.push_back({first+static_cast<unsigned>(a-1),first+static_cast<unsigned>(b-1),
                                     first+static_cast<unsigned>(c-1)});
                tags.push_back(face.tag);
            }
            TopLoc_Location originalLocation;
            const auto initial=BRep_Tool::Triangulation(face.original,originalLocation);
            require(initial == face.initialTriangulation &&
                    (initial.IsNull() || (initial->NbNodes() == face.initialNodes &&
                                         initial->NbTriangles() == face.initialTriangles)),
                    "authoritative source triangulation changed");
        }
        require(cad.model->identity() == identity && tags.size() == triangles.size(),
                "native model identity/output accounting changed");
        std::ostringstream json;
        json << std::setprecision(17)
             << "{\"schema\":\"paacs-native-cad-display-v1\",\"display_geometry_only\":true"
             << ",\"source_model_unmodified\":true,\"face_tags_preserved\":true"
             << ",\"face_locations_and_orientation_applied\":true"
             << ",\"linear_deflection\":" << deflection << ",\"angular_deflection_radians\":0.4"
             << ",\"native_solids\":" << cad.solids << ",\"native_faces\":" << cad.faces
             << ",\"native_curved_faces\":" << cad.curvedFaces << ",\"source_file_bytes\":" << cad.fileBytes
             << ",\"units\":" << quoted(cad.millimetresPerModelUnit == 1 ? "mm" : "unspecified BREP coordinates")
             << ",\"scope\":\"OCCT display tessellation of a deep shape copy. This is not a PAACS surface/volume mesh or a geometric-error certificate. Output budgets do not bound OCCT parser/tessellator allocations.\""
             << ",\"points\":[";
        for (std::size_t i=0; i<points.size(); ++i) {
            if (i) json << ',';
            json << '[' << points[i][0] << ',' << points[i][1] << ',' << points[i][2] << ']';
        }
        json << "],\"triangles\":[";
        for (std::size_t i=0; i<triangles.size(); ++i) {
            if (i) json << ',';
            json << '[' << triangles[i][0] << ',' << triangles[i][1] << ',' << triangles[i][2] << ']';
        }
        json << "],\"face_tags\":[";
        for (std::size_t i=0; i<tags.size(); ++i) {
            if (i) json << ',';
            json << tags[i];
        }
        json << "]}\n";
        const auto bytes=json.str();
        require(bytes.size() <= maxJsonBytes,"display JSON output byte budget exceeded");
        publish(output,bytes);
        std::cout << "{\"display_geometry_only\":true,\"published\":true,\"source_model_unmodified\":true"
                  << ",\"points\":" << points.size() << ",\"triangles\":" << triangles.size()
                  << ",\"native_faces\":" << cad.faces << ",\"json_bytes\":" << bytes.size() << "}\n";
        return 0;
    } catch (const Standard_Failure &error) {
        std::cerr << "native CAD display: " << error.GetMessageString() << '\n';
    } catch (const std::exception &error) {
        std::cerr << "native CAD display: " << error.what() << '\n';
    }
    return 1;
}
