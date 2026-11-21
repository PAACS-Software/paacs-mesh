// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Optional OCCT display exporter. No PAACS surface/volume mesh is generated.
#include "paacs/geometry/OcctPreview.hpp"
#include "NumberArgument.hpp"
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
        geo::OcctPreviewOptions displayOptions;
        displayOptions.linearDeflection=deflection;
        displayOptions.maximumPoints=maxPoints;
        displayOptions.maximumTriangles=maxTriangles;
        displayOptions.maximumOutputBytes=maxJsonBytes;
        geo::OcctPreview display;
        const auto displayed=geo::previewOcctCad(*cad.model,displayOptions,display);
        require(displayed.completed,geo::occtPreviewReasonName(displayed.reason));
        const auto& points=display.points;
        const auto& triangles=display.triangles;
        const auto& tags=display.faceTags;
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
            json << '[' << points[i].x << ',' << points[i].y << ',' << points[i].z << ']';
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
    } catch (const std::exception &error) {
        std::cerr << "native CAD display: " << error.what() << '\n';
    }
    return 1;
}
