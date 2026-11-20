// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <BRepBndLib.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <Interface_Static.hxx>
#include <NCollection_Sequence.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <Standard_Version.hxx>
#include <TCollection_AsciiString.hxx>
#include <TopExp_Explorer.hxx>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
void checkNear(double actual, double expected, const char *message) {
    require(std::isfinite(actual) && std::abs(actual - expected) < 1e-6, message);
}
} // namespace

// SDK-only check: PAACS geometry does not participate in the measurements.
int main(int argc, char **argv) {
    const auto directory = argc == 2
        ? std::filesystem::path(argv[1])
        : std::filesystem::temp_directory_path() /
              ("paacs-occt-sdk-" + std::to_string(
                  std::chrono::steady_clock::now().time_since_epoch().count()));
    bool ownsDirectory = false;
    try {
        if (argc == 2)
            std::filesystem::create_directories(directory);
        else {
            ownsDirectory = std::filesystem::create_directory(directory);
            require(ownsDirectory, "exclusive fixture directory");
        }
        const auto original = BRepPrimAPI_MakeBox(2., 3., 4.).Shape();
        require(BRepCheck_Analyzer(original).IsValid(), "invalid constructed box");
        for (const char *unit : {"MM", "M"}) {
            STEPControl_Writer writer;
            require(Interface_Static::SetCVal("write.step.unit", unit), "cannot set STEP units");
            require(writer.Transfer(original, STEPControl_AsIs) == IFSelect_RetDone,
                    "STEP transfer failed");
            const auto path = directory / (std::string("box-") + unit + ".step");
            require(writer.Write(path.string().c_str()) == IFSelect_RetDone,
                    "STEP write failed");
            STEPControl_Reader reader;
            require(reader.ReadFile(path.string().c_str()) == IFSelect_RetDone,
                    "STEP read failed");
            NCollection_Sequence<TCollection_AsciiString> lengths, angles, solidAngles;
            reader.FileUnits(lengths, angles, solidAngles);
            require(lengths.Length() == 1, "expected one declared length unit");
            const std::string declared = lengths(1).ToCString();
            require(declared == (std::string(unit) == "MM" ? "millimetre" : "metre"),
                    "incorrect STEP length-unit declaration");
            // OCCT expresses this factor in millimetres; 1 fixes the output units to mm.
            reader.SetSystemLengthUnit(1.);
            checkNear(reader.SystemLengthUnit(), 1., "incorrect transfer units");
            require(reader.TransferRoots() == 1, "expected one transferred root");
            const auto restored = reader.OneShape();
            require(!restored.IsNull() && BRepCheck_Analyzer(restored).IsValid(),
                    "invalid read-back shape");
            int solids = 0;
            for (TopExp_Explorer ex(restored, TopAbs_SOLID); ex.More(); ex.Next())
                ++solids;
            require(solids == 1, "expected one solid after read-back");
            Bnd_Box bounds;
            BRepBndLib::Add(restored, bounds, false);
            double xmin, ymin, zmin, xmax, ymax, zmax;
            bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);
            checkNear(xmax - xmin, 2., "incorrect X dimension");
            checkNear(ymax - ymin, 3., "incorrect Y dimension");
            checkNear(zmax - zmin, 4., "incorrect Z dimension");
            GProp_GProps properties;
            BRepGProp::VolumeProperties(restored, properties);
            checkNear(properties.Mass(), 24., "incorrect solid volume");
            std::cout << "{\"occt\":\"" << OCC_VERSION_COMPLETE
                      << "\",\"declared_unit\":\"" << declared
                      << "\",\"model_unit\":\"mm\",\"dimensions\":["
                      << xmax - xmin << ',' << ymax - ymin << ',' << zmax - zmin
                      << "],\"volume\":" << properties.Mass() << ",\"solids\":"
                      << solids << "}\n";
        }
        if (argc != 2)
            std::filesystem::remove_all(directory);
        return 0;
    } catch (const std::exception &error) {
        if (ownsDirectory) {
            std::error_code ignored;
            std::filesystem::remove_all(directory, ignored);
        }
        std::cerr << "OCCT SDK smoke: " << error.what() << '\n';
        return 1;
    } catch (...) {
        std::cerr << "OCCT SDK smoke: kernel exception\n";
        return 1;
    }
}
