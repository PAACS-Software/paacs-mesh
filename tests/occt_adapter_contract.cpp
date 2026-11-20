// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/OcctImport.hpp"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools.hxx>
#include <Interface_Static.hxx>
#include <STEPControl_Writer.hxx>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <stdexcept>

namespace {
void require(bool valid, const char *message) {
    if (!valid) throw std::runtime_error(message);
}
bool close(double actual, double expected) {
    return std::abs(actual - expected) <= 1e-8 * std::max(1.0, std::abs(expected));
}
struct UnitSetting {
    std::string previous;
    UnitSetting() : previous(Interface_Static::CVal("xstep.cascade.unit")) {
        require(Interface_Static::SetCVal("xstep.cascade.unit", "M"), "set hostile global units");
    }
    ~UnitSetting() { Interface_Static::SetCVal("xstep.cascade.unit", previous.c_str()); }
};
void inspect(const paacs::geo::CadImportOutcome &cad, double scale) {
    using namespace paacs;
    require(cad.imported(), geo::cadImportReasonName(cad.reason));
    require(cad.solids == 1 && cad.faces == 6 && cad.edges == 12 && cad.vertices == 8,
            "box topology counts");
    require(close(cad.cadVolume, 24 * scale * scale * scale), "independent CAD volume");
    require(close(cad.polyhedralVolume, cad.cadVolume), "baked volume agreement");
    const auto *body = cad.model->solids().at(0);
    Vec3 lo, hi;
    body->bbox(lo, hi);
    require(close(hi.x - lo.x, 2 * scale) && close(hi.y - lo.y, 3 * scale) &&
                close(hi.z - lo.z, 4 * scale), "dimensions and unit conversion");
    const auto center = (lo + hi) * .5;
    require(body->isInside(center) && !body->isInside(hi + Vec3{scale, scale, scale}),
            "inside/outside classification");
    std::map<int, int> owners;
    for (const auto *face : cad.model->faces()) {
        Vec2 lower, upper;
        face->paramBounds(lower, upper);
        const Vec2 uv{(lower.u + upper.u) * .5, (lower.v + upper.v) * .5};
        const auto point = face->pointAt(uv);
        const auto normal = face->normalAt(uv);
        require(close(normal.norm2(), 1), "unit face normal");
        require((point - center).dot(normal) > 0, "outward face orientation");
        Vec2 projectedUv;
        Vec3 projected;
        require(face->projectPoint(point + normal * scale, projectedUv, projected),
                "projection capability");
        require((projected - point).norm2() < 1e-14 * scale * scale,
                "face projection geometric error");
        require(face->containsParam(projectedUv), "projection lies in face trim");
        for (const auto &loop : face->boundaryLoops())
            for (int tag : loop) {
                require(cad.model->edge(tag) != nullptr, "shared edge identity resolves");
                ++owners[tag];
            }
    }
    require(owners.size() == 12, "unique shared box edges");
    for (const auto &[tag, count] : owners) {
        require(count == 2, "two adjacent faces share each edge identity");
        const auto *edge = cad.model->edge(tag);
        const auto *start = cad.model->vertex(edge->startVertexTag());
        const auto *end = cad.model->vertex(edge->endVertexTag());
        require(start && end, "edge endpoint identities");
        require((edge->pointAt(0) - start->position()).norm2() < 1e-14 * scale * scale &&
                    (edge->pointAt(1) - end->position()).norm2() < 1e-14 * scale * scale,
                "edge evaluation endpoint orientation");
    }
}
} // namespace
int main() {
    const auto directory = std::filesystem::temp_directory_path() /
        ("paacs-occt-contract-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    bool ownsDirectory = false;
    try {
        ownsDirectory = std::filesystem::create_directory(directory);
        require(ownsDirectory, "exclusive fixture directory");
        const auto shape = BRepPrimAPI_MakeBox(2., 3., 4.).Shape();
        const auto mmPath = directory / "box-mm.step";
        STEPControl_Writer writer;
        require(writer.Transfer(shape, STEPControl_AsIs) == IFSelect_RetDone, "STEP transfer");
        require(writer.Write(mmPath.string().c_str()) == IFSelect_RetDone, "STEP write");
        std::ifstream input(mmPath);
        std::string text{std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
        input.close();
        const std::string declaration = "SI_UNIT(.MILLI.,.METRE.)";
        const auto at = text.find(declaration);
        require(at != std::string::npos, "fixture declares millimetres");
        text.replace(at, declaration.size(), "SI_UNIT($,.METRE.)");
        const auto mPath = directory / "box-m.step";
        std::ofstream(mPath) << text;
        // A GUI may already have changed OCCT's global transfer setting.
        UnitSetting hostileUnits;
        auto mm = paacs::geo::importPlanarCad(mmPath);
        inspect(mm, 1);
        require(mm.millimetresPerModelUnit == 1 && !mm.sourceLengthUnits.empty(), "STEP unit receipt");
        auto metres = paacs::geo::importPlanarCad(mPath);
        inspect(metres, 1000);
        require(metres.millimetresPerModelUnit == 1 && !metres.sourceLengthUnits.empty(), "metre unit receipt");
        const auto brepPath = directory / "box.brep";
        require(BRepTools::Write(shape, brepPath.string().c_str()), "BREP write");
        paacs::geo::CadImportOptions options;
        options.format = paacs::geo::CadFileFormat::Brep;
        const auto brep = paacs::geo::importPlanarCad(brepPath, options);
        inspect(brep, 1);
        require(brep.millimetresPerModelUnit == 0 && brep.sourceLengthUnits.empty(), "BREP unit unknown");
        std::filesystem::remove_all(directory);
        std::cout << "CAD adapter dimensions, units, orientation, projection, classification and shared identities passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::error_code ignored;
        if (ownsDirectory)
            std::filesystem::remove_all(directory, ignored);
        std::cerr << e.what() << '\n';
        return 1;
    }
}
