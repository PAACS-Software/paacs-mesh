// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Native preview regression: link to PAACS::occt; src/geometry is
// private to this test for constructing an already-triangulated owned snapshot.
#include "paacs/geometry/OcctPreview.hpp"
#include "OcctModelInternal.hpp"
#include <BRep_Builder.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRep_Tool.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <Poly_Triangulation.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
struct SourceFace {
    TopoDS_Face face;
    TopLoc_Location location;
    occ::handle<Poly_Triangulation> mesh;
    std::vector<gp_Pnt> points;
    std::vector<std::array<int, 3>> triangles;
    double deflection = 0;
};
std::vector<SourceFace> snapshot(const TopoDS_Shape& shape) {
    std::vector<SourceFace> result;
    for (TopExp_Explorer explorer(shape, TopAbs_FACE); explorer.More(); explorer.Next()) {
        SourceFace entry;
        entry.face = TopoDS::Face(explorer.Current());
        entry.mesh = BRep_Tool::Triangulation(entry.face, entry.location);
        require(!entry.mesh.IsNull(), "test fixture has source display triangles");
        entry.deflection = entry.mesh->Deflection();
        for (int n = 1; n <= entry.mesh->NbNodes(); ++n) entry.points.push_back(entry.mesh->Node(n));
        for (int t = 1; t <= entry.mesh->NbTriangles(); ++t) {
            int a, b, c;
            entry.mesh->Triangle(t).Get(a, b, c);
            entry.triangles.push_back({a, b, c});
        }
        result.push_back(std::move(entry));
    }
    return result;
}
void unchanged(const std::vector<SourceFace>& source) {
    for (const auto& entry : source) {
        TopLoc_Location location;
        const auto mesh = BRep_Tool::Triangulation(entry.face, location);
        require(mesh == entry.mesh && location == entry.location && mesh->Deflection() == entry.deflection &&
                mesh->NbNodes() == static_cast<int>(entry.points.size()) &&
                mesh->NbTriangles() == static_cast<int>(entry.triangles.size()), "source mesh metadata unchanged");
        for (int n = 1; n <= mesh->NbNodes(); ++n)
            require(mesh->Node(n).IsEqual(entry.points[n - 1], 0), "source point coordinates unchanged");
        for (int t = 1; t <= mesh->NbTriangles(); ++t) {
            int a, b, c;
            mesh->Triangle(t).Get(a, b, c);
            require(std::array<int, 3>{a, b, c} == entry.triangles[t - 1], "source triangle connectivity unchanged");
        }
    }
}
void referenceWinding(const paacs::geo::OcctModel& model, const paacs::geo::OcctPreview& display) {
    for (std::size_t t = 0; t < display.triangles.size(); ++t) {
        const auto* face = model.face(display.faceTags[t]);
        require(face, "display face reference resolves");
        const auto& tri = display.triangles[t];
        const auto &a = display.points[tri[0]], &b = display.points[tri[1]], &c = display.points[tri[2]];
        const paacs::Vec3 centre{(a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3, (a.z + b.z + c.z) / 3};
        const auto uv = face->paramAt(centre);
        const auto expected = face->normalAt(uv);
        const double ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
        const double vx = c.x - a.x, vy = c.y - a.y, vz = c.z - a.z;
        require((uy * vz - uz * vy) * expected.x + (uz * vx - ux * vz) * expected.y +
                (ux * vy - uy * vx) * expected.z > 0, "display winding agrees with transformed face reference");
    }
}
}

int main() {
    using namespace paacs;
    try {
        // Two uses of the same topology at distinct locations catch accidental
        // placement-free mapping. The second placement includes a rotation.
        auto cylinder = BRepPrimAPI_MakeCylinder(1., 2.).Shape();
        BRepMesh_IncrementalMesh prepared(cylinder, .08, false, .4, false);
        require(prepared.IsDone(), "prepare own cylinder source mesh");
        gp_Trsf placed;
        placed.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), .43);
        placed.SetTranslationPart(gp_Vec(5, 0, 0));
        TopoDS_Compound assembly;
        BRep_Builder builder;
        builder.MakeCompound(assembly);
        builder.Add(assembly, cylinder);
        builder.Add(assembly, cylinder.Located(TopLoc_Location(placed)));
        const auto source = snapshot(assembly);
        auto imported = geo::detail::OcctModelAccess::fromShape(assembly);
        require(imported.imported() && imported.solids == 2 && imported.faces == 6,
                "import own placed assembly without changing body identities");
        geo::OcctPreviewOptions options;
        options.linearDeflection = .02;
        geo::OcctPreview display;
        const auto result = geo::previewOcctCad(*imported.model, options, display);
        require(result.completed && result.reason == geo::OcctPreviewReason::None, "placed CAD preview succeeds");
        unchanged(source);
        referenceWinding(*imported.model, display);
        require(display.triangles.size() == display.faceTags.size(), "display face classifications align");
        std::vector<int> tags = display.faceTags;
        std::sort(tags.begin(), tags.end());
        tags.erase(std::unique(tags.begin(), tags.end()), tags.end());
        require(tags.size() == 6, "all placed face references remain distinct");
        for (std::size_t i = 0; i < display.triangles.size(); ++i) {
            const auto tag = display.faceTags[i];
            require(imported.model->face(tag), "triangle references an imported face");
            const auto owners = imported.model->faceOwners(tag);
            require(owners.size() == 1, "placed face has its own body identity");
            for (const auto n : display.triangles[i]) {
                require(n < display.points.size(), "display connectivity in range");
                const auto& p = display.points[n];
                require(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z), "display points finite");
                require(owners.front().solidTag == 1 ? std::abs(p.x) <= 1.000001 : p.x >= 3.999999 && p.x <= 6.000001,
                        "display placement agrees with body ownership");
            }
        }
        // Failed admission must keep a previously published output intact.
        const auto pointsBefore = display.points;
        const auto trianglesBefore = display.triangles;
        const auto tagsBefore = display.faceTags;
        for (int budget = 0; budget < 3; ++budget) {
            auto limited = options;
            if (budget == 0) limited.maximumPoints = 1;
            if (budget == 1) limited.maximumTriangles = 1;
            if (budget == 2) limited.maximumOutputBytes = 1;
            const auto refused = geo::previewOcctCad(*imported.model, limited, display);
            require(!refused.completed && refused.reason == geo::OcctPreviewReason::OutputBudget, "output budget refusal");
            require(display.triangles == trianglesBefore && display.faceTags == tagsBefore &&
                    display.points.size() == pointsBefore.size(), "refusal preserves published output");
            for (std::size_t n = 0; n < display.points.size(); ++n)
                require(display.points[n].x == pointsBefore[n].x && display.points[n].y == pointsBefore[n].y &&
                        display.points[n].z == pointsBefore[n].z, "refusal preserves published coordinates");
            unchanged(source);
        }
        options.linearDeflection = 0;
        const auto invalid = geo::previewOcctCad(*imported.model, options, display);
        require(!invalid.completed && invalid.reason == geo::OcctPreviewReason::InvalidRequest, "invalid request refusal");
        unchanged(source);
        // TopLoc in the pinned kernel permits a negative transform when the
        // optional TopoDS location validation is not requested. Admission of
        // such a solid remains the native importer's responsibility. Whenever
        // admitted, preview winding must agree with GFace::normalAt, which
        // crosses transformed derivatives and already contains the sign.
        gp_Trsf reflection;
        reflection.SetMirror(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)));
        const auto mirrored = cylinder.Located(TopLoc_Location(reflection));
        const auto mirroredSource = snapshot(mirrored);
        auto reflected = geo::detail::OcctModelAccess::fromShape(mirrored);
        if (reflected.imported()) {
            geo::OcctPreview reflectedDisplay;
            options.linearDeflection = .02;
            const auto reflectedResult = geo::previewOcctCad(*reflected.model, options, reflectedDisplay);
            require(reflectedResult.completed, "admitted reflected CAD previews successfully");
            referenceWinding(*reflected.model, reflectedDisplay);
            unchanged(mirroredSource);
        }
        std::cout << "negative placement admitted: " << (reflected.imported() ? "yes" : "no") << '\n';
        std::cout << "CAD display copy, placement, source immutability and output budgets passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
