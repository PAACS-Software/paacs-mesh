// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/OcctPreview.hpp"
#include "OcctModelInternal.hpp"
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <Standard_Failure.hxx>
#include <Standard_OutOfMemory.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <cmath>
#include <new>
#include <utility>

namespace paacs::geo {
namespace {
struct Failure { OcctPreviewReason reason; };
void require(bool condition, OcctPreviewReason reason = OcctPreviewReason::KernelFailure) {
    if (!condition) throw Failure{reason};
}
struct FaceCopy {
    int tag;
    TopoDS_Face original, copied;
    occ::handle<Poly_Triangulation> initial, display;
    TopLoc_Location initialLocation, displayLocation;
    int nodes = 0, triangles = 0;
    double deflection = 0;
};
// These checks are diagnostics in addition to the deep-copy construction;
// identity() alone is only an import-instance token, not a content hash.
bool sourceUnchanged(const OcctModel& model, std::uint64_t identity,
                     const std::vector<FaceCopy>& faces) noexcept {
    try {
        if (model.identity() != identity) return false;
        for (const auto& face : faces) {
            const auto current = detail::OcctModelAccess::face(model, face.tag);
            if (!current.IsEqual(face.original)) return false;
            TopLoc_Location location;
            const auto mesh = BRep_Tool::Triangulation(current, location);
            if (mesh != face.initial || location != face.initialLocation ||
                (!mesh.IsNull() && (mesh->NbNodes() != face.nodes ||
                                   mesh->NbTriangles() != face.triangles ||
                                   mesh->Deflection() != face.deflection))) return false;
        }
        return true;
    } catch (...) {
        return false;
    }
}
}

const char* occtPreviewReasonName(OcctPreviewReason reason) noexcept {
    switch (reason) {
    case OcctPreviewReason::None: return "none";
    case OcctPreviewReason::InvalidRequest: return "invalid_request";
    case OcctPreviewReason::OutputBudget: return "output_budget";
    case OcctPreviewReason::KernelFailure: return "kernel_failure";
    case OcctPreviewReason::AllocationFailed: return "allocation_failed";
    case OcctPreviewReason::SourceChanged: return "source_changed";
    }
    return "unknown";
}

OcctPreviewOutcome previewOcctCad(const OcctModel& model, const OcctPreviewOptions& options,
                                  OcctPreview& output) noexcept {
    OcctPreviewOutcome outcome;
    const auto identity = model.identity();
    std::vector<FaceCopy> faces;
    OcctPreview result;
    bool observedSource = false;
    try {
        require(std::isfinite(options.linearDeflection) && options.linearDeflection > 0 &&
                std::isfinite(options.angularDeflection) && options.angularDeflection > 0 &&
                options.angularDeflection <= 3.141592653589793 && options.maximumPoints &&
                options.maximumTriangles && options.maximumOutputBytes, OcctPreviewReason::InvalidRequest);
        // Capture the source before invoking either copy or tessellation.
        for (const auto* entity : model.faces()) {
            FaceCopy face;
            face.tag = entity->tag();
            face.original = detail::OcctModelAccess::face(model, face.tag);
            face.initial = BRep_Tool::Triangulation(face.original, face.initialLocation);
            if (!face.initial.IsNull()) {
                face.nodes = face.initial->NbNodes();
                face.triangles = face.initial->NbTriangles();
                face.deflection = face.initial->Deflection();
            }
            faces.push_back(std::move(face));
        }
        observedSource = true;
        BRepBuilderAPI_Copy copy(detail::OcctModelAccess::shape(model), true, false);
        require(copy.IsDone());
        auto root = copy.Shape();
        require(!root.IsNull() && root.TShape() != detail::OcctModelAccess::shape(model).TShape());
        for (auto& face : faces) {
            face.copied = TopoDS::Face(copy.ModifiedShape(face.original));
            require(!face.copied.IsNull() && face.copied.TShape() != face.original.TShape() &&
                    face.copied.Location() == face.original.Location());
            face.copied.Orientation(face.original.Orientation());
        }
        BRepMesh_IncrementalMesh mesher(root, options.linearDeflection, false,
                                       options.angularDeflection, false);
        require(mesher.IsDone());
        // Admit the whole output before allocating point/connectivity arrays.
        // Reserve exactly once instead of overallocating during per-face appends.
        std::uint64_t bytes = 0, totalPoints = 0, totalTriangles = 0;
        for (auto& face : faces) {
            face.display = BRep_Tool::Triangulation(face.copied, face.displayLocation);
            const auto& mesh = face.display;
            require(!mesh.IsNull() && mesh->NbNodes() > 0 && mesh->NbTriangles() > 0);
            const auto nodes = static_cast<std::uint64_t>(mesh->NbNodes());
            const auto triangles = static_cast<std::uint64_t>(mesh->NbTriangles());
            require(nodes <= options.maximumPoints - totalPoints &&
                    triangles <= options.maximumTriangles - totalTriangles, OcctPreviewReason::OutputBudget);
            constexpr auto pointBytes = sizeof(Vec3);
            constexpr auto triangleBytes = sizeof(std::array<std::uint64_t, 3>) + sizeof(int);
            require(nodes <= (options.maximumOutputBytes - bytes) / pointBytes, OcctPreviewReason::OutputBudget);
            bytes += nodes * pointBytes;
            require(triangles <= (options.maximumOutputBytes - bytes) / triangleBytes, OcctPreviewReason::OutputBudget);
            bytes += triangles * triangleBytes;
            totalPoints += nodes;
            totalTriangles += triangles;
        }
        require(totalPoints && totalTriangles && totalPoints <= result.points.max_size() &&
                totalTriangles <= result.triangles.max_size() && totalTriangles <= result.faceTags.max_size(),
                OcctPreviewReason::OutputBudget);
        result.points.reserve(static_cast<std::size_t>(totalPoints));
        result.triangles.reserve(static_cast<std::size_t>(totalTriangles));
        result.faceTags.reserve(static_cast<std::size_t>(totalTriangles));
        for (const auto& face : faces) {
            const auto& mesh = face.display;
            const auto first = static_cast<std::uint64_t>(result.points.size());
            const auto& transform = face.displayLocation.Transformation();
            for (std::int64_t n = 1; n <= mesh->NbNodes(); ++n) {
                const auto p = mesh->Node(static_cast<int>(n)).Transformed(transform);
                require(std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()));
                result.points.push_back({p.X(), p.Y(), p.Z()});
            }
            for (std::int64_t t = 1; t <= mesh->NbTriangles(); ++t) {
                int a, b, c;
                mesh->Triangle(static_cast<int>(t)).Get(a, b, c);
                require(a >= 1 && b >= 1 && c >= 1 && a <= mesh->NbNodes() && b <= mesh->NbNodes() &&
                        c <= mesh->NbNodes() && a != b && b != c && c != a);
                // Transforming the vertices already applies the determinant
                // sign to their cross product, just as GFace::normalAt crosses
                // the two transformed surface derivatives. Match that face
                // reference; do not add another reflection correction here.
                if (face.copied.Orientation() == TopAbs_REVERSED) std::swap(b, c);
                result.triangles.push_back({first + static_cast<unsigned>(a - 1), first + static_cast<unsigned>(b - 1),
                                           first + static_cast<unsigned>(c - 1)});
                result.faceTags.push_back(face.tag);
            }
        }
        require(result.faceTags.size() == result.triangles.size());
        outcome.completed = true;
        outcome.reason = OcctPreviewReason::None;
    } catch (const Failure& error) {
        outcome.reason = error.reason;
    } catch (const std::bad_alloc&) {
        outcome.reason = OcctPreviewReason::AllocationFailed;
    } catch (const Standard_OutOfMemory&) {
        outcome.reason = OcctPreviewReason::AllocationFailed;
    } catch (const Standard_Failure&) {
        outcome.reason = OcctPreviewReason::KernelFailure;
    } catch (...) {
        outcome.reason = OcctPreviewReason::KernelFailure;
    }
    // Check failure paths too: a refused preview must not prepare the source
    // model for later meshing or leak a partial display into the caller.
    if (observedSource && !sourceUnchanged(model, identity, faces)) {
        outcome.completed = false;
        outcome.reason = OcctPreviewReason::SourceChanged;
    }
    if (outcome.completed) output = std::move(result);
    return outcome;
}
} // namespace paacs::geo
