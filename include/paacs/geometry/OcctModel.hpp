// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/geometry/OcctImport.hpp"
#include <span>

namespace paacs::geo {
namespace detail {
struct OcctModelAccess;
}
enum class CadOrientation : std::uint8_t { Forward, Reversed, Internal, External };
struct OcctFaceUse {
    int faceTag = -1;
    CadOrientation orientation = CadOrientation::Forward;
    // Relative to the shared GFace normal: +1 same, -1 opposite, 0 internal/external.
    int sense = 0;
};
struct OcctFaceOwner {
    int solidTag = -1;
    int sense = 0;
};
struct OcctEdgeUse {
    std::uint64_t useId = 0;
    int edgeTag = -1, startVertexTag = -1, endVertexTag = -1;
    CadOrientation orientation = CadOrientation::Forward;
    // Relative to increasing native 3D curve parameter, including closed curves.
    int sense = 0;
    bool seam = false, degenerated = false;
    double firstParameter = 0, lastParameter = 0, tolerance = 0;
};
struct OcctWireUse {
    bool outer = false;
    std::vector<OcctEdgeUse> edges;
};
struct OcctSurfaceDescriptor {
    SurfaceKind kind = SurfaceKind::Unknown;
    Vec3 origin{}, axis{}, xDirection{};
    double radius = 0, majorRadius = 0, semiAngle = 0;
    Vec2 parameterLower{}, parameterUpper{};
    bool periodicU = false, periodicV = false;
    CadOrientation orientation = CadOrientation::Forward;
};
struct NativeCadImportOptions {
    CadFileFormat format = CadFileFormat::Step;
    std::uint64_t maximumFileBytes = 64ULL * 1024 * 1024;
    std::uint64_t maximumSolids = 10000, maximumFaces = 100000;
    std::uint64_t maximumEdges = 200000, maximumVertices = 200000;
    std::uint64_t maximumUses = 1000000;
    // Explicitly retain valid standalone document faces/edges/vertices in the
    // immutable model. Only solid-owned geometry participates in volume meshes.
    bool retainAuxiliaryGeometry = false;
};
struct NativeCadImportOutcome;

// Optional PAACS::occt native BRep snapshot. Handles/locations and all entity
// tags remain stable for its lifetime. identity() is an import-instance token;
// neither tags nor token assert correspondence across different STEP files.
// Queries use local OCCT adaptors/classifiers; no mutable OCCT object escapes.
class OcctModel final : public Model {
  public:
    ~OcctModel() override;
    std::vector<const GVertex *> vertices() const override;
    std::vector<const GEdge *> edges() const override;
    std::vector<const GFace *> faces() const override;
    std::vector<const GSolid *> solids() const override;
    const GVertex *vertex(int tag) const override;
    const GEdge *edge(int tag) const override;
    const GFace *face(int tag) const override;
    const GSolid *solid(int tag) const override;
    std::uint64_t identity() const noexcept;
    double tolerance(int dimension, int tag) const noexcept;
    double volume(int solidTag) const noexcept;
    std::span<const OcctFaceUse> faceUses(int solidTag) const noexcept;
    std::span<const OcctFaceOwner> faceOwners(int faceTag) const noexcept;
    std::span<const OcctWireUse> boundaryUses(int faceTag) const noexcept;
    const OcctSurfaceDescriptor *surfaceDescriptor(int faceTag) const noexcept;
    // Use-specific pcurve evaluation. Each affine-range result is checked
    // against its native 3D curve within source tolerances. On failure
    // outputs remain unchanged; allocation failures propagate.
    bool parameterOnUse(std::uint64_t useId, double nativeCurveParameter, Vec2 &uv) const;

  private:
    struct Impl;
    explicit OcctModel(std::unique_ptr<Impl> impl);
    std::unique_ptr<Impl> impl_;
    friend struct detail::OcctModelAccess;
    friend NativeCadImportOutcome importOcctCad(const std::filesystem::path &,
                                                const NativeCadImportOptions &) noexcept;
};
struct NativeCadImportOutcome {
    std::unique_ptr<OcctModel> model;
    CadImportReason reason = CadImportReason::KernelFailure;
    std::uint64_t fileBytes = 0, solids = 0, faces = 0, edges = 0, vertices = 0;
    std::uint64_t curvedFaces = 0, curvedEdges = 0;
    double cadVolume = 0, millimetresPerModelUnit = 0;
    std::vector<std::string> sourceLengthUnits;
    std::uint64_t auxiliaryFaces = 0, auxiliaryEdges = 0, auxiliaryVertices = 0;
    bool auxiliaryGeometryRetained = false;
    bool imported() const noexcept { return model != nullptr; }
};
// Accepts valid solid BReps, including disconnected bodies and curved/periodic
// faces. Unbounded material is refused. Extra free geometry is retained only
// with the explicit option, counted separately and excluded from volume filling.
// STEP transfers explicitly
// to mm; BREP retains unspecified coordinate units. No healing is requested.
// Bounds apply to the input/topology receipt, not OCCT parser allocations.
NativeCadImportOutcome importOcctCad(const std::filesystem::path &path,
                                     const NativeCadImportOptions &options = {}) noexcept;
} // namespace paacs::geo
