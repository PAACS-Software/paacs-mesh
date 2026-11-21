// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/geometry/OcctImport.hpp"
#include <array>
#include <optional>
#include <span>
#include <string_view>

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
// Inspection categories are independent of SurfaceKind, whose existing values
// remain the mesher's analytic capability contract.
enum class OcctSurfaceType : std::uint8_t {
    Unknown, Plane, Cylinder, Cone, Sphere, Torus, Bezier, BSpline,
    Revolution, Extrusion, Offset, Other
};
enum class OcctPropertyReason : std::uint8_t { None, KernelFailure, NonfiniteResult };
struct OcctBounds {
    Vec3 lower{}, upper{};
};
struct OcctFaceProperties {
    int faceTag = -1;
    OcctSurfaceType surfaceType = OcctSurfaceType::Unknown;
    OcctPropertyReason reason = OcctPropertyReason::None;
    std::optional<OcctBounds> bounds;
    std::optional<double> area, areaRelativeErrorEstimate;
    // Area centroids need not lie on the face or inside its trimmed region.
    std::optional<Vec3> areaCentroid;
    // Constant plane normal in the existing shared GFace orientation.
    // Multiply by a nonzero faceOwners() sense for a selected body's normal.
    std::optional<Vec3> planeReferenceNormal;
};
struct OcctSolidProperties {
    int solidTag = -1;
    OcctBounds bounds;
    double volume = 0;
    std::optional<Vec3> volumeCentroid;
    OcctPropertyReason reason = OcctPropertyReason::None;
};
enum class OcctSourceNameReason : std::uint8_t {
    Disabled, None, UnavailableSource, TransferUnavailable, BudgetExceeded, KernelFailure
};
struct OcctSourceNameSummary {
    bool requested = false, complete = false;
    OcctSourceNameReason reason = OcctSourceNameReason::Disabled;
    std::uint64_t mappingsVisited = 0, namedMappings = 0, unmatchedNamedMappings = 0;
    // Vertex/edge/face/solid source names that could not be matched to final
    // placed topology. Generic shell/compound/geometry names are separate so
    // an unrelated unsupported name does not invalidate a face-name query.
    std::array<std::uint64_t, 4> unmatchedByDimension{};
    std::uint64_t unsupportedNamedMappings = 0;
    std::uint64_t candidateNames = 0, candidateBytes = 0;
    std::uint64_t namedEntities = 0, ambiguousEntities = 0;
};
const char* occtKernelVersion() noexcept;
const char* occtSurfaceTypeName(OcctSurfaceType type) noexcept;
const char* occtPropertyReasonName(OcctPropertyReason reason) noexcept;
const char* occtSourceNameReasonName(OcctSourceNameReason reason) noexcept;
struct NativeCadImportOptions {
    CadFileFormat format = CadFileFormat::Step;
    std::uint64_t maximumFileBytes = 64ULL * 1024 * 1024;
    std::uint64_t maximumSolids = 10000, maximumFaces = 100000;
    std::uint64_t maximumEdges = 200000, maximumVertices = 200000;
    std::uint64_t maximumUses = 1000000;
    // Explicitly retain valid standalone document faces/edges/vertices in the
    // immutable model. Only solid-owned geometry participates in volume meshes.
    bool retainAuxiliaryGeometry = false;
    // Optional names of directly shape-backed STEP representation items.
    // Product names, colors and assembly paths are outside this scope.
    bool retainRepresentationNames = false;
    // Bounds on visited transfer mappings and unique candidate UTF-8 bytes;
    // these do not cap parser allocations or total process memory.
    std::uint64_t maximumNameMappings = 2'000'000, maximumNamesPerEntity = 64;
    std::uint64_t maximumNameBytes = 4ULL * 1024 * 1024;
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
    // Invalid tags return nullopt. Integration uses the trimmed CAD surface,
    // not a mesh triangulation; unavailable metrics carry a reason. Allocation
    // failures propagate, as with the other native geometry queries.
    std::optional<OcctFaceProperties> faceProperties(int faceTag,
                                                   double relativeAreaTolerance = 1e-8) const;
    std::optional<OcctSolidProperties> solidProperties(int solidTag) const;
    // Candidate names are copied from exact placed transfer-shape matches.
    // Conflicts remain multiple candidates. Entity::name() is populated only
    // when extraction is complete and this entity has one unique candidate.
    // The span remains valid for the model lifetime. "complete" means that
    // this supported transfer-name scan finished, not complete XDE metadata.
    std::span<const std::string> sourceNames(int dimension, int tag) const noexcept;
    OcctSourceNameSummary sourceNameSummary() const noexcept;
    // Sorted unique face tags sharing an edge, cached during topology build.
    std::span<const int> edgeFaces(int edgeTag) const noexcept;
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
    friend NativeCadImportOutcome importOcctStepData(std::span<const char>, std::string_view,
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
// Parse an immutable caller-owned STEP snapshot. Only the synchronous call
// borrows the bytes; no source-buffer pointers escape into the returned model.
// sourceName is a nonempty UTF-8 diagnostic/context name and need not exist.
// The source byte count obeys maximumFileBytes. As for file imports, the bound
// does not cap OCCT parser allocations. Geometry admission and units match the
// path importer; this overload does not request healing or repair.
NativeCadImportOutcome importOcctStepData(std::span<const char> data,
                                         std::string_view sourceName = "memory.step",
                                         const NativeCadImportOptions &options = {}) noexcept;
} // namespace paacs::geo
