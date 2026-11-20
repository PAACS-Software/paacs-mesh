// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/sizing/CappedSizeProfile.hpp"
#include <array>
#include <atomic>
#include <cstdint>
#include <numbers>

namespace paacs::algo {
enum class PrimitiveKind : std::uint8_t { Cylinder, Cone, Frustum, Tube, Torus, Elbow };
enum class PrimitiveElements : std::uint8_t { HexDominant, Tet4, Tet10 };
enum class PrimitiveConeStrategy : std::uint8_t { Automatic, Swept, Graded };
enum class PrimitiveCylinderStrategy : std::uint8_t { Automatic, Swept, Uniform };
enum class PrimitiveCircularStrategy : std::uint8_t { Automatic, Swept, Uniform };
enum class PrimitiveProfileStrategy : std::uint8_t { Sweep, Volume };
enum class PrimitiveQualityPolicy : std::uint8_t { Warn, Reject };
struct PrimitiveQualityIssue {
    Index cell = kInvalidIndex;
    double quality = 0;
    bool touchesSharpFeature = false;
};
struct PrimitiveQualityReport {
    std::uint64_t lowQualityCells = 0, sharpFeatureLowQualityCells = 0;
    // [0,.1), [.1,.2), ... [.9,1]. Measured per volume cell.
    std::array<std::uint64_t, 10> histogram{};
    // First affected cell ordinals, in deterministic order. Counts are complete.
    std::array<PrimitiveQualityIssue, 32> examples{};
    std::uint32_t exampleCount = 0;
    bool examplesTruncated = false;
    Index worstCell = kInvalidIndex;
    double worstQuality = 1;
};
struct PrimitiveShape {
    PrimitiveKind kind = PrimitiveKind::Cylinder;
    Vec3 origin{}, axis{0, 0, 1}, xDirection{1, 0, 0};
    double radius = 1, height = 2;
    double endRadius = .5;                    // Frustum only; Cone has a zero-radius end.
    double innerRadius = 0;                   // Tube requires >0; Torus/Elbow optionally hollow.
    double majorRadius = 2;                   // Torus/Elbow: strictly larger than radius.
    double sweepAngle = std::numbers::pi / 2; // Elbow only, in (0, 2*pi).
};
struct PrimitiveMeshingOptions {
    double targetSize = .2;
    // Absolute conservative chord-error bound. Zero selects .0025*radius.
    double maximumChordError = 0;
    double maximumRelativeVolumeError = .01;
    double minimumQuality = 0; // Minimum signed sub-tet mean ratio.
    PrimitiveElements elements = PrimitiveElements::HexDominant;
    // Graded applies to tetrahedral cones only; Automatic selects it for cones.
    PrimitiveConeStrategy coneStrategy = PrimitiveConeStrategy::Automatic;
    // Automatic uses Uniform for cylinder Tet4/Tet10; Swept retains the O-grid.
    // Uniform balances radial and axial spacing, refining both for curvature.
    PrimitiveCylinderStrategy cylinderStrategy = PrimitiveCylinderStrategy::Automatic;
    // A warning never permits nonpositive or unverified Jacobians.
    double warningQuality = .1;
    PrimitiveQualityPolicy qualityPolicy = PrimitiveQualityPolicy::Warn;
    int threads = 1;
    std::uint64_t maximumNodes = 5'000'000, maximumCells = 5'000'000;
    std::uint64_t maximumOutputBytes = 512ULL * 1024 * 1024;
    // Conservative sum of auxiliary array payloads (not process RSS).
    std::uint64_t maximumWorkingBytes = 512ULL * 1024 * 1024;
    mesh::LinearVolumeAuditLimits auditLimits{};
    const std::atomic<bool> *cancellation = nullptr;
    // Automatic balances solid torus/elbow Tet4/Tet10 sections. Hollow sweeps
    // retain the annular O-grid. Explicit Uniform requires a solid section.
    PrimitiveCircularStrategy circularStrategy = PrimitiveCircularStrategy::Automatic;
    // Optional physical axial profile for Cylinder/Tube. Keep the immutable
    // profile alive during this call. The section uses the finest target;
    // axial widths grow to the global cap. Defaults retain uniform behavior.
    const sizing::CappedSizeProfile *sizeProfile = nullptr;
    sizing::ProfileDistributionLimits profileLimits{};
    // Optional reusable preplanned distribution. Mutually exclusive with
    // sizeProfile; must span exactly 0 ... height with positive widths.
    // Caller-owned storage; targetSize controls the transverse section.
    std::span<const double> axialPositions{};
    double maximumAxialGrowth = 1.3;
    // Volume varies the interior section as well as axial spacing, for solid
    // cylinders. Tet4/Tet10 use bounded two-plane Delaunay slabs; Hex8 uses
    // a conforming four-hex subdivision of each tetrahedron.
    PrimitiveProfileStrategy profileStrategy = PrimitiveProfileStrategy::Sweep;
    std::uint64_t maximumProfileSectionNodes = 2048, maximumProfileSlabWork = 20'000'000;
    // Disable to compare the bounded slab builder against uncached generation.
    bool reuseProfileSlabs = true;
};
enum class PrimitiveReason : std::uint8_t {
    None,
    InvalidRequest,
    OutputBudget,
    WorkingBudget,
    Cancelled,
    NumericallyUnresolved,
    InvalidCell,
    BoundaryAuditFailed,
    GeometryTolerance,
    QualityRequirement,
    AllocationFailed,
    InternalError
};
struct PrimitiveMeshingResult {
    bool published = false, destinationUnchanged = true;
    PrimitiveReason reason = PrimitiveReason::InternalError;
    std::uint64_t nodes = 0, hexes = 0, pyramids = 0, boundaryQuads = 0, boundaryTriangles = 0;
    std::uint64_t tet4s = 0, tet10s = 0, certifiedTet10s = 0;
    std::uint32_t coneRefinementLevels = 0;
    bool gradedCone = false, qualityWarning = false, quadraticConnectivityVerified = false;
    bool uniformCylinder = false, cylinderGeometryLimited = false;
    // Nominal spacings for Uniform cylinders, not maximum tetrahedron edges.
    double cylinderRadialSpacing = 0, cylinderAxialSpacing = 0;
    PrimitiveQualityReport quality{};
    std::uint64_t outputBytes = 0, plannedWorkingBytes = 0;
    std::uint32_t segmentsPerQuarter = 0, radialLayers = 0, sweepLayers = 0;
    double analyticVolume = 0, meshVolume = 0, relativeVolumeError = 0;
    double chordErrorBound = 0, minimumQuality = 0, minimumCellVolume = 0;
    double generationSeconds = 0, validationSeconds = 0, auditSeconds = 0;
    mesh::LinearVolumeAudit audit{};
    bool uniformCircular = false, circularGeometryLimited = false;
    // Radial spacing and chord lengths along the innermost/outermost sweep
    // curves. These are nominal directions, not maximum tetrahedron edges.
    double circularRadialSpacing = 0, circularInnerSpacing = 0, circularOuterSpacing = 0;
    bool profiledSizing = false;
    double minimumAxialSpacing = 0, maximumAxialSpacing = 0, maximumAxialGrowth = 1;
    double transverseTarget = 0, profilePlanningSeconds = 0;
    std::uint64_t profileWorkingBytes = 0, plateauCellsRemoved = 0;
    bool profiledVolume = false, convertedHex = false;
    std::uint64_t minimumSectionNodes = 0, maximumSectionNodes = 0, profileSlabWork = 0;
    std::uint64_t profileDelaunayBytes = 0, profileDelaunaySlabs = 0, profileReusedSlabs = 0;
    double profileConversionSeconds = 0;
};
// HexDominant uses O-grid sections: Hex8 plus Pyramid5 at a cone's one apex.
// Tet4/Tet10 use a globally consistent pulling triangulation. Automatic cones
// use projected graded refinement; Automatic cylinders use concentric triangular
// sections and balanced axial spacing. Swept keeps the previous O-grid routes.
// Automatic solid torus/elbow tetrahedra use the same sections with a shared
// periodic seam for tori. Tight bends compress inner cells; quality policy
// remains authoritative. Hollow torus/elbow routes retain their annular grid.
// Tet10 has shared, straight midside nodes and Tri6 boundary faces. It does NOT
// project midsides onto CAD curves. All actual Tet10 Jacobians are checked over
// the complete reference cell with interval bounds (Bernstein fallback).
// Hex8 use the conservative full-reference Bernstein test. Invalid/unresolved
// cells always fail, regardless of quality policy. The audit member records
// linear corner incidence; quadraticConnectivityVerified additionally checks
// every actual high-order edge binding. Quality uses corner/sub-tet mean ratio.
// Warning details are bounded to 32 examples; counts/histogram cover all cells.
// Semantic classifications: volume 1; outer face 1, start cap 2, end cap 3,
// inner face 4; outer rims 101/102, inner rims 103/104; cone apex vertex 201.
// These are analytic primitive tags, not an imported CAD model's entity IDs.
// Strong guarantee: destination is exchanged only after ALL checks succeed.
// Thread pool concurrency contract is the same as core::parallelForWithThreads.
PrimitiveMeshingResult meshPrimitive(const PrimitiveShape &shape,
                                     const PrimitiveMeshingOptions &options,
                                     mesh::Mesh &destination) noexcept;
const char *primitiveKindName(PrimitiveKind kind) noexcept;
const char *primitiveReasonName(PrimitiveReason reason) noexcept;
const char *primitiveElementsName(PrimitiveElements elements) noexcept;
} // namespace paacs::algo
