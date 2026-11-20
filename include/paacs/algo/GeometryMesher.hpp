// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/algo/ConvexToTetMesher.hpp"
#include "paacs/algo/MeshSizeAudit.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/sizing/CappedSizeProfile.hpp"
#include "paacs/bnd/TetEdgeGeneralCavity.hpp"
#include <atomic>
#include <cstdint>
#include <optional>
#include <vector>

namespace paacs::algo {
enum class GeometryStrategy : std::uint8_t {
    Auto,
    MappedHex8,
    LinearPrismSweep,
    BoundaryConformingTet4,
    Unsupported
};
enum class GeometryElements : std::uint8_t { Native, Tet4 };
enum class GeometryReason : std::uint8_t {
    None,
    InvalidRequest,
    InvalidTopology,
    NonfiniteGeometry,
    UnsupportedChart,
    NoMatchingStrategy,
    UnresolvedBodyInterface,
    AnalysisBudget,
    MeshingBudget,
    OutputBudget,
    SourceMeshingFailed,
    VolumeMeshingFailed,
    BoundaryAuditFailed,
    QualityRequirementNotMet,
    Cancelled,
    AllocationFailed,
    InternalError,
    ProfileRejected,
    ConversionFailed,
    SpatialFieldRejected,
    SizeRequirementNotMet
};
struct GeometryTetRefinement {
    bool enabled = false;
    double maximumEdgeToTarget = 1.5, maximumRadiusEdge = 3;
    std::uint64_t maximumInsertions = 64, maximumWork = 500'000;
};
struct GeometryBodyControl {
    int solidTag = -1;
    // Zero inherits the global target size; strategy, when set, overrides
    // the global choice for this body. Unknown/duplicate tags are rejected.
    double targetSize = 0.0;
    std::optional<GeometryStrategy> strategy;
    // Selects the starting CAD face, including reverse sweep directions.
    // -1 retains automatic recognition. Only verified blocks/sweeps qualify.
    int sourceFaceTag = -1;
    // Zero derives the layer count from targetSize. A positive value sets
    // the exact count along the recognized translation, independently of
    // the source-face resolution. Applies to mapped Hex8 and Prism6 routes.
    std::uint32_t layers = 0;
    double endToStartLayerWidthRatio = 1.0;
    double maximumAdjacentLayerRatio = 1.3;
    // Optional caller-owned physical-distance profile for a recognized
    // translation. Select sourceFaceTag to fix its origin/direction. The
    // profile length must match the translation length. Tangential sizing
    // remains targetSize. Cannot combine with layers or total-ratio grading.
    const sizing::CappedSizeProfile *sizeProfile = nullptr;
    sizing::ProfileDistributionLimits profileLimits{};
    // Native = mapped Hex8 or Prism6 (or constrained Tet4). Tet4 converts
    // convex planar-faced mapped/swept cells without adding nodes.
    GeometryElements elements = GeometryElements::Native;
    // Caller-owned immutable field; nullptr inherits the global spatialField.
    // Auto chooses constrained Tet4. Explicit mapped/sweep strategies and
    // axial-profile/layer controls are rejected rather than ignoring the field.
    const sizing::BoundedSizeField *spatialField = nullptr;
    // Optional bounded interior refinement. Protected boundaries cannot move;
    // unresolved boundary refinement returns a diagnostic, not a partial mesh.
    GeometryTetRefinement tetRefinement{};
};
struct GeometryMeshingOptions {
    double targetSize = 0.25;
    double relativeRecognitionTolerance = 1.0e-10;
    GeometryStrategy strategy = GeometryStrategy::Auto;
    int threads = 1;
    bool fullBoundaryAudit = true;
    // Zero disables the extra O(cells log cells) distribution audit. A
    // positive value requires the minimum SubTetMeanRatio in [0,1] before
    // publication; geometric validity alone does not imply solver quality.
    double minimumElementQuality = 0.0;
    std::uint64_t maximumNodes = 1000000;
    std::uint64_t maximumCells = 5000000;
    // Logical output storage only, not a process RSS limit. Source meshing,
    // retained results and audits have separate work/size bounds.
    std::uint64_t maximumOutputBytes = 512ULL * 1024 * 1024;
    std::uint64_t maximumSourcePoints = 100000;
    std::uint64_t maximumAnalysisWork = 1000000;
    // Planar source-recovery and optional tet refinement/Tet10 phase bounds,
    // plus a per-body structured construction floor of one operation per
    // output node and cell (mapped/sweep). Baseline constrained Tet4
    // construction has no work counter. This does not bound import, analysis,
    // audits, or total process work.
    std::uint64_t maximumMeshingWork = 50000000;
    std::uint64_t maximumAuditFaceUses = 12000000;
    std::uint64_t maximumAuditTemporaryBytes = 512ULL * 1024 * 1024;
    const std::atomic<bool> *cancellation = nullptr;
    std::vector<GeometryBodyControl> bodyControls;
    // Aggregate profile planner/retained array payload across bodies.
    std::uint64_t maximumProfileWorkingBytes = 16ULL * 1024 * 1024;
    // Converter temporary-array payload; excludes native/output mesh arrays.
    std::uint64_t maximumConversionWorkingBytes = 512ULL * 1024 * 1024;
    // Zero disables warning-only quality assessment. Any positive threshold
    // reports qualityWarning without rejecting; minimumElementQuality is strict.
    double warningElementQuality = 0;
    const sizing::BoundedSizeField *spatialField = nullptr;
    double warningElementSizeRatio = 0, maximumElementSizeRatio = 0;
    std::uint64_t maximumSizingAuditWork = 100'000'000;
    // Explicit bounded recovery extensions. Defaults preserve prior routing.
    bool usePlanarSurfaceRecovery = false;
    bnd::TetEdgeGeneralRecoveryPolicy edgeGeneralRecovery{};
};
struct GeometryBodyReport {
    int solidTag = -1;
    GeometryStrategy strategy = GeometryStrategy::Unsupported;
    GeometryReason reason = GeometryReason::None;
    std::uint64_t vertices = 0, edges = 0, faces = 0, holes = 0, curvedFaces = 0;
    int sourceFaceTag = -1, targetFaceTag = -1;
    Vec3 sweep{};
    double targetSize = 0.0;
    std::uint32_t layers = 0;
    double endToStartLayerWidthRatio = 1.0;
    std::uint64_t geometryHash = 0;
    std::uint64_t estimatedNodes = 0, estimatedCells = 0, estimatedOutputBytes = 0;
    bool exactCountEstimate = false;
    std::uint64_t outputNodes = 0, outputCells = 0, outputBytes = 0;
    std::uint32_t engineStatus = 0, engineReason = 0, enginePhase = 0;
    std::uint32_t recoveryKind = 0, recoveryDetail = 0;
    // Elapsed generation attempt, also recorded when this body fails/refuses.
    double generationSeconds = 0.0;
    GeometryElements elements = GeometryElements::Native;
    bool usedSizeProfile = false;
    sizing::ProfileReason profileReason = sizing::ProfileReason::None;
    double minimumLayerWidth = 0, maximumLayerWidth = 0, maximumAdjacentLayerGrowth = 1;
    double maximumLayerWidthToTarget = 0;
    std::uint64_t profileWorkingBytes = 0, plateauCellsRemoved = 0;
    ConvexToTetResult conversion{};
    // Compact detail from the general constrained-volume stage, when reached.
    // Numeric values use BoundaryConformingV3Status/Reason/Phase enums.
    bool hasBoundaryDiagnostic = false;
    std::uint32_t boundaryStatus = 0, boundaryReason = 0, boundaryPhase = 0;
    // Nested surface construction receipt, including the first stopped face.
    // Numeric values use SurfaceSetStatus/SurfaceReason. recoveryKind/Detail
    // use ProductionSurfaceRecoveryKind and the corresponding Steiner enum.
    bool hasSurfaceDiagnostic = false;
    std::uint32_t surfaceStatus = 0, surfaceReason = 0;
    std::uint64_t surfaceFailureFacePersistentId = 0;
    std::uint64_t surfaceProcessedFaces = 0, surfaceUnattemptedFaces = 0;
    bool usedSpatialField = false;
    double spatialMinimumBound = 0;
    MeshSizeAudit sizeAudit{};
    bool hasRefinementDiagnostic = false;
    std::uint32_t refinementStatus = 0, refinementReason = 0, refinementPhase = 0;
    std::uint64_t refinementInsertions = 0;
};
struct GeometryAnalysis {
    bool complete = false;
    GeometryReason reason = GeometryReason::InternalError;
    std::uint64_t work = 0;
    std::vector<GeometryBodyReport> bodies;
};
struct GeometryMeshingOutcome {
    bool published = false;
    bool destinationUnchanged = true;
    GeometryReason reason = GeometryReason::InternalError;
    GeometryAnalysis analysis;
    mesh::LinearVolumeAudit audit;
    bool qualityAssessed = false;
    bool qualityWarning = false;
    MixedMeshQualityReport quality;
    std::uint64_t outputBytes = 0;
    // Generation sums completed and failed body attempts; it does not imply publication.
    double analysisSeconds = 0.0, generationSeconds = 0.0, auditSeconds = 0.0;
    double qualitySeconds = 0.0;
};
// Recognizes verified affine blocks and straight polygonal sweeps from
// topology + backend analytic types. Otherwise proposes the existing bounded
// constrained Tet route. It does not decompose arbitrary solids, heal CAD, or
// guess a conformal mesh across shared/touching body interfaces.
GeometryAnalysis analyzeGeometryStrategies(const geo::Model &model,
                                           const GeometryMeshingOptions &options = {}) noexcept;

// Re-analyzes the immutable input for each call, privately builds each body,
// checks bounds/audits and publishes by swap. No caller-editable capability
// receipt is trusted. Separated bodies may use different strategies. Bodies
// with overlapping/touching bounding boxes require an explicit hybrid plan.
GeometryMeshingOutcome meshGeometryModel(const geo::Model &model,
                                         const GeometryMeshingOptions &options,
                                         mesh::Mesh &destination) noexcept;
const char *geometryStrategyName(GeometryStrategy strategy) noexcept;
const char *geometryReasonName(GeometryReason reason) noexcept;
} // namespace paacs::algo
