// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// CK-2 bounded production surface integration.
//
// The structured entry consumes ordered GFace boundary loops: one outer
// loop and optional holes on a regular, non-periodic chart. Embedded CAD
// curves, seams and periodic charts are deliberately not inferred.
// Candidate construction and constraint recovery are local; the destination
// changes only in SurfaceMesherAdapter's audited set-level append tail.

#include <cstddef>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <vector>

#include "paacs/algo/Delaunay2DScalable.hpp"
#include "paacs/algo/PlanarConstraintRecovery.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/sizing/SizeField.hpp"

namespace paacs::algo {

struct SurfaceSetResult;

enum class ProductionSurfaceStatus {
    Published,
    Rejected,
    Unsupported,
    Failed
};

enum class ProductionSurfaceReason {
    None,
    EmptyFaceSet,
    InvalidSizingPolicy,
    SizingEvaluationFailed,
    InvalidFaceIdentity,
    InvalidBoundaryContract,
    MissingBoundaryDiscretization,
    AmbiguousBoundaryDiscretization,
    DuplicateCandidateCoordinate,
    UnsupportedPeriodicChart,
    UnsupportedSeamChart,
    UnsupportedHoles,
    UnsupportedEmbeddedCurves,
    GeometryProviderFailure,
    CandidateConstructionFailure,
    CandidateBoundaryRecoveryStopped,
    SurfaceAdapterStopped,
    AllocationFailed,
    InternalError
};

enum class ProductionSurfacePhase {
    None,
    Preflight,
    BoundaryMapping,
    Spacing,
    ScalableConstruction,
    CandidateBoundaryRecovery,
    SurfaceBuild,
    Publication
};

struct ProductionSurfaceFaceMetrics {
    std::uint64_t facePersistentId = 0;
    std::uint64_t boundaryNodes = 0;
    std::uint64_t interiorCandidates = 0;
    std::uint64_t spacingQueries = 0;
    std::uint64_t spacingBucketsVisited = 0;
    std::uint64_t spacingDistanceTests = 0;
    std::uint64_t spacingIndexHighWater = 0;
    std::uint64_t constraintPointTests = 0;
    std::uint64_t protectedSegments = 0;
    std::uint64_t recoveredProtectedSegments = 0;
    Delaunay2DScalableMetrics construction{};
};

enum class ProductionSurfaceRecoveryKind {
    None,
    Rejected,
    Refused,
    Failed
};

struct ProductionSurfaceOutcome {
    ProductionSurfaceStatus status = ProductionSurfaceStatus::Failed;
    ProductionSurfaceReason reason = ProductionSurfaceReason::InternalError;
    ProductionSurfacePhase phase = ProductionSurfacePhase::None;
    bool destinationUnchanged = true;
    bool reportComplete = false;
    std::size_t processedFaces = 0;
    std::size_t unattemptedFaces = 0;
    ProductionSurfaceRecoveryKind recoveryKind =
        ProductionSurfaceRecoveryKind::None;
    std::uint32_t recoveryDetail = 0;
    PlanarRecoveryOutcome planarRecovery;
    std::vector<ProductionSurfaceFaceMetrics> faceMetrics;
    std::unique_ptr<SurfaceSetResult> retainedArtifact;

    ProductionSurfaceOutcome();
    ~ProductionSurfaceOutcome();
    ProductionSurfaceOutcome(ProductionSurfaceOutcome&&) noexcept;
    ProductionSurfaceOutcome& operator=(ProductionSurfaceOutcome&&) noexcept;
    ProductionSurfaceOutcome(const ProductionSurfaceOutcome&) = delete;
    ProductionSurfaceOutcome& operator=(const ProductionSurfaceOutcome&) =
        delete;

    bool published() const noexcept {
        return status == ProductionSurfaceStatus::Published;
    }
};

static_assert(std::is_nothrow_move_constructible_v<ProductionSurfaceFaceMetrics>);
static_assert(std::is_nothrow_move_constructible_v<ProductionSurfaceOutcome>);

// The destination must already contain the Line2 curve discretization (the
// normal pipeline obtains it from meshEdges()). On every non-published
// outcome, node/cell contents and ordering are unchanged; capacity residue
// from a failed publication reserve is permitted by Mesh::reserveAppend.
ProductionSurfaceOutcome meshFacesConstrained(
    const geo::Model& model, const sizing::MeshSpec& spec,
    mesh::Mesh& destination) noexcept;

struct ProductionSurfaceOptions {
    // Insert literal constraints in the dense candidate before entering the
    // existing audited publication adapter. No new boundary vertices.
    bool recoverPlanarConstraints = false;
    std::uint64_t maximumConstraintWork = 50000000;
};
ProductionSurfaceOutcome meshFacesConstrained(
    const geo::Model& model, const sizing::MeshSpec& spec,
    mesh::Mesh& destination, const ProductionSurfaceOptions& options) noexcept;

}  // namespace paacs::algo
