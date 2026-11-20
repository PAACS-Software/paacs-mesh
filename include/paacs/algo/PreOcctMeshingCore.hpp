// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/algo/ModelGeometryProviders.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/algo/Tet10CurvingMesher.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/mesh/Quality.hpp"

#include <cstdint>
#include <memory>

namespace paacs::algo {

inline constexpr const char* kPreOcctV010RouteId =
    "pre-occt-production-core-v010";
inline constexpr std::uint32_t kPreOcctV010RouteVersion = 1;

enum class PreOcctTetOutput : std::uint8_t {
    BoundaryConformingTet4,
    RefinedTet4,
    QualityTet4,
    CertifiedTet10
};

enum class PreOcctMeshingStatus : std::uint8_t {
    Published,
    Incomplete,
    Rejected,
    Unsupported,
    Failed
};

enum class PreOcctMeshingPhase : std::uint8_t {
    Preflight,
    BoundaryConforming,
    Refinement,
    Quality,
    Tet10,
    FinalAudit,
    Publication
};

enum class PreOcctMeshingReason : std::uint16_t {
    None,
    InvalidRequest,
    MissingTet10PlacementProvider,
    BoundaryConformingStopped,
    RefinementStopped,
    QualityStopped,
    Tet10Stopped,
    FinalAuditFailed,
    AllocationFailed,
    InternalError
};

struct PreOcctResourceLimits {
    std::uint64_t maximumNodes = 1000000;
    std::uint64_t maximumTets = 5000000;
    std::uint64_t maximumWorkUnits = 50000000;
    std::uint64_t maximumReportRows = 250000;
};

// One nested request drives every accepted Tet route.  The lower stages use
// the corresponding prefix: conformingRequest, source.source, source, or the
// complete T4B request.  This avoids independent, drifting production APIs.
struct PreOcctTetRequest {
    PreOcctTetOutput output = PreOcctTetOutput::QualityTet4;
    T4BRequest pipeline{};
    // Deterministic smart-Laplacian movement of volume-classified nodes after
    // the accepted constrained/refined Tet4 route. It never moves boundary or
    // interface nodes and the final exact boundary audit remains mandatory.
    int interiorSmoothingIterations = 0;
    bool requireExactBoundaryAudit = true;
};

PreOcctTetRequest makePreOcctTetRequest(
    PreOcctTetOutput output,
    VolumeBoundaryProviderStamp ownershipStamp,
    Tet10PlacementStamp placementStamp = {},
    const PreOcctResourceLimits& limits = {}) noexcept;

// Enables the bounded accepted T4-A local cavity improvement stage.  The
// factory defaults to report-only because residual low-quality cells are a
// structured incomplete outcome rather than a reason to hide an otherwise
// valid conforming mesh.
void enableBoundedTetQuality(
    PreOcctTetRequest& request,
    double minimumRadiusRatio = 0.15,
    std::uint64_t maximumOperations = 16) noexcept;

struct PreOcctMeshingOutcome {
    PreOcctMeshingStatus status = PreOcctMeshingStatus::Failed;
    PreOcctMeshingReason reason = PreOcctMeshingReason::InternalError;
    PreOcctMeshingPhase phase = PreOcctMeshingPhase::Preflight;
    PreOcctTetOutput output = PreOcctTetOutput::QualityTet4;
    bool destinationUnchanged = true;
    std::uint64_t publishedNodes = 0;
    std::uint64_t publishedCells = 0;
    std::uint64_t strictHash = 0;
    mesh::SemanticHash semanticHash{};
    mesh::AuditResult audit{};
    mesh::QualityReport tetQualityBeforePost{};
    mesh::QualityReport tetQuality{};
    MixedMeshQualityReport mixedQuality{};

    std::unique_ptr<BoundaryConformingV3Outcome> conforming;
    std::unique_ptr<T3AOutcome> refined;
    std::unique_ptr<T4AOutcome> quality;
    std::unique_ptr<T4BOutcome> tet10;

    bool published() const noexcept {
        return status == PreOcctMeshingStatus::Published;
    }
};

// Coherent V0.10 Tet production entry.  Every lower-level engine writes into
// a private Mesh.  The final classified result is audited and then published
// through one noexcept ownership swap, so any failure leaves destination
// unchanged.
PreOcctMeshingOutcome meshPreparedTetModel(
    const geo::Model& model,
    const sizing::MeshSpec& spec,
    const PreOcctTetRequest& request,
    mesh::Mesh& destination,
    const VolumeBoundaryOwnershipProvider& ownershipProvider,
    const Tet10PlacementProvider* placementProvider = nullptr) noexcept;

}  // namespace paacs::algo
