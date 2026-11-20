// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// CK-2 constrained-v2 full-pipeline entry.
//
// legacy-v1 remains algo::meshModel()/meshFaces(). This structured route
// constructs edges, the qualified constrained surface, and the prototype
// volume entirely in a private Mesh. The caller destination is exchanged
// only after every stage succeeds.

#include <cstdint>
#include <type_traits>

#include "paacs/algo/ProductionSurfaceMesher.hpp"
#include "paacs/algo/TetCandidate.hpp"

namespace paacs::algo {

inline constexpr const char* kLegacyV1RouteId = "legacy-v1";
inline constexpr std::uint32_t kLegacyV1RouteVersion = 1;
inline constexpr const char* kConstrainedV2RouteId = "constrained-v2";
inline constexpr std::uint32_t kConstrainedV2RouteVersion = 2;

enum class ProductionModelStatus {
    Published,
    Rejected,
    Unsupported,
    Failed
};

enum class ProductionModelReason {
    None,
    InvalidSizingPolicy,
    SizingEvaluationFailed,
    EdgeConstructionFailed,
    SurfaceStopped,
    VolumeConstructionFailed,
    OrderElevationFailed,
    FinalAuditFailed,
    AllocationFailed,
    InternalError
};

enum class ProductionModelPhase {
    None,
    Preflight,
    EdgeConstruction,
    SurfaceConstruction,
    VolumeConstruction,
    OrderElevation,
    FinalAudit,
    Publication
};

struct ProductionModelAuditSummary {
    Index inverted = 0;
    Index invalidTet10 = 0;
    Index nonFinite = 0;
    Index malformed = 0;
    Index deferred = 0;
    Index prescribed = 0;
    Index matched = 0;
    Index orphan = 0;
    Index unexpected = 0;
    Index duplicate = 0;
    Index nonManifold = 0;
    int connectedComponents = 0;
};

static_assert(std::is_standard_layout_v<ProductionModelAuditSummary>);
static_assert(std::is_trivially_copyable_v<ProductionModelAuditSummary>);

struct ProductionModelOutcome {
    ProductionModelStatus status = ProductionModelStatus::Failed;
    ProductionModelReason reason = ProductionModelReason::InternalError;
    ProductionModelPhase phase = ProductionModelPhase::None;
    std::uint32_t routeVersion = kConstrainedV2RouteVersion;
    bool destinationUnchanged = true;
    Index publishedNodes = 0;
    Index publishedCells = 0;
    ProductionModelAuditSummary audit;
    ProductionSurfaceOutcome surface;
    TetHandoffMetrics tetHandoff;

    bool published() const noexcept {
        return status == ProductionModelStatus::Published;
    }
};

static_assert(std::is_nothrow_move_constructible_v<ProductionModelOutcome>);

// Success replaces destination with the complete constrained-v2 artifact.
// Every non-success leaves its complete observable contents and ordering
// unchanged. Capacity residue is impossible because publication is a
// noexcept Mesh ownership swap after a separate target allocation succeeds.
ProductionModelOutcome meshModelConstrained(
    const geo::Model& model, const sizing::MeshSpec& spec,
    mesh::Mesh& destination) noexcept;

}  // namespace paacs::algo
