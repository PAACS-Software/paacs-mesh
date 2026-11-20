// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/algo/ProductionSurfaceMesher.hpp"
#include "paacs/algo/SurfaceMesherAdapter.hpp"

#include <span>

namespace paacs::algo::t2c_detail {
struct PrivateCurveNodeIdentity;
struct PrivateFacePoint;
struct DyadicEdgeSplit;
class SurfaceFaceCacheSession;
}

namespace paacs::algo::detail {

// Internal construction choice for the production surface candidate. Both
// routes create the same scalable exact-Delaunay triangulation. The direct-
// constraint route may initialize already-present constraint edges while the
// transient workspace is assembled; otherwise it falls back to the unchanged
// exact recovery path without changing the public outcome contract.
enum class SurfaceCandidateConstructionMode : std::uint8_t {
    ScalableDelaunay,
    ScalableDelaunayDirectConstraints,
    ScalablePlanarRecovery
};

struct CandidateBoundaryRecoveryResult {
    bool complete = false;
    std::uint64_t recoveredSegments = 0;
    bnd::SteinerReject reject = bnd::SteinerReject::None;
    bnd::SteinerRefuse refuse = bnd::SteinerRefuse::None;
    bnd::SteinerFail fail = bnd::SteinerFail::None;
};

// Additive CK-2 bridge implemented beside the already-qualified surface/BND
// composition. It mutates only the passed transient candidate.
CandidateBoundaryRecoveryResult recoverCandidateBoundaryForCk2(
    const SurfaceFaceRequest& request, SurfaceCandidate& candidate,
    const SurfaceGeometryProvider& geometry,
    bool allowDirectConstraintInitialization = false);

struct ProductionSurfaceSetStage {
    SurfaceSetOutcome surface;
    ProductionSurfaceReason reason = ProductionSurfaceReason::None;
    ProductionSurfacePhase phase = ProductionSurfacePhase::None;
    std::vector<ProductionSurfaceFaceMetrics> metrics;
};

// T2-B additive composition seam. It uses the same qualified production
// adapters as constrained-v2, but preserves the complete nested surface
// outcome and copies the caller-owned set authority into the request.
ProductionSurfaceSetStage buildProductionSurfaceSetForV3(
    const geo::Model& model, const sizing::MeshSpec& spec,
    ArtifactStamp setStamp, mesh::Mesh& destination,
    std::span<const t2c_detail::PrivateCurveNodeIdentity>
        curveNodeIdentities = {},
    std::span<const t2c_detail::PrivateFacePoint> facePoints = {},
    t2c_detail::SurfaceFaceCacheSession* cacheSession = nullptr,
    std::span<const t2c_detail::DyadicEdgeSplit>
        proposedEdgeSplits = {},
    std::span<const t2c_detail::PrivateFacePoint>
        proposedFacePoints = {},
    SurfaceCandidateConstructionMode constructionMode =
        SurfaceCandidateConstructionMode::ScalableDelaunay) noexcept;

}  // namespace paacs::algo::detail
