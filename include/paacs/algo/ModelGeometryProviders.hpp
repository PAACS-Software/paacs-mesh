// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/algo/BoundaryConformingMesher.hpp"
#include "paacs/algo/Tet10CurvingMesher.hpp"
#include "paacs/geometry/Kernel.hpp"

#include <cstdint>
#include <map>
#include <span>
#include <vector>

namespace paacs::algo {

// Geometry-kernel-independent adapter for the prepared PAACS Model.  It
// derives volume ownership from each solid's explicit boundaryFaceTags().
// Real CAD backends may supply a richer provider, but the Tet algorithms do
// not need to change when OCCT is introduced.
class ModelBoundaryOwnershipProvider final
    : public VolumeBoundaryOwnershipProvider {
public:
    enum class BuildReason : std::uint8_t {
        None,
        EmptyModel,
        MissingFace,
        InvalidFaceGeometry,
        AmbiguousFaceOrientation,
        InvalidSharedFace,
        AllocationFailed,
        InternalError
    };

    struct BuildOutcome;

    VolumeBoundaryProviderStamp stamp() const noexcept override {
        return stamp_;
    }
    VolumeFaceOwnership ownership(FacePersistentId face) const override;
    ProviderRecoveryEnvelope recoveryEnvelope(
        std::span<const Vec3> constructionCoordinates) const override;

    bool valid() const noexcept { return valid_; }
    BuildReason buildReason() const noexcept { return reason_; }
    std::size_t faceCount() const noexcept { return owners_.size(); }

private:
    friend BuildOutcome buildModelBoundaryOwnershipProvider(
        const geo::Model&, VolumeBoundaryProviderStamp) noexcept;

    VolumeBoundaryProviderStamp stamp_{};
    std::map<FacePersistentId, std::vector<bnd::VolumeRegionUse>> owners_;
    bool valid_ = false;
    BuildReason reason_ = BuildReason::InternalError;
};

struct ModelBoundaryOwnershipProvider::BuildOutcome {
    bool built = false;
    BuildReason reason = BuildReason::InternalError;
    ModelBoundaryOwnershipProvider provider{};
};

ModelBoundaryOwnershipProvider::BuildOutcome
buildModelBoundaryOwnershipProvider(
    const geo::Model& model,
    VolumeBoundaryProviderStamp stamp = {0x5041414353563031ULL, 1}) noexcept;

// Default Tet10 placement adapter over the existing abstract PAACS geometry
// model.  Curve nodes project through GEdge::closestPoint(); surface nodes use
// GFace::paramAt()/pointAt().  No OCCT type or assumption appears here.
class ModelTet10PlacementProvider final : public Tet10PlacementProvider {
public:
    explicit ModelTet10PlacementProvider(
        const geo::Model& model,
        Tet10PlacementStamp stamp = {0x5041414353543130ULL, 1}) noexcept
        : model_(&model), stamp_(stamp) {}

    Tet10PlacementStamp stamp() const override { return stamp_; }
    Tet10PlacementResult place(
        const Tet10PlacementRequest& request) const override;

private:
    const geo::Model* model_ = nullptr;
    Tet10PlacementStamp stamp_{};
};

}  // namespace paacs::algo
