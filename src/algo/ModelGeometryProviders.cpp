// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/ModelGeometryProviders.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>

namespace paacs::algo {
namespace {

FacePersistentId persistentFaceId(int tag) noexcept {
    return static_cast<FacePersistentId>(
               static_cast<std::uint32_t>(tag)) +
           1U;
}

bool finite(const Vec2& p) noexcept {
    return std::isfinite(p.u) && std::isfinite(p.v);
}

bool finite(const Vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) &&
           std::isfinite(p.z);
}

int faceSenseForSolid(const geo::GFace& face,
                      const geo::GSolid& solid) noexcept {
    Vec2 lo{}, hi{};
    face.paramBounds(lo, hi);
    if (!finite(lo) || !finite(hi)) return 0;
    const Vec2 uv{(lo.u + hi.u) * 0.5, (lo.v + hi.v) * 0.5};
    const Vec3 point = face.pointAt(uv);
    Vec3 normal = face.normalAt(uv);
    if (!finite(point) || !finite(normal)) return 0;
    const double normalNorm = normal.norm();
    if (!(normalNorm > 0.0) || !std::isfinite(normalNorm)) return 0;
    normal = normal / normalNorm;

    Vec3 boxLo{}, boxHi{};
    solid.bbox(boxLo, boxHi);
    if (!finite(boxLo) || !finite(boxHi)) return 0;
    const double scale = (boxHi - boxLo).norm();
    if (!(scale > 0.0) || !std::isfinite(scale)) return 0;

    // Prefer the solid's own point-in-volume classification.  The bbox-centre
    // dot product remains a deterministic fallback for analytic providers that
    // classify boundary-near samples conservatively.
    const double offset = std::max(scale * 1.0e-7,
                                   std::numeric_limits<double>::epsilon() *
                                       scale * 256.0);
    const bool minusInside = solid.isInside(point - normal * offset);
    const bool plusInside = solid.isInside(point + normal * offset);
    if (minusInside != plusInside) return minusInside ? 1 : -1;

    const Vec3 centre = (boxLo + boxHi) * 0.5;
    const double dot = normal.dot(point - centre);
    if (!std::isfinite(dot) || dot == 0.0) return 0;
    return dot > 0.0 ? 1 : -1;
}

}  // namespace

VolumeFaceOwnership ModelBoundaryOwnershipProvider::ownership(
        FacePersistentId face) const {
    VolumeFaceOwnership out;
    out.providerStamp = stamp_;
    out.facePersistentId = face;
    const auto found = owners_.find(face);
    if (found == owners_.end()) {
        out.availability =
            VolumeProviderAvailability::MissingAuthoritativeOrientation;
        return out;
    }
    out.owners = found->second;
    return out;
}

ProviderRecoveryEnvelope ModelBoundaryOwnershipProvider::recoveryEnvelope(
        std::span<const Vec3>) const {
    ProviderRecoveryEnvelope out;
    out.providerStamp = stamp_;
    // The normal V0.10 request uses DerivedFiniteBox.  A CAD adapter may
    // override this provider and supply a tighter certified convex envelope.
    out.availability =
        VolumeProviderAvailability::RecoveryEnvelopeUnavailable;
    return out;
}

ModelBoundaryOwnershipProvider::BuildOutcome
buildModelBoundaryOwnershipProvider(
        const geo::Model& model,
        VolumeBoundaryProviderStamp stamp) noexcept {
    ModelBoundaryOwnershipProvider::BuildOutcome out;
    out.provider.stamp_ = stamp;
    try {
        auto solids = model.solids();
        if (solids.empty()) {
            out.reason =
                ModelBoundaryOwnershipProvider::BuildReason::EmptyModel;
            out.provider.reason_ = out.reason;
            return out;
        }
        std::sort(solids.begin(), solids.end(),
                  [](const geo::GSolid* a, const geo::GSolid* b) {
                      return a->tag() < b->tag();
                  });
        for (const geo::GSolid* solid : solids) {
            if (!solid) {
                out.reason = ModelBoundaryOwnershipProvider::BuildReason::
                    InternalError;
                out.provider.reason_ = out.reason;
                return out;
            }
            auto faceTags = solid->boundaryFaceTags();
            std::sort(faceTags.begin(), faceTags.end());
            faceTags.erase(std::unique(faceTags.begin(), faceTags.end()),
                           faceTags.end());
            for (const int faceTag : faceTags) {
                const geo::GFace* face = model.face(faceTag);
                if (!face) {
                    out.reason =
                        ModelBoundaryOwnershipProvider::BuildReason::
                            MissingFace;
                    out.provider.reason_ = out.reason;
                    return out;
                }
                const int sense = faceSenseForSolid(*face, *solid);
                if (sense == 0) {
                    out.reason =
                        ModelBoundaryOwnershipProvider::BuildReason::
                            AmbiguousFaceOrientation;
                    out.provider.reason_ = out.reason;
                    return out;
                }
                out.provider.owners_[persistentFaceId(faceTag)].push_back(
                    {solid->tag(), sense});
            }
        }

        for (auto& [face, owners] : out.provider.owners_) {
            (void)face;
            std::sort(owners.begin(), owners.end());
            owners.erase(std::unique(owners.begin(), owners.end()),
                         owners.end());
            if (owners.empty() || owners.size() > 2U ||
                (owners.size() == 2U &&
                 (owners[0].regionTag == owners[1].regionTag ||
                  owners[0].directedSense == owners[1].directedSense))) {
                out.reason =
                    ModelBoundaryOwnershipProvider::BuildReason::
                        InvalidSharedFace;
                out.provider.reason_ = out.reason;
                return out;
            }
        }
        out.provider.valid_ = true;
        out.provider.reason_ =
            ModelBoundaryOwnershipProvider::BuildReason::None;
        out.built = true;
        out.reason = ModelBoundaryOwnershipProvider::BuildReason::None;
        return out;
    } catch (const std::bad_alloc&) {
        out.reason =
            ModelBoundaryOwnershipProvider::BuildReason::AllocationFailed;
    } catch (...) {
        out.reason =
            ModelBoundaryOwnershipProvider::BuildReason::InternalError;
    }
    out.provider.reason_ = out.reason;
    return out;
}

Tet10PlacementResult ModelTet10PlacementProvider::place(
        const Tet10PlacementRequest& request) const {
    Tet10PlacementResult out;
    out.observedStamp = stamp_;
    if (!model_) {
        out.status = Tet10PlacementStatus::Failed;
        return out;
    }
    const Vec3 midpoint =
        (request.endpointCoordinates[0] + request.endpointCoordinates[1]) *
        0.5;
    if (!finite(midpoint)) {
        out.status = Tet10PlacementStatus::Rejected;
        out.providerReason = 1;
        return out;
    }
    try {
        if (request.kind == Tet10PlacementKind::Curve ||
            request.classification.dim == 1) {
            const geo::GEdge* edge = model_->edge(request.classification.tag);
            if (!edge) {
                out.status = Tet10PlacementStatus::Unsupported;
                out.providerReason = 2;
                return out;
            }
            out.coordinate = edge->closestPoint(midpoint);
        } else if (request.kind == Tet10PlacementKind::Surface ||
                   request.classification.dim == 2) {
            const geo::GFace* face = model_->face(request.classification.tag);
            if (!face) {
                out.status = Tet10PlacementStatus::Unsupported;
                out.providerReason = 3;
                return out;
            }
            const Vec2 uv = face->paramAt(midpoint);
            if (!finite(uv)) {
                out.status = Tet10PlacementStatus::Rejected;
                out.providerReason = 4;
                return out;
            }
            out.coordinate = face->pointAt(uv);
        } else {
            out.coordinate = midpoint;
        }
        if (!finite(out.coordinate)) {
            out.status = Tet10PlacementStatus::Rejected;
            out.providerReason = 5;
            return out;
        }
        out.status = Tet10PlacementStatus::Placed;
        return out;
    } catch (...) {
        out.status = Tet10PlacementStatus::Failed;
        out.providerReason = 6;
        return out;
    }
}

}  // namespace paacs::algo
