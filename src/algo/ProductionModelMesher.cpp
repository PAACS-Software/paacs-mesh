// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/ProductionModelMesher.hpp"

#include <algorithm>
#include <cmath>
#include <memory>

#include "CheckedSizing.hpp"
#include "paacs/algo/Meshers.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/mesh/MeshAudit.hpp"

namespace paacs::algo {
namespace {

bool validSizingPolicy(const sizing::MeshSpec& spec) noexcept {
    return (spec.order == 1 || spec.order == 2) &&
           std::isfinite(spec.orderFactor) && spec.orderFactor > 0.0 &&
           std::isfinite(spec.gradation) && spec.gradation >= 1.0 &&
           spec.elementsPerCircle >= 0 && spec.threads >= 0;
}

void configureDedupTolerance(const geo::Model& model, mesh::Mesh& mesh) {
    Vec3 lo{1e300, 1e300, 1e300};
    Vec3 hi{-1e300, -1e300, -1e300};
    bool haveSolid = false;
    for (const geo::GSolid* solid : model.solids()) {
        Vec3 a, b;
        solid->bbox(a, b);
        lo.x = std::min(lo.x, a.x);
        lo.y = std::min(lo.y, a.y);
        lo.z = std::min(lo.z, a.z);
        hi.x = std::max(hi.x, b.x);
        hi.y = std::max(hi.y, b.y);
        hi.z = std::max(hi.z, b.z);
        haveSolid = true;
    }
    if (!haveSolid) return;
    const double diagonal = (hi - lo).norm();
    if (diagonal > 0.0 && diagonal < 1e300)
        mesh.setDedupTolerance(1e-9 * diagonal);
}

ProductionModelStatus modelStatus(ProductionSurfaceStatus status) noexcept {
    switch (status) {
        case ProductionSurfaceStatus::Published:
            return ProductionModelStatus::Published;
        case ProductionSurfaceStatus::Rejected:
            return ProductionModelStatus::Rejected;
        case ProductionSurfaceStatus::Unsupported:
            return ProductionModelStatus::Unsupported;
        case ProductionSurfaceStatus::Failed:
            return ProductionModelStatus::Failed;
    }
    return ProductionModelStatus::Failed;
}

}  // namespace

ProductionModelOutcome meshModelConstrained(
    const geo::Model& model, const sizing::MeshSpec& spec,
    mesh::Mesh& destination) noexcept {
    ProductionModelOutcome out;
    if (!validSizingPolicy(spec)) {
        out.status = ProductionModelStatus::Rejected;
        out.reason = ProductionModelReason::InvalidSizingPolicy;
        out.phase = ProductionModelPhase::Preflight;
        return out;
    }

    std::unique_ptr<mesh::Mesh> local;
    sizing::MeshSpec checkedSpec;
    try {
        out.phase = ProductionModelPhase::EdgeConstruction;
        // MSVC's unordered_map default construction may allocate. Keep the
        // complete scratch object inside the caught boundary; constructing a
        // named Mesh before this try would turn allocation #1 into a
        // throw-through-noexcept fast-fail.
        local = std::make_unique<mesh::Mesh>();
        checkedSpec = detail::makeCheckedSpec(spec);
        core::setThreadCount(checkedSpec.threads);
        configureDedupTolerance(model, *local);
        meshEdges(model, checkedSpec, *local);
    } catch (const detail::SizingEvaluationError&) {
        out.reason = ProductionModelReason::SizingEvaluationFailed;
        return out;
    } catch (const std::bad_alloc&) {
        out.reason = ProductionModelReason::AllocationFailed;
        return out;
    } catch (...) {
        out.reason = ProductionModelReason::EdgeConstructionFailed;
        return out;
    }

    out.phase = ProductionModelPhase::SurfaceConstruction;
    out.surface = meshFacesConstrained(model, checkedSpec, *local);
    if (!out.surface.published()) {
        out.status = modelStatus(out.surface.status);
        if (out.surface.reason == ProductionSurfaceReason::AllocationFailed)
            out.reason = ProductionModelReason::AllocationFailed;
        else if (out.surface.reason ==
                 ProductionSurfaceReason::SizingEvaluationFailed)
            out.reason = ProductionModelReason::SizingEvaluationFailed;
        else
            out.reason = ProductionModelReason::SurfaceStopped;
        return out;
    }

    try {
        out.phase = ProductionModelPhase::VolumeConstruction;
        meshVolumesAuthoritative(model, checkedSpec, *local, out.tetHandoff);
    } catch (const detail::SizingEvaluationError&) {
        out.reason = ProductionModelReason::SizingEvaluationFailed;
        return out;
    } catch (const std::bad_alloc&) {
        out.reason = ProductionModelReason::AllocationFailed;
        return out;
    } catch (...) {
        out.reason = ProductionModelReason::VolumeConstructionFailed;
        return out;
    }

    if (spec.order == 2) {
        try {
            out.phase = ProductionModelPhase::OrderElevation;
            const auto elevation = elevateToQuadratic(model, *local);
            if (elevation.invalidAfter != 0) {
                out.reason = ProductionModelReason::OrderElevationFailed;
                return out;
            }
        } catch (const std::bad_alloc&) {
            out.reason = ProductionModelReason::AllocationFailed;
            return out;
        } catch (...) {
            out.reason = ProductionModelReason::OrderElevationFailed;
            return out;
        }
    }

    try {
        out.phase = ProductionModelPhase::FinalAudit;
        mesh::AuditOptions options;
        options.boundaryMode = mesh::BoundaryComparisonMode::Required;
        options.boundaryOrphanBudget = 0;
        options.boundaryUnexpectedBudget = 0;
        options.duplicateFacetsAreErrors = true;
        options.unusedNodesAreErrors = true;
        const auto solids = model.solids();
        if (!solids.empty())
            options.expectedComponents = static_cast<int>(solids.size());

        const auto audit = mesh::auditMesh(*local, options);
        out.audit.inverted = audit.nInverted;
        out.audit.invalidTet10 = audit.nInvalidTet10;
        out.audit.nonFinite = audit.nNonFinite;
        out.audit.malformed = audit.malformedCells;
        out.audit.deferred = audit.deferredAuditCells;
        out.audit.prescribed = audit.boundaryPrescribed;
        out.audit.matched = audit.boundaryMatched;
        out.audit.orphan = audit.boundaryOrphans;
        out.audit.unexpected = audit.boundaryUnexpected;
        out.audit.duplicate = audit.boundaryDuplicates;
        out.audit.nonManifold = audit.nonManifoldFaces;
        out.audit.connectedComponents = audit.connectedComponents;
        const bool exact = audit.passed() && audit.nInverted == 0 &&
                           audit.nInvalidTet10 == 0 &&
                           audit.nNonFinite == 0 &&
                           audit.malformedCells == 0 &&
                           audit.deferredAuditCells == 0 &&
                           audit.boundaryDuplicates == 0 &&
                           audit.boundaryOrphans == 0 &&
                           audit.boundaryUnexpected == 0 &&
                           audit.nonManifoldFaces == 0 &&
                           audit.boundaryPrescribed == audit.boundaryMatched;
        if (!exact) {
            out.reason = ProductionModelReason::FinalAuditFailed;
            return out;
        }
    } catch (const std::bad_alloc&) {
        out.reason = ProductionModelReason::AllocationFailed;
        return out;
    } catch (...) {
        out.reason = ProductionModelReason::FinalAuditFailed;
        return out;
    }

    try {
        out.phase = ProductionModelPhase::Publication;
        // Destination preparation is deliberately separate. If this one
        // allocation fails, the completed local result and caller destination
        // have not been touched. Everything after it is a noexcept swap.
        auto prepared = std::make_unique<mesh::Mesh>();
        prepared->swap(*local);
        out.publishedNodes = prepared->nodeCount();
        out.publishedCells = prepared->cellCount();
        destination.swap(*prepared);
    } catch (const std::bad_alloc&) {
        out.reason = ProductionModelReason::AllocationFailed;
        return out;
    } catch (...) {
        out.reason = ProductionModelReason::InternalError;
        return out;
    }

    out.status = ProductionModelStatus::Published;
    out.reason = ProductionModelReason::None;
    out.destinationUnchanged = false;
    return out;
}

}  // namespace paacs::algo
