// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MeshingPlan.hpp"

#include <algorithm>
#include <limits>
#include <new>

namespace paacs::algo {
namespace {

bool checkedAdd(std::uint64_t& total, std::uint64_t value) noexcept {
    if (value > (std::numeric_limits<std::uint64_t>::max)() - total)
        return false;
    total += value;
    return true;
}
void hashValue(std::uint64_t& hash, std::uint64_t value) noexcept {
    for (unsigned i = 0; i < 8; ++i) {
        hash ^= static_cast<unsigned char>(value >> (i * 8U));
        hash *= 1099511628211ULL;
    }
}
bool needsReceipt(sizing::VolumeEngineKind engine) noexcept {
    return engine == sizing::VolumeEngineKind::MappedHex8 ||
           engine == sizing::VolumeEngineKind::StructuredPolyhedral ||
           engine == sizing::VolumeEngineKind::LayeredPrism6;
}
bool supports(const BodyMeshingPlanInput& body,
              sizing::VolumeEngineKind engine) noexcept {
    switch (engine) {
        case sizing::VolumeEngineKind::Auto: return true;
        case sizing::VolumeEngineKind::Tet4:
            return hasCapability(body.capabilities,
                                 BodyMeshingCapability::GeneralTet4);
        case sizing::VolumeEngineKind::MappedHex8:
            return hasCapability(body.capabilities,
                                 BodyMeshingCapability::MappedBlock);
        case sizing::VolumeEngineKind::StructuredPolyhedral:
            return hasCapability(body.capabilities,
                                 BodyMeshingCapability::StructuredPolyhedral);
        case sizing::VolumeEngineKind::LayeredPrism6:
            return hasCapability(body.capabilities,
                                 BodyMeshingCapability::PreparedSweep) ||
                   hasCapability(body.capabilities,
                                 BodyMeshingCapability::PrismLayer);
    }
    return false;
}
sizing::VolumeEngineKind autoEngine(const BodyMeshingPlanInput& body) noexcept {
    if (hasCapability(body.capabilities, BodyMeshingCapability::MappedBlock))
        return sizing::VolumeEngineKind::MappedHex8;
    if (hasCapability(body.capabilities, BodyMeshingCapability::PreparedSweep))
        return sizing::VolumeEngineKind::LayeredPrism6;
    if (hasCapability(body.capabilities, BodyMeshingCapability::GeneralTet4))
        return sizing::VolumeEngineKind::Tet4;
    if (hasCapability(body.capabilities,
                      BodyMeshingCapability::StructuredPolyhedral))
        return sizing::VolumeEngineKind::StructuredPolyhedral;
    return sizing::VolumeEngineKind::Auto;
}
bool isQuadBased(sizing::VolumeEngineKind engine) noexcept {
    return engine == sizing::VolumeEngineKind::MappedHex8 ||
           engine == sizing::VolumeEngineKind::StructuredPolyhedral;
}
}  // namespace

const PlannedBody* MeshingPlan::body(int solidTag) const noexcept {
    const auto it = std::lower_bound(
        bodies_.begin(), bodies_.end(), solidTag,
        [](const PlannedBody& body, int tag) { return body.solidTag < tag; });
    return it != bodies_.end() && it->solidTag == solidTag ? &*it : nullptr;
}

MeshingPlanOutcome compileMeshingPlan(
        const MeshingPlanRequest& request) noexcept {
    MeshingPlanOutcome out;
    if (request.bodies.empty()) {
        out.reason = MeshingPlanReason::EmptyModel;
        return out;
    }
    try {
        auto bodies = request.bodies;
        std::sort(bodies.begin(), bodies.end(),
                  [](const auto& a, const auto& b) {
                      return a.solidTag < b.solidTag;
                  });
        out.plan.bodies_.reserve(bodies.size());
        for (std::size_t i = 0; i < bodies.size(); ++i) {
            const auto& input = bodies[i];
            if (input.solidTag < 0 || input.estimatedCells == 0 ||
                input.estimatedNodes == 0 ||
                input.capabilities == BodyMeshingCapability::None) {
                out.reason = MeshingPlanReason::InvalidBody;
                return out;
            }
            if (i && bodies[i-1].solidTag == input.solidTag) {
                out.reason = MeshingPlanReason::DuplicateBody;
                return out;
            }
            sizing::VolumeEngineKind engine = input.requestedEngine;
            if (engine == sizing::VolumeEngineKind::Auto)
                engine = autoEngine(input);
            if (engine == sizing::VolumeEngineKind::Auto ||
                !supports(input, engine)) {
                if (input.allowAutoFallbackToTet &&
                    supports(input, sizing::VolumeEngineKind::Tet4))
                    engine = sizing::VolumeEngineKind::Tet4;
                else {
                    out.reason = MeshingPlanReason::UnsupportedRequestedEngine;
                    return out;
                }
            }
            if (needsReceipt(engine) && input.capabilityReceiptHash == 0) {
                out.reason = MeshingPlanReason::MissingCapabilityReceipt;
                return out;
            }
            if (!checkedAdd(out.plan.totalNodes_, input.estimatedNodes) ||
                !checkedAdd(out.plan.totalCells_, input.estimatedCells)) {
                out.reason = MeshingPlanReason::CountOverflow;
                return out;
            }
            out.plan.bodies_.push_back({
                input.solidTag, engine, input.capabilities,
                input.capabilityReceiptHash, input.estimatedNodes,
                input.estimatedCells});
            switch (engine) {
                case sizing::VolumeEngineKind::Tet4: ++out.metrics.tetBodies; break;
                case sizing::VolumeEngineKind::MappedHex8: ++out.metrics.hexBodies; break;
                case sizing::VolumeEngineKind::StructuredPolyhedral:
                    ++out.metrics.polyhedralBodies; break;
                case sizing::VolumeEngineKind::LayeredPrism6:
                    ++out.metrics.sweptBodies; break;
                case sizing::VolumeEngineKind::Auto: break;
            }
        }
        if (out.plan.totalNodes_ > request.maximumNodes ||
            out.plan.totalCells_ > request.maximumCells) {
            out.reason = MeshingPlanReason::BudgetExceeded;
            return out;
        }

        auto portals = request.portals;
        std::sort(portals.begin(), portals.end(), [](const auto& a,
                                                     const auto& b) {
            if (a.faceTag != b.faceTag) return a.faceTag < b.faceTag;
            if (a.firstSolidTag != b.firstSolidTag)
                return a.firstSolidTag < b.firstSolidTag;
            return a.secondSolidTag < b.secondSolidTag;
        });
        out.plan.portals_.reserve(portals.size());
        for (std::size_t i = 0; i < portals.size(); ++i) {
            const auto& portal = portals[i];
            if (portal.faceTag < 0 || portal.firstSolidTag < 0 ||
                portal.secondSolidTag < 0 ||
                portal.firstSolidTag == portal.secondSolidTag ||
                portal.nodeCount == 0 || portal.faceCount == 0 ||
                portal.identityHash == 0) {
                out.reason = MeshingPlanReason::InvalidPortal;
                return out;
            }
            if (i && portals[i-1].faceTag == portal.faceTag) {
                out.reason = MeshingPlanReason::DuplicatePortal;
                return out;
            }
            const PlannedBody* first = out.plan.body(portal.firstSolidTag);
            const PlannedBody* second = out.plan.body(portal.secondSolidTag);
            if (!first || !second) {
                out.reason = MeshingPlanReason::MissingPortalBody;
                return out;
            }
            PortalStrategy strategy = PortalStrategy::SharedTri3;
            if (portal.topology == PortalTopology::Tri3) {
                strategy = PortalStrategy::SharedTri3;
            } else if (portal.topology == PortalTopology::MixedTriQuad) {
                strategy = PortalStrategy::SharedMixed;
            } else {
                const bool firstQuad = isQuadBased(first->selectedEngine);
                const bool secondQuad = isQuadBased(second->selectedEngine);
                if (firstQuad && secondQuad)
                    strategy = PortalStrategy::SharedQuad4;
                else if (firstQuad != secondQuad &&
                         portal.allowPyramidTransition &&
                         hasCapability(first->capabilities,
                             BodyMeshingCapability::PyramidTransition) &&
                         hasCapability(second->capabilities,
                             BodyMeshingCapability::PyramidTransition)) {
                    strategy = PortalStrategy::Pyramid5Bridge;
                    ++out.metrics.pyramidPortals;
                } else {
                    out.reason = MeshingPlanReason::IncompatiblePortal;
                    return out;
                }
            }
            out.plan.portals_.push_back({
                portal.faceTag, portal.firstSolidTag, portal.secondSolidTag,
                portal.topology, strategy, portal.nodeCount, portal.faceCount,
                portal.identityHash});
        }

        out.plan.controlHash_ = request.compiledControlHash;
        out.plan.fieldHash_ = request.compiledSpatialFieldHash;
        std::uint64_t hash = 1469598103934665603ULL;
        hashValue(hash, 1);
        hashValue(hash, out.plan.controlHash_);
        hashValue(hash, out.plan.fieldHash_);
        for (const auto& body : out.plan.bodies_) {
            hashValue(hash, static_cast<std::uint64_t>(body.solidTag));
            hashValue(hash, static_cast<std::uint64_t>(body.selectedEngine));
            hashValue(hash, static_cast<std::uint64_t>(body.capabilities));
            hashValue(hash, body.capabilityReceiptHash);
            hashValue(hash, body.estimatedNodes);
            hashValue(hash, body.estimatedCells);
        }
        for (const auto& portal : out.plan.portals_) {
            hashValue(hash, static_cast<std::uint64_t>(portal.faceTag));
            hashValue(hash, static_cast<std::uint64_t>(portal.firstSolidTag));
            hashValue(hash, static_cast<std::uint64_t>(portal.secondSolidTag));
            hashValue(hash, static_cast<std::uint64_t>(portal.topology));
            hashValue(hash, static_cast<std::uint64_t>(portal.strategy));
            hashValue(hash, portal.nodeCount);
            hashValue(hash, portal.faceCount);
            hashValue(hash, portal.identityHash);
        }
        out.plan.hash_ = hash;
        out.metrics.bodies = out.plan.bodies_.size();
        out.metrics.portals = out.plan.portals_.size();
        out.metrics.estimatedNodes = out.plan.totalNodes_;
        out.metrics.estimatedCells = out.plan.totalCells_;
        out.compiled = true;
        out.reason = MeshingPlanReason::None;
    } catch (const std::bad_alloc&) {
        out.reason = MeshingPlanReason::AllocationFailed;
    } catch (...) {
        out.reason = MeshingPlanReason::InternalError;
    }
    return out;
}

}  // namespace paacs::algo
