// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/sizing/MeshControls.hpp"

#include <cstdint>
#include <limits>
#include <vector>

namespace paacs::algo {

enum class BodyMeshingCapability : std::uint32_t {
    None = 0,
    GeneralTet4 = 1U << 0,
    MappedBlock = 1U << 1,
    PreparedSweep = 1U << 2,
    StructuredPolyhedral = 1U << 3,
    PrismLayer = 1U << 4,
    PyramidTransition = 1U << 5
};
constexpr BodyMeshingCapability operator|(BodyMeshingCapability a,
                                          BodyMeshingCapability b) noexcept {
    return static_cast<BodyMeshingCapability>(
        static_cast<std::uint32_t>(a) | static_cast<std::uint32_t>(b));
}
constexpr BodyMeshingCapability operator&(BodyMeshingCapability a,
                                          BodyMeshingCapability b) noexcept {
    return static_cast<BodyMeshingCapability>(
        static_cast<std::uint32_t>(a) & static_cast<std::uint32_t>(b));
}
constexpr bool hasCapability(BodyMeshingCapability value,
                             BodyMeshingCapability required) noexcept {
    return (value & required) == required;
}

enum class PortalTopology : std::uint8_t { Tri3, Quad4, MixedTriQuad };
enum class PortalStrategy : std::uint8_t {
    SharedTri3,
    SharedQuad4,
    SharedMixed,
    Pyramid5Bridge
};

struct BodyMeshingPlanInput {
    int solidTag = -1;
    sizing::VolumeEngineKind requestedEngine =
        sizing::VolumeEngineKind::Auto;
    BodyMeshingCapability capabilities = BodyMeshingCapability::None;
    std::uint64_t capabilityReceiptHash = 0;
    std::uint64_t estimatedNodes = 0;
    std::uint64_t estimatedCells = 0;
    bool allowAutoFallbackToTet = true;
};
struct InterfacePortalInput {
    int faceTag = -1;
    int firstSolidTag = -1;
    int secondSolidTag = -1;
    PortalTopology topology = PortalTopology::Tri3;
    std::uint64_t nodeCount = 0;
    std::uint64_t faceCount = 0;
    std::uint64_t identityHash = 0;
    bool allowPyramidTransition = true;
};
struct MeshingPlanRequest {
    std::vector<BodyMeshingPlanInput> bodies;
    std::vector<InterfacePortalInput> portals;
    std::uint64_t compiledControlHash = 0;
    std::uint64_t compiledSpatialFieldHash = 0;
    std::uint64_t maximumNodes =
        (std::numeric_limits<std::uint64_t>::max)();
    std::uint64_t maximumCells =
        (std::numeric_limits<std::uint64_t>::max)();
};
struct PlannedBody {
    int solidTag = -1;
    sizing::VolumeEngineKind selectedEngine = sizing::VolumeEngineKind::Tet4;
    BodyMeshingCapability capabilities = BodyMeshingCapability::None;
    std::uint64_t capabilityReceiptHash = 0;
    std::uint64_t estimatedNodes = 0;
    std::uint64_t estimatedCells = 0;
};
struct PlannedPortal {
    int faceTag = -1;
    int firstSolidTag = -1;
    int secondSolidTag = -1;
    PortalTopology topology = PortalTopology::Tri3;
    PortalStrategy strategy = PortalStrategy::SharedTri3;
    std::uint64_t nodeCount = 0;
    std::uint64_t faceCount = 0;
    std::uint64_t identityHash = 0;
};

struct MeshingPlanOutcome;
class MeshingPlan {
public:
    const std::vector<PlannedBody>& bodies() const noexcept { return bodies_; }
    const std::vector<PlannedPortal>& portals() const noexcept { return portals_; }
    const PlannedBody* body(int solidTag) const noexcept;
    std::uint64_t totalEstimatedNodes() const noexcept { return totalNodes_; }
    std::uint64_t totalEstimatedCells() const noexcept { return totalCells_; }
    std::uint64_t deterministicHash() const noexcept { return hash_; }
    std::uint64_t compiledControlHash() const noexcept { return controlHash_; }
    std::uint64_t compiledSpatialFieldHash() const noexcept { return fieldHash_; }
private:
    friend MeshingPlanOutcome compileMeshingPlan(
        const MeshingPlanRequest&) noexcept;
    std::vector<PlannedBody> bodies_;
    std::vector<PlannedPortal> portals_;
    std::uint64_t totalNodes_ = 0;
    std::uint64_t totalCells_ = 0;
    std::uint64_t controlHash_ = 0;
    std::uint64_t fieldHash_ = 0;
    std::uint64_t hash_ = 0;
};

enum class MeshingPlanReason : std::uint8_t {
    None,
    EmptyModel,
    InvalidBody,
    DuplicateBody,
    UnsupportedRequestedEngine,
    MissingCapabilityReceipt,
    InvalidPortal,
    DuplicatePortal,
    MissingPortalBody,
    IncompatiblePortal,
    CountOverflow,
    BudgetExceeded,
    AllocationFailed,
    InternalError
};
struct MeshingPlanMetrics {
    std::uint64_t bodies = 0;
    std::uint64_t portals = 0;
    std::uint64_t estimatedNodes = 0;
    std::uint64_t estimatedCells = 0;
    std::uint64_t tetBodies = 0;
    std::uint64_t hexBodies = 0;
    std::uint64_t polyhedralBodies = 0;
    std::uint64_t sweptBodies = 0;
    std::uint64_t pyramidPortals = 0;
};
struct MeshingPlanOutcome {
    bool compiled = false;
    MeshingPlanReason reason = MeshingPlanReason::InternalError;
    MeshingPlan plan{};
    MeshingPlanMetrics metrics{};
};

MeshingPlanOutcome compileMeshingPlan(
    const MeshingPlanRequest& request) noexcept;

}  // namespace paacs::algo
