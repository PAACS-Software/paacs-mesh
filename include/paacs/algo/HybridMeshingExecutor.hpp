// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/algo/MeshingPlan.hpp"
#include "paacs/algo/PreparedSweepMesher.hpp"
#include "paacs/algo/StructuredPolyhedralMesher.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/mesh/PolyhedralMesh.hpp"

#include <cstdint>
#include <functional>
#include <limits>
#include <vector>

namespace paacs::algo {

enum class HybridBodySourceKind : std::uint8_t {
    PreparedFixedMesh,
    MappedHexBlock,
    PreparedMixedSweep,
    StructuredPolyhedralBlock,
    PreparedPolyhedralMesh
};

struct HybridBodyInput {
    int solidTag = -1;
    HybridBodySourceKind sourceKind = HybridBodySourceKind::PreparedFixedMesh;
    const mesh::Mesh* preparedFixedMesh = nullptr;
    const mesh::PolyhedralMesh* preparedPolyhedralMesh = nullptr;
    MappedHexBlockRequest mappedHex{};
    const mesh::Mesh* sweepSourceSurface = nullptr;
    PreparedMixedSweepRequest sweep{};
    StructuredPolyBlockRequest structuredPolyhedral{};
};

// Explicit identity binding across two body-local meshes. This prevents
// coordinate-only guessing at interfaces and lets the executor publish one
// shared global node for conformal fixed-cell portals.
struct HybridPortalNodeBinding {
    int faceTag = -1;
    int firstSolidTag = -1;
    Index firstLocalNode = kInvalidIndex;
    int secondSolidTag = -1;
    Index secondLocalNode = kInvalidIndex;
};

enum class HybridExecutionStage : std::uint8_t {
    Preflight, BodyGeneration, PortalReconciliation, Assembly, Audit, Publication
};
using HybridProgressCallback = std::function<bool(
    HybridExecutionStage, std::uint64_t, std::uint64_t)>;

struct HybridExecutionRequest {
    const MeshingPlan* plan = nullptr;
    std::vector<HybridBodyInput> bodies;
    std::vector<HybridPortalNodeBinding> portalNodeBindings;
    double portalCoordinateTolerance = 1.0e-11;
    bool suppressDuplicateBoundaryFaces = true;
    std::uint64_t maximumOutputNodes =
        (std::numeric_limits<std::uint64_t>::max)();
    std::uint64_t maximumOutputCells =
        (std::numeric_limits<std::uint64_t>::max)();
    HybridProgressCallback progress;
};

struct HybridBodyRecord {
    int solidTag = -1;
    sizing::VolumeEngineKind engine = sizing::VolumeEngineKind::Tet4;
    bool polyhedral = false;
    Index firstCell = 0;
    Index cellCount = 0;
    Index generatedNodes = 0;
    Index generatedCells = 0;
    std::uint64_t strictHash = 0;
};
struct HybridPortalRecord {
    int faceTag = -1;
    int firstSolidTag = -1;
    int secondSolidTag = -1;
    PortalStrategy strategy = PortalStrategy::SharedTri3;
    std::uint64_t sharedNodeCount = 0;
    std::uint64_t identityHash = 0;
};

class HybridMesh {
public:
    mesh::Mesh fixed;
    mesh::PolyhedralMesh polyhedral;
    std::vector<HybridBodyRecord> bodies;
    std::vector<HybridPortalRecord> portals;
    void swap(HybridMesh& other) noexcept;
};

enum class HybridExecutionReason : std::uint8_t {
    None,
    MissingPlan,
    PlanBodyMismatch,
    DuplicateBody,
    MissingBody,
    EngineSourceMismatch,
    BodyGenerationFailed,
    InvalidPreparedMesh,
    InvalidPortalBinding,
    PortalCoordinateMismatch,
    PortalClassificationMismatch,
    UnsupportedPolyhedralPortal,
    CountOverflow,
    BudgetExceeded,
    Cancelled,
    AuditFailed,
    AllocationFailed,
    InternalError
};
struct HybridExecutionMetrics {
    std::uint64_t bodiesGenerated = 0;
    std::uint64_t fixedBodies = 0;
    std::uint64_t polyhedralBodies = 0;
    std::uint64_t portalBindings = 0;
    std::uint64_t mergedPortalNodes = 0;
    std::uint64_t duplicateBoundaryFacesSuppressed = 0;
    std::uint64_t outputNodes = 0;
    std::uint64_t outputCells = 0;
    std::uint64_t outputPolyhedralCells = 0;
};
struct HybridExecutionOutcome {
    bool published = false;
    HybridExecutionReason reason = HybridExecutionReason::InternalError;
    HybridExecutionStage stage = HybridExecutionStage::Preflight;
    bool destinationUnchanged = true;
    HybridExecutionMetrics metrics{};
};

HybridExecutionOutcome executeHybridMeshingPlan(
    const HybridExecutionRequest& request,
    HybridMesh& destination) noexcept;

}  // namespace paacs::algo
