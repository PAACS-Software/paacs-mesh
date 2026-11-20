// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// BND-2A Increment 5: explicitly bounded SYNTHETIC FaceMesher adapter.
//
// This is not the normal meshFaces() route. It accepts one hand-built planar,
// non-periodic chart with one already-discretized outer loop, builds and
// retains a chart-local TopologyWorkspace, runs qualified BND recovery, audits
// the complete local result, and only then permits atomic publication.

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

#include "paacs/bnd/SteinerRecovery.hpp"
#include "paacs/mesh/Mesh.hpp"

namespace paacs::algo {

struct SyntheticFaceNode {
    Vec2 uv{};
    mesh::Classification classification{};
    // Existing global identity. Boundary-loop nodes MUST provide this.
    std::optional<Index> existingGlobalIndex;
    // Only declared constraint vertices receive an id. An existing interior
    // waypoint may intentionally omit it so BND can reuse it as StoredParam.
    std::optional<bnd::ConstraintVertexId> constraintVertex;
};

struct SyntheticFaceTriangle {
    std::array<std::size_t, 3> nodes{};
};

struct SyntheticFaceInput {
    bnd::ChartId chart = 0;
    bool periodic = false;
    mesh::Classification faceClass{2, -1};
    std::vector<SyntheticFaceNode> nodes;
    std::vector<SyntheticFaceTriangle> triangles;
    // Exactly one simple loop, ordered; implicit closing edge last->first.
    std::vector<std::size_t> outerLoop;
    // Includes one already-present segment for EVERY outer-loop edge plus any
    // interior prescribed constraints. Boundary Steiner creation is thereby
    // prohibited structurally: every boundary segment routes AlreadyPresent.
    std::vector<bnd::ConstraintSegment> constraints;
};

struct FaceMeshCell {
    std::array<std::size_t, 3> nodes{};
    mesh::Classification classification{};
};

struct FaceRecoveredConstraint {
    bnd::ConstraintEntityId entity = 0;
    mesh::Classification classification{};
    bnd::SteinerStatus status = bnd::SteinerStatus::Failed;
    std::vector<std::size_t> chain;
};

struct FaceMeshAudit {
    bool workspaceReady = false;
    bool tri3Manifold = false;
    bool oneOuterLoop = false;
    bool orientationConsistent = false;
    bool classificationsPreserved = false;
    bool mappingsInjective = false;
    bool constraintsRecovered = false;
    bool boundarySteinerFree = false;

    bool passed() const noexcept {
        return workspaceReady && tri3Manifold && oneOuterLoop &&
               orientationConsistent && classificationsPreserved &&
               mappingsInjective && constraintsRecovered &&
               boundarySteinerFree;
    }
};

struct FaceMeshResult {
    bnd::ChartId chart = 0;
    mesh::Classification faceClass{};
    std::vector<topo::NodePool::NodeData> localNodes;
    std::vector<Vec2> localParams;
    std::vector<topo::NodeHandle> localHandles;
    std::vector<FaceMeshCell> localCells;
    // local result index -> pre-existing global mesh node identity.
    std::unordered_map<std::size_t, Index> existingGlobalIndex;
    std::vector<std::size_t> outerLoop;
    std::vector<FaceRecoveredConstraint> recoveredConstraints;
    FaceMeshAudit audit;
    // Retained local topology is the authority behind the exported values.
    std::unique_ptr<topo::TopologyWorkspace> localWorkspace;
};

enum class FaceAdapterStatus { Ready, Refused, Failed, Rejected };
enum class FaceAdapterReason {
    None,
    UnsupportedPeriodicChart,
    InvalidInput,
    BoundarySteinerProhibited,
    RecoveryStopped,
    LocalAuditFailed,
    AllocationFailed,
    InternalError
};

struct FaceAdapterOutcome {
    FaceAdapterStatus status = FaceAdapterStatus::Failed;
    FaceAdapterReason reason = FaceAdapterReason::InternalError;
    bnd::SteinerReject reject = bnd::SteinerReject::None;
    bnd::SteinerRefuse refuse = bnd::SteinerRefuse::None;
    bnd::SteinerFail fail = bnd::SteinerFail::None;
    bnd::SteinerSegmentsReport recovery;
    std::unique_ptr<FaceMeshResult> result;
};

// Builds local topology only. `global` is const and can never be mutated by a
// rejection/failure. No exception escapes this five-outcome adapter.
FaceAdapterOutcome buildSyntheticFaceMesh(
    const SyntheticFaceInput& input, const mesh::Mesh& global,
    const bnd::GeometryProvider& provider,
    const bnd::GeometryProviderStamp& expectedProvider,
    const bnd::SteinerPolicy& policy = {});

enum class FacePublishStatus { Published, Failed, Rejected };
enum class FacePublishReason {
    None,
    InvalidResult,
    ExistingIdentityMismatch,
    AllocationFailed,
    InternalError
};

struct FacePublication {
    FacePublishStatus status = FacePublishStatus::Failed;
    FacePublishReason reason = FacePublishReason::InternalError;
    std::vector<Index> localToGlobal;
};

// Audits again against the destination mesh, prepares all mappings and
// capacities off-mesh, then publishes with Mesh::addNode() (never
// addNodeUnique()). Any reported failure/rejection leaves logical global mesh
// content unchanged.
FacePublication publishFaceMeshResult(mesh::Mesh& global,
                                      const FaceMeshResult& result) noexcept;

}  // namespace paacs::algo
