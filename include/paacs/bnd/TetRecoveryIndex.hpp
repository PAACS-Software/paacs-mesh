// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// T2-B: one-time audited, deterministic Tet recovery index.

#include "paacs/topo/TetCavity.hpp"
#include "paacs/bnd/TetFacetGeometry.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <limits>
#include <span>
#include <type_traits>
#include <vector>

namespace paacs::bnd {

using DenseTetNodeId = std::uint32_t;
inline constexpr DenseTetNodeId kInvalidDenseTetNode =
    std::numeric_limits<DenseTetNodeId>::max();

struct TetAabb {
    Vec3 lo{}, hi{};
};

struct DenseTetEdgeKey {
    std::array<DenseTetNodeId, 2> nodes{};
    static DenseTetEdgeKey of(DenseTetNodeId a, DenseTetNodeId b) noexcept {
        if (b < a) std::swap(a, b);
        return {{a, b}};
    }
    friend bool operator==(const DenseTetEdgeKey&,
                           const DenseTetEdgeKey&) = default;
    friend bool operator<(const DenseTetEdgeKey& a,
                          const DenseTetEdgeKey& b) noexcept {
        return a.nodes < b.nodes;
    }
};

struct DenseTetFaceKey {
    std::array<DenseTetNodeId, 3> nodes{};
    static DenseTetFaceKey of(std::array<DenseTetNodeId, 3> value) noexcept {
        std::sort(value.begin(), value.end());
        return {value};
    }
    friend bool operator==(const DenseTetFaceKey&,
                           const DenseTetFaceKey&) = default;
    friend bool operator<(const DenseTetFaceKey& a,
                          const DenseTetFaceKey& b) noexcept {
        return a.nodes < b.nodes;
    }
};

struct CurrentTetRecord {
    topo::CellHandle handle{};
    std::uint64_t sequence = 0;
    std::array<DenseTetNodeId, 4> nodes{};
    mesh::Classification classification{};
    TetAabb bounds{};
};

struct EdgeResolverRow {
    topo::EdgeHandle handle{};
    std::vector<topo::CellHandle> completeStar;
    mesh::Classification classification{};
    bool required = false;
};

struct FaceResolverIncidence {
    topo::CellHandle cell{};
    std::int8_t sense = 0;
};

struct FaceResolverRow {
    topo::FaceHandle handle{};
    std::vector<FaceResolverIncidence> owners;
    mesh::Classification classification{};
};

class DynamicTetAabbIndex {
public:
    struct Node {
        TetAabb bounds{};
        std::int32_t left = -1, right = -1;
        topo::CellHandle cell{};
        std::uint64_t minimumSequence = 0;
        bool leaf = false;
    };

    void build(const std::map<topo::CellHandle, CurrentTetRecord>& cells);
    std::vector<topo::CellHandle> query(
        const TetAabb& bounds, std::uint64_t* nodesVisited = nullptr) const;
    bool covers(const std::map<topo::CellHandle, CurrentTetRecord>& cells)
        const noexcept;
    std::size_t accountedBytes() const noexcept {
        return nodes_.capacity() * sizeof(Node) +
               sortedLeafHandles_.capacity() * sizeof(topo::CellHandle);
    }
    const std::vector<Node>& nodes() const noexcept { return nodes_; }

private:
    std::vector<Node> nodes_;
    // Built once with the BVH and kept in canonical handle order.  This makes
    // covers() exact and allocation-free even on standard libraries whose
    // empty ordered containers allocate a sentinel.
    std::vector<topo::CellHandle> sortedLeafHandles_;
};

enum class TetRecoveryContextReason : std::uint8_t {
    None,
    StaleWorkspace,
    InvalidClassification,
    UnsupportedCell,
    NonfiniteNode,
    DuplicateNodeSequence,
    InvalidConnectivity,
    NonPositiveTet,
    DuplicateTet,
    MissingEdgeIncidence,
    MissingFaceIncidence,
    NonManifoldFace,
    InvalidFaceSense,
    OrphanNode,
    IndexCoverageMismatch,
    AllocationFailed,
    InternalError
};

struct TetRecoveryContextMetrics {
    std::uint64_t nodes = 0, cells = 0, edges = 0, faces = 0;
    std::uint64_t auditVisits = 0, indexNodes = 0,
                  segmentQueries = 0, facetQueries = 0, boundsQueries = 0,
                  indexNodesVisited = 0, returnedCandidates = 0,
                  indexRemovals = 0, indexInsertions = 0,
                  indexRebuilds = 0, deltaHighWater = 0;
    std::size_t indexBytes = 0;
    bool runComplete = false;
};

struct TetRecoveryContext {
    topo::WorkspaceStamp stamp{};
    std::uint64_t contextRevision = 0;
    mesh::Classification region{};
    std::vector<topo::NodeHandle> denseToNode;
    std::vector<TetFacetVertex> plainVertices;
    std::map<topo::NodeHandle, DenseTetNodeId> nodeToDense;
    std::map<topo::CellHandle, CurrentTetRecord> cells;
    std::map<DenseTetEdgeKey, EdgeResolverRow> edges;
    std::map<DenseTetFaceKey, FaceResolverRow> faces;
    std::vector<std::uint64_t> nodeUseCount;
    std::vector<std::vector<topo::CellHandle>> nodeStars;
    std::vector<TetFacetEdgeKey> protectedEdges;
    std::vector<TetFacetProtectedFace> protectedFaces;
    DynamicTetAabbIndex index;
    std::vector<topo::CellHandle> deltaCells;
    mutable TetRecoveryContextMetrics metrics{};
    bool poisoned = false;

    std::vector<topo::CellHandle> querySegment(const Vec3& a,
                                               const Vec3& b) const;
    std::vector<topo::CellHandle> queryTriangle(
        const std::array<Vec3, 3>& triangle) const;
    std::vector<topo::CellHandle> queryBounds(const TetAabb& bounds) const;
    bool currentFor(const topo::TopologyWorkspace& workspace) const noexcept;
};

struct TetRecoveryContextUpdateOutcome {
    TetRecoveryContextReason reason = TetRecoveryContextReason::InternalError;
    bool rebuilt = false;
    bool updated() const noexcept {
        return reason == TetRecoveryContextReason::None;
    }
};

struct TetRecoveryContextBuildOutcome {
    TetRecoveryContextReason reason = TetRecoveryContextReason::InternalError;
    std::unique_ptr<TetRecoveryContext> context;
    bool built() const noexcept {
        return reason == TetRecoveryContextReason::None && context != nullptr;
    }
};

static_assert(std::is_standard_layout_v<TetAabb>);
static_assert(std::is_trivially_copyable_v<TetAabb>);
static_assert(std::is_standard_layout_v<DenseTetEdgeKey>);
static_assert(std::is_trivially_copyable_v<DenseTetEdgeKey>);
static_assert(std::is_standard_layout_v<DenseTetFaceKey>);
static_assert(std::is_trivially_copyable_v<DenseTetFaceKey>);
static_assert(std::is_standard_layout_v<TetRecoveryContextMetrics>);
static_assert(std::is_trivially_copyable_v<TetRecoveryContextMetrics>);

TetRecoveryContextBuildOutcome buildTetRecoveryContext(
    const topo::TopologyWorkspace& workspace,
    mesh::Classification temporaryRegion) noexcept;

TetRecoveryContextUpdateOutcome updateTetRecoveryContext(
    TetRecoveryContext& context, const topo::TopologyWorkspace& workspace,
    std::span<const topo::CellHandle> retired,
    std::span<const topo::CellHandle> created) noexcept;

}  // namespace paacs::bnd
