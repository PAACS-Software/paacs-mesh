// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// T2-A -- atomic Tet4 cavity replacement.
//
// This additive transaction is the implementation boundary accepted by
// T2A_DESIGN.md revision 2.  It owns the complete cell/edge/face net delta:
// Phase 1 validates and canonicalizes using existing node creation sequence,
// Phase 2 prepares every destination slot/table/incidence buffer, and Phase 3
// consumes only prepared state.  The three target nodes already exist; this
// primitive never creates nodes.

#include "paacs/predicates/Predicates.hpp"
#include "paacs/topo/SteinerCavity.hpp"  // WorkspaceStamp + workspace surface

#include <algorithm>
#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace paacs::topo {

struct TetCavityCell {
    std::array<NodeHandle, 4> nodes{};
    mesh::Classification cls{};
};

struct TetCavitySpec {
    WorkspaceStamp stamp{};
    std::vector<CellHandle> removed;
    std::vector<TetCavityCell> added;
    std::array<NodeHandle, 3> target{};
    mesh::Classification targetClassification{};
    mesh::Classification region{};
};

struct TetEdgeCavitySpec {
    WorkspaceStamp stamp{};
    std::vector<CellHandle> removed;
    std::vector<TetCavityCell> added;
    std::array<NodeHandle, 2> target{};
    mesh::Classification targetClassification{};
    bool requiredConstraint = true;
    mesh::Classification region{};
};

struct TetEdgeMetadataSpec {
    WorkspaceStamp stamp{};
    EdgeHandle target{};
    mesh::Classification targetClassification{};
    bool requiredConstraint = true;
};

// Plain values only.  No validator can observe a workspace or a half-commit.
struct TetCavityPreview {
    std::vector<CellHandle> removedCells;
    std::vector<std::array<NodeHandle, 4>> addedConnectivity;
    std::vector<mesh::Classification> addedClassifications;
    std::array<NodeHandle, 3> target{};
    std::size_t cellsToCreate = 0;
    std::size_t cellsToRetire = 0;
    std::size_t edgesToCreate = 0;
    std::size_t edgesToRetain = 0;
    std::size_t edgesToRetire = 0;
    std::size_t facesToCreate = 0;
    std::size_t facesToRetain = 0;
    std::size_t facesToRetire = 0;
};

struct TetEdgeCavityPreview {
    std::vector<CellHandle> removedCells;
    std::vector<std::array<NodeHandle, 4>> addedConnectivity;
    std::vector<mesh::Classification> addedClassifications;
    std::array<NodeHandle, 2> target{};
    std::size_t cellsToCreate = 0;
    std::size_t cellsToRetire = 0;
    std::size_t edgesToCreate = 0;
    std::size_t edgesToRetain = 0;
    std::size_t edgesToRetire = 0;
    std::size_t facesToCreate = 0;
    std::size_t facesToRetain = 0;
    std::size_t facesToRetire = 0;
};

struct TetCavityCommit {
    std::vector<CellHandle> resultCells;
    FaceHandle targetFace{};
    std::uint64_t revision = 0;
};

struct TetEdgeCavityCommit {
    std::vector<CellHandle> resultCells;
    EdgeHandle targetEdge{};
    std::uint64_t revision = 0;
};

static_assert(std::is_nothrow_move_constructible_v<TetCavityCommit>);
static_assert(std::is_nothrow_move_constructible_v<TetEdgeCavityCommit>);

class TetCavityTransaction {
public:
    TetCavityTransaction(const TetCavityTransaction&) = delete;
    TetCavityTransaction& operator=(const TetCavityTransaction&) = delete;
    TetCavityTransaction(TetCavityTransaction&& other) noexcept
        : ws_(other.ws_), resume_(other.resume_), spec_(std::move(other.spec_)),
          edgeSpec_(std::move(other.edgeSpec_)),
          edgeMetadataSpec_(std::move(other.edgeMetadataSpec_)),
          preview_(std::move(other.preview_)),
          edgePreview_(std::move(other.edgePreview_)),
          plan_(std::move(other.plan_)),
          edgeMetadataPlan_(std::move(other.edgeMetadataPlan_)),
          resultCells_(std::move(other.resultCells_)),
          predictedTarget_(other.predictedTarget_),
          predictedTargetEdge_(other.predictedTargetEdge_),
          poisoned_(other.poisoned_) {
        other.ws_ = nullptr;
        other.plan_.prepared = false;
    }
    TetCavityTransaction& operator=(TetCavityTransaction&&) = delete;
    ~TetCavityTransaction() {
        if (ws_) abandonImpl();
    }

    const TetCavityPreview& stage(const TetCavitySpec& input) {
        requireUsable("stage");
        if (spec_ || edgeSpec_ || edgeMetadataSpec_)
            throw std::logic_error(
                "TetCavityTransaction: exactly one cavity may be staged");
        validateStamp(input.stamp);
        if (!workspaceTopologyReady(ws_->nodes_, ws_->cells_, ws_->edges_,
                                    ws_->faces_))
            throw std::invalid_argument(
                "TetCavityTransaction: workspace topology is not ready");
        if (input.removed.empty() || input.added.empty())
            throw std::invalid_argument(
                "TetCavityTransaction: cavity must be nonempty");
        if (input.targetClassification.dim != 2 ||
            input.targetClassification.tag < 0 || input.region.dim != 3 ||
            input.region.tag < 0)
            throw std::invalid_argument(
                "TetCavityTransaction: invalid classification domain");

        TetCavitySpec staged = input;
        auto cellLess = [&](CellHandle a, CellHandle b) {
            return ws_->cells_.seqOf(a) < ws_->cells_.seqOf(b);
        };
        std::sort(staged.removed.begin(), staged.removed.end(), cellLess);
        if (std::adjacent_find(staged.removed.begin(), staged.removed.end()) !=
            staged.removed.end())
            throw std::invalid_argument(
                "TetCavityTransaction: duplicate removed cell");

        auto nodeKey = [&](const TetCavityCell& cell) {
            std::array<std::uint64_t, 4> key{};
            for (std::size_t i = 0; i < 4; ++i)
                key[i] = ws_->nodes_.seqOf(cell.nodes[i]);
            std::sort(key.begin(), key.end());
            return key;
        };
        std::sort(staged.added.begin(), staged.added.end(),
                  [&](const TetCavityCell& a, const TetCavityCell& b) {
                      return nodeKey(a) < nodeKey(b);
                  });

        for (const auto removed : staged.removed) {
            if (!ws_->cells_.isCurrent(removed))
                throw std::invalid_argument(
                    "TetCavityTransaction: stale removed cell");
            const auto cell = ws_->cells_.get(removed);
            if (cell.type != ElementType::Tet4 ||
                !sameClass(cell.cls, staged.region))
                throw std::invalid_argument(
                    "TetCavityTransaction: removed cell region mismatch");
        }
        for (const auto& cell : staged.added) {
            if (!sameClass(cell.cls, staged.region))
                throw std::invalid_argument(
                    "TetCavityTransaction: candidate region mismatch");
            for (const auto node : cell.nodes)
                if (!ws_->nodes_.isCurrent(node))
                    throw std::invalid_argument(
                        "TetCavityTransaction: stale candidate node");
            const auto a = ws_->nodes_.get(cell.nodes[0]).coord;
            const auto b = ws_->nodes_.get(cell.nodes[1]).coord;
            const auto c = ws_->nodes_.get(cell.nodes[2]).coord;
            const auto d = ws_->nodes_.get(cell.nodes[3]).coord;
            if (predicates::orient3d(a, b, c, d) !=
                predicates::Sign::Positive)
                throw std::invalid_argument(
                    "TetCavityTransaction: candidate Tet is not exact-positive");
        }
        for (const auto targetNode : staged.target)
            if (!ws_->nodes_.isCurrent(targetNode))
                throw std::invalid_argument(
                    "TetCavityTransaction: stale target node");
        if (staged.target[0] == staged.target[1] ||
            staged.target[1] == staged.target[2] ||
            staged.target[0] == staged.target[2])
            throw std::invalid_argument(
                "TetCavityTransaction: degenerate target identity");
        {
            const auto a = ws_->nodes_.get(staged.target[0]).coord;
            const auto b = ws_->nodes_.get(staged.target[1]).coord;
            const auto c = ws_->nodes_.get(staged.target[2]).coord;
            if (predicates::orient2d({a.x, a.y}, {b.x, b.y}, {c.x, c.y}) ==
                    predicates::Sign::Zero &&
                predicates::orient2d({a.y, a.z}, {b.y, b.z}, {c.y, c.z}) ==
                    predicates::Sign::Zero &&
                predicates::orient2d({a.z, a.x}, {b.z, b.x}, {c.z, c.x}) ==
                    predicates::Sign::Zero)
                throw std::invalid_argument(
                    "TetCavityTransaction: target nodes are collinear");
        }

        std::vector<std::array<NodeHandle, 4>> connectivity;
        std::vector<mesh::Classification> classifications;
        connectivity.reserve(staged.added.size());
        classifications.reserve(staged.added.size());
        for (const auto& cell : staged.added) {
            connectivity.push_back(cell.nodes);
            classifications.push_back(cell.cls);
        }

        auto planned = detail::TransactionAccess::planTetCavity(
            ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_,
            staged.removed, connectivity, classifications, staged.target,
            staged.targetClassification, staged.region);

        TetCavityPreview preview;
        preview.removedCells = staged.removed;
        preview.addedConnectivity = std::move(connectivity);
        preview.addedClassifications = std::move(classifications);
        preview.target = staged.target;
        preview.cellsToCreate = planned.adds.size();
        preview.cellsToRetire = planned.removes.size();
        for (const auto& edge : planned.edges) {
            preview.edgesToCreate += edge.create ? 1U : 0U;
            preview.edgesToRetire += edge.retire ? 1U : 0U;
            preview.edgesToRetain +=
                (!edge.create && !edge.retire && edge.finalUse > 0) ? 1U : 0U;
        }
        for (const auto& face : planned.faces) {
            preview.facesToCreate += face.create ? 1U : 0U;
            preview.facesToRetire += face.retire ? 1U : 0U;
            preview.facesToRetain +=
                (!face.create && !face.retire && face.finalUse > 0) ? 1U : 0U;
        }
        spec_.emplace(std::move(staged));
        preview_.emplace(std::move(preview));
        plan_ = std::move(planned);
        return *preview_;
    }

    const TetEdgeCavityPreview& stage(const TetEdgeCavitySpec& input) {
        requireUsable("stage(edge)");
        if (spec_ || edgeSpec_ || edgeMetadataSpec_)
            throw std::logic_error(
                "TetCavityTransaction: exactly one cavity may be staged");
        validateStamp(input.stamp);
        if (!workspaceTopologyReady(ws_->nodes_, ws_->cells_, ws_->edges_,
                                    ws_->faces_))
            throw std::invalid_argument(
                "TetCavityTransaction: workspace topology is not ready");
        if (input.removed.empty() || input.added.empty())
            throw std::invalid_argument(
                "TetCavityTransaction: cavity must be nonempty");
        if ((input.targetClassification.dim != 1 &&
             input.targetClassification.dim != 2) ||
            input.targetClassification.tag < 0 ||
            !input.requiredConstraint || input.region.dim != 3 ||
            input.region.tag < 0)
            throw std::invalid_argument(
                "TetCavityTransaction: invalid edge classification domain");

        TetEdgeCavitySpec staged = input;
        auto cellLess = [&](CellHandle a, CellHandle b) {
            return ws_->cells_.seqOf(a) < ws_->cells_.seqOf(b);
        };
        std::sort(staged.removed.begin(), staged.removed.end(), cellLess);
        if (std::adjacent_find(staged.removed.begin(), staged.removed.end()) !=
            staged.removed.end())
            throw std::invalid_argument(
                "TetCavityTransaction: duplicate removed cell");
        auto nodeKey = [&](const TetCavityCell& cell) {
            std::array<std::uint64_t, 4> key{};
            for (std::size_t i = 0; i < 4; ++i)
                key[i] = ws_->nodes_.seqOf(cell.nodes[i]);
            std::sort(key.begin(), key.end());
            return key;
        };
        std::sort(staged.added.begin(), staged.added.end(),
                  [&](const TetCavityCell& a, const TetCavityCell& b) {
                      return nodeKey(a) < nodeKey(b);
                  });
        for (const auto removed : staged.removed) {
            if (!ws_->cells_.isCurrent(removed))
                throw std::invalid_argument(
                    "TetCavityTransaction: stale removed cell");
            const auto cell = ws_->cells_.get(removed);
            if (cell.type != ElementType::Tet4 ||
                !sameClass(cell.cls, staged.region))
                throw std::invalid_argument(
                    "TetCavityTransaction: removed cell region mismatch");
        }
        for (const auto& cell : staged.added) {
            if (!sameClass(cell.cls, staged.region))
                throw std::invalid_argument(
                    "TetCavityTransaction: candidate region mismatch");
            for (const auto node : cell.nodes)
                if (!ws_->nodes_.isCurrent(node))
                    throw std::invalid_argument(
                        "TetCavityTransaction: stale candidate node");
            if (predicates::orient3d(
                    ws_->nodes_.get(cell.nodes[0]).coord,
                    ws_->nodes_.get(cell.nodes[1]).coord,
                    ws_->nodes_.get(cell.nodes[2]).coord,
                    ws_->nodes_.get(cell.nodes[3]).coord) !=
                predicates::Sign::Positive)
                throw std::invalid_argument(
                    "TetCavityTransaction: candidate Tet is not exact-positive");
        }
        if (!ws_->nodes_.isCurrent(staged.target[0]) ||
            !ws_->nodes_.isCurrent(staged.target[1]) ||
            staged.target[0] == staged.target[1])
            throw std::invalid_argument(
                "TetCavityTransaction: invalid target edge identity");
        const Vec3 a = ws_->nodes_.get(staged.target[0]).coord;
        const Vec3 b = ws_->nodes_.get(staged.target[1]).coord;
        if (a.x == b.x && a.y == b.y && a.z == b.z)
            throw std::invalid_argument(
                "TetCavityTransaction: zero-length target edge");

        std::vector<std::array<NodeHandle, 4>> connectivity;
        std::vector<mesh::Classification> classifications;
        connectivity.reserve(staged.added.size());
        classifications.reserve(staged.added.size());
        for (const auto& cell : staged.added) {
            connectivity.push_back(cell.nodes);
            classifications.push_back(cell.cls);
        }
        auto planned = detail::TransactionAccess::planTetEdgeCavity(
            ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_,
            staged.removed, connectivity, classifications, staged.target,
            staged.targetClassification, staged.region);
        TetEdgeCavityPreview preview;
        preview.removedCells = staged.removed;
        preview.addedConnectivity = std::move(connectivity);
        preview.addedClassifications = std::move(classifications);
        preview.target = staged.target;
        preview.cellsToCreate = planned.adds.size();
        preview.cellsToRetire = planned.removes.size();
        for (const auto& edge : planned.edges) {
            preview.edgesToCreate += edge.create ? 1U : 0U;
            preview.edgesToRetire += edge.retire ? 1U : 0U;
            preview.edgesToRetain +=
                (!edge.create && !edge.retire && edge.finalUse > 0) ? 1U : 0U;
        }
        for (const auto& face : planned.faces) {
            preview.facesToCreate += face.create ? 1U : 0U;
            preview.facesToRetire += face.retire ? 1U : 0U;
            preview.facesToRetain +=
                (!face.create && !face.retire && face.finalUse > 0) ? 1U : 0U;
        }
        edgeSpec_.emplace(std::move(staged));
        edgePreview_.emplace(std::move(preview));
        plan_ = std::move(planned);
        return *edgePreview_;
    }

    const TetEdgeCavityPreview& stage(const TetEdgeMetadataSpec& input) {
        requireUsable("stage(edge metadata)");
        if (spec_ || edgeSpec_ || edgeMetadataSpec_)
            throw std::logic_error(
                "TetCavityTransaction: exactly one cavity may be staged");
        validateStamp(input.stamp);
        if (!workspaceTopologyReady(ws_->nodes_, ws_->cells_, ws_->edges_,
                                    ws_->faces_))
            throw std::invalid_argument(
                "TetCavityTransaction: workspace topology is not ready");
        auto planned = detail::TransactionAccess::planTetEdgeMetadata(
            ws_->edges_,input.target,input.targetClassification,
            input.requiredConstraint);
        TetEdgeCavityPreview preview;
        const auto data = ws_->edges_.get(input.target);
        preview.target = {data.a,data.b};
        preview.edgesToRetain = 1;
        edgeMetadataSpec_.emplace(input);
        edgePreview_.emplace(std::move(preview));
        edgeMetadataPlan_ = std::move(planned);
        return *edgePreview_;
    }

    const TetCavityPreview& preview() const {
        requireActive("preview");
        if (!preview_)
            throw std::logic_error("TetCavityTransaction: nothing staged");
        return *preview_;
    }

    const TetEdgeCavityPreview& edgePreview() const {
        requireActive("edgePreview");
        if (!edgePreview_)
            throw std::logic_error("TetCavityTransaction: no edge staged");
        return *edgePreview_;
    }

    void prepare() {
        requireUsable("prepare");
        if (!spec_ && !edgeSpec_ && !edgeMetadataSpec_)
            throw std::logic_error("TetCavityTransaction: nothing staged");
        if (edgeMetadataSpec_) {
            if (edgeMetadataPlan_.prepared) return;
            try {
                detail::TransactionAccess::prepareTetEdgeMetadata(
                    edgeMetadataPlan_);
                predictedTargetEdge_ = edgeMetadataSpec_->target;
            } catch (...) {
                poisoned_ = true;
                throw;
            }
            ws_->state_ = WorkspaceState::TransactionPrepared;
            return;
        }
        if (plan_.prepared) return;
        try {
            // Result append is part of Phase 3 and therefore reserved before
            // the first pool placeholder is placed.
            resultCells_.reserve(plan_.adds.size());
            detail::TransactionAccess::prepareTetCavity(
                ws_->cells_, ws_->edges_, ws_->faces_, plan_);
            if (edgeSpec_)
                predictedTargetEdge_ =
                    detail::TransactionAccess::predictedTetTargetEdge(
                        ws_->edges_, plan_);
            else
                predictedTarget_ =
                    detail::TransactionAccess::predictedTetTargetFace(
                        ws_->faces_, plan_);
        } catch (...) {
            poisoned_ = true;
            detail::TransactionAccess::abandonTetCavity(
                ws_->edges_, ws_->faces_, plan_);
            throw;
        }
        ws_->state_ = WorkspaceState::TransactionPrepared;
    }

    bool prepared() const noexcept {
        return edgeMetadataSpec_ ? edgeMetadataPlan_.prepared : plan_.prepared;
    }
    FaceHandle predictedTargetFace() const {
        requireActive("predictedTargetFace");
        if (edgeSpec_ || edgeMetadataSpec_ || !plan_.prepared)
            throw std::logic_error(
                "TetCavityTransaction: target is predicted after preparation");
        return predictedTarget_;
    }

    EdgeHandle predictedTargetEdge() const {
        requireActive("predictedTargetEdge");
        if ((!edgeSpec_ && !edgeMetadataSpec_) || !prepared())
            throw std::logic_error(
                "TetCavityTransaction: target edge is predicted after edge preparation");
        return predictedTargetEdge_;
    }

    TetCavityCommit commit() {
        requireUsable("commit");
        if (!spec_ || edgeSpec_ || edgeMetadataSpec_)
            throw std::logic_error("TetCavityTransaction: nothing staged");
        if (!plan_.prepared) prepare();
        resultCells_.clear();
        FaceHandle materializedTarget{};
        detail::TransactionAccess::TetCommitBuffers buffers;
        buffers.resultCells = &resultCells_;
        buffers.targetFace = &materializedTarget;
#if defined(PAACS_TOPO_TESTING)
        detail::phase3Active() = true;
#endif
        detail::TransactionAccess::commitTetCavity(
            ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_, plan_,
            buffers);
#if defined(PAACS_TOPO_TESTING)
        detail::phase3Active() = false;
        if (materializedTarget != predictedTarget_) std::abort();
#endif
        ws_->revision_ += 1;
        ws_->state_ = WorkspaceState::CommittedAdjacencyStale;
        TetCavityCommit result;
        result.resultCells = std::move(resultCells_);
        result.targetFace = materializedTarget;
        result.revision = ws_->revision_;
        ws_ = nullptr;
        return result;
    }

    TetEdgeCavityCommit commitEdge() {
        requireUsable("commitEdge");
        if ((!edgeSpec_ && !edgeMetadataSpec_) || spec_)
            throw std::logic_error("TetCavityTransaction: no edge staged");
        if (!prepared()) prepare();
        if (edgeMetadataSpec_) {
#if defined(PAACS_TOPO_TESTING)
            detail::phase3Active() = true;
#endif
            detail::TransactionAccess::commitTetEdgeMetadata(
                ws_->edges_,edgeMetadataPlan_);
#if defined(PAACS_TOPO_TESTING)
            detail::phase3Active() = false;
#endif
            ws_->revision_ += 1;
            ws_->state_ = WorkspaceState::CommittedAdjacencyStale;
            TetEdgeCavityCommit result;
            result.targetEdge = predictedTargetEdge_;
            result.revision = ws_->revision_;
            ws_ = nullptr;
            return result;
        }
        resultCells_.clear();
        EdgeHandle materializedTarget{};
        detail::TransactionAccess::TetCommitBuffers buffers;
        buffers.resultCells = &resultCells_;
        buffers.targetEdge = &materializedTarget;
#if defined(PAACS_TOPO_TESTING)
        detail::phase3Active() = true;
#endif
        detail::TransactionAccess::commitTetCavity(
            ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_, plan_,
            buffers);
#if defined(PAACS_TOPO_TESTING)
        detail::phase3Active() = false;
        if (materializedTarget != predictedTargetEdge_) std::abort();
#endif
        ws_->revision_ += 1;
        ws_->state_ = WorkspaceState::CommittedAdjacencyStale;
        TetEdgeCavityCommit result;
        result.resultCells = std::move(resultCells_);
        result.targetEdge = materializedTarget;
        result.revision = ws_->revision_;
        ws_ = nullptr;
        return result;
    }

    void abandon() {
        requireActive("abandon");
        abandonImpl();
    }

private:
    friend class TopologyWorkspace;
    explicit TetCavityTransaction(TopologyWorkspace& ws)
        : ws_(&ws), resume_(ws.state_) {
        ws_->state_ = WorkspaceState::TransactionStaging;
    }

    static bool sameClass(const mesh::Classification& a,
                          const mesh::Classification& b) noexcept {
        return a.dim == b.dim && a.tag == b.tag;
    }
    void validateStamp(const WorkspaceStamp& stamp) const {
        const WorkspaceStamp current{ws_->nodes_.poolId(),
                                     ws_->cells_.poolId(),
                                     ws_->edges_.poolId(),
                                     ws_->faces_.poolId(), ws_->revision_};
        if (!(stamp == current))
            throw std::invalid_argument(
                "TetCavityTransaction: stale or foreign workspace stamp");
    }
    void requireActive(const char* operation) const {
        if (!ws_)
            throw std::logic_error(std::string("TetCavityTransaction: '") +
                                   operation + "' on a finished transaction");
    }
    void requireUsable(const char* operation) const {
        requireActive(operation);
        if (poisoned_)
            throw std::logic_error(std::string("TetCavityTransaction: '") +
                                   operation +
                                   "' on a prepare-poisoned transaction");
    }
    void abandonImpl() noexcept {
        detail::TransactionAccess::abandonTetCavity(ws_->edges_, ws_->faces_,
                                                    plan_);
        edgeMetadataPlan_.prepared = false;
        ws_->state_ = resume_;
        ws_ = nullptr;
    }

    TopologyWorkspace* ws_ = nullptr;
    WorkspaceState resume_ = WorkspaceState::Ready;
    std::optional<TetCavitySpec> spec_;
    std::optional<TetEdgeCavitySpec> edgeSpec_;
    std::optional<TetEdgeMetadataSpec> edgeMetadataSpec_;
    std::optional<TetCavityPreview> preview_;
    std::optional<TetEdgeCavityPreview> edgePreview_;
    detail::TransactionAccess::TetPlan plan_;
    detail::TransactionAccess::TetEdgeMetadataPlan edgeMetadataPlan_;
    std::vector<CellHandle> resultCells_;
    FaceHandle predictedTarget_{};
    EdgeHandle predictedTargetEdge_{};
    bool poisoned_ = false;
};

inline TetCavityTransaction TopologyWorkspace::beginTetCavityTransaction() {
    if (!readable())
        throw std::logic_error(
            "TopologyWorkspace::beginTetCavityTransaction: a transaction is "
            "already active");
    return TetCavityTransaction(*this);
}

}  // namespace paacs::topo
