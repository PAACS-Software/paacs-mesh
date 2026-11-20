// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// T3-A -- atomic staged-node Tet4 refinement replacement.
//
// Phase 1 owns only immutable validation/canonicalization. Phase 2 predicts and
// prepares one interior node plus the complete Tet cell/edge/face net delta.
// Phase 3 materializes both through prepared pool state, advances the workspace
// revision once, and performs no allocation or recoverable operation.

#include "paacs/predicates/Predicates.hpp"
#include "paacs/topo/SteinerCavity.hpp" // WorkspaceStamp + node-batch contact

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace paacs::topo {

enum class TetRefinementNodeRefKind : std::uint8_t { Existing, Staged };

struct TetRefinementNodeRef {
    TetRefinementNodeRefKind kind = TetRefinementNodeRefKind::Existing;
    NodeHandle existing{};
    std::uint32_t stagedLocalId = 0;

    static TetRefinementNodeRef existingNode(NodeHandle h) noexcept {
        return {TetRefinementNodeRefKind::Existing, h, 0};
    }
    static TetRefinementNodeRef stagedNode() noexcept {
        return {TetRefinementNodeRefKind::Staged, NodeHandle{}, 0};
    }
};

struct TetRefinementStagedNode {
    std::uint32_t localId = 0;
    Vec3 coordinate{};
    mesh::Classification classification{};
};

struct TetRefinementCellSpec {
    std::array<TetRefinementNodeRef, 4> vertices{};
    mesh::Classification classification{};
};

struct TetRefinementCavitySpec {
    WorkspaceStamp stamp{};
    TetRefinementStagedNode stagedNode{};
    std::vector<CellHandle> removed;
    std::vector<TetRefinementCellSpec> added;
    bool recoverTargetEdge = false;
    std::array<NodeHandle, 2> targetEdge{};
    mesh::Classification targetClassification{};
};

struct TetRefinementPreview {
    WorkspaceStamp sourceStamp{};
    std::uint32_t stagedLocalId = 0;
    mesh::Classification region{};
    std::uint64_t removedCellCount = 0;
    std::uint64_t addedCellCount = 0;
    std::uint64_t boundaryFaceCount = 0;
    std::uint64_t cavityDigest = 0;
    bool recoversTargetEdge = false;
    bool valid = false;
};

enum class TetRefinementTransition : std::uint8_t {
    Retained,
    Created,
    Retired
};

struct TetRefinementCellTransition {
    TetRefinementTransition transition = TetRefinementTransition::Retained;
    std::array<TetRefinementNodeRef, 4> key{};
    CellHandle predicted{};
    CellHandle actual{};
    mesh::Classification classification{};
};

struct TetRefinementEdgeTransition {
    TetRefinementTransition transition = TetRefinementTransition::Retained;
    std::array<TetRefinementNodeRef, 2> key{};
    EdgeHandle predicted{};
    EdgeHandle actual{};
    mesh::Classification classification{};
    bool required = false;
    std::uint32_t initialUseCount = 0;
    std::uint32_t finalUseCount = 0;
};

struct TetRefinementFaceTransition {
    TetRefinementTransition transition = TetRefinementTransition::Retained;
    std::array<TetRefinementNodeRef, 3> key{};
    FaceHandle predicted{};
    FaceHandle actual{};
    mesh::Classification classification{};
    std::array<CellHandle, 2> initialOwners{};
    std::array<std::int8_t, 2> initialSenses{};
    std::array<CellHandle, 2> finalOwners{};
    std::array<std::int8_t, 2> finalSenses{};
    std::uint8_t initialOwnerCount = 0;
    std::uint8_t finalOwnerCount = 0;
};

struct TetRefinementPrediction {
    NodeHandle stagedNode{};
    EdgeHandle targetEdge{};
    WorkspaceStamp postCommitStamp{};
    std::vector<TetRefinementCellTransition> cells;
    std::vector<TetRefinementEdgeTransition> edges;
    std::vector<TetRefinementFaceTransition> faces;
};

struct TetRefinementCommit {
    NodeHandle predictedNode{};
    NodeHandle actualNode{};
    EdgeHandle predictedTargetEdge{};
    EdgeHandle actualTargetEdge{};
    WorkspaceStamp stamp{};
    std::vector<TetRefinementCellTransition> cells;
    std::vector<TetRefinementEdgeTransition> edges;
    std::vector<TetRefinementFaceTransition> faces;
};

static_assert(std::is_nothrow_move_constructible_v<TetRefinementPrediction>);
static_assert(std::is_nothrow_move_constructible_v<TetRefinementCommit>);

class TetRefinementTransaction {
public:
    TetRefinementTransaction(const TetRefinementTransaction&) = delete;
    TetRefinementTransaction& operator=(const TetRefinementTransaction&) =
        delete;
    TetRefinementTransaction(TetRefinementTransaction&& other) noexcept
        : ws_(other.ws_), resume_(other.resume_), spec_(std::move(other.spec_)),
          preview_(other.preview_), concreteAdded_(std::move(other.concreteAdded_)),
          addedClasses_(std::move(other.addedClasses_)),
          nodeSpecs_(std::move(other.nodeSpecs_)),
          nodePrediction_(std::move(other.nodePrediction_)),
          plan_(std::move(other.plan_)), prediction_(std::move(other.prediction_)),
          resultCells_(std::move(other.resultCells_)),
          createdNodes_(std::move(other.createdNodes_)),
          poisoned_(other.poisoned_) {
        other.ws_ = nullptr;
        other.plan_.prepared = false;
    }
    TetRefinementTransaction& operator=(TetRefinementTransaction&&) = delete;
    ~TetRefinementTransaction() {
        if (ws_) abandonImpl();
    }

    const TetRefinementPreview& stage(const TetRefinementCavitySpec& input) {
        requireUsable("stage");
        if (spec_)
            throw std::logic_error(
                "TetRefinementTransaction: exactly one cavity may be staged");
        validateStamp(input.stamp);
        // The coordinator performed the complete immutable workspace audit
        // once and preserves it through audited local deltas. Repeating the
        // global pool audit here would make every insertion O(total Tets).
        // The plan below still validates every removed cell and all affected
        // edge/face incidence before any preparation.
        if (input.removed.empty() || input.added.empty())
            throw std::invalid_argument(
                "TetRefinementTransaction: cavity must be nonempty");
        if (input.stagedNode.localId != 0 ||
            input.stagedNode.classification.dim != 3 ||
            input.stagedNode.classification.tag < 0 ||
            !finite(input.stagedNode.coordinate))
            throw std::invalid_argument(
                "TetRefinementTransaction: invalid staged node");

        TetRefinementCavitySpec staged = input;
        const auto region = staged.stagedNode.classification;
        std::sort(staged.removed.begin(), staged.removed.end(),
                  [&](CellHandle a, CellHandle b) {
                      return ws_->cells_.seqOf(a) < ws_->cells_.seqOf(b);
                  });
        if (std::adjacent_find(staged.removed.begin(), staged.removed.end()) !=
            staged.removed.end())
            throw std::invalid_argument(
                "TetRefinementTransaction: duplicate removed cell");
        for (const auto cell : staged.removed) {
            if (!ws_->cells_.isCurrent(cell))
                throw std::invalid_argument(
                    "TetRefinementTransaction: stale removed cell");
            const auto snap = ws_->cells_.get(cell);
            if (snap.type != ElementType::Tet4 || !sameClass(snap.cls, region))
                throw std::invalid_argument(
                    "TetRefinementTransaction: removed region mismatch");
        }

        const auto symbolicRank = [&](const TetRefinementNodeRef& ref) {
            return ref.kind == TetRefinementNodeRefKind::Staged
                       ? (std::numeric_limits<std::uint64_t>::max)()
                       : ws_->nodes_.seqOf(ref.existing);
        };
        const auto coordinate = [&](const TetRefinementNodeRef& ref) {
            return ref.kind == TetRefinementNodeRefKind::Staged
                       ? staged.stagedNode.coordinate
                       : ws_->nodes_.get(ref.existing).coord;
        };
        const auto cellKey = [&](const TetRefinementCellSpec& cell) {
            std::array<std::uint64_t, 4> key{};
            for (std::size_t i = 0; i < 4; ++i)
                key[i] = symbolicRank(cell.vertices[i]);
            std::sort(key.begin(), key.end());
            return key;
        };

        bool stagedNodeUsed = false;
        for (auto& cell : staged.added) {
            if (!sameClass(cell.classification, region))
                throw std::invalid_argument(
                    "TetRefinementTransaction: candidate region mismatch");
            std::size_t stagedUses = 0;
            for (const auto& ref : cell.vertices) {
                if (ref.kind == TetRefinementNodeRefKind::Staged) {
                    if (ref.existing.isSpecified() || ref.stagedLocalId != 0)
                        throw std::invalid_argument(
                            "TetRefinementTransaction: malformed staged reference");
                    ++stagedUses;
                } else {
                    if (ref.stagedLocalId != 0 ||
                        !ws_->nodes_.isCurrent(ref.existing))
                        throw std::invalid_argument(
                            "TetRefinementTransaction: stale existing reference");
                }
            }
            if (stagedUses > 1 ||
                (!staged.recoverTargetEdge && stagedUses != 1))
                throw std::invalid_argument(
                    "TetRefinementTransaction: invalid staged node use");
            stagedNodeUsed = stagedNodeUsed || stagedUses == 1;
            const auto a = coordinate(cell.vertices[0]);
            const auto b = coordinate(cell.vertices[1]);
            const auto c = coordinate(cell.vertices[2]);
            const auto d = coordinate(cell.vertices[3]);
            const auto sign = predicates::orient3d(a, b, c, d);
            if (sign == predicates::Sign::Zero)
                throw std::invalid_argument(
                    "TetRefinementTransaction: degenerate candidate Tet");
            if (sign == predicates::Sign::Negative)
                std::swap(cell.vertices[0], cell.vertices[1]);
        }
        if (!stagedNodeUsed)
            throw std::invalid_argument(
                "TetRefinementTransaction: staged node would be unused");
        if (staged.recoverTargetEdge) {
            if (!ws_->nodes_.isCurrent(staged.targetEdge[0]) ||
                !ws_->nodes_.isCurrent(staged.targetEdge[1]) ||
                staged.targetEdge[0] == staged.targetEdge[1] ||
                (staged.targetClassification.dim != 1 &&
                 staged.targetClassification.dim != 2) ||
                staged.targetClassification.tag < 0)
                throw std::invalid_argument(
                    "TetRefinementTransaction: invalid target edge");
        }
        std::sort(staged.added.begin(), staged.added.end(),
                  [&](const auto& a, const auto& b) {
                      return cellKey(a) < cellKey(b);
                  });
        for (std::size_t i = 1; i < staged.added.size(); ++i)
            if (cellKey(staged.added[i - 1]) == cellKey(staged.added[i]))
                throw std::invalid_argument(
                    "TetRefinementTransaction: duplicate candidate Tet");

        preview_ = {};
        preview_.sourceStamp = staged.stamp;
        preview_.stagedLocalId = 0;
        preview_.region = region;
        preview_.removedCellCount = staged.removed.size();
        preview_.addedCellCount = staged.added.size();
        preview_.recoversTargetEdge = staged.recoverTargetEdge;
        preview_.cavityDigest = digest(staged, symbolicRank);
        preview_.valid = true;
        spec_.emplace(std::move(staged));
        return preview_;
    }

    const TetRefinementPrediction& prepare() {
        requireUsable("prepare");
        if (!spec_)
            throw std::logic_error("TetRefinementTransaction: nothing staged");
        if (plan_.prepared) return prediction_;
        try {
            using TA = detail::TransactionAccess;
            if (!TA::nodeBatchFits(ws_->nodes_, 1))
                throw std::length_error(
                    "TetRefinementTransaction: node allocator exhausted");
            nodePrediction_ = TA::predictNodeBatch(ws_->nodes_, 1);
            const auto& np = nodePrediction_.front();
            const NodeHandle stagedHandle{np.slot, np.gen, ws_->nodes_.poolId()};

            concreteAdded_.clear();
            addedClasses_.clear();
            concreteAdded_.reserve(spec_->added.size());
            addedClasses_.reserve(spec_->added.size());
            for (const auto& cell : spec_->added) {
                std::array<NodeHandle, 4> conn{};
                for (std::size_t i = 0; i < 4; ++i)
                    conn[i] = cell.vertices[i].kind ==
                                      TetRefinementNodeRefKind::Staged
                                  ? stagedHandle
                                  : cell.vertices[i].existing;
                concreteAdded_.push_back(conn);
                addedClasses_.push_back(cell.classification);
            }
            plan_ = spec_->recoverTargetEdge
                ? TA::planTetEdgeRefinement(
                      ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_,
                      spec_->removed, concreteAdded_, addedClasses_,
                      stagedHandle, np.seq, spec_->targetEdge,
                      spec_->targetClassification,
                      spec_->stagedNode.classification)
                : TA::planTetRefinement(
                      ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_,
                      spec_->removed, concreteAdded_, addedClasses_,
                      stagedHandle, np.seq,
                      spec_->stagedNode.classification);

            nodeSpecs_.clear();
            nodeSpecs_.push_back(TA::StagedNodeSpec{
                spec_->stagedNode.coordinate, spec_->stagedNode.classification});
            createdNodes_.reserve(1);
            resultCells_.reserve(plan_.adds.size());
            TA::prepareTetRefinement(ws_->nodes_, ws_->cells_, ws_->edges_,
                                     ws_->faces_, nodePrediction_, plan_);
            buildPrediction(stagedHandle);
        } catch (...) {
            poisoned_ = true;
            detail::TransactionAccess::abandonTetRefinementPrepared(
                ws_->edges_, ws_->faces_, plan_);
            throw;
        }
        ws_->state_ = WorkspaceState::TransactionPrepared;
        return prediction_;
    }

    const TetRefinementPrediction& prediction() const noexcept {
        if (!ws_ || !plan_.prepared) std::abort();
        return prediction_;
    }

    TetRefinementCommit commit() noexcept {
        if (!ws_ || !plan_.prepared || poisoned_) std::abort();
        TetRefinementCommit result;
        result.predictedNode = prediction_.stagedNode;
        result.predictedTargetEdge = prediction_.targetEdge;
        result.stamp = prediction_.postCommitStamp;
        result.cells = std::move(prediction_.cells);
        result.edges = std::move(prediction_.edges);
        result.faces = std::move(prediction_.faces);
        createdNodes_.clear();
        resultCells_.clear();
        detail::TransactionAccess::TetCommitBuffers buffers;
        buffers.resultCells = &resultCells_;
        EdgeHandle committedTargetEdge{};
        if (spec_->recoverTargetEdge)
            buffers.targetEdge = &committedTargetEdge;
#if defined(PAACS_TOPO_TESTING)
        detail::phase3Active() = true;
#endif
        detail::TransactionAccess::commitTetRefinementPrepared(
            ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_, nodeSpecs_,
            nodePrediction_, createdNodes_, plan_, buffers);
#if defined(PAACS_TOPO_TESTING)
        detail::phase3Active() = false;
#endif
        result.actualNode = createdNodes_.front();
        result.actualTargetEdge = committedTargetEdge;
        for (std::size_t i = 0; i < plan_.removes.size(); ++i)
            result.cells[i].actual = result.cells[i].predicted;
        for (std::size_t i = 0; i < plan_.adds.size(); ++i)
            result.cells[plan_.removes.size() + i].actual = resultCells_[i];
        for (std::size_t i = 0; i < plan_.edges.size(); ++i)
            result.edges[i].actual = plan_.edges[i].retire
                                         ? result.edges[i].predicted
                                         : plan_.edges[i].finalHandle;
        for (std::size_t i = 0; i < plan_.faces.size(); ++i)
            result.faces[i].actual = plan_.faces[i].retire
                                         ? result.faces[i].predicted
                                         : plan_.faces[i].finalHandle;
#if defined(PAACS_TOPO_TESTING)
        if (result.actualNode != result.predictedNode) {
            std::abort();
        }
        for (const auto& row : result.cells)
            if (row.actual != row.predicted) {
                std::abort();
            }
        for (const auto& row : result.edges)
            if (row.actual != row.predicted) {
                std::abort();
            }
        for (const auto& row : result.faces)
            if (row.actual != row.predicted) {
                std::abort();
            }
#endif
        ws_->revision_ += 1;
        ws_->state_ = WorkspaceState::CommittedAdjacencyStale;
        result.stamp.revision = ws_->revision_;
        ws_ = nullptr;
        return result;
    }

    void abandon() noexcept {
        if (ws_) abandonImpl();
    }

private:
    friend class TopologyWorkspace;
    explicit TetRefinementTransaction(TopologyWorkspace& ws)
        : ws_(&ws), resume_(ws.state_) {
        ws_->state_ = WorkspaceState::TransactionStaging;
    }

    static bool finite(const Vec3& p) noexcept {
        return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
    }
    static bool sameClass(const mesh::Classification& a,
                          const mesh::Classification& b) noexcept {
        return a.dim == b.dim && a.tag == b.tag;
    }
    void validateStamp(const WorkspaceStamp& stamp) const {
        const WorkspaceStamp current{ws_->nodes_.poolId(),
                                     ws_->cells_.poolId(),
                                     ws_->edges_.poolId(),
                                     ws_->faces_.poolId(),ws_->revision_};
        if (!(stamp == current))
            throw std::invalid_argument(
                "TetRefinementTransaction: stale or foreign workspace stamp");
    }
    void requireActive(const char* operation) const {
        if (!ws_)
            throw std::logic_error(std::string("TetRefinementTransaction: '") +
                                   operation + "' on a finished transaction");
    }
    void requireUsable(const char* operation) const {
        requireActive(operation);
        if (poisoned_)
            throw std::logic_error(std::string("TetRefinementTransaction: '") +
                                   operation + "' on a poisoned transaction");
    }

    template <class RankFn>
    static std::uint64_t digest(const TetRefinementCavitySpec& spec,
                                const RankFn& rank) noexcept {
        std::uint64_t h = 1469598103934665603ull;
        auto add = [&](std::uint64_t value) {
            h ^= value;
            h *= 1099511628211ull;
        };
        add(spec.removed.size());
        add(spec.added.size());
        if (spec.recoverTargetEdge) {
            add(1U);
            add(rank(TetRefinementNodeRef::existingNode(
                spec.targetEdge[0])));
            add(rank(TetRefinementNodeRef::existingNode(
                spec.targetEdge[1])));
            add(static_cast<std::uint64_t>(
                static_cast<std::int64_t>(
                    spec.targetClassification.dim)));
            add(static_cast<std::uint64_t>(
                static_cast<std::int64_t>(
                    spec.targetClassification.tag)));
        }
        for (const auto& cell : spec.added)
            for (const auto& ref : cell.vertices) add(rank(ref));
        return h;
    }

    TetRefinementNodeRef refFor(NodeHandle handle,
                                NodeHandle staged) const noexcept {
        return handle == staged ? TetRefinementNodeRef::stagedNode()
                                : TetRefinementNodeRef::existingNode(handle);
    }

    template <std::size_t N>
    void sortRefs(std::array<TetRefinementNodeRef, N>& refs,
                  NodeHandle staged) const {
        const auto rank = [&](const TetRefinementNodeRef& ref) {
            return ref.kind == TetRefinementNodeRefKind::Staged
                       ? nodePrediction_.front().seq
                       : ws_->nodes_.seqOf(ref.existing);
        };
        (void)staged;
        std::sort(refs.begin(), refs.end(),
                  [&](const auto& a, const auto& b) { return rank(a) < rank(b); });
    }

    void buildPrediction(NodeHandle staged) {
        prediction_ = {};
        prediction_.stagedNode = staged;
        if (spec_->recoverTargetEdge)
            prediction_.targetEdge =
                detail::TransactionAccess::predictedTetTargetEdge(
                    ws_->edges_, plan_);
        prediction_.postCommitStamp = WorkspaceStamp{
            ws_->nodes_.poolId(),ws_->cells_.poolId(),ws_->edges_.poolId(),
            ws_->faces_.poolId(),ws_->revision_};
        prediction_.postCommitStamp.revision += 1;
        prediction_.cells.reserve(plan_.removes.size() + plan_.adds.size());
        prediction_.edges.reserve(plan_.edges.size());
        prediction_.faces.reserve(plan_.faces.size());

        for (const auto& remove : plan_.removes) {
            TetRefinementCellTransition row;
            row.transition = TetRefinementTransition::Retired;
            row.predicted = remove.cell;
            const auto cell = ws_->cells_.get(remove.cell);
            row.classification = cell.cls;
            for (std::size_t i = 0; i < 4; ++i)
                row.key[i] = TetRefinementNodeRef::existingNode(cell.nodes[i]);
            sortRefs(row.key, staged);
            prediction_.cells.push_back(row);
        }
        for (const auto& add : plan_.adds) {
            TetRefinementCellTransition row;
            row.transition = TetRefinementTransition::Created;
            row.predicted = CellHandle{add.predictedSlot, add.predictedGen,
                                       ws_->cells_.poolId()};
            row.classification = add.cls;
            for (std::size_t i = 0; i < 4; ++i)
                row.key[i] = refFor(add.nodes[i], staged);
            sortRefs(row.key, staged);
            prediction_.cells.push_back(row);
        }

        for (const auto& edge : plan_.edges) {
            TetRefinementEdgeTransition row;
            row.transition = edge.create
                                 ? TetRefinementTransition::Created
                                 : (edge.retire ? TetRefinementTransition::Retired
                                                : TetRefinementTransition::Retained);
            row.key = {refFor(edge.a, staged), refFor(edge.b, staged)};
            sortRefs(row.key, staged);
            row.predicted = edge.create
                                ? EdgeHandle{edge.predictedSlot, edge.predictedGen,
                                             ws_->edges_.poolId()}
                                : edge.existing;
            row.initialUseCount = edge.initialUse;
            row.finalUseCount = edge.finalUse;
            if (edge.existing.isSpecified()) {
                const auto data = ws_->edges_.get(edge.existing);
                row.classification = data.cls;
                row.required = data.requiredConstraint;
            }
            prediction_.edges.push_back(row);
        }

        for (std::size_t fi = 0; fi < plan_.faces.size(); ++fi) {
            const auto& face = plan_.faces[fi];
            TetRefinementFaceTransition row;
            row.transition = face.create
                                 ? TetRefinementTransition::Created
                                 : (face.retire ? TetRefinementTransition::Retired
                                                : TetRefinementTransition::Retained);
            for (std::size_t i = 0; i < 3; ++i)
                row.key[i] = refFor(face.staged.outward[i], staged);
            sortRefs(row.key, staged);
            row.predicted = face.create
                                ? FaceHandle{face.predictedSlot, face.predictedGen,
                                             ws_->faces_.poolId()}
                                : face.existing;
            if (face.existing.isSpecified()) {
                const auto data = ws_->faces_.get(face.existing);
                row.classification = data.cls;
                row.initialOwnerCount = static_cast<std::uint8_t>(data.incident.size());
                for (std::size_t i = 0; i < data.incident.size(); ++i) {
                    row.initialOwners[i] = CellHandle{
                        data.incident[i].cellSlot, data.incident[i].cellGen,
                        ws_->cells_.poolId()};
                    row.initialSenses[i] = data.incident[i].sense;
                }
            }
            std::size_t owner = 0;
            for (const auto& surviving : face.surviving) {
                row.finalOwners[owner] = CellHandle{
                    surviving.cellSlot, surviving.cellGen, ws_->cells_.poolId()};
                row.finalSenses[owner++] = surviving.sense;
            }
            for (std::size_t ai = 0; ai < plan_.adds.size(); ++ai)
                for (std::size_t lf = 0; lf < 4; ++lf)
                    if (plan_.adds[ai].faces[lf] == fi) {
                        row.finalOwners[owner] = CellHandle{
                            plan_.adds[ai].predictedSlot,
                            plan_.adds[ai].predictedGen, ws_->cells_.poolId()};
                        row.finalSenses[owner++] = plan_.adds[ai].faceSenses[lf];
                    }
            row.finalOwnerCount = static_cast<std::uint8_t>(owner);
            prediction_.faces.push_back(row);
        }
    }

    void abandonImpl() noexcept {
        detail::TransactionAccess::abandonTetRefinementPrepared(
            ws_->edges_, ws_->faces_, plan_);
        ws_->state_ = resume_;
        ws_ = nullptr;
    }

    TopologyWorkspace* ws_ = nullptr;
    WorkspaceState resume_ = WorkspaceState::Ready;
    std::optional<TetRefinementCavitySpec> spec_;
    TetRefinementPreview preview_{};
    std::vector<std::array<NodeHandle, 4>> concreteAdded_;
    std::vector<mesh::Classification> addedClasses_;
    std::vector<detail::TransactionAccess::StagedNodeSpec> nodeSpecs_;
    std::vector<detail::TransactionAccess::StagedNodePredict> nodePrediction_;
    detail::TransactionAccess::TetPlan plan_;
    TetRefinementPrediction prediction_;
    std::vector<CellHandle> resultCells_;
    std::vector<NodeHandle> createdNodes_;
    bool poisoned_ = false;
};

inline TetRefinementTransaction
TopologyWorkspace::beginTetRefinementTransaction() {
    if (!readable())
        throw std::logic_error(
            "TopologyWorkspace::beginTetRefinementTransaction: a transaction "
            "is already active");
    return TetRefinementTransaction(*this);
}

} // namespace paacs::topo
