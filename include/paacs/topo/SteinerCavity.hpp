// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// BND-2A / ASCR - Atomic Steiner Cavity Replacement (Increment 2 primitive).
// Design authority: docs/PMBND002_DESIGN.md rev 3 (accepted).
//
// This NEW, unfrozen header adds a SECOND cavity transaction that differs
// from the frozen M4-B.8 CavityTransaction in exactly one way: it may CREATE
// a batch of interior Steiner nodes and commit them together with the
// replacement triangles ATOMICALLY. The frozen cavity never creates a node
// and enforces nodes(M) subset of nodes(R) with |M| == |R|; a Steiner
// replacement has nodes(M) = nodes(R) union {staged} and |M| > |R|, so it
// cannot reuse the frozen buildPreview body. Instead it replicates EVERY
// combinatorial protection over the additive rank-based cavaudit helpers
// (disk AND node-link for BOTH R and M, boundary preservation with inside-
// sense match, orientation-insensitive duplicate rejection, surviving-cell
// duplicate rejection, closure-only ensure-only requirements, required-edge
// preservation, orphan-free node accounting) and drives the additive
// TransactionAccess batched-node + staged-cavity helpers. It reuses the
// FROZEN commitCavityPrepared unchanged (which already supports
// |added| != |removed| through its claims->frees->reuse ordering); by the
// time that runs, the staged nodes have been materialized in the SAME
// noexcept Phase 3, so their handles are live.
//
// IDENTITY MODEL (rev-3, review round 2). Staged nodes carry an explicit
// structural localId; every reference (cell corner, requirement endpoint) is
// a tagged Existing|Staged(localId), NEVER a container index. Phase-1
// identity is SYMBOLIC and allocator-independent: closure nodes are mapped to
// dense order-preserving COMPACT ranks (existing nodes by live sequence,
// staged nodes after them by localId), so validation cannot depend on caller
// container order and cannot construct wrapped ranks. Real slot/gen/sequence
// prediction happens ONLY in Phase 2 through the single predictSlots batch
// path, after the full node batch demand (sequence AND slot space) is
// verified. The commit returns an explicit localId -> NodeHandle mapping.
//
// SCOPE (Increment 2): the atomic primitive only. It does NOT choose Steiner
// positions, recover PSLG segments, integrate FaceMesher, publish to
// mesh::Mesh, or create boundary Steiner nodes - those are later increments.

#include "paacs/topo/Cavity.hpp"   // Workspace, cavaudit, TransactionAccess,
                                   // value-type deps

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace paacs::topo {

// Topology-layer ownership stamp (review blocker 1). Deliberately independent
// of paacs/bnd: all four pool identities plus the workspace revision. The
// caller captures it while the workspace is readable; the transaction
// validates it (against the private pools/revision) BEFORE any candidate
// processing, so a candidate built before a committed metadata/revision
// change is refused even while its handles remain current.
struct WorkspaceStamp {
    std::uint64_t nodePool = 0, cellPool = 0, edgePool = 0, facePool = 0;
    std::uint64_t revision = 0;
    static WorkspaceStamp of(const TopologyWorkspace& ws) {
        return {ws.nodes().poolId(), ws.cells().poolId(), ws.edges().poolId(),
                ws.faces().poolId(), ws.topologyRevision()};
    }
    bool operator==(const WorkspaceStamp& o) const {
        return nodePool == o.nodePool && cellPool == o.cellPool &&
               edgePool == o.edgePool && facePool == o.facePool &&
               revision == o.revision;
    }
};

// ------------------------------------------------------------ value types ---
// Tagged node reference: an EXISTING current node, or a STAGED node addressed
// by its structural localId (NOT a vector position).
struct SteinerNodeRef {
    bool isStaged = false;
    NodeHandle existing{};       // when !isStaged
    std::uint32_t localId = 0;   // when isStaged
    static SteinerNodeRef Existing(const NodeHandle& h) {
        SteinerNodeRef r;
        r.isStaged = false;
        r.existing = h;
        return r;
    }
    static SteinerNodeRef Staged(std::uint32_t id) {
        SteinerNodeRef r;
        r.isStaged = true;
        r.localId = id;
        return r;
    }
};

// An interior Steiner node to be created. localId is caller-assigned and must
// be unique within the spec; it is the ONLY way cells/requirements address it.
struct SteinerNodeSpec {
    std::uint32_t localId = 0;
    Vec3 coord{0, 0, 0};
    mesh::Classification cls{};
};

struct SteinerCavityCell {
    std::array<SteinerNodeRef, 3> nodes{};
    mesh::Classification cls{};
};

// Ensure-only edge requirement (design section 3), staged-aware.
struct SteinerEdgeRequirement {
    SteinerNodeRef a, b;
    mesh::Classification cls{};
};

struct SteinerCavitySpec {
    WorkspaceStamp stamp;                         // rev-3 ownership stamp
    std::vector<CellHandle> removed;              // Tri3 disk, existing cells
    std::vector<SteinerNodeSpec> steinerNodes;    // >= 1 interior additions
    std::vector<SteinerCavityCell> added;         // retriangulation of R
    std::vector<SteinerEdgeRequirement> requirements;
};

// PLAIN VALUES ONLY (the frozen validator rule): a validator never sees the
// workspace or pools. Ranks are COMPACT symbolic ids (Phase-1 domain).
struct SteinerCavityPreview {
    std::vector<cavaudit::Rank> boundaryRanks;    // canonical compact ring
    std::size_t steinerCount = 0;
    std::vector<std::uint32_t> steinerLocalIds;   // sorted, canonical
    std::vector<std::array<cavaudit::Rank, 3>> removedRanks;
    std::vector<std::array<cavaudit::Rank, 3>> addedRanks;   // canonical
    std::size_t cellsRemoved = 0, cellsAdded = 0;
    std::size_t edgesToCreate = 0, edgesToRetain = 0, edgesToRetire = 0;
};

struct SteinerCavityCommit {
    // Explicit structural mapping: localId -> materialized handle, ascending.
    std::vector<std::pair<std::uint32_t, NodeHandle>> createdNodes;
    std::vector<CellHandle> resultCells;    // new triangles
    std::vector<EdgeHandle> requiredEdges;  // exactly the requested edges
    std::uint64_t revision = 0;
    NodeHandle nodeOf(std::uint32_t localId) const {
        for (const auto& [id, h] : createdNodes)
            if (id == localId) return h;
        return NodeHandle{};
    }
};

// Typed Phase-2 preparation failures so a caller can attribute the EXACT §9
// tier: staged-node reservation vs cell/edge/incidence preparation. Both
// derive from std::bad_alloc, so existing catch(std::bad_alloc) sites (probes)
// still treat them as allocation failures; a caller that wants the tier
// catches the specific type first.
struct StagedNodeReservationError : std::bad_alloc {
    const char* what() const noexcept override {
        return "SteinerCavity: staged-node reservation failed";
    }
};
struct CavityPreparationError : std::bad_alloc {
    const char* what() const noexcept override {
        return "SteinerCavity: cell/edge/incidence preparation failed";
    }
};

// ------------------------------------------------------------ transaction ---
class SteinerCavityTransaction {
public:
    SteinerCavityTransaction(const SteinerCavityTransaction&) = delete;
    SteinerCavityTransaction& operator=(const SteinerCavityTransaction&) =
        delete;
    SteinerCavityTransaction(SteinerCavityTransaction&& o) noexcept
        : ws_(o.ws_), resume_(o.resume_), staged_(std::move(o.staged_)),
          prep_(std::move(o.prep_)),
          resultCellsBuf_(std::move(o.resultCellsBuf_)),
          requiredEdgesBuf_(std::move(o.requiredEdgesBuf_)),
          createdNodesBuf_(std::move(o.createdNodesBuf_)),
          createdMapBuf_(std::move(o.createdMapBuf_)),
          rejected_(o.rejected_), poisoned_(o.poisoned_) {
        o.ws_ = nullptr;
        o.prep_.prepared = false;
    }
    SteinerCavityTransaction& operator=(SteinerCavityTransaction&&) = delete;
    ~SteinerCavityTransaction() {
        if (ws_) abandonImpl();
    }

    // Phase 1: stamp validation + full combinatorial validation +
    // canonicalization over SYMBOLIC compact ranks. Const reads only; pools
    // untouched (no sequence consumed, no allocator prediction).
    const SteinerCavityPreview& stage(const SteinerCavitySpec& spec) {
        requireUsable("stage");
        if (staged_)
            throw std::logic_error(
                "SteinerCavityTransaction: exactly one cavity may be staged "
                "per transaction");
        buildStaged(spec);
        return staged_->preview;
    }

    bool validate(const std::function<bool(const SteinerCavityPreview&)>& fn) {
        requireUsable("validate");
        if (ws_->state_ != WorkspaceState::TransactionStaging)
            throw std::logic_error(
                "SteinerCavityTransaction::validate: validators run in "
                "Phase 1 only");
        if (!staged_)
            throw std::logic_error(
                "SteinerCavityTransaction::validate: nothing staged");
        const bool ok = fn(staged_->preview);
        if (!ok) rejected_ = true;
        return ok;
    }
    bool rejected() const noexcept { return rejected_; }
    const SteinerCavityPreview& preview() const {
        requireActive("preview");
        if (!staged_)
            throw std::logic_error(
                "SteinerCavityTransaction::preview: nothing staged");
        return staged_->preview;
    }

    void abandon() {
        requireActive("abandon");
        abandonImpl();
    }

    // Phase 2: full node-batch demand check, real slot/gen/sequence
    // prediction (single predictSlots batch path), pool plans, staged cavity
    // preparation. No handle minted, no sequence consumed, no revision/
    // binding touched.
    //
    // PREPARE-FAILURE CONTRACT (review): a throw from prepare() POISONS the
    // transaction. It is NOT re-preparable in place; the caller must abandon
    // (or destroy) it. Any subsequent stage/validate/prepare/commit throws.
    void prepare() {
        requireUsable("prepare");
        if (!staged_)
            throw std::logic_error(
                "SteinerCavityTransaction::prepare: nothing staged");
        if (rejected_)
            throw std::logic_error(
                "SteinerCavityTransaction::prepare: transaction was "
                "rejected by a validator");
        if (prep_.prepared) return;
        try {
            prepareImpl();
        } catch (...) {
            poisoned_ = true;
            detail::TransactionAccess::abandonCavityPrepared(
                ws_->cells_, ws_->edges_, prep_);
            throw;
        }
        ws_->state_ = WorkspaceState::TransactionPrepared;
    }
    bool prepared() const noexcept { return prep_.prepared; }

    // Predicted staged-node handles (valid AFTER prepare()): the exact
    // slot/gen/pool each staged node WILL have at commit (commitNodeBatch
    // asserts actual == predicted). Lets a caller build a metadata shadow
    // keyed by the final handle BEFORE the commit, so post-commit publication
    // is bounded and non-throwing. Ordered ascending by staged localId.
    std::vector<NodeHandle> predictedStagedNodes() const {
        std::vector<NodeHandle> out;
        if (!staged_ || !prep_.prepared) return out;
        const std::uint64_t pool = ws_->nodes_.poolId();
        out.reserve(staged_->nodePred.size());
        for (const auto& p : staged_->nodePred)
            out.push_back(NodeHandle{p.slot, p.gen, pool});
        return out;
    }

    // Phase 3: bounded, noexcept, allocation-free. Materialize the staged
    // nodes, then apply the prepared cavity (both no-throw). Atomic: nothing
    // is visible until the single revision bump.
    SteinerCavityCommit commit() {
        requireUsable("commit");
        if (!staged_)
            throw std::logic_error(
                "SteinerCavityTransaction::commit: nothing staged");
        if (rejected_)
            throw std::logic_error(
                "SteinerCavityTransaction::commit: transaction was "
                "rejected by a validator");
        using TA = detail::TransactionAccess;
        if (!prep_.prepared) prepare();      // may throw+poison (auto-prep)
        St& s = *staged_;

        // ---- build CavityCommitInput (all throwing work BEFORE Phase 3) --
        TA::CavityCommitInput in;
        in.removedCells = s.removed;
        in.removedExhausted.resize(in.removedCells.size(), 0);
        {
            const auto& a = TA::cellAllocC(ws_->cells_);
            for (std::size_t k = 0; k < in.removedCells.size(); ++k)
                if (TA::slotGen(a, in.removedCells[k].slot) ==
                    std::numeric_limits<std::uint32_t>::max())
                    in.removedExhausted[k] = 1;
        }
        for (const AddedTri& t : s.addedTri) {
            TA::CavityCommitInput::AddedCell ac;
            for (int i = 0; i < 3; ++i) {
                const NodeHandle n = resolveHandle(s, t.ref[std::size_t(i)]);
                ac.nodeSlot[std::size_t(i)] = n.slot;
                ac.nodeGen[std::size_t(i)] = n.gen;
            }
            ac.cls = t.cls;
            in.added.push_back(ac);
        }
        in.cellSlots = s.cellPlan.predictedAllocationOrder;
        in.cellGens = s.cellPlan.predictedGenerations;
        in.cellPhase2Claims = s.cellPlan.phase2Claims;
        in.edgeSlots = s.edgePlan.predictedAllocationOrder;
        in.edgeGens = s.edgePlan.predictedGenerations;
        in.edgePhase2Claims = s.edgePlan.phase2Claims;

        std::map<std::pair<cavaudit::Rank, cavaudit::Rank>,
                 TA::CavityCommitInput::ERef>
            byKey;
        std::size_t createIdx = 0;
        for (const EdgeRow& r : s.edges) {
            const auto key = std::make_pair(r.ra, r.rb);
            if (r.initialUse > 0 && r.finalUse > 0) {          // retain
                TA::CavityCommitInput::Retain rt;
                rt.h = r.liveHandle;
                rt.finalUse = r.finalUse;
                rt.makeRequired = r.willBeRequired ? 1 : 0;
                rt.setCls = (r.finalCls.dim != r.initialCls.dim ||
                             r.finalCls.tag != r.initialCls.tag)
                                ? 1
                                : 0;
                rt.cls = r.finalCls;
                byKey[key] = TA::CavityCommitInput::ERef{0, rt.h, 0};
                in.retains.push_back(rt);
            } else if (r.initialUse > 0 && r.finalUse == 0) {  // retire
                in.retires.push_back(r.liveHandle);
            } else if (r.initialUse == 0 && r.finalUse > 0) {  // create
                TA::CavityCommitInput::Create cs;
                cs.a = resolveHandle(s, r.aRef);   // lo real-seq endpoint
                cs.b = resolveHandle(s, r.bRef);   // hi real-seq endpoint
                cs.finalUse = r.finalUse;
                cs.required = r.willBeRequired ? 1 : 0;
                cs.cls = r.finalCls;
                byKey[key] = TA::CavityCommitInput::ERef{
                    1, EdgeHandle{}, static_cast<std::uint32_t>(createIdx++)};
                in.creates.push_back(cs);
            }
        }
        for (const AddedTri& t : s.addedTri) {
            std::array<TA::CavityCommitInput::ERef, 3> refs;
            for (int i = 0; i < 3; ++i) {
                const cavaudit::Rank rx = t.rank[std::size_t(i)];
                const cavaudit::Rank ry = t.rank[std::size_t((i + 1) % 3)];
                refs[std::size_t(i)] =
                    byKey.at(rx < ry ? std::make_pair(rx, ry)
                                     : std::make_pair(ry, rx));
            }
            in.addedEdges.push_back(refs);
        }
        for (const EdgeRow& r : s.edges)
            if (r.requested)
                in.requested.push_back(byKey.at(std::make_pair(r.ra, r.rb)));

        std::vector<EdgeHandle> createdScratch(in.creates.size());
        createdNodesBuf_.clear();
        resultCellsBuf_.clear();
        requiredEdgesBuf_.clear();
        TA::CavityCommitBuffers bufs;
        bufs.resultCells = &resultCellsBuf_;
        bufs.requiredEdges = &requiredEdgesBuf_;
        bufs.createdScratch = &createdScratch;

        // ---- Phase 3: composite, bounded, noexcept, allocation-free ----
#if defined(PAACS_TOPO_TESTING)
        detail::phase3Active() = true;
#endif
        TA::commitNodeBatch(ws_->nodes_, s.nodeSpecs, s.nodePred,
                            createdNodesBuf_);
        TA::commitCavityPrepared(ws_->cells_, ws_->edges_, in, prep_, bufs);
#if defined(PAACS_TOPO_TESTING)
        detail::phase3Active() = false;
#endif
        // Build the localId -> handle mapping into a PRE-RESERVED buffer
        // (reserved in prepare()) so this cannot allocate: an allocation
        // AFTER the revision bump would throw post-commit and break
        // atomicity (a bug the commit-zone allocation sweep caught).
        // createdNodesBuf_ is parallel to s.nodePred / s.nodeSpecs, i.e. in
        // ascending-localId order (the canonicalization).
        createdMapBuf_.clear();
        for (std::size_t i = 0; i < s.sortedLocalIds.size(); ++i)
            createdMapBuf_.emplace_back(s.sortedLocalIds[i],
                                        createdNodesBuf_[i]);   // reserved
        ws_->revision_ += 1;                 // exactly ONCE
        ws_->state_ = WorkspaceState::CommittedAdjacencyStale;
        // Only no-throw moves remain past this point.
        SteinerCavityCommit result;
        result.createdNodes = std::move(createdMapBuf_);
        result.resultCells = std::move(resultCellsBuf_);
        result.requiredEdges = std::move(requiredEdgesBuf_);
        result.revision = ws_->revision_;
        ws_ = nullptr;
        return result;
    }

private:
    friend class TopologyWorkspace;
    explicit SteinerCavityTransaction(TopologyWorkspace& ws)
        : ws_(&ws), resume_(ws.state_) {
        ws_->state_ = WorkspaceState::TransactionStaging;
    }
    void requireActive(const char* what) const {
        if (!ws_)
            throw std::logic_error(
                std::string("SteinerCavityTransaction: '") + what +
                "' on a finished/abandoned transaction");
    }
    void requireUsable(const char* what) const {
        requireActive(what);
        if (poisoned_)
            throw std::logic_error(
                std::string("SteinerCavityTransaction: '") + what +
                "' on a poisoned transaction (prepare() failed; abandon it)");
    }
    void abandonImpl() noexcept {
        detail::TransactionAccess::abandonCavityPrepared(ws_->cells_,
                                                         ws_->edges_, prep_);
        ws_->state_ = resume_;
        ws_ = nullptr;
    }
    static void rej(const std::string& why) {
        throw std::invalid_argument("SteinerCavityTransaction::stage: " + why);
    }

    // ---- internal staged model ----
    struct AddedTri {
        std::array<SteinerNodeRef, 3> ref{};    // canonical rotation
        std::array<cavaudit::Rank, 3> rank{};   // compact, canonical rotation
        mesh::Classification cls{};
    };
    struct EdgeRow {
        SteinerNodeRef aRef{}, bRef{};   // lo-compact, hi-compact endpoints
        cavaudit::Rank ra = 0, rb = 0;   // compact ranks (ra < rb)
        std::uint32_t initialUse = 0, finalUse = 0;
        std::uint32_t localRemovedUse = 0, localAddedUse = 0;
        int removedSenseSum = 0, addedSenseSum = 0;
        std::int8_t removedSense = 0, addedSense = 0;
        bool live = false, wasRequired = false, requested = false,
             willBeRequired = false;
        EdgeHandle liveHandle{};
        mesh::Classification initialCls{}, finalCls{};
    };
    struct St {
        SteinerCavityPreview preview;
        std::vector<CellHandle> removed;                 // canonical
        // staged nodes, canonicalized ascending by localId
        std::vector<std::uint32_t> sortedLocalIds;
        std::map<std::uint32_t, std::uint32_t> localToSorted;  // id -> index
        std::vector<detail::TransactionAccess::StagedNodeSpec> nodeSpecs;
        std::vector<detail::TransactionAccess::StagedNodePredict> nodePred;
        std::map<std::uint32_t, NodeHandle> stagedHandle;      // Phase 2
        std::map<std::uint32_t, std::uint64_t> stagedSeq;      // Phase 2
        std::vector<AddedTri> addedTri;                  // canonical
        std::vector<EdgeRow> edges;                      // compact-key order
        std::vector<EdgeKey> createKeys;                 // creates order (P2)
        std::size_t retireCount = 0;
        std::size_t requestedCount = 0;
        CavityPoolPlan cellPlan, edgePlan;
    };

    static NodeHandle resolveHandle(const St& s, const SteinerNodeRef& r) {
        return r.isStaged ? s.stagedHandle.at(r.localId) : r.existing;
    }

    // Fully general per-pool plan (mirrors CavityTransaction::fillPoolPlan;
    // makes NO balance assumption - demand may exceed retirements).
    template <class Alloc>
    static void fillPlan(CavityPoolPlan& q, const Alloc& a,
                         const std::vector<std::int32_t>& freedOrder) {
        using TA = detail::TransactionAccess;
        q.reusableThisCommit = freedOrder.size();
        q.inCommitReuse = std::min(q.demand, q.reusableThisCommit);
        q.phase2Claims = q.demand - q.inCommitReuse;
        q.freeListAvailable = a.reusableFreeSlots();
        const std::size_t materialized =
            static_cast<std::size_t>(a.materializedSlots());
        q.sequenceDemand = q.demand;
        q.remainingSequences = TA::remainingSequences(a);
        const auto claims = TA::predictSlots(a, q.phase2Claims);
        q.preExistingConsumed =
            std::min(q.phase2Claims, q.freeListAvailable);
        q.freshSlots = q.phase2Claims - q.preExistingConsumed;
        q.predictedAllocationOrder = claims;
        for (std::size_t i = 0; i < q.inCommitReuse; ++i)
            q.predictedAllocationOrder.push_back(
                freedOrder[freedOrder.size() - 1 - i]);
        for (std::size_t i = 0; i < q.predictedAllocationOrder.size(); ++i) {
            const auto sl = q.predictedAllocationOrder[i];
            if (i >= claims.size())
                q.predictedGenerations.push_back(TA::slotGen(a, sl) + 1);
            else if (static_cast<std::size_t>(sl) < materialized)
                q.predictedGenerations.push_back(TA::slotGen(a, sl));
            else
                q.predictedGenerations.push_back(0);
        }
        q.finalFreeListSize =
            (q.freeListAvailable - q.preExistingConsumed) +
            (q.reusableThisCommit - q.inCommitReuse);
    }

    // ================================================= Phase 1 (symbolic) ===
    void buildStaged(const SteinerCavitySpec& spec) {
        const NodePool& np = ws_->nodes_;
        const CellPool& cp = ws_->cells_;
        const EdgePool& ep = ws_->edges_;

        // (0) ownership stamp - BEFORE any candidate processing.
        if (!(spec.stamp.nodePool == np.poolId() &&
              spec.stamp.cellPool == cp.poolId() &&
              spec.stamp.edgePool == ep.poolId() &&
              spec.stamp.facePool == ws_->faces_.poolId() &&
              spec.stamp.revision == ws_->revision_))
            rej("workspace ownership stamp is stale or foreign (pool ids / "
                "revision differ from the current workspace)");

        if (!workspaceTopologyReady(np, cp, ep, ws_->faces_))
            throw std::logic_error(
                "SteinerCavityTransaction::stage: workspace is not "
                "topology-ready (composite gate failed)");

        St s;
        // (1) staged nodes: >= 1, unique localIds, canonicalize by localId.
        if (spec.steinerNodes.empty())
            rej("no Steiner nodes (use the frozen cavity for a pure "
                "retriangulation without new nodes)");
        {
            std::vector<std::uint32_t> ids;
            for (const SteinerNodeSpec& n : spec.steinerNodes) {
                if (!(std::isfinite(n.coord.x) && std::isfinite(n.coord.y) &&
                      std::isfinite(n.coord.z)))
                    rej("Steiner node coordinate is not finite");
                ids.push_back(n.localId);
            }
            std::vector<std::uint32_t> sorted = ids;
            std::sort(sorted.begin(), sorted.end());
            for (std::size_t i = 1; i < sorted.size(); ++i)
                if (sorted[i - 1] == sorted[i])
                    rej("duplicate Steiner node localId");
            s.sortedLocalIds = sorted;
            for (std::uint32_t i = 0; i < sorted.size(); ++i)
                s.localToSorted[sorted[i]] = i;
            // nodeSpecs in canonical (sorted-localId) order
            std::map<std::uint32_t, const SteinerNodeSpec*> byId;
            for (const SteinerNodeSpec& n : spec.steinerNodes)
                byId[n.localId] = &n;
            for (std::uint32_t id : sorted)
                s.nodeSpecs.push_back({byId[id]->coord, byId[id]->cls});
        }

        // (2) removed: non-empty, current, distinct, Tri3; canonical.
        if (spec.removed.empty()) rej("removed set is empty");
        std::vector<CellHandle> removed = spec.removed;
        for (const CellHandle& h : removed) {
            if (!cp.isCurrent(h))
                rej("removed cell is not current (stale or foreign)");
            if (cp.get(h).type != ElementType::Tri3)
                rej("BND-2A qualifies Tri3 cavities only");
        }
        std::sort(removed.begin(), removed.end(),
                  [&](const CellHandle& x, const CellHandle& y) {
                      return cp.seqOf(x) < cp.seqOf(y);
                  });
        for (std::size_t i = 1; i < removed.size(); ++i)
            if (removed[i - 1] == removed[i])
                rej("removed set contains duplicates");
        s.removed = removed;

        // (3) gather R connectivity + COMPACT ranks (symbolic, order-
        //     preserving): existing nodes by live sequence -> 0..ne-1;
        //     staged nodes after them by sorted localId -> ne + idx.
        std::vector<std::array<NodeHandle, 3>> rConn;
        std::set<std::uint64_t> rSeqSet;
        for (const CellHandle& h : removed) {
            const auto sn = cp.get(h);
            rConn.push_back({sn.nodes[0], sn.nodes[1], sn.nodes[2]});
            for (int i = 0; i < 3; ++i) rSeqSet.insert(np.seqOf(sn.nodes[i]));
        }
        std::vector<std::uint64_t> rSeqs(rSeqSet.begin(), rSeqSet.end());
        std::map<std::uint64_t, cavaudit::Rank> existingCompact;
        for (std::size_t i = 0; i < rSeqs.size(); ++i)
            existingCompact[rSeqs[i]] = static_cast<cavaudit::Rank>(i);
        const cavaudit::Rank ne = static_cast<cavaudit::Rank>(rSeqs.size());
        std::map<cavaudit::Rank, NodeHandle> compactToNode;
        for (const auto& [seq, c] : existingCompact) {
            (void)seq;
        }
        // build compact->handle for existing (needed for boundary preview)
        for (const CellHandle& h : removed) {
            const auto sn = cp.get(h);
            for (int i = 0; i < 3; ++i)
                compactToNode[existingCompact[np.seqOf(sn.nodes[i])]] =
                    sn.nodes[i];
        }
        auto stagedCompact = [&](std::uint32_t localId) {
            return ne + static_cast<cavaudit::Rank>(s.localToSorted.at(localId));
        };
        auto refCompact = [&](const SteinerNodeRef& r) -> cavaudit::Rank {
            if (r.isStaged) {
                if (!s.localToSorted.count(r.localId))
                    rej("reference names an unknown staged localId");
                return stagedCompact(r.localId);
            }
            if (!np.isCurrent(r.existing))
                rej("added/requirement references a non-current node");
            const std::uint64_t sq = np.seqOf(r.existing);
            const auto it = existingCompact.find(sq);
            if (it == existingCompact.end())
                rej("added triangle pulls in a node that is neither in the "
                    "removed cavity nor a staged Steiner node");
            return it->second;
        };

        // (4) added: resolve refs -> compact ranks, canonical rotation.
        std::set<cavaudit::Rank> usedStaged, mRanks;
        std::vector<AddedTri> aTri;
        for (const SteinerCavityCell& c : spec.added) {
            std::array<SteinerNodeRef, 3> ref = c.nodes;
            std::array<cavaudit::Rank, 3> r{};
            for (int i = 0; i < 3; ++i) {
                r[std::size_t(i)] = refCompact(ref[std::size_t(i)]);
                if (ref[std::size_t(i)].isStaged)
                    usedStaged.insert(r[std::size_t(i)]);
            }
            if (r[0] == r[1] || r[1] == r[2] || r[0] == r[2])
                rej("added triangle repeats a node");
            for (int i = 0; i < 3; ++i) mRanks.insert(r[std::size_t(i)]);
            int lo = 0;
            for (int i = 1; i < 3; ++i)
                if (r[std::size_t(i)] < r[std::size_t(lo)]) lo = i;
            AddedTri t;
            for (int i = 0; i < 3; ++i) {
                t.ref[std::size_t(i)] = ref[std::size_t((lo + i) % 3)];
                t.rank[std::size_t(i)] = r[std::size_t((lo + i) % 3)];
            }
            t.cls = c.cls;
            aTri.push_back(t);
        }
        std::sort(aTri.begin(), aTri.end(),
                  [](const AddedTri& x, const AddedTri& y) {
                      return x.rank < y.rank;
                  });
        {
            std::set<std::array<cavaudit::Rank, 3>> seen;
            for (const AddedTri& t : aTri) {
                std::array<cavaudit::Rank, 3> k = t.rank;
                std::sort(k.begin(), k.end());
                if (!seen.insert(k).second)
                    rej("added set contains duplicate triangles");
            }
        }

        // (5) every staged node used; every R node survives.
        for (std::uint32_t id : s.sortedLocalIds)
            if (!usedStaged.count(stagedCompact(id)))
                rej("a staged Steiner node is never used by any added "
                    "triangle");
        for (cavaudit::Rank r = 0; r < ne; ++r)
            if (!mRanks.count(r))
                rej("a removed-cavity node would be orphaned (every node of "
                    "R must appear in M)");

        // (6) disk + node-link for BOTH R and M; boundary preserved.
        std::vector<cavaudit::RankTri> Rc, Mc;
        for (const auto& tri : rConn) {
            std::array<cavaudit::Rank, 3> tr{
                existingCompact[np.seqOf(tri[0])],
                existingCompact[np.seqOf(tri[1])],
                existingCompact[np.seqOf(tri[2])]};
            Rc.push_back({tr[0], tr[1], tr[2]});
        }
        for (const AddedTri& t : aTri)
            Mc.push_back({t.rank[0], t.rank[1], t.rank[2]});
        if (!cavaudit::auditCavityDiskCombinatorial(Rc))
            rej("removed set R is not a combinatorial disk");
        if (!cavaudit::auditCavityNodeLink(Rc))
            rej("removed set R fails the node-link audit");
        if (!cavaudit::auditCavityDiskCombinatorial(Mc))
            rej("candidate set M is not a combinatorial disk");
        if (!cavaudit::auditCavityNodeLink(Mc))
            rej("candidate set M fails the node-link audit (pinched)");
        std::vector<cavaudit::Rank> cycR, cycM;
        if (!cavaudit::boundaryCycle(Rc, cycR))
            rej("removed set R has no single boundary cycle");
        if (!cavaudit::boundaryCycle(Mc, cycM))
            rej("candidate set M has no single boundary cycle");
        if (cycR != cycM)
            rej("boundary(M) != boundary(R) (a Steiner node must be interior; "
                "the boundary ring must be preserved exactly)");

        // (7) edge closure over compact ranks (roles + senses), keeping the
        //     endpoint REFERENCES so Phase 2 can resolve real keys.
        std::map<std::pair<cavaudit::Rank, cavaudit::Rank>, EdgeRow> em;
        auto touch = [&](const SteinerNodeRef& xr, cavaudit::Rank rx,
                         const SteinerNodeRef& yr,
                         cavaudit::Rank ry) -> EdgeRow& {
            const auto key = rx < ry ? std::make_pair(rx, ry)
                                     : std::make_pair(ry, rx);
            auto it = em.find(key);
            if (it == em.end()) {
                EdgeRow e;
                e.ra = key.first;
                e.rb = key.second;
                if (rx < ry) {
                    e.aRef = xr;
                    e.bRef = yr;
                } else {
                    e.aRef = yr;
                    e.bRef = xr;
                }
                if (!xr.isStaged && !yr.isStaged) {
                    const auto eh = ep.find(np, xr.existing, yr.existing);
                    if (eh) {
                        e.live = true;
                        e.liveHandle = *eh;
                        const auto d = ep.get(*eh);
                        e.initialUse = d.useCount;
                        e.wasRequired = d.requiredConstraint;
                        e.initialCls = d.cls;
                    }
                }
                it = em.emplace(key, e).first;
            }
            return it->second;
        };
        auto senseByRank = [](cavaudit::Rank a, cavaudit::Rank b) {
            return a < b ? std::int8_t(+1) : std::int8_t(-1);
        };
        for (std::size_t ti = 0; ti < rConn.size(); ++ti) {
            std::array<SteinerNodeRef, 3> ref{
                SteinerNodeRef::Existing(rConn[ti][0]),
                SteinerNodeRef::Existing(rConn[ti][1]),
                SteinerNodeRef::Existing(rConn[ti][2])};
            std::array<cavaudit::Rank, 3> rk{
                existingCompact[np.seqOf(rConn[ti][0])],
                existingCompact[np.seqOf(rConn[ti][1])],
                existingCompact[np.seqOf(rConn[ti][2])]};
            for (int i = 0; i < 3; ++i) {
                const int j = (i + 1) % 3;
                EdgeRow& e = touch(ref[std::size_t(i)], rk[std::size_t(i)],
                                   ref[std::size_t(j)], rk[std::size_t(j)]);
                e.localRemovedUse += 1;
                e.removedSense = senseByRank(rk[std::size_t(i)],
                                             rk[std::size_t(j)]);
                e.removedSenseSum += e.removedSense;
            }
        }
        for (const AddedTri& t : aTri)
            for (int i = 0; i < 3; ++i) {
                const int j = (i + 1) % 3;
                EdgeRow& e = touch(t.ref[std::size_t(i)], t.rank[std::size_t(i)],
                                   t.ref[std::size_t(j)], t.rank[std::size_t(j)]);
                e.localAddedUse += 1;
                e.addedSense = senseByRank(t.rank[std::size_t(i)],
                                           t.rank[std::size_t(j)]);
                e.addedSenseSum += e.addedSense;
            }

        // (8) global-manifold precondition + per-edge rules.
        for (auto& [k, e] : em) {
            (void)k;
            if (e.live && e.initialUse > 2)
                rej("source edge in the cavity closure is already "
                    "nonmanifold (initialUse > 2)");
            if (e.localRemovedUse > 2 || e.localAddedUse > 2)
                rej("a cavity edge is used more than twice within one "
                    "complex (nonmanifold candidate)");
            if (e.localRemovedUse > e.initialUse)
                rej("removed cells use an edge more often than its live use "
                    "count (missing or corrupted incidence)");
            const std::uint32_t outsideUse = e.initialUse - e.localRemovedUse;
            if (e.localRemovedUse == 2) {
                if (outsideUse != 0) rej("interior(R) edge has outside use");
                if (e.removedSenseSum != 0)
                    rej("interior edge of R has non-opposite senses");
            }
            if (e.localAddedUse == 2) {
                if (outsideUse != 0)
                    rej("candidate interior edge is already used by a "
                        "surviving outside triangle (outsideUse != 0)");
                if (e.addedSenseSum != 0)
                    rej("interior edge of M has non-opposite senses");
            }
        }
        for (auto& [k, e] : em) {
            (void)k;
            if (e.localRemovedUse == 1 && e.localAddedUse == 1 &&
                e.removedSense != e.addedSense)
                rej("a boundary edge's inside sense differs between R and M");
        }

        // (9) no candidate duplicates a surviving cell.
        {
            std::set<std::array<cavaudit::Rank, 3>> mset;
            for (const AddedTri& t : aTri) {
                std::array<cavaudit::Rank, 3> k = t.rank;
                std::sort(k.begin(), k.end());
                mset.insert(k);
            }
            std::set<std::int32_t> removedSlots;
            for (const CellHandle& h : removed) removedSlots.insert(h.slot);
            for (const CellHandle& h : cp.activeBySequence()) {
                if (removedSlots.count(h.slot)) continue;
                const auto sn = cp.get(h);
                if (sn.type != ElementType::Tri3) continue;
                std::array<cavaudit::Rank, 3> k{};
                bool allR = true;
                for (int i = 0; i < 3; ++i) {
                    const auto it = existingCompact.find(np.seqOf(sn.nodes[i]));
                    if (it == existingCompact.end()) { allR = false; break; }
                    k[std::size_t(i)] = it->second;
                }
                if (!allR) continue;
                std::sort(k.begin(), k.end());
                if (mset.count(k))
                    rej("an added triangle duplicates a surviving cell");
            }
        }

        // (10) requirements: staged-aware, canonical, ensure-only.
        std::map<std::pair<cavaudit::Rank, cavaudit::Rank>,
                 mesh::Classification>
            reqs;
        for (const SteinerEdgeRequirement& r : spec.requirements) {
            const cavaudit::Rank ra = refCompact(r.a);
            const cavaudit::Rank rb = refCompact(r.b);
            if (ra == rb) rej("requirement endpoints must be distinct");
            const auto key = ra < rb ? std::make_pair(ra, rb)
                                     : std::make_pair(rb, ra);
            const auto it = reqs.find(key);
            if (it != reqs.end()) {
                if (it->second.dim != r.cls.dim || it->second.tag != r.cls.tag)
                    rej("duplicate contradictory requirements for one edge");
                continue;
            }
            reqs.emplace(key, r.cls);
        }
        const mesh::Classification defCls{};
        auto isDefault = [&](const mesh::Classification& c) {
            return c.dim == defCls.dim && c.tag == defCls.tag;
        };
        for (const auto& [k, cls] : reqs) {
            const auto it = em.find(k);
            if (it == em.end())
                rej("requirement names an edge outside the cavity closure");
            EdgeRow& e = it->second;
            const bool inFinal =
                (e.initialUse - e.localRemovedUse + e.localAddedUse) > 0;
            if (!inFinal)
                rej("requirement names an edge absent from the final "
                    "topology (it retires in this cavity)");
            if (e.live && e.wasRequired) {
                if (!(e.initialCls.dim == cls.dim && e.initialCls.tag == cls.tag))
                    rej("conflicting classification on an already "
                        "constrained edge");
            } else if (e.live && !isDefault(e.initialCls) &&
                       !(e.initialCls.dim == cls.dim &&
                         e.initialCls.tag == cls.tag)) {
                rej("conflicting non-default classification on an existing "
                    "edge");
            }
        }

        // (11) required-edge preservation.
        for (auto& [k, e] : em) {
            (void)k;
            const std::uint32_t finalUse =
                e.initialUse - e.localRemovedUse + e.localAddedUse;
            if (e.live && e.wasRequired && finalUse == 0)
                rej("an initially required edge would be retired by this "
                    "cavity");
        }

        // (12) finalize edge rows, create order + counts, and the preview.
        for (auto& [k, e] : em) {
            e.finalUse = e.initialUse - e.localRemovedUse + e.localAddedUse;
            const auto rit = reqs.find(k);
            e.willBeRequired = e.wasRequired || rit != reqs.end();
            e.requested = rit != reqs.end();
            e.finalCls = rit != reqs.end() ? rit->second : e.initialCls;
            if (e.requested) ++s.requestedCount;
            s.edges.push_back(e);
        }
        SteinerCavityPreview& pv = s.preview;
        for (const EdgeRow& e : s.edges) {
            if (e.initialUse == 0 && e.finalUse > 0)
                ++pv.edgesToCreate;
            else if (e.initialUse > 0 && e.finalUse > 0)
                ++pv.edgesToRetain;
            else if (e.initialUse > 0 && e.finalUse == 0)
                ++s.retireCount, ++pv.edgesToRetire;
        }
        s.addedTri = std::move(aTri);
        pv.boundaryRanks = cycM;
        pv.steinerCount = s.nodeSpecs.size();
        pv.steinerLocalIds = s.sortedLocalIds;
        for (const auto& tri : Rc) pv.removedRanks.push_back(tri);
        for (const AddedTri& t : s.addedTri)
            pv.addedRanks.push_back({t.rank[0], t.rank[1], t.rank[2]});
        pv.cellsRemoved = s.removed.size();
        pv.cellsAdded = s.addedTri.size();
        (void)compactToNode;

        staged_.emplace(std::move(s));
    }

    // ================================================ Phase 2 (allocator) ===
    void prepareImpl() {
        using TA = detail::TransactionAccess;
        St& s = *staged_;

        // ===== NODE ZONE: any failure here is StagedNodeReservationError =====
        try {
            // (A) full node-batch demand FIRST (sequence AND slot space), so
            //     no wrapped rank can ever be produced.
            if (!TA::nodeBatchFits(ws_->nodes_, s.nodeSpecs.size()))
                throw StagedNodeReservationError{};
            // (B) real slot/gen/sequence prediction via the single batch path.
            s.nodePred = TA::predictNodeBatch(ws_->nodes_, s.nodeSpecs.size());
            const std::uint64_t poolId = ws_->nodes_.poolId();
            for (std::size_t i = 0; i < s.sortedLocalIds.size(); ++i) {
                const std::uint32_t id = s.sortedLocalIds[i];
                s.stagedHandle[id] =
                    NodeHandle{s.nodePred[i].slot, s.nodePred[i].gen, poolId};
                s.stagedSeq[id] = s.nodePred[i].seq;
            }
            // (E) node-batch capacity (own pool).
            TA::reserveNodeBatch(ws_->nodes_, s.nodePred);
        } catch (const StagedNodeReservationError&) {
            throw;
        } catch (const std::bad_alloc&) {
            throw StagedNodeReservationError{};
        }

        auto refSeq = [&](const SteinerNodeRef& r) -> std::uint64_t {
            return r.isStaged ? s.stagedSeq.at(r.localId)
                              : ws_->nodes_.seqOf(r.existing);
        };

        // ===== CAVITY ZONE: any failure here is CavityPreparationError =====
        try {
            // (C) pool plans (fully general; demand may exceed retirements).
            {
                const auto& a = TA::cellAllocC(ws_->cells_);
                std::vector<std::int32_t> freed;
                for (const CellHandle& h : s.removed)
                    if (TA::slotGen(a, h.slot) ==
                        std::numeric_limits<std::uint32_t>::max())
                        ++s.cellPlan.exhaustedRetirements;
                    else
                        freed.push_back(h.slot);
                s.cellPlan.demand = s.addedTri.size();
                s.cellPlan.plannedRetirements = s.removed.size();
                fillPlan(s.cellPlan, a, freed);
            }
            {
                const auto& a = TA::edgeAllocC(ws_->edges_);
                std::vector<std::int32_t> freed;
                for (const EdgeRow& r : s.edges) {
                    if (!(r.live && r.finalUse == 0)) continue;
                    if (TA::slotGen(a, r.liveHandle.slot) ==
                        std::numeric_limits<std::uint32_t>::max())
                        ++s.edgePlan.exhaustedRetirements;
                    else
                        freed.push_back(r.liveHandle.slot);
                }
                s.edgePlan.demand = 0;
                for (const EdgeRow& r : s.edges)
                    if (r.initialUse == 0 && r.finalUse > 0)
                        ++s.edgePlan.demand;
                s.edgePlan.plannedRetirements = s.retireCount;
                fillPlan(s.edgePlan, a, freed);
            }
            if (s.cellPlan.remainingSequences < s.cellPlan.sequenceDemand)
                throw CavityPreparationError{};
            if (s.edgePlan.remainingSequences < s.edgePlan.sequenceDemand)
                throw CavityPreparationError{};

            // (D) Phase-3 result buffers reserved BEFORE any placeholder, so
            //     the post-Phase-3 tail (the localId->handle mapping) is
            //     allocation-free and cannot throw after the revision bump.
            createdNodesBuf_.reserve(s.nodeSpecs.size());
            createdMapBuf_.reserve(s.nodeSpecs.size());
            resultCellsBuf_.reserve(s.addedTri.size());
            requiredEdgesBuf_.reserve(s.requestedCount);

            // (F) REAL create keys through the staged-edge helper for EVERY
            //     create (Existing-Existing, Existing-Staged, Staged-Staged).
            s.createKeys.clear();
            for (const EdgeRow& r : s.edges) {
                if (!(r.initialUse == 0 && r.finalUse > 0)) continue;
                EdgeKey k;
                const bool ok = TA::stagedEdgeKey(
                    ws_->nodes_, r.aRef.isStaged, refSeq(r.aRef),
                    resolveHandle(s, r.aRef), r.bRef.isStaged, refSeq(r.bRef),
                    resolveHandle(s, r.bRef), k);
                if (!ok)
                    throw std::logic_error(
                        "SteinerCavityTransaction::prepare: created edge "
                        "resolved to a self-edge (sequence collision)");
                s.createKeys.push_back(k);
            }

            // (G) staged cavity preparation.
            TA::CavityStagedPrepareInput in;
            in.cellAllocOrder = s.cellPlan.predictedAllocationOrder;
            in.edgeAllocOrder = s.edgePlan.predictedAllocationOrder;
            in.cellFreePushes = s.cellPlan.reusableThisCommit;
            in.edgeFreePushes = s.edgePlan.reusableThisCommit;
            in.cellFresh = s.cellPlan.freshSlots;
            in.edgeFresh = s.edgePlan.freshSlots;
            in.newEdgeKeys = s.createKeys;
            TA::prepareStagedCavity(ws_->cells_, ws_->edges_, in, prep_);
        } catch (const CavityPreparationError&) {
            throw;
        } catch (const std::bad_alloc&) {
            throw CavityPreparationError{};
        }
    }

    TopologyWorkspace* ws_ = nullptr;
    WorkspaceState resume_ = WorkspaceState::Ready;
    std::optional<St> staged_;
    detail::TransactionAccess::CavityPrepared prep_;
    std::vector<CellHandle> resultCellsBuf_;
    std::vector<EdgeHandle> requiredEdgesBuf_;
    std::vector<NodeHandle> createdNodesBuf_;
    std::vector<std::pair<std::uint32_t, NodeHandle>> createdMapBuf_;
    bool rejected_ = false;
    bool poisoned_ = false;
};

inline SteinerCavityTransaction
TopologyWorkspace::beginSteinerCavityTransaction() {
    if (!readable())
        throw std::logic_error(
            "TopologyWorkspace::beginSteinerCavityTransaction: a "
            "transaction is already active (exactly one may be active)");
    return SteinerCavityTransaction(*this);
}

}  // namespace paacs::topo
