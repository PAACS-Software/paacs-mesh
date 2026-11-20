// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// paacs::topo::detail - M4-B.6 prepared-commit access layer (stage 3).
//
// STRUCTURAL access control ("detail does not mean private"): every
// member is PRIVATE; the ONLY friends are TopologyTransaction and the
// PAACS_TOPO_TESTING probe (production-absent, seam-checked). All pool-
// internal touches for the transaction machinery live HERE, behind the
// friend wall - the frozen pools gained only one friend declaration each.
//
// Phase discipline (docs/M4B6_DESIGN.md): planAddCell/planRemoveCell =
// Phase 1 (throwing; pools untouched; all validation + structured
// exhaustion prechecks). prepareAddCell/prepareRemoveCell = Phase 2
// (throwing; ONLY capacity-class residue plus erasable placeholder map
// nodes; NO production handles, NO sequence consumption, materialized
// slot count unchanged). commitAddCell/commitRemoveCell = Phase 3
// (noexcept, validation-free, allocation-free; every step's throw
// source and preparation is a design-note evidence-table row).
// abandonAddPrepared = noexcept placeholder erase.
//
#include <cstdlib>
#include <array>
#include <limits>
#include <map>
#include <set>
#include <span>
#include <vector>

#include "paacs/topo/FaceTopology.hpp"

namespace paacs::topo {

// Transaction policy (design decision; lives here so Phase-1 planning
// can consume it without a circular include with Workspace.hpp).
enum class NonManifoldPolicy { Reject, AllowAndDiagnose };

struct TopologyPolicy {
    NonManifoldPolicy nonManifold = NonManifoldPolicy::Reject;
    bool requireOppositeInteriorOrientation = true;
    bool permitOrphanNodes = false;
};

// M4-B.8 STAGE A (disclosed frozen-file additive edit #2): the cavity
// primitive in the unfrozen Cavity.hpp will use the prepared no-alloc
// primitives in Stages B/C; the friend is granted now so both access
// edits freeze together with Stage A's review.
class CavityTransaction;
// BND-2A / ASCR: the staged-node atomic cavity transaction lives in the
// (unfrozen) SteinerCavity.hpp and drives the additive batched helpers below;
// friended here (rev-3 §0.5.2) so it can reach the prepared no-alloc pool
// primitives, exactly as CavityTransaction does.
class SteinerCavityTransaction;
class TetCavityTransaction;  // T2-A accepted revision-2 additive contact
class TetRefinementTransaction; // T3-A accepted revision-2 additive contact

} // namespace paacs::topo

namespace paacs::topo::detail {

#if defined(PAACS_TOPO_TESTING)
struct TransactionAccessTestProbe;   // test-only friend (seam-checked)
// Set around Phase 3 so a global operator-new probe can assert zero
// allocations during commit (tests/probe_transaction.cpp).
inline bool& phase3Active() { static bool v = false; return v; }
#endif

class TransactionAccess {
    friend class paacs::topo::TopologyTransaction;
    friend class paacs::topo::CavityTransaction;
    friend class paacs::topo::SteinerCavityTransaction;   // BND-2A / ASCR
    friend class paacs::topo::TetCavityTransaction;       // T2-A
    friend class paacs::topo::TetRefinementTransaction;   // T3-A
#if defined(PAACS_TOPO_TESTING)
    friend struct TransactionAccessTestProbe;
#endif

private:
    // ---- Stage-1 primitives (contract unchanged) ---------------------------
    template <class Tag>
    static Handle<Tag> allocPrepared(HandleAllocator<Tag>& a) noexcept {
#if defined(PAACS_TOPO_TESTING)
        if (a.poolId() == 0) std::abort();
        if (a.nextSeq_ == std::numeric_limits<std::uint64_t>::max())
            std::abort();                    // sequence exhaustion
        if (a.freeList_.empty() &&
            HandleAllocator<Tag>::atSlotCapacity(a.slots_.size()))
            std::abort();                    // slot-count exhaustion
        if (a.freeList_.empty() &&
            a.slots_.capacity() <= a.slots_.size())
            std::abort();                    // unprepared fresh capacity
#endif
        std::int32_t s;
        if (!a.freeList_.empty()) {
            s = a.freeList_.back();
            a.freeList_.pop_back();
            a.slots_[static_cast<std::size_t>(s)].alive = true;
        } else {
            a.slots_.push_back(typename HandleAllocator<Tag>::Slot{});
            s = static_cast<std::int32_t>(a.slots_.size() - 1);
            a.slots_[static_cast<std::size_t>(s)].alive = true;
        }
        auto& sl = a.slots_[static_cast<std::size_t>(s)];
        sl.seq = a.nextSeq_++;
        ++a.live_;
        return Handle<Tag>{s, sl.gen, a.poolId()};
    }
    template <class Tag>
    static void freePrepared(HandleAllocator<Tag>& a,
                             const Handle<Tag>& h) noexcept {
        auto& sl = a.slots_[static_cast<std::size_t>(h.slot)];
        if (sl.gen == std::numeric_limits<std::uint32_t>::max()) {
            sl.alive = false;
            sl.retired = true;
            ++a.retired_;
            --a.live_;
            return;
        }
#if defined(PAACS_TOPO_TESTING)
        if (a.freeList_.capacity() <= a.freeList_.size()) std::abort();
#endif
        a.freeList_.push_back(h.slot);
        sl.alive = false;
        ++sl.gen;
        --a.live_;
    }

    // ---- Phase-1/2 read helpers (no mutation) ------------------------------
    // Predict the slots the NEXT k allocations will use, without allocating
    // or consuming sequences (single-threaded, single active transaction;
    // commit asserts each actual slot equals its prediction under testing).
    template <class Tag>
    static std::vector<std::int32_t> predictSlots(
        const HandleAllocator<Tag>& a, std::size_t k) {
        std::vector<std::int32_t> out;
        out.reserve(k);
        std::size_t freeIdx = a.freeList_.size();
        std::size_t nextFresh = a.slots_.size();
        for (std::size_t i = 0; i < k; ++i) {
            if (freeIdx > 0)
                out.push_back(a.freeList_[--freeIdx]);   // pops from back
            else
                out.push_back(static_cast<std::int32_t>(nextFresh++));
        }
        return out;
    }
    template <class Tag>
    static std::uint32_t slotGen(const HandleAllocator<Tag>& a,
                                 std::int32_t slot) {
        return static_cast<std::size_t>(slot) < a.slots_.size()
                   ? a.slots_[static_cast<std::size_t>(slot)].gen
                   : 0;                                   // fresh slot
    }
    template <class Tag>
    static bool sequencesExhausted(const HandleAllocator<Tag>& a) {
        return a.nextSeq_ == std::numeric_limits<std::uint64_t>::max();
    }
    static std::int8_t faceSense(const FacePool& fp, const NodePool& nodes,
                                 const FaceHandle& fh,
                                 const FacePool::Staged& st) {
        const auto& r =
            fp.rec_[static_cast<std::size_t>(raw::slotUnchecked(fh))];
        return fp.orientationSense(nodes, r.data, st);
    }
    static std::size_t faceIncidenceCount(const FacePool& fp,
                                          const FaceHandle& fh) {
        return fp.rec_[static_cast<std::size_t>(raw::slotUnchecked(fh))]
            .data.incident.size();
    }
    static std::int8_t faceSoleSense(const FacePool& fp,
                                     const FaceHandle& fh) {
        const auto& inc =
            fp.rec_[static_cast<std::size_t>(raw::slotUnchecked(fh))]
                .data.incident;
        return inc.empty() ? std::int8_t{0} : inc[0].sense;
    }

    // Exact per-pool Phase-2 accounting (review: every prepared capacity
    // claim states the exact reusable slots consumed and fresh slots
    // materialized; becomes essential once replacement transactions
    // stage additions and removals together).
    struct PoolAccount {
        std::size_t requested = 0;      // handles this op will allocate
        std::size_t reusable = 0;       // satisfied from the free list
        std::size_t fresh = 0;          // materialized slots required
        std::size_t retirements = 0;    // handles this op will free
        std::size_t slotCapacity = 0;   // reserveMaterializedCapacity target
        std::size_t freeCapacity = 0;   // reserveFreeCapacity amount
    };
    template <class Tag>
    static PoolAccount account(const HandleAllocator<Tag>& a,
                               std::size_t requested,
                               std::size_t retirements) {
        PoolAccount acc;
        acc.requested = requested;
        acc.reusable = std::min(requested, a.freeList_.size());
        acc.fresh = requested - acc.reusable;
        acc.retirements = retirements;
        acc.slotCapacity = a.materializedSlots() + acc.fresh;
        acc.freeCapacity = retirements;
        return acc;
    }
    // ==== M4-B.8 STAGE B (disclosed additive frozen-file edit #3, ====
    // ==== reviewer-approved AMENDED OPTION 1, 2026-07-18)          ====
    // An op-shaped prepared-cavity surface: exact placeholders and
    // capacity, NOT merely counts. It must not mint handles, consume
    // sequences, change bindings/revision, or expose topology mutation;
    // placeholders are invisible behind the transaction visibility
    // contract and are fully removed by abandonCavityPrepared.
    template <class Tag>
    static std::uint64_t remainingSequences(const HandleAllocator<Tag>& a) {
        return std::numeric_limits<std::uint64_t>::max() - a.nextSeq_;
    }
    struct CavityPrepareInput {
        std::vector<std::int32_t> cellAllocOrder;   // predicted, full
        std::vector<std::int32_t> edgeAllocOrder;   // (created edges)
        std::vector<std::pair<NodeHandle, NodeHandle>> newEdgeEndpoints;
        std::size_t cellFreePushes = 0;   // non-exhausted cell frees
        std::size_t edgeFreePushes = 0;
        std::size_t cellFresh = 0;        // fresh slots to materialize
        std::size_t edgeFresh = 0;
    };
    struct CavityPrepared {
        bool prepared = false;
        using TabIt = decltype(std::declval<EdgePool&>().table_.begin());
        std::vector<TabIt> tabIts;        // new-edge table placeholders
        std::vector<char> tabPlaced;      // parallel; cleanup flags
        std::vector<std::int32_t> freshIncidenceSlots;  // emplaced entries
    };
    static void prepareCavity(const NodePool& nodes, CellPool& cells,
                              EdgePool& edges, const CavityPrepareInput& in,
                              CavityPrepared& out) {
#if defined(PAACS_TOPO_TESTING)
        maybeFault();                              // seam: entry
#endif
        // Allocator capacity (public frozen contracts; capacity only).
        cells.alloc_.reserveMaterializedCapacity(
            cells.alloc_.materializedSlots() + in.cellFresh);
        cells.alloc_.reserveFreeCapacity(in.cellFreePushes);
        edges.alloc_.reserveMaterializedCapacity(
            edges.alloc_.materializedSlots() + in.edgeFresh);
        edges.alloc_.reserveFreeCapacity(in.edgeFreePushes);
        // Record arrays sized through the HIGHEST predicted slot, so no
        // Phase-3 write can touch rec_ beyond size.
        {
            std::size_t need = cells.rec_.size();
            for (const auto s : in.cellAllocOrder)
                need = std::max(need, static_cast<std::size_t>(s) + 1);
            if (need > cells.rec_.size()) {
                detail::reserveFor(cells.rec_, need);
                cells.rec_.resize(need);
            }
        }
        {
            std::size_t need = edges.rec_.size();
            for (const auto s : in.edgeAllocOrder)
                need = std::max(need, static_cast<std::size_t>(s) + 1);
            if (need > edges.rec_.size()) {
                detail::reserveFor(edges.rec_, need);
                edges.rec_.resize(need);
            }
        }
#if defined(PAACS_TOPO_TESTING)
        maybeFault();                              // seam: post-capacity
#endif
        // EXCEPTION-SAFETY ORDER (Stage-B review blocker 1):
        // 1. table BUCKET reserve for current size + new edges - a later
        //    insertion may otherwise REHASH and invalidate the earlier
        //    iterators recorded for cleanup (undefined behaviour);
        // 2. pre-reserve ALL cleanup-tracking vectors, so every
        //    record-a-placement push below is non-allocating;
        // 3. insert each node, RECORD it immediately, and only then
        //    reserve its embedded vector;
        // 4. a failure seam after each recorded placement.
        edges.table_.reserve(edges.table_.size() +
                             in.newEdgeEndpoints.size());
        out.tabIts.reserve(in.newEdgeEndpoints.size());
        out.tabPlaced.reserve(in.newEdgeEndpoints.size());
        out.freshIncidenceSlots.reserve(in.cellAllocOrder.size());
#if defined(PAACS_TOPO_TESTING)
        maybeFault();                              // seam: post-reserve
#endif
        for (const auto& [x, y] : in.newEdgeEndpoints) {
            EdgeKey k;
            if (!edges.makeKey(nodes, x, y, k))
                throw std::logic_error(
                    "TransactionAccess::prepareCavity: new-edge endpoint "
                    "is not current");
            auto [it, inserted] = edges.table_.emplace(k, EdgeHandle{});
            if (!inserted)
                throw std::logic_error(
                    "TransactionAccess::prepareCavity: a 'new' edge is "
                    "already in the table (staging inconsistency)");
            out.tabIts.push_back(it);              // non-allocating
            out.tabPlaced.push_back(1);            // non-allocating
#if defined(PAACS_TOPO_TESTING)
            maybeFault();                          // seam: per table node
#endif
        }
        for (const auto s : in.cellAllocOrder) {
            auto it = edges.incident_.find(s);
            if (it == edges.incident_.end()) {
                auto& entry = edges.incident_[s];  // insert the map node
                out.freshIncidenceSlots.push_back(s);   // record FIRST,
                entry.edges.reserve(3);            // reserve AFTER - no
                                                   // untracked ghost if
                                                   // the reserve throws
            } else {
                it->second.edges.reserve(3);
            }
#if defined(PAACS_TOPO_TESTING)
            maybeFault();                          // seam: per incidence
#endif
        }
        out.prepared = true;
    }
    static void abandonCavityPrepared(CellPool& cells, EdgePool& edges,
                                      CavityPrepared& p) noexcept {
        (void)cells;                      // capacity residue stays (legal)
        for (std::size_t i = 0; i < p.tabIts.size(); ++i)
            if (p.tabPlaced[i]) edges.table_.erase(p.tabIts[i]);
        for (const auto s : p.freshIncidenceSlots) edges.incident_.erase(s);
        p.tabIts.clear();
        p.tabPlaced.clear();
        p.freshIncidenceSlots.clear();
        p.prepared = false;
    }
    // ---- 10c APPROVED (amended) commit-side counterpart -----------------
    // Mutates PREPARED pool state only: it does NOT touch revision_,
    // workspace state, or adjacency fields - CavityTransaction alone
    // increments the revision once and sets CommittedAdjacencyStale
    // after this returns. All handle validation, key construction and
    // other potentially throwing work happens BEFORE Phase 3; only
    // no-throw operations remain after this helper begins.
    struct CavityCommitInput {
        // cells (canonical order throughout)
        std::vector<CellHandle> removedCells;
        std::vector<std::uint8_t> removedExhausted;   // parallel
        struct AddedCell {
            std::array<std::int32_t, 3> nodeSlot;
            std::array<std::uint32_t, 3> nodeGen;
            mesh::Classification cls;
        };
        std::vector<AddedCell> added;
        std::vector<std::int32_t> cellSlots;          // Stage-B predicted
        std::vector<std::uint32_t> cellGens;
        std::size_t cellPhase2Claims = 0;
        // edges
        struct Retain {
            EdgeHandle h;
            std::uint32_t finalUse = 0;
            std::uint8_t makeRequired = 0;            // ensure-only
            std::uint8_t setCls = 0;
            mesh::Classification cls;
        };
        struct Create {
            NodeHandle a, b;                          // lo->hi by seq
            std::uint32_t finalUse = 0;
            std::uint8_t required = 0;
            mesh::Classification cls;
        };
        std::vector<Retain> retains;
        std::vector<EdgeHandle> retires;              // canonical order
        std::vector<Create> creates;                  // == tabIts order
        std::vector<std::int32_t> edgeSlots;          // Stage-B predicted
        std::vector<std::uint32_t> edgeGens;
        std::size_t edgePhase2Claims = 0;
        // per-added-cell edge references (created -> index into creates,
        // else the live retained handle), triangle-edge order
        struct ERef {
            std::uint8_t created = 0;
            EdgeHandle h{};
            std::uint32_t idx = 0;
        };
        std::vector<std::array<ERef, 3>> addedEdges;  // parallel to added
        std::vector<ERef> requested;                  // canonical req order
    };
    struct CavityCommitBuffers {                      // ALL preallocated
        std::vector<CellHandle>* resultCells = nullptr;
        std::vector<EdgeHandle>* requiredEdges = nullptr;
        std::vector<EdgeHandle>* createdScratch = nullptr;  // size ==
                                                            // creates
    };
    static void commitCavityPrepared(CellPool& cells, EdgePool& edges,
                                     const CavityCommitInput& in,
                                     CavityPrepared& prep,
                                     CavityCommitBuffers& out) noexcept {
        // ---------------- CELLS: claims -> frees -> in-commit reuse ----
        auto writeCell = [&](std::size_t j) {
            const CellHandle ch = allocPrepared(cells.alloc_);
#if defined(PAACS_TOPO_TESTING)
            if (ch.slot != in.cellSlots[j] || ch.gen != in.cellGens[j])
                std::abort();          // slot AND generation vs Stage B
#endif
            const auto& a = in.added[j];
            auto& r = cells.rec_[static_cast<std::size_t>(ch.slot)];
            r.type = ElementType::Tri3;
            r.count = 3;
            for (int i = 0; i < 3; ++i)
                r.conn[static_cast<std::size_t>(i)] =
                    CellPool::LocalNodeRef{a.nodeSlot[std::size_t(i)],
                                           a.nodeGen[std::size_t(i)]};
            r.cls = a.cls;
            r.seq =
                cells.alloc_.slots_[static_cast<std::size_t>(ch.slot)].seq;
            r.gen = ch.gen;
            cells.linkTail(ch.slot);
            out.resultCells->push_back(ch);           // reserved
        };
        for (std::size_t j = 0; j < in.cellPhase2Claims; ++j) writeCell(j);
        for (std::size_t k = 0; k < in.removedCells.size(); ++k) {
            const CellHandle& h = in.removedCells[k];
            cells.unlink(h.slot);
            freePrepared(cells.alloc_, h);
        }
        for (std::size_t j = in.cellPhase2Claims; j < in.added.size(); ++j)
            writeCell(j);
        cells.bumpPeaks();
        // ---------------- EDGES: claims -> retains -> retires -> reuse -
        auto writeCreate = [&](std::size_t i) {
            const EdgeHandle eh = allocPrepared(edges.alloc_);
#if defined(PAACS_TOPO_TESTING)
            if (eh.slot != in.edgeSlots[i] || eh.gen != in.edgeGens[i])
                std::abort();
#endif
            const auto& cs = in.creates[i];
            auto& r = edges.rec_[static_cast<std::size_t>(eh.slot)];
            r.data.a = cs.a;
            r.data.b = cs.b;
            r.data.cls = cs.cls;                      // ensure-only,
            r.data.requiredConstraint = cs.required != 0;   // BEFORE
            r.data.useCount = cs.finalUse;                  // exposure
            r.key = prep.tabIts[i]->first;
            r.seq =
                edges.alloc_.slots_[static_cast<std::size_t>(eh.slot)].seq;
            r.gen = eh.gen;
            edges.linkTail(eh.slot);
            prep.tabIts[i]->second = eh;              // consume the
            prep.tabPlaced[i] = 0;                    // placeholder
            (*out.createdScratch)[i] = eh;
        };
        for (std::size_t i = 0; i < in.edgePhase2Claims; ++i)
            writeCreate(i);
        for (const auto& rt : in.retains) {           // NET final counts:
            auto& r = edges.rec_[static_cast<std::size_t>(
                raw::slotUnchecked(rt.h))];           // NEVER transit zero
            r.data.useCount = rt.finalUse;
            if (rt.makeRequired) r.data.requiredConstraint = true;
            if (rt.setCls) r.data.cls = rt.cls;
        }
        for (const EdgeHandle& h : in.retires) {
            auto& r =
                edges.rec_[static_cast<std::size_t>(raw::slotUnchecked(h))];
            edges.table_.erase(r.key);                // key erase: no
            edges.unlink(h.slot);                     // growth, no throw
            freePrepared(edges.alloc_, h);
        }
        for (std::size_t i = in.edgePhase2Claims; i < in.creates.size();
             ++i)
            writeCreate(i);
        edges.bumpPeaks();
        // ---------------- incidence: consume/overwrite prepared entries
        auto resolve = [&](const CavityCommitInput::ERef& e) {
            return e.created ? (*out.createdScratch)[e.idx] : e.h;
        };
        for (std::size_t j = 0; j < in.added.size(); ++j) {
            const CellHandle& ch = (*out.resultCells)[j];
            auto it = edges.incident_.find(ch.slot);  // entry EXISTS:
            it->second.gen = ch.gen;                  // reused or Phase-2
            it->second.edges.clear();                 // placeholder
            for (const auto& e : in.addedEdges[j])
                it->second.edges.push_back(resolve(e));   // reserved
        }
        for (std::size_t k = 0; k < in.removedCells.size(); ++k)
            if (in.removedExhausted[k])               // never reallocated:
                edges.incident_.erase(in.removedCells[k].slot);
        // ---------------- requested requirements, canonical order ------
        for (const auto& e : in.requested)
            out.requiredEdges->push_back(resolve(e)); // reserved
        // ---------------- consume ALL placeholder tracking --------------
        prep.tabIts.clear();
        prep.tabPlaced.clear();
        prep.freshIncidenceSlots.clear();
        prep.prepared = false;   // completion is NOT later abandonment
    }

    // ==== BND-2A / ASCR additive helpers (rev-3 §0.5.2) ==================
    // Atomic Steiner Cavity Replacement composes THREE frozen-compatible
    // steps in ONE noexcept Phase 3: (1) materialize a batch of interior
    // Steiner nodes, then (2) apply a prepared cavity whose added triangles
    // reference those now-live nodes. Nothing below changes any existing
    // helper's behaviour: prepareStagedCavity duplicates prepareCavity's
    // body EXACTLY but takes pre-computed staged-aware edge keys (the only
    // point prepareCavity needs live nodes is makeKey, which cannot resolve
    // a not-yet-materialized Steiner node); commitCavityPrepared is reused
    // UNCHANGED (it already supports |added| != |removed| via the
    // claims->frees->reuse ordering, and by Phase 3 the staged endpoints are
    // live). The frozen cavity itself still never creates a node - node
    // creation is this explicit, separately-prepared batch.

    // Phase-2 prediction of the next k node allocations: slot, generation
    // AND sequence. Uses the SINGLE frozen batch path predictSlots(alloc, k)
    // (free-list LIFO, then fresh) - NOT a hand-rolled duplicate - and adds
    // the predicted generation (reused slot: its current gen; fresh slot: 0)
    // and the predicted allocator sequence (nextSeq_ + i). It neither
    // allocates nor consumes a sequence. CALLERS MUST FIRST verify the batch
    // fits (nodeBatchFits): with the space checked, seq = nextSeq_ + i is
    // wrap-free for every i < k. commitNodeBatch asserts each actual value
    // equals its prediction under testing.
    struct StagedNodePredict {
        std::int32_t slot = -1;
        std::uint32_t gen = 0;
        std::uint64_t seq = 0;
    };
    // Full-batch demand check (review blocker 3): both the creation-sequence
    // space AND the slot space must admit k additions BEFORE any prediction
    // narrows to concrete ranks, so Phase 1 can never expose wrapped ranks.
    static bool nodeBatchFits(const NodePool& nodes, std::size_t k) {
        const auto& a = nodes.alloc_;
        if (remainingSequences(a) < k) return false;         // sequence space
        const std::size_t reusable = a.reusableFreeSlots();
        const std::size_t fresh = k > reusable ? k - reusable : 0;
        // fresh slots materialize at materializedSlots()+1 .. +fresh; the
        // last must not exceed the allocator's slot-space bound.
        for (std::size_t i = 1; i <= fresh; ++i)
            if (HandleAllocator<NodeTag>::atSlotCapacity(
                    a.materializedSlots() + i - 1))
                return false;                                // slot space
        return true;
    }
    static std::vector<StagedNodePredict> predictNodeBatch(
        const NodePool& nodes, std::size_t k) {
        const auto& a = nodes.alloc_;
        const auto slots = predictSlots(a, k);   // THE single batch path
        const std::size_t materialized =
            static_cast<std::size_t>(a.materializedSlots());
        const std::uint64_t base = a.nextSeq_;
        std::vector<StagedNodePredict> out;
        out.reserve(k);
        for (std::size_t i = 0; i < k; ++i) {
            StagedNodePredict p;
            p.slot = slots[i];
            p.gen = static_cast<std::size_t>(p.slot) < materialized
                        ? a.slots_[static_cast<std::size_t>(p.slot)].gen
                        : 0;
            p.seq = base + static_cast<std::uint64_t>(i);   // wrap-free
            out.push_back(p);
        }
        return out;
    }
    // Phase-2: reserve node allocator capacity for the predicted batch and
    // grow rec_ through the highest predicted slot, so Phase-3
    // materialization is allocation-free. Capacity only; no handle minted,
    // no sequence consumed. No node retirement occurs (Steiner adds only),
    // so no free-list capacity is reserved.
    static void reserveNodeBatch(NodePool& nodes,
                                 const std::vector<StagedNodePredict>& pred) {
        auto& a = nodes.alloc_;
        std::size_t fresh = 0, maxSlot = nodes.rec_.size();
        for (const auto& p : pred) {
            const std::size_t s = static_cast<std::size_t>(p.slot);
            maxSlot = std::max(maxSlot, s + 1);
            if (s >= a.materializedSlots()) ++fresh;   // fresh vs reused
        }
        a.reserveMaterializedCapacity(a.materializedSlots() + fresh);
        if (maxSlot > nodes.rec_.size()) {
            detail::reserveFor(nodes.rec_, maxSlot);
            nodes.rec_.resize(maxSlot);
        }
    }
    struct StagedNodeSpec {
        Vec3 coord{0, 0, 0};
        mesh::Classification cls{};
    };
    // Phase-3 (noexcept, allocation-free): materialize the predicted batch,
    // writing each node's payload, linking it live, and appending the live
    // handle to `out` (caller pre-reserved). Asserts actual slot/gen/seq ==
    // prediction (single active transaction).
    static void commitNodeBatch(NodePool& nodes,
                                const std::vector<StagedNodeSpec>& specs,
                                const std::vector<StagedNodePredict>& pred,
                                std::vector<NodeHandle>& out) noexcept {
#if !defined(PAACS_TOPO_TESTING)
        (void)pred;
#endif
        for (std::size_t i = 0; i < specs.size(); ++i) {
            const NodeHandle h = allocPrepared(nodes.alloc_);
#if defined(PAACS_TOPO_TESTING)
            if (h.slot != pred[i].slot || h.gen != pred[i].gen)
                std::abort();
            if (nodes.alloc_.slots_[static_cast<std::size_t>(h.slot)].seq !=
                pred[i].seq)
                std::abort();
#endif
            auto& r = nodes.rec_[static_cast<std::size_t>(h.slot)];
            r.data.coord = specs[i].coord;
            r.data.cls = specs[i].cls;
            r.seq =
                nodes.alloc_.slots_[static_cast<std::size_t>(h.slot)].seq;
            r.gen = h.gen;
            nodes.linkTail(h.slot);
            out.push_back(h);                          // reserved
        }
        nodes.bumpPeaks();
    }

    // Canonical EdgeKey for an edge whose endpoints may be staged (not-yet-
    // materialized) Steiner nodes: existing endpoints use their live
    // sequence; staged endpoints use the Phase-1 predicted sequence. Returns
    // false on a self-edge (equal sequences), exactly like EdgePool::makeKey.
    static bool stagedEdgeKey(const NodePool& nodes, bool aStaged,
                              std::uint64_t aSeq, const NodeHandle& aH,
                              bool bStaged, std::uint64_t bSeq,
                              const NodeHandle& bH, EdgeKey& out) {
        const std::uint64_t sa = aStaged ? aSeq : nodes.seqOf(aH);
        const std::uint64_t sb = bStaged ? bSeq : nodes.seqOf(bH);
        if (sa == sb) return false;
        out.lo = std::min(sa, sb);
        out.hi = std::max(sa, sb);
        return true;
    }

    // Staged-aware cavity preparation. IDENTICAL to prepareCavity except the
    // new-edge table placeholders are keyed by pre-computed staged-aware keys
    // (in.newEdgeKeys) instead of makeKey(nodes,x,y) - a Steiner endpoint is
    // not current yet, so makeKey would (correctly) refuse it. Every other
    // step - allocator capacity, rec_ sizing, bucket reserve, cleanup-vector
    // reserve, incidence placeholders, the failure seams - is byte-for-byte
    // the frozen prepareCavity behaviour.
    struct CavityStagedPrepareInput {
        std::vector<std::int32_t> cellAllocOrder;
        std::vector<std::int32_t> edgeAllocOrder;
        std::vector<EdgeKey> newEdgeKeys;      // staged-aware, canonical
        std::size_t cellFreePushes = 0;
        std::size_t edgeFreePushes = 0;
        std::size_t cellFresh = 0;
        std::size_t edgeFresh = 0;
    };
    static void prepareStagedCavity(CellPool& cells, EdgePool& edges,
                                    const CavityStagedPrepareInput& in,
                                    CavityPrepared& out) {
#if defined(PAACS_TOPO_TESTING)
        maybeFault();                              // seam: entry
#endif
        cells.alloc_.reserveMaterializedCapacity(
            cells.alloc_.materializedSlots() + in.cellFresh);
        cells.alloc_.reserveFreeCapacity(in.cellFreePushes);
        edges.alloc_.reserveMaterializedCapacity(
            edges.alloc_.materializedSlots() + in.edgeFresh);
        edges.alloc_.reserveFreeCapacity(in.edgeFreePushes);
        {
            std::size_t need = cells.rec_.size();
            for (const auto s : in.cellAllocOrder)
                need = std::max(need, static_cast<std::size_t>(s) + 1);
            if (need > cells.rec_.size()) {
                detail::reserveFor(cells.rec_, need);
                cells.rec_.resize(need);
            }
        }
        {
            std::size_t need = edges.rec_.size();
            for (const auto s : in.edgeAllocOrder)
                need = std::max(need, static_cast<std::size_t>(s) + 1);
            if (need > edges.rec_.size()) {
                detail::reserveFor(edges.rec_, need);
                edges.rec_.resize(need);
            }
        }
#if defined(PAACS_TOPO_TESTING)
        maybeFault();                              // seam: post-capacity
#endif
        edges.table_.reserve(edges.table_.size() + in.newEdgeKeys.size());
        out.tabIts.reserve(in.newEdgeKeys.size());
        out.tabPlaced.reserve(in.newEdgeKeys.size());
        out.freshIncidenceSlots.reserve(in.cellAllocOrder.size());
#if defined(PAACS_TOPO_TESTING)
        maybeFault();                              // seam: post-reserve
#endif
        for (const EdgeKey& k : in.newEdgeKeys) {
            auto [it, inserted] = edges.table_.emplace(k, EdgeHandle{});
            if (!inserted)
                throw std::logic_error(
                    "TransactionAccess::prepareStagedCavity: a 'new' edge "
                    "key is already in the table (staging inconsistency)");
            out.tabIts.push_back(it);
            out.tabPlaced.push_back(1);
#if defined(PAACS_TOPO_TESTING)
            maybeFault();                          // seam: per table node
#endif
        }
        for (const auto s : in.cellAllocOrder) {
            auto it = edges.incident_.find(s);
            if (it == edges.incident_.end()) {
                auto& entry = edges.incident_[s];
                out.freshIncidenceSlots.push_back(s);
                entry.edges.reserve(3);
            } else {
                it->second.edges.reserve(3);
            }
#if defined(PAACS_TOPO_TESTING)
            maybeFault();                          // seam: per incidence
#endif
        }
        out.prepared = true;
    }

#if defined(PAACS_TOPO_TESTING)
    // Fixture seam DISCLOSED BEYOND the three approved names (flagged
    // for the reviewer): the Stage-B review requires genuine
    // exhausted-retirement fixtures, and a LIVE cell with an exhausted
    // generation needs BOTH the allocator and the record generation set
    // (record generations rebuild handles in activeBySequence). Testing
    // builds only; no production surface.
    // Testing-only remaining-sequence setter (Stage-B evidence
    // correction: the nonzero off-by-one boundary needs remaining ==
    // demand and remaining == demand - 1 with demand > 0).
    template <class Tag>
    static void testSetRemainingSequences(HandleAllocator<Tag>& a,
                                          std::uint64_t remaining) {
        a.nextSeq_ = std::numeric_limits<std::uint64_t>::max() - remaining;
    }
    static CellHandle testExhaustLiveCell(CellPool& cells,
                                          const CellHandle& h) {
        const auto s = static_cast<std::size_t>(cells.alloc_.deref(h));
        cells.alloc_.testForceGeneration(
            h, std::numeric_limits<std::uint32_t>::max());
        cells.rec_[s].gen = std::numeric_limits<std::uint32_t>::max();
        return CellHandle{static_cast<std::int32_t>(s),
                          std::numeric_limits<std::uint32_t>::max(),
                          cells.alloc_.poolId()};
    }
#endif

    // Probe reach (testing builds only use these; harmless privates).
    template <class Tag>
    static std::size_t slotCapacityOf(const HandleAllocator<Tag>& a) {
        return a.slots_.capacity();
    }
    // Mutable node allocator accessor (BND-2A additive; symmetric completion
    // of the cell/edge/face mutable accessors that already exist - the node
    // one is needed ONLY by the sequence-boundary TEST fixtures, so unlike
    // its siblings it is compiled out of the production surface entirely.
#if defined(PAACS_TOPO_TESTING)
    static HandleAllocator<NodeTag>& nodeAlloc(NodePool& n) { return n.alloc_; }
#endif
    static HandleAllocator<CellTag>& cellAlloc(CellPool& c) { return c.alloc_; }
    static HandleAllocator<EdgeTag>& edgeAlloc(EdgePool& e) { return e.alloc_; }
    static HandleAllocator<FaceTag>& faceAlloc(FacePool& f) { return f.alloc_; }
    static const HandleAllocator<NodeTag>& nodeAllocC(const NodePool& n) {
        return n.alloc_;
    }
    static const HandleAllocator<CellTag>& cellAllocC(const CellPool& c) {
        return c.alloc_;
    }
    static const HandleAllocator<EdgeTag>& edgeAllocC(const EdgePool& e) {
        return e.alloc_;
    }
    static const HandleAllocator<FaceTag>& faceAllocC(const FacePool& f) {
        return f.alloc_;
    }

    struct PoolBindings {
        std::uint64_t cellNode = 0, edgeNode = 0, edgeCell = 0,
                      faceNode = 0, faceCell = 0;
    };
    static PoolBindings debugBindings(const CellPool& cells,
                                      const EdgePool& edges,
                                      const FacePool& faces) {
        return PoolBindings{cells.boundNodePool_, edges.boundNodePool_,
                            edges.boundCellPool_, faces.boundNodePool_,
                            faces.boundCellPool_};
    }

    // ---- Plans (Phase-1 output; owned by TopologyTransaction) --------------
    struct EdgeOp {
        bool create = false;
        EdgeKey key{};
        NodeHandle a{}, b{};                 // canonical lo->hi by seq
        EdgeHandle existing{};               // reuse only
        decltype(EdgePool::table_)::iterator tabIt{};   // create placeholder
        bool tabPlaced = false;              // placeholder actually emplaced
        std::int32_t predictedSlot = -1;     // create only
    };
    struct FaceOp {
        bool create = false;
        FacePool::Staged staged{};           // key/shape/n/outward
        FaceHandle existing{};               // reuse only
        std::int8_t sense = +1;
        decltype(FacePool::table_)::iterator tabIt{};   // create placeholder
        bool tabPlaced = false;              // placeholder actually emplaced
        std::int32_t predictedSlot = -1;     // create only
    };
    struct AddPlan {
        ElementType type{};
        std::uint8_t nodeCount = 0;
        std::array<NodeHandle, 10> nodes{};
        mesh::Classification cls{};
        std::vector<EdgeOp> edges;           // traits local-edge order
        std::vector<FaceOp> faces;           // traits local-face order
        std::int32_t cellSlot = -1;          // predicted
        std::uint32_t cellGen = 0;           // predicted
        // First-time pool bindings deferred to Phase 3 (review blocker:
        // binding in Phase 2 left BEHAVIORAL state behind on failure,
        // invisible to both hashes). Computed in Phase 1, written by the
        // no-throw commit as plain integer assignments.
        bool setCellBinding = false;
        bool setEdgeBinding = false;
        bool setFaceBinding = false;
        // Incidence placeholder ITERATORS captured at Phase-2 emplacement
        // and dereferenced in Phase 3 (review: no topology lookup may
        // begin after commit starts - find() during commit is forbidden
        // even when it cannot throw).
        decltype(EdgePool::incident_)::iterator edgeIncIt{};
        decltype(FacePool::cellFaces_)::iterator faceCfIt{};
        PoolAccount cellsAcct{}, edgesAcct{}, facesAcct{};
        bool prepared = false;
        bool cellIncidencePlaced = false;    // edge incident_ placeholder
        bool cellFacesPlaced = false;        // face cellFaces_ placeholder
    };
    enum class RelKind : std::uint8_t { Retained, AddReuse, AddCreate };
    struct REdgeOp {
        RelKind kind = RelKind::AddCreate;
        EdgeKey key{};
        NodeHandle a{}, b{};                 // canonical lo->hi (creates)
        EdgeHandle existing{};               // Retained / AddReuse
        decltype(EdgePool::table_)::iterator tabIt{};
        bool tabPlaced = false;
        std::int32_t predictedSlot = -1;
    };
    struct RFaceOp {
        RelKind kind = RelKind::AddCreate;
        FacePool::Staged staged{};
        FaceHandle existing{};               // Retained / AddReuse
        std::int8_t sense = +1;
        std::size_t retainedIncidIdx = 0;    // Retained only
        decltype(FacePool::table_)::iterator tabIt{};
        bool tabPlaced = false;
        std::int32_t predictedSlot = -1;
    };
    struct RemoveEdge {
        EdgeHandle h{};
        bool retire = false;
        decltype(EdgePool::table_)::iterator tabIt{};   // retire only
    };
    struct RemoveFace {
        FaceHandle h{};
        bool retire = false;
        std::size_t incidIdx = 0;            // this cell's Incidence index
        decltype(FacePool::table_)::iterator tabIt{};   // retire only
    };
    struct RemovePlan {
        CellHandle cell{};
        std::vector<RemoveEdge> edges;
        std::vector<RemoveFace> faces;
        decltype(EdgePool::incident_)::iterator incidIt{};
        bool hasIncid = false;
        decltype(FacePool::cellFaces_)::iterator cfIt{};
        bool hasCf = false;
        PoolAccount cellsAcct{}, edgesAcct{}, facesAcct{};
        bool prepared = false;
    };

    struct ReplacePlan {
        CellHandle cell{};                   // PRESERVED: same pool/slot/
                                             // gen/sequence, new connectivity
        ElementType type{};                  // unchanged (v1)
        std::uint8_t nodeCount = 0;
        std::array<NodeHandle, 10> newNodes{};
        std::vector<REdgeOp> edges;          // NEW incidence, traits order
        std::vector<RemoveEdge> removedEdges;    // obsolete incidence
        std::vector<RFaceOp> faces;
        std::vector<RemoveFace> removedFaces;
        decltype(EdgePool::incident_)::iterator incidIt{};
        decltype(FacePool::cellFaces_)::iterator cfIt{};
        bool hasCf = false;
        PoolAccount cellsAcct{}, edgesAcct{}, facesAcct{};
        bool prepared = false;
    };

    // ---- Phase 1: planAddCell (throwing; pools untouched) ------------------
    static AddPlan planAddCell(const NodePool& nodes, const CellPool& cells,
                               const EdgePool& edges, const FacePool& faces,
                               ElementType type,
                               std::span<const NodeHandle> conn,
                               const mesh::Classification& cls,
                               const TopologyPolicy& policy) {
        AddPlan p;
        const ElementTraits& tr = checkedTraits(type);
        if (static_cast<int>(conn.size()) != tr.nodeCount)
            throw std::invalid_argument(
                "transaction addCell: connectivity count does not match "
                "element type");
        for (const NodeHandle& nh : conn)
            if (nh.pool != nodes.poolId() || !nodes.isCurrent(nh))
                throw std::invalid_argument(
                    "transaction addCell: stale or foreign node handle");
        if (!validateHandleConnectivity(type, conn,
                [&](const NodeHandle& nh) { return nodes.isCurrent(nh); }))
            throw std::invalid_argument(
                "transaction addCell: invalid connectivity (repeated or "
                "degenerate nodes)");
        // Structured CELL sequence precheck (always: even a reused slot
        // consumes a sequence). Edge/face sequence checks are CONDITIONAL
        // and happen after delta derivation below (review: an operation
        // that reuses all existing edges/faces allocates no handles there
        // and must not be rejected for unrelated exhausted spaces).
        if (sequencesExhausted(cells.alloc_))
            throw std::length_error(
                "transaction addCell: cell handle sequence space exhausted");
        // Desired pool-ownership validation + first-time-binding plan
        // (Phase 1 reads; Phase 3 writes).
        if ((cells.boundNodePool_ != 0 &&
             cells.boundNodePool_ != nodes.poolId()) ||
            (edges.boundNodePool_ != 0 &&
             (edges.boundNodePool_ != nodes.poolId() ||
              edges.boundCellPool_ != cells.poolId())) ||
            (faces.boundNodePool_ != 0 &&
             (faces.boundNodePool_ != nodes.poolId() ||
              faces.boundCellPool_ != cells.poolId())))
            throw std::invalid_argument(
                "transaction addCell: pools are bound to foreign pools");
        p.setCellBinding = (cells.boundNodePool_ == 0);
        p.setEdgeBinding = (edges.boundNodePool_ == 0);
        p.setFaceBinding = (faces.boundNodePool_ == 0);
        p.type = type;
        p.nodeCount = static_cast<std::uint8_t>(conn.size());
        for (std::size_t i = 0; i < conn.size(); ++i) p.nodes[i] = conn[i];
        p.cls = cls;
        // Edge deltas, traits local order.
        for (int e = 0; e < tr.edgeCount; ++e) {
            EdgeOp op;
            NodeHandle x = conn[static_cast<std::size_t>(tr.edges[e].a)];
            NodeHandle y = conn[static_cast<std::size_t>(tr.edges[e].b)];
            if (nodes.seqOf(x) > nodes.seqOf(y)) std::swap(x, y);
            op.a = x;
            op.b = y;
            op.key.lo = nodes.seqOf(x);
            op.key.hi = nodes.seqOf(y);
            if (auto f = edges.find(nodes, x, y)) {
                op.create = false;
                op.existing = *f;
            } else {
                op.create = true;
            }
            p.edges.push_back(op);
        }
        // Face deltas + orientation + policy, traits local order.
        for (int f = 0; f < tr.faceCount; ++f) {
            const LocalFace& lf = tr.faces[f];
            FaceOp op;
            op.staged.shape = lf.shape;
            op.staged.n = static_cast<std::uint8_t>(lf.n);
            for (int i = 0; i < lf.n; ++i)
                op.staged.outward[static_cast<std::size_t>(i)] =
                    conn[static_cast<std::size_t>(
                        lf.v[static_cast<std::size_t>(i)])];
            std::array<std::uint64_t, kMaxFaceNodes> ss{};
            for (int i = 0; i < lf.n; ++i)
                ss[static_cast<std::size_t>(i)] = nodes.seqOf(
                    op.staged.outward[static_cast<std::size_t>(i)]);
            std::sort(ss.begin(), ss.begin() + lf.n);
            op.staged.key.n = static_cast<std::uint8_t>(lf.n);
            for (int i = 0; i < lf.n; ++i)
                op.staged.key.s[static_cast<std::size_t>(i)] =
                    ss[static_cast<std::size_t>(i)];
            const auto tit = faces.table_.find(op.staged.key);
            if (tit != faces.table_.end()) {
                op.create = false;
                op.existing = tit->second;
                op.sense = faceSense(faces, nodes, op.existing, op.staged);
                if (op.sense == 0)
                    throw std::invalid_argument(
                        "transaction addCell: face corner set matches but "
                        "winding is not a rotation or reflection");
                const std::size_t inc =
                    faceIncidenceCount(faces, op.existing);
                if (inc + 1 >= 3 &&
                    policy.nonManifold == NonManifoldPolicy::Reject)
                    throw std::invalid_argument(
                        "transaction addCell: would create nonmanifold "
                        "face incidence under NonManifoldPolicy::Reject");
                // PRODUCT form (fix of a latent stage-3 defect, found
                // while designing replacement: the old check assumed the
                // existing incidence has sense +1, which is false when
                // the face's creating cell was later removed - the
                // survivor's sense can be -1, and the new cell must then
                // be +1).
                if (inc == 1 && policy.requireOppositeInteriorOrientation &&
                    op.sense * faceSoleSense(faces, op.existing) != -1)
                    throw std::invalid_argument(
                        "transaction addCell: interior face would not have "
                        "opposite orientation senses");
            } else {
                op.create = true;
                op.sense = +1;
            }
            p.faces.push_back(op);
        }
        {
            std::size_t newEdges = 0, newFaces = 0;
            for (const auto& op : p.edges) newEdges += op.create ? 1 : 0;
            for (const auto& op : p.faces) newFaces += op.create ? 1 : 0;
            if (newEdges > 0 && sequencesExhausted(edges.alloc_))
                throw std::length_error(
                    "transaction addCell: edge handle sequence space "
                    "exhausted");
            if (newFaces > 0 && sequencesExhausted(faces.alloc_))
                throw std::length_error(
                    "transaction addCell: face handle sequence space "
                    "exhausted");
        }
        return p;
    }

    // ---- Phase 2: prepareAddCell (throwing; capacity + placeholders) -------
    static void prepareAddCell(NodePool& nodes, CellPool& cells,
                               EdgePool& edges, FacePool& faces,
                               AddPlan& p) {
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        (void)nodes;
        // Bindings are NOT touched here (review blocker): Phase 2 may
        // leave only capacity-class residue on failure; binding is
        // behavioral state and is written by the no-throw Phase 3.
        // Cell: slot capacity + payload sizing + predicted identity.
        const auto cellPred = predictSlots(cells.alloc_, 1);
        p.cellSlot = cellPred[0];
        p.cellGen = slotGen(cells.alloc_, p.cellSlot);
        p.cellsAcct = account(cells.alloc_, 1, 0);
        cells.alloc_.reserveMaterializedCapacity(p.cellsAcct.slotCapacity);
        if (static_cast<std::size_t>(p.cellSlot) >= cells.rec_.size()) {
            detail::reserveFor(cells.rec_,
                               static_cast<std::size_t>(p.cellSlot) + 1);
            cells.rec_.resize(static_cast<std::size_t>(p.cellSlot) + 1);
        }
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        // Edges: fresh slot capacity, payload sizing, table placeholders.
        std::size_t newEdges = 0;
        for (const auto& op : p.edges) newEdges += op.create ? 1 : 0;
        p.edgesAcct = account(edges.alloc_, newEdges, 0);
        edges.alloc_.reserveMaterializedCapacity(p.edgesAcct.slotCapacity);
        const auto edgePred = predictSlots(edges.alloc_, newEdges);
        std::size_t need = edges.rec_.size();
        for (auto s : edgePred)
            need = std::max(need, static_cast<std::size_t>(s) + 1);
        if (need > edges.rec_.size()) {
            detail::reserveFor(edges.rec_, need);
            edges.rec_.resize(need);
        }
        edges.table_.reserve(edges.table_.size() + newEdges);
        {
            std::size_t pi = 0;
            for (auto& op : p.edges)
                if (op.create) {
                    op.predictedSlot = edgePred[pi++];
                    op.tabIt =
                        edges.table_.emplace(op.key, EdgeHandle{}).first;
                    op.tabPlaced = true;
                }
        }
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        // Edge incidence placeholder for the predicted cell slot; reserve
        // the list so Phase-3 pushes cannot allocate.
        {
            auto [it, inserted] = edges.incident_.emplace(
                p.cellSlot, EdgePool::IncidenceEntry{});
            it->second.gen = p.cellGen;
            it->second.edges.reserve(p.edges.size());
            p.edgeIncIt = it;                // dereferenced by Phase 3
            p.cellIncidencePlaced = inserted;
        }
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        // Faces: fresh slot capacity, payload sizing (incl. incident-vector
        // reserves at PREDICTED slots for creates, +1 for reuses), table
        // placeholders.
        std::size_t newFaces = 0;
        for (const auto& op : p.faces) newFaces += op.create ? 1 : 0;
        p.facesAcct = account(faces.alloc_, newFaces, 0);
        faces.alloc_.reserveMaterializedCapacity(p.facesAcct.slotCapacity);
        const auto facePred = predictSlots(faces.alloc_, newFaces);
        need = faces.rec_.size();
        for (auto s : facePred)
            need = std::max(need, static_cast<std::size_t>(s) + 1);
        if (need > faces.rec_.size()) {
            detail::reserveFor(faces.rec_, need);
            faces.rec_.resize(need);
        }
        faces.table_.reserve(faces.table_.size() + newFaces);
        {
            std::size_t pi = 0;
            for (auto& op : p.faces) {
                if (op.create) {
                    op.predictedSlot = facePred[pi++];
                    faces.rec_[static_cast<std::size_t>(op.predictedSlot)]
                        .data.incident.reserve(2);
                    op.tabIt = faces.table_
                                   .emplace(op.staged.key, FaceHandle{})
                                   .first;
                    op.tabPlaced = true;
                } else {
                    auto& fr = faces.rec_[static_cast<std::size_t>(
                        raw::slotUnchecked(op.existing))];
                    fr.data.incident.reserve(fr.data.incident.size() + 1);
                }
            }
        }
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        // Face cellFaces placeholder for the predicted cell slot.
        {
            auto [it, inserted] = faces.cellFaces_.emplace(
                p.cellSlot, FacePool::CellFaceEntry{});
            it->second.gen = p.cellGen;
            it->second.faces.reserve(p.faces.size());
            p.faceCfIt = it;                 // dereferenced by Phase 3
            p.cellFacesPlaced = inserted;
        }
        p.prepared = true;
    }

    // ---- abandon after prepare: erase placeholders (noexcept) --------------
    // unordered_map::erase(iterator) is noexcept; the int32-keyed erases
    // hash with std::hash<int32> (noexcept).
    static void abandonAddPrepared(EdgePool& edges, FacePool& faces,
                                   AddPlan& p) noexcept {
        // Erase EXACTLY the placeholders that were emplaced - a fault can
        // interrupt prepare mid-way, so per-op flags (not p.prepared,
        // which is set only at the END of prepare) decide what exists
        // (review-loop lesson: a mid-prepare fault previously left ghost
        // placeholder table entries behind, failing the audits).
        for (auto& op : p.edges)
            if (op.tabPlaced) { edges.table_.erase(op.tabIt);
                                op.tabPlaced = false; }
        for (auto& op : p.faces)
            if (op.tabPlaced) { faces.table_.erase(op.tabIt);
                                op.tabPlaced = false; }
        if (p.cellIncidencePlaced) edges.incident_.erase(p.edgeIncIt);
        if (p.cellFacesPlaced) faces.cellFaces_.erase(p.faceCfIt);
        p.prepared = false;
        p.cellIncidencePlaced = p.cellFacesPlaced = false;
    }

    // ---- Phase 3: commitAddCell (noexcept, validation- & alloc-free) -------
    static CellHandle commitAddCell(NodePool& nodes, CellPool& cells,
                                    EdgePool& edges, FacePool& faces,
                                    AddPlan& p) noexcept {
        // First-time pool bindings: plain integer writes, planned in
        // Phase 1, applied only by the commit that cannot fail.
        if (p.setCellBinding) cells.boundNodePool_ = nodes.poolId();
        if (p.setEdgeBinding) {
            edges.boundNodePool_ = nodes.poolId();
            edges.boundCellPool_ = cells.poolId();
        }
        if (p.setFaceBinding) {
            faces.boundNodePool_ = nodes.poolId();
            faces.boundCellPool_ = cells.poolId();
        }
        const CellHandle ch = allocPrepared(cells.alloc_);
#if defined(PAACS_TOPO_TESTING)
        if (ch.slot != p.cellSlot || ch.gen != p.cellGen) std::abort();
#endif
        {
            auto& r = cells.rec_[static_cast<std::size_t>(ch.slot)];
            r.type = p.type;
            r.count = p.nodeCount;
            for (std::uint8_t i = 0; i < p.nodeCount; ++i)
                r.conn[i] = CellPool::LocalNodeRef{p.nodes[i].slot,
                                                   p.nodes[i].gen};
            r.cls = p.cls;
            r.seq =
                cells.alloc_.slots_[static_cast<std::size_t>(ch.slot)].seq;
            r.gen = ch.gen;
        }
        cells.linkTail(ch.slot);
        cells.bumpPeaks();
        // Edge deltas + incidence: the PREPARED iterator is dereferenced;
        // no hash lookup occurs after commit begins (review contract).
        auto incIt = p.edgeIncIt;
        incIt->second.gen = ch.gen;
        for (auto& op : p.edges) {
            EdgeHandle eh;
            if (op.create) {
                eh = allocPrepared(edges.alloc_);
#if defined(PAACS_TOPO_TESTING)
                if (eh.slot != op.predictedSlot) std::abort();
#endif
                auto& r = edges.rec_[static_cast<std::size_t>(eh.slot)];
                r.data.a = op.a;
                r.data.b = op.b;
                r.data.cls = mesh::Classification{};
                r.data.requiredConstraint = false;
                r.data.useCount = 1;
                r.key = op.key;
                r.seq = edges.alloc_
                            .slots_[static_cast<std::size_t>(eh.slot)]
                            .seq;
                r.gen = eh.gen;
                edges.linkTail(eh.slot);
                op.tabIt->second = eh;       // overwrite placeholder value
                op.tabPlaced = false;        // consumed
            } else {
                eh = op.existing;
                ++edges
                      .rec_[static_cast<std::size_t>(raw::slotUnchecked(eh))]
                      .data.useCount;
            }
            incIt->second.edges.push_back(eh);       // reserved: noexcept
        }
        edges.bumpPeaks();
        // Face deltas + incidence: prepared iterator, no lookup.
        auto cfIt = p.faceCfIt;
        cfIt->second.gen = ch.gen;
        for (auto& op : p.faces) {
            FaceHandle fh;
            if (op.create) {
                fh = allocPrepared(faces.alloc_);
#if defined(PAACS_TOPO_TESTING)
                if (fh.slot != op.predictedSlot) std::abort();
#endif
                auto& r = faces.rec_[static_cast<std::size_t>(fh.slot)];
                r.data.incident.clear();     // keeps Phase-2 reserve
                r.data.shape = op.staged.shape;
                r.data.cornerCount = op.staged.n;
                for (std::uint8_t i = 0; i < op.staged.n; ++i)
                    r.data.corners[i] = op.staged.outward[i];
                r.data.cls = mesh::Classification{};
                r.key = op.staged.key;
                r.seq = faces.alloc_
                            .slots_[static_cast<std::size_t>(fh.slot)]
                            .seq;
                r.gen = fh.gen;
                faces.linkTail(fh.slot);
                op.tabIt->second = fh;
                op.tabPlaced = false;        // consumed
            } else {
                fh = op.existing;
            }
            faces.rec_[static_cast<std::size_t>(raw::slotUnchecked(fh))]
                .data.incident.push_back(FacePool::Incidence{
                    p.cellSlot, ch.gen, op.sense});  // reserved: noexcept
            cfIt->second.faces.push_back(fh);        // reserved: noexcept
        }
        faces.bumpPeaks();
        p.prepared = false;                  // consumed
        p.cellIncidencePlaced = p.cellFacesPlaced = false;
        return ch;
    }

    // ---- M4-B.6C node lifecycle -----------------------------------------
    // 6C-1 creation: one node, no implicit cells/edges/faces. A brand-new
    // node is BY DEFINITION isolated, so under the default
    // permitOrphanNodes == false it is REJECTED (structured) - callers
    // creating nodes opt in with permitOrphanNodes = true (v1's
    // one-operation transactions cannot pair the node with a referencing
    // cell; multi-operation cavity transactions will).
    struct AddNodePlan {
        Vec3 coord{};
        mesh::Classification cls{};
        std::int32_t predictedSlot = -1;
        std::uint32_t predictedGen = 0;
        PoolAccount nodesAcct{};
        bool prepared = false;
    };
    static AddNodePlan planAddNode(const NodePool& nodes, const Vec3& p0,
                                   const mesh::Classification& c,
                                   const TopologyPolicy& policy) {
        // RECORDED DECISION (M4-B.6C closure review): the transaction
        // layer rejects nonfinite coordinates in Phase 1 - NaN/Inf
        // poison later hashing, orientation and geometry computations,
        // so they are a structured plan-time rejection rather than a
        // deferred validator concern. (The low-level NodePool::add
        // contract is unchanged.)
        if (!std::isfinite(p0.x) || !std::isfinite(p0.y) ||
            !std::isfinite(p0.z))
            throw std::invalid_argument(
                "transaction addNode: coordinates must be finite");
        if (!policy.permitOrphanNodes)
            throw std::invalid_argument(
                "transaction addNode: a new node is necessarily isolated; "
                "set TopologyPolicy::permitOrphanNodes to create nodes in "
                "a single-operation transaction");
        if (sequencesExhausted(nodes.alloc_))
            throw std::length_error(
                "transaction addNode: node handle sequence space "
                "exhausted");
        AddNodePlan p;
        p.coord = p0;
        p.cls = c;
        return p;
    }
    static void prepareAddNode(NodePool& nodes, AddNodePlan& p) {
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        p.nodesAcct = account(nodes.alloc_, 1, 0);
        const auto pred = predictSlots(nodes.alloc_, 1);
        p.predictedSlot = pred[0];
        p.predictedGen = slotGen(nodes.alloc_, p.predictedSlot);
        nodes.alloc_.reserveMaterializedCapacity(p.nodesAcct.slotCapacity);
        if (static_cast<std::size_t>(p.predictedSlot) >= nodes.rec_.size()) {
            detail::reserveFor(nodes.rec_,
                               static_cast<std::size_t>(p.predictedSlot) + 1);
            nodes.rec_.resize(static_cast<std::size_t>(p.predictedSlot) + 1);
        }
        p.prepared = true;
    }
    static NodeHandle commitAddNode(NodePool& nodes,
                                    AddNodePlan& p) noexcept {
        const NodeHandle nh = allocPrepared(nodes.alloc_);
#if defined(PAACS_TOPO_TESTING)
        if (nh.slot != p.predictedSlot || nh.gen != p.predictedGen)
            std::abort();
#endif
        auto& r = nodes.rec_[static_cast<std::size_t>(nh.slot)];
        r.data.coord = p.coord;
        r.data.cls = p.cls;
        r.seq = nodes.alloc_.slots_[static_cast<std::size_t>(nh.slot)].seq;
        r.gen = nh.gen;
        nodes.linkTail(nh.slot);
        nodes.bumpPeaks();
        p.prepared = false;
        return nh;
    }

    // 6C-2 removal: only a node with NO surviving reference anywhere -
    // current cells, live edge endpoints, live face corners. The plan
    // requires the composite workspace state (caller-checked) and then
    // INDEPENDENTLY confirms zero references; no cascade deletion (that
    // is future multi-operation cavity work).
    struct RemoveNodePlan {
        NodeHandle node{};
        PoolAccount nodesAcct{};
        bool prepared = false;
    };
    static RemoveNodePlan planRemoveNode(const NodePool& nodes,
                                         const CellPool& cells,
                                         const EdgePool& edges,
                                         const FacePool& faces,
                                         const NodeHandle& h) {
        if (h.pool != nodes.poolId() || !nodes.isCurrent(h))
            throw std::invalid_argument(
                "transaction removeNode: stale or foreign node handle");
        for (const CellHandle& ch : cells.activeBySequence()) {
            const auto snap = cells.get(ch);
            for (std::uint8_t i = 0; i < snap.count; ++i)
                if (snap.nodes[i] == h)
                    throw std::invalid_argument(
                        "transaction removeNode: node is referenced by a "
                        "current cell");
        }
        for (const EdgeHandle& eh : edges.activeBySequence()) {
            const auto d = edges.get(eh);
            if (d.a == h || d.b == h)
                throw std::invalid_argument(
                    "transaction removeNode: node is an endpoint of a "
                    "live edge");
        }
        for (const FaceHandle& fh : faces.activeBySequence()) {
            const auto d = faces.get(fh);
            for (std::uint8_t i = 0; i < d.cornerCount; ++i)
                if (d.corners[i] == h)
                    throw std::invalid_argument(
                        "transaction removeNode: node is a corner of a "
                        "live face");
        }
        RemoveNodePlan p;
        p.node = h;
        return p;
    }
    static void prepareRemoveNode(NodePool& nodes, RemoveNodePlan& p) {
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        p.nodesAcct = account(nodes.alloc_, 0, 1);
        nodes.alloc_.reserveFreeCapacity(1);
        p.prepared = true;
    }
    static void commitRemoveNode(NodePool& nodes,
                                 RemoveNodePlan& p) noexcept {
        nodes.unlink(p.node.slot);
        freePrepared(nodes.alloc_, p.node);
        p.prepared = false;
    }

    // ---- M4-B.6B metadata (existing compact classification IDs ONLY;
    // NO implicit propagation: a transaction modifies exactly the entity
    // it names - the GroupId registry remains PM-SOLV-001) ---------------
    enum class MetaKind : std::uint8_t {
        NodeCls, CellCls, EdgeCls, FaceCls, EdgeConstraint
    };
    struct MetadataPlan {
        MetaKind kind = MetaKind::NodeCls;
        NodeHandle node{};
        CellHandle cell{};
        EdgeHandle edge{};
        FaceHandle face{};
        mesh::Classification newCls{}, oldCls{};
        bool newConstraint = false, oldConstraint = false;
        bool noOp = false;                   // new value == current value
        bool prepared = false;
    };
    static bool clsEq(const mesh::Classification& a,
                      const mesh::Classification& b) {
        return a.dim == b.dim && a.tag == b.tag;
    }
    // Phase 1 planners (throwing; pools untouched).
    static MetadataPlan planSetNodeCls(const NodePool& nodes,
                                       const NodeHandle& h,
                                       const mesh::Classification& c) {
        if (h.pool != nodes.poolId() || !nodes.isCurrent(h))
            throw std::invalid_argument(
                "transaction setNodeClassification: stale or foreign "
                "node handle");
        MetadataPlan p;
        p.kind = MetaKind::NodeCls;
        p.node = h;
        p.newCls = c;
        p.oldCls = nodes.get(h).cls;
        p.noOp = clsEq(p.newCls, p.oldCls);
        return p;
    }
    static MetadataPlan planSetCellCls(const CellPool& cells,
                                       const CellHandle& h,
                                       const mesh::Classification& c) {
        if (h.pool != cells.poolId() || !cells.isCurrent(h))
            throw std::invalid_argument(
                "transaction setCellClassification: stale or foreign "
                "cell handle");
        MetadataPlan p;
        p.kind = MetaKind::CellCls;
        p.cell = h;
        p.newCls = c;
        p.oldCls = cells.get(h).cls;
        p.noOp = clsEq(p.newCls, p.oldCls);
        return p;
    }
    static MetadataPlan planSetEdgeCls(const EdgePool& edges,
                                       const EdgeHandle& h,
                                       const mesh::Classification& c) {
        if (h.pool != edges.poolId() || !edges.isCurrent(h))
            throw std::invalid_argument(
                "transaction setEdgeClassification: stale or foreign "
                "edge handle");
        MetadataPlan p;
        p.kind = MetaKind::EdgeCls;
        p.edge = h;
        p.newCls = c;
        p.oldCls = edges.get(h).cls;
        p.noOp = clsEq(p.newCls, p.oldCls);
        return p;
    }
    static MetadataPlan planSetFaceCls(const FacePool& faces,
                                       const FaceHandle& h,
                                       const mesh::Classification& c) {
        if (h.pool != faces.poolId() || !faces.isCurrent(h))
            throw std::invalid_argument(
                "transaction setFaceClassification: stale or foreign "
                "face handle");
        MetadataPlan p;
        p.kind = MetaKind::FaceCls;
        p.face = h;
        p.newCls = c;
        p.oldCls = faces.get(h).cls;
        p.noOp = clsEq(p.newCls, p.oldCls);
        return p;
    }
    static MetadataPlan planSetEdgeConstraint(const EdgePool& edges,
                                              const EdgeHandle& h,
                                              bool required) {
        if (h.pool != edges.poolId() || !edges.isCurrent(h))
            throw std::invalid_argument(
                "transaction setEdgeRequiredConstraint: stale or foreign "
                "edge handle");
        MetadataPlan p;
        p.kind = MetaKind::EdgeConstraint;
        p.edge = h;
        p.newConstraint = required;
        p.oldConstraint = edges.get(h).requiredConstraint;
        p.noOp = (p.newConstraint == p.oldConstraint);
        return p;
    }
    // Phase 2: a field overwrite needs NO allocation preparation; the
    // seam gives the failure-path tests their injection point.
    static void prepareMetadata(MetadataPlan& p) {
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        p.prepared = true;
    }
    // Phase 3: one in-place field write (evidence row: no throw source).
    static void commitMetadata(NodePool& nodes, CellPool& cells,
                               EdgePool& edges, FacePool& faces,
                               MetadataPlan& p) noexcept {
        switch (p.kind) {
            case MetaKind::NodeCls:
                nodes.rec_[static_cast<std::size_t>(
                               raw::slotUnchecked(p.node))]
                    .data.cls = p.newCls;
                break;
            case MetaKind::CellCls:
                cells.rec_[static_cast<std::size_t>(
                               raw::slotUnchecked(p.cell))]
                    .cls = p.newCls;
                break;
            case MetaKind::EdgeCls:
                edges.rec_[static_cast<std::size_t>(
                               raw::slotUnchecked(p.edge))]
                    .data.cls = p.newCls;
                break;
            case MetaKind::FaceCls:
                faces.rec_[static_cast<std::size_t>(
                               raw::slotUnchecked(p.face))]
                    .data.cls = p.newCls;
                break;
            case MetaKind::EdgeConstraint:
                edges.rec_[static_cast<std::size_t>(
                               raw::slotUnchecked(p.edge))]
                    .data.requiredConstraint = p.newConstraint;
                break;
        }
        p.prepared = false;
    }

    // ---- Phase 1: planReplaceCell (throwing; NET deltas; pools untouched)
    static ReplacePlan planReplaceCell(const NodePool& nodes,
                                       CellPool& cells, EdgePool& edges,
                                       FacePool& faces,
                                       const CellHandle& cell,
                                       std::span<const NodeHandle> conn,
                                       const TopologyPolicy& policy) {
        ReplacePlan p;
        if (cell.pool != cells.poolId() || !cells.isCurrent(cell))
            throw std::invalid_argument(
                "transaction replaceCellConnectivity: stale or foreign "
                "cell handle");
        const auto snap = cells.get(cell);
        const ElementTraits& tr = checkedTraits(snap.type);
        if (!edges.hasCellEdges(cells, cell) ||
            (tr.faceCount > 0 && !faces.hasCellFaces(cells, cell)))
            throw std::invalid_argument(
                "transaction replaceCellConnectivity: cell incidence is "
                "incomplete");
        if (static_cast<int>(conn.size()) != tr.nodeCount)
            throw std::invalid_argument(
                "transaction replaceCellConnectivity: connectivity count "
                "does not match the cell's element type");
        for (const NodeHandle& nh : conn)
            if (nh.pool != nodes.poolId() || !nodes.isCurrent(nh))
                throw std::invalid_argument(
                    "transaction replaceCellConnectivity: stale or foreign "
                    "node handle");
        if (!validateHandleConnectivity(snap.type, conn,
                [&](const NodeHandle& nh) { return nodes.isCurrent(nh); }))
            throw std::invalid_argument(
                "transaction replaceCellConnectivity: invalid replacement "
                "connectivity");
        p.cell = cell;
        p.type = snap.type;
        p.nodeCount = static_cast<std::uint8_t>(conn.size());
        for (std::size_t i = 0; i < conn.size(); ++i)
            p.newNodes[i] = conn[i];
        // Orphan policy: old nodes DROPPED by the replacement.
        if (!policy.permitOrphanNodes) {
            for (std::uint8_t i = 0; i < snap.count; ++i) {
                const NodeHandle nh = snap.nodes[i];
                bool kept = false;
                for (const NodeHandle& nn : conn)
                    if (nn == nh) { kept = true; break; }
                if (kept) continue;
                bool referenced = false;
                for (const CellHandle& other : cells.activeBySequence()) {
                    if (other.slot == cell.slot) continue;
                    const auto os = cells.get(other);
                    for (std::uint8_t j = 0; j < os.count; ++j)
                        if (os.nodes[j] == nh) { referenced = true; break; }
                    if (referenced) break;
                }
                if (!referenced)
                    throw std::invalid_argument(
                        "transaction replaceCellConnectivity: would orphan "
                        "a node (TopologyPolicy::permitOrphanNodes is "
                        "false)");
            }
        }
        // OLD edge incidence -> key map (NET deltas, not remove+re-add).
        p.incidIt = edges.incident_.find(cell.slot);
        std::vector<std::pair<EdgeKey, EdgeHandle>> oldEdges;
        std::vector<bool> oldEdgeRetained;
        for (const EdgeHandle& eh : p.incidIt->second.edges) {
            oldEdges.emplace_back(
                edges.rec_[static_cast<std::size_t>(raw::slotUnchecked(eh))]
                    .key,
                eh);
            oldEdgeRetained.push_back(false);
        }
        for (int e = 0; e < tr.edgeCount; ++e) {
            REdgeOp op;
            NodeHandle x = conn[static_cast<std::size_t>(tr.edges[e].a)];
            NodeHandle y = conn[static_cast<std::size_t>(tr.edges[e].b)];
            if (nodes.seqOf(x) > nodes.seqOf(y)) std::swap(x, y);
            op.a = x; op.b = y;
            op.key.lo = nodes.seqOf(x);
            op.key.hi = nodes.seqOf(y);
            bool retained = false;
            for (std::size_t k = 0; k < oldEdges.size(); ++k)
                if (oldEdges[k].first.lo == op.key.lo &&
                    oldEdges[k].first.hi == op.key.hi) {
                    op.kind = RelKind::Retained;
                    op.existing = oldEdges[k].second;
                    oldEdgeRetained[k] = true;
                    retained = true;
                    break;
                }
            if (!retained) {
                if (auto f = edges.find(nodes, x, y)) {
                    op.kind = RelKind::AddReuse;
                    op.existing = *f;
                } else {
                    op.kind = RelKind::AddCreate;
                }
            }
            p.edges.push_back(op);
        }
        for (std::size_t k = 0; k < oldEdges.size(); ++k) {
            if (oldEdgeRetained[k]) continue;
            RemoveEdge re;
            re.h = oldEdges[k].second;
            const auto& r = edges.rec_[static_cast<std::size_t>(
                raw::slotUnchecked(re.h))];
            re.retire = (r.data.useCount == 1);
            if (re.retire) re.tabIt = edges.table_.find(r.key);
            p.removedEdges.push_back(re);
        }
        // OLD face incidence -> key map.
        p.hasCf = (tr.faceCount > 0);
        std::vector<FaceHandle> oldFaces;
        std::vector<std::size_t> oldFaceIdx;
        std::vector<bool> oldFaceRetained;
        if (p.hasCf) {
            p.cfIt = faces.cellFaces_.find(cell.slot);
            for (const FaceHandle& fh : p.cfIt->second.faces) {
                const auto& r = faces.rec_[static_cast<std::size_t>(
                    raw::slotUnchecked(fh))];
                std::size_t idx = 0;
                for (std::size_t k = 0; k < r.data.incident.size(); ++k)
                    if (r.data.incident[k].cellSlot == cell.slot &&
                        r.data.incident[k].cellGen == cell.gen) {
                        idx = k;
                        break;
                    }
                oldFaces.push_back(fh);
                oldFaceIdx.push_back(idx);
                oldFaceRetained.push_back(false);
            }
        }
        for (int f = 0; f < tr.faceCount; ++f) {
            const LocalFace& lf = tr.faces[f];
            RFaceOp op;
            op.staged.shape = lf.shape;
            op.staged.n = static_cast<std::uint8_t>(lf.n);
            for (int i = 0; i < lf.n; ++i)
                op.staged.outward[static_cast<std::size_t>(i)] =
                    conn[static_cast<std::size_t>(
                        lf.v[static_cast<std::size_t>(i)])];
            std::array<std::uint64_t, kMaxFaceNodes> ss{};
            for (int i = 0; i < lf.n; ++i)
                ss[static_cast<std::size_t>(i)] = nodes.seqOf(
                    op.staged.outward[static_cast<std::size_t>(i)]);
            std::sort(ss.begin(), ss.begin() + lf.n);
            op.staged.key.n = static_cast<std::uint8_t>(lf.n);
            for (int i = 0; i < lf.n; ++i)
                op.staged.key.s[static_cast<std::size_t>(i)] =
                    ss[static_cast<std::size_t>(i)];
            bool retained = false;
            for (std::size_t k = 0; k < oldFaces.size(); ++k) {
                const auto& r = faces.rec_[static_cast<std::size_t>(
                    raw::slotUnchecked(oldFaces[k]))];
                if (r.key == op.staged.key) {
                    op.kind = RelKind::Retained;
                    op.existing = oldFaces[k];
                    op.retainedIncidIdx = oldFaceIdx[k];
                    op.sense = faceSense(faces, nodes, op.existing,
                                         op.staged);
                    if (op.sense == 0)
                        throw std::invalid_argument(
                            "transaction replaceCellConnectivity: retained "
                            "face winding is not a rotation or reflection");
                    // Interior retained face: the UPDATED sense must still
                    // oppose the other cell's sense.
                    if (policy.requireOppositeInteriorOrientation &&
                        r.data.incident.size() == 2) {
                        const auto& other = r.data.incident
                            [oldFaceIdx[k] == 0 ? 1 : 0];
                        if (op.sense * other.sense != -1)
                            throw std::invalid_argument(
                                "transaction replaceCellConnectivity: "
                                "retained interior face would not have "
                                "opposite orientation senses");
                    }
                    oldFaceRetained[k] = true;
                    retained = true;
                    break;
                }
            }
            if (!retained) {
                const auto tit = faces.table_.find(op.staged.key);
                if (tit != faces.table_.end()) {
                    op.kind = RelKind::AddReuse;
                    op.existing = tit->second;
                    op.sense = faceSense(faces, nodes, op.existing,
                                         op.staged);
                    if (op.sense == 0)
                        throw std::invalid_argument(
                            "transaction replaceCellConnectivity: face "
                            "corner set matches but winding is not a "
                            "rotation or reflection");
                    const std::size_t inc =
                        faceIncidenceCount(faces, op.existing);
                    if (inc + 1 >= 3 &&
                        policy.nonManifold == NonManifoldPolicy::Reject)
                        throw std::invalid_argument(
                            "transaction replaceCellConnectivity: would "
                            "create nonmanifold face incidence under "
                            "NonManifoldPolicy::Reject");
                    if (inc == 1 &&
                        policy.requireOppositeInteriorOrientation &&
                        op.sense * faceSoleSense(faces, op.existing) != -1)
                        throw std::invalid_argument(
                            "transaction replaceCellConnectivity: interior "
                            "face would not have opposite orientation "
                            "senses");
                } else {
                    op.kind = RelKind::AddCreate;
                    op.sense = +1;
                }
            }
            p.faces.push_back(op);
        }
        for (std::size_t k = 0; k < oldFaces.size(); ++k) {
            if (oldFaceRetained[k]) continue;
            RemoveFace rf;
            rf.h = oldFaces[k];
            rf.incidIdx = oldFaceIdx[k];
            const auto& r = faces.rec_[static_cast<std::size_t>(
                raw::slotUnchecked(rf.h))];
            rf.retire = (r.data.incident.size() == 1);
            if (rf.retire) rf.tabIt = faces.table_.find(r.key);
            p.removedFaces.push_back(rf);
        }
        // Conditional sequence prechecks: replacement allocates NO cell
        // handle (identity preserved); edges/faces only when creating.
        {
            std::size_t ne = 0, nf = 0;
            for (const auto& op : p.edges)
                ne += op.kind == RelKind::AddCreate ? 1 : 0;
            for (const auto& op : p.faces)
                nf += op.kind == RelKind::AddCreate ? 1 : 0;
            if (ne > 0 && sequencesExhausted(edges.alloc_))
                throw std::length_error(
                    "transaction replaceCellConnectivity: edge handle "
                    "sequence space exhausted");
            if (nf > 0 && sequencesExhausted(faces.alloc_))
                throw std::length_error(
                    "transaction replaceCellConnectivity: face handle "
                    "sequence space exhausted");
        }
        return p;
    }

    // ---- Phase 2: prepareReplaceCell ----------------------------------------
    static void prepareReplaceCell(NodePool& nodes, CellPool& cells,
                                   EdgePool& edges, FacePool& faces,
                                   ReplacePlan& p) {
        (void)nodes;
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        std::size_t ne = 0, er = 0;
        for (const auto& op : p.edges)
            ne += op.kind == RelKind::AddCreate ? 1 : 0;
        for (const auto& re : p.removedEdges) er += re.retire ? 1 : 0;
        p.cellsAcct = account(cells.alloc_, 0, 0);   // handle PRESERVED
        p.edgesAcct = account(edges.alloc_, ne, er);
        edges.alloc_.reserveMaterializedCapacity(p.edgesAcct.slotCapacity);
        edges.alloc_.reserveFreeCapacity(er);
        const auto edgePred = predictSlots(edges.alloc_, ne);
        std::size_t need = edges.rec_.size();
        for (auto sl : edgePred)
            need = std::max(need, static_cast<std::size_t>(sl) + 1);
        if (need > edges.rec_.size()) {
            detail::reserveFor(edges.rec_, need);
            edges.rec_.resize(need);
        }
        edges.table_.reserve(edges.table_.size() + ne);
        {
            std::size_t pi = 0;
            for (auto& op : p.edges)
                if (op.kind == RelKind::AddCreate) {
                    op.predictedSlot = edgePred[pi++];
                    op.tabIt =
                        edges.table_.emplace(op.key, EdgeHandle{}).first;
                    op.tabPlaced = true;
                }
        }
        p.incidIt->second.edges.reserve(p.edges.size());
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        std::size_t nf = 0, fr = 0;
        for (const auto& op : p.faces)
            nf += op.kind == RelKind::AddCreate ? 1 : 0;
        for (const auto& rf : p.removedFaces) fr += rf.retire ? 1 : 0;
        p.facesAcct = account(faces.alloc_, nf, fr);
        faces.alloc_.reserveMaterializedCapacity(p.facesAcct.slotCapacity);
        faces.alloc_.reserveFreeCapacity(fr);
        const auto facePred = predictSlots(faces.alloc_, nf);
        need = faces.rec_.size();
        for (auto sl : facePred)
            need = std::max(need, static_cast<std::size_t>(sl) + 1);
        if (need > faces.rec_.size()) {
            detail::reserveFor(faces.rec_, need);
            faces.rec_.resize(need);
        }
        faces.table_.reserve(faces.table_.size() + nf);
        {
            std::size_t pi = 0;
            for (auto& op : p.faces) {
                if (op.kind == RelKind::AddCreate) {
                    op.predictedSlot = facePred[pi++];
                    faces.rec_[static_cast<std::size_t>(op.predictedSlot)]
                        .data.incident.reserve(2);
                    op.tabIt = faces.table_
                                   .emplace(op.staged.key, FaceHandle{})
                                   .first;
                    op.tabPlaced = true;
                } else if (op.kind == RelKind::AddReuse) {
                    auto& fr2 = faces.rec_[static_cast<std::size_t>(
                        raw::slotUnchecked(op.existing))];
                    fr2.data.incident.reserve(fr2.data.incident.size() + 1);
                }
            }
        }
        if (p.hasCf) p.cfIt->second.faces.reserve(p.faces.size());
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        cells.alloc_.reserveFreeCapacity(0);          // explicit: none
        p.prepared = true;
    }

    static void abandonReplacePrepared(EdgePool& edges, FacePool& faces,
                                       ReplacePlan& p) noexcept {
        for (auto& op : p.edges)
            if (op.tabPlaced) { edges.table_.erase(op.tabIt);
                                op.tabPlaced = false; }
        for (auto& op : p.faces)
            if (op.tabPlaced) { faces.table_.erase(op.tabIt);
                                op.tabPlaced = false; }
        p.prepared = false;
    }

    // ---- Phase 3: commitReplaceCell (noexcept; reviewer's order) ------------
    static void commitReplaceCell(NodePool& nodes, CellPool& cells,
                                  EdgePool& edges, FacePool& faces,
                                  ReplacePlan& p) noexcept {
        (void)nodes;
        // 1. Remove obsolete face incidence (retained faces untouched, so
        //    their captured indices stay valid).
        for (auto& rf : p.removedFaces) {
            auto& r = faces.rec_[static_cast<std::size_t>(
                raw::slotUnchecked(rf.h))];
            r.data.incident.erase(
                r.data.incident.begin() +
                static_cast<std::ptrdiff_t>(rf.incidIdx));
        }
        // 2. Remove obsolete edge incidence (use counts).
        for (auto& re : p.removedEdges) {
            auto& r = edges.rec_[static_cast<std::size_t>(
                raw::slotUnchecked(re.h))];
            if (r.data.useCount > 0) --r.data.useCount;
        }
        // 3. Overwrite cell connectivity IN PLACE: same pool, slot,
        //    generation and creation sequence - the handle is preserved.
        {
            auto& r = cells.rec_[static_cast<std::size_t>(p.cell.slot)];
            for (std::uint8_t i = 0; i < p.nodeCount; ++i)
                r.conn[i] = CellPool::LocalNodeRef{p.newNodes[i].slot,
                                                   p.newNodes[i].gen};
        }
        // 4. Apply retained/new edge incidence (cleared list keeps its
        //    capacity; every push is reserved).
        p.incidIt->second.edges.clear();
        for (auto& op : p.edges) {
            EdgeHandle eh;
            if (op.kind == RelKind::AddCreate) {
                eh = allocPrepared(edges.alloc_);
#if defined(PAACS_TOPO_TESTING)
                if (eh.slot != op.predictedSlot) std::abort();
#endif
                auto& r = edges.rec_[static_cast<std::size_t>(eh.slot)];
                r.data.a = op.a;
                r.data.b = op.b;
                r.data.cls = mesh::Classification{};
                r.data.requiredConstraint = false;
                r.data.useCount = 1;
                r.key = op.key;
                r.seq = edges.alloc_
                            .slots_[static_cast<std::size_t>(eh.slot)]
                            .seq;
                r.gen = eh.gen;
                edges.linkTail(eh.slot);
                op.tabIt->second = eh;
                op.tabPlaced = false;
            } else if (op.kind == RelKind::AddReuse) {
                eh = op.existing;
                ++edges
                      .rec_[static_cast<std::size_t>(raw::slotUnchecked(eh))]
                      .data.useCount;
            } else {                          // Retained
                eh = op.existing;
            }
            p.incidIt->second.edges.push_back(eh);
        }
        edges.bumpPeaks();
        // 5. Apply retained/new face incidence.
        if (p.hasCf) {
            p.cfIt->second.faces.clear();
            for (auto& op : p.faces) {
                FaceHandle fh;
                if (op.kind == RelKind::AddCreate) {
                    fh = allocPrepared(faces.alloc_);
#if defined(PAACS_TOPO_TESTING)
                    if (fh.slot != op.predictedSlot) std::abort();
#endif
                    auto& r =
                        faces.rec_[static_cast<std::size_t>(fh.slot)];
                    r.data.incident.clear();
                    r.data.shape = op.staged.shape;
                    r.data.cornerCount = op.staged.n;
                    for (std::uint8_t i = 0; i < op.staged.n; ++i)
                        r.data.corners[i] = op.staged.outward[i];
                    r.data.cls = mesh::Classification{};
                    r.key = op.staged.key;
                    r.seq = faces.alloc_
                                .slots_[static_cast<std::size_t>(fh.slot)]
                                .seq;
                    r.gen = fh.gen;
                    faces.linkTail(fh.slot);
                    op.tabIt->second = fh;
                    op.tabPlaced = false;
                    r.data.incident.push_back(FacePool::Incidence{
                        p.cell.slot, p.cell.gen, op.sense});
                } else if (op.kind == RelKind::AddReuse) {
                    fh = op.existing;
                    faces.rec_[static_cast<std::size_t>(
                                   raw::slotUnchecked(fh))]
                        .data.incident.push_back(FacePool::Incidence{
                            p.cell.slot, p.cell.gen, op.sense});
                } else {                      // Retained: update the sense
                    fh = op.existing;
                    faces.rec_[static_cast<std::size_t>(
                                   raw::slotUnchecked(fh))]
                        .data.incident[op.retainedIncidIdx]
                        .sense = op.sense;
                }
                p.cfIt->second.faces.push_back(fh);
            }
            faces.bumpPeaks();
        }
        // 6. Retire obsolete records (AFTER creates: retiring slots are
        //    deliberately NOT counted reusable for this commit's creates).
        for (auto& rf : p.removedFaces)
            if (rf.retire) {
                faces.table_.erase(rf.tabIt);
                faces.unlink(rf.h.slot);
                freePrepared(faces.alloc_, rf.h);
            }
        for (auto& re : p.removedEdges)
            if (re.retire) {
                edges.table_.erase(re.tabIt);
                edges.unlink(re.h.slot);
                freePrepared(edges.alloc_, re.h);
            }
        p.prepared = false;
    }

    // ---- Phase 1: planRemoveCell (throwing; pools untouched) ---------------
    static RemovePlan planRemoveCell(const NodePool& nodes, CellPool& cells,
                                     EdgePool& edges, FacePool& faces,
                                     const CellHandle& cell,
                                     const TopologyPolicy& policy) {
        (void)nodes;
        RemovePlan p;
        if (cell.pool != cells.poolId() || !cells.isCurrent(cell))
            throw std::invalid_argument(
                "transaction removeCell: stale or foreign cell handle");
        const auto snap = cells.get(cell);
        if (!edges.hasCellEdges(cells, cell) ||
            (checkedTraits(snap.type).faceCount > 0 &&
             !faces.hasCellFaces(cells, cell)))
            throw std::invalid_argument(
                "transaction removeCell: cell incidence is incomplete");
        p.cell = cell;
        // Orphan-node policy (Phase 1 may be expensive; Phase 3 may not).
        if (!policy.permitOrphanNodes) {
            for (std::uint8_t i = 0; i < snap.count; ++i) {
                const NodeHandle nh = snap.nodes[i];
                bool referenced = false;
                for (const CellHandle& other : cells.activeBySequence()) {
                    if (other.slot == cell.slot) continue;
                    const auto os = cells.get(other);
                    for (std::uint8_t j = 0; j < os.count; ++j)
                        if (os.nodes[j] == nh) {
                            referenced = true;
                            break;
                        }
                    if (referenced) break;
                }
                if (!referenced)
                    throw std::invalid_argument(
                        "transaction removeCell: would orphan a node "
                        "(TopologyPolicy::permitOrphanNodes is false)");
            }
        }
        // Edge consequences; iterators captured while everything is live.
        p.incidIt = edges.incident_.find(cell.slot);
        p.hasIncid = p.incidIt != edges.incident_.end();
        if (p.hasIncid)
            for (const EdgeHandle& eh : p.incidIt->second.edges) {
                RemoveEdge re;
                re.h = eh;
                const auto& r = edges.rec_[static_cast<std::size_t>(
                    raw::slotUnchecked(eh))];
                re.retire = (r.data.useCount == 1);
                if (re.retire) re.tabIt = edges.table_.find(r.key);
                p.edges.push_back(re);
            }
        // Face consequences.
        p.cfIt = faces.cellFaces_.find(cell.slot);
        p.hasCf = p.cfIt != faces.cellFaces_.end();
        if (p.hasCf)
            for (const FaceHandle& fh : p.cfIt->second.faces) {
                RemoveFace rf;
                rf.h = fh;
                const auto& r = faces.rec_[static_cast<std::size_t>(
                    raw::slotUnchecked(fh))];
                rf.retire = (r.data.incident.size() == 1);
                for (std::size_t k = 0; k < r.data.incident.size(); ++k)
                    if (r.data.incident[k].cellSlot == cell.slot &&
                        r.data.incident[k].cellGen == cell.gen) {
                        rf.incidIdx = k;
                        break;
                    }
                if (rf.retire) rf.tabIt = faces.table_.find(r.key);
                p.faces.push_back(rf);
            }
        return p;
    }

    // ---- Phase 2: prepareRemoveCell (free-list capacity only) --------------
    static void prepareRemoveCell(CellPool& cells, EdgePool& edges,
                                  FacePool& faces, RemovePlan& p) {
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        std::size_t er = 0, fr = 0;
        for (const auto& e : p.edges) er += e.retire ? 1 : 0;
        for (const auto& f : p.faces) fr += f.retire ? 1 : 0;
        p.cellsAcct = account(cells.alloc_, 0, 1);
        p.edgesAcct = account(edges.alloc_, 0, er);
        p.facesAcct = account(faces.alloc_, 0, fr);
        cells.alloc_.reserveFreeCapacity(1);
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        edges.alloc_.reserveFreeCapacity(er);
        faces.alloc_.reserveFreeCapacity(fr);
        p.prepared = true;
    }

    // ---- Phase 3: commitRemoveCell (noexcept; reviewer's order) ------------
    static void commitRemoveCell(CellPool& cells, EdgePool& edges,
                                 FacePool& faces, RemovePlan& p) noexcept {
        // 1. Face incidence first; everything was captured in Phase 1 -
        //    no query or traits recomputation happens here.
        for (auto& rf : p.faces) {
            auto& r = faces.rec_[static_cast<std::size_t>(
                raw::slotUnchecked(rf.h))];
            r.data.incident.erase(
                r.data.incident.begin() +
                static_cast<std::ptrdiff_t>(rf.incidIdx));   // memmove
            if (rf.retire) {
                faces.table_.erase(rf.tabIt);    // iterator erase: noexcept
                faces.unlink(rf.h.slot);
                freePrepared(faces.alloc_, rf.h);
            }
        }
        if (p.hasCf) faces.cellFaces_.erase(p.cfIt);
        // 2. Edge incidence.
        for (auto& re : p.edges) {
            auto& r = edges.rec_[static_cast<std::size_t>(
                raw::slotUnchecked(re.h))];
            if (r.data.useCount > 0) --r.data.useCount;
            if (re.retire) {
                edges.table_.erase(re.tabIt);
                edges.unlink(re.h.slot);
                freePrepared(edges.alloc_, re.h);
            }
        }
        if (p.hasIncid) edges.incident_.erase(p.incidIt);
        // 3. Cell retirement.
        cells.unlink(p.cell.slot);
        freePrepared(cells.alloc_, p.cell);
        p.prepared = false;
    }

    // ------------------------------------------------------------------
    // T2-A accepted revision-2 contact: complete Tet4 cavity net delta.
    // The unfrozen TetCavityTransaction is the only caller. Phase 1 builds
    // canonical transition rows, Phase 2 places all hash/incidence
    // placeholders without minting handles, and Phase 3 consumes captured
    // slots/iterators only.
    struct TetEdgeRankKey {
        std::uint64_t lo = 0, hi = 0;
        bool operator<(const TetEdgeRankKey& o) const {
            return std::tie(lo, hi) < std::tie(o.lo, o.hi);
        }
    };
    struct TetFaceRankKey {
        std::array<std::uint64_t, 3> ranks{};
        bool operator<(const TetFaceRankKey& o) const {
            return ranks < o.ranks;
        }
    };
    struct TetEdgeRow {
        EdgeKey key{};
        NodeHandle a{}, b{};
        EdgeHandle existing{};
        EdgeHandle finalHandle{};
        std::uint32_t initialUse = 0, removedUse = 0, addedUse = 0,
                      finalUse = 0;
        bool create = false, retire = false, target = false;
        std::int32_t predictedSlot = -1;
        std::uint32_t predictedGen = 0;
        decltype(EdgePool::table_)::iterator tabIt{};
        bool tabPlaced = false;
    };
    struct TetFaceRow {
        FaceKey key{};
        FacePool::Staged staged{};
        FaceHandle existing{};
        FaceHandle finalHandle{};
        std::vector<FacePool::Incidence> surviving;
        std::uint32_t initialUse = 0, removedUse = 0, addedUse = 0,
                      finalUse = 0;
        bool create = false, retire = false, target = false;
        std::int32_t predictedSlot = -1;
        std::uint32_t predictedGen = 0;
        decltype(FacePool::table_)::iterator tabIt{};
        bool tabPlaced = false;
    };
    struct TetAddRow {
        std::array<NodeHandle, 4> nodes{};
        mesh::Classification cls{};
        std::array<std::size_t, 6> edges{};
        std::array<std::size_t, 4> faces{};
        std::array<std::int8_t, 4> faceSenses{};
        std::int32_t predictedSlot = -1;
        std::uint32_t predictedGen = 0;
        CellHandle materialized{};
        decltype(EdgePool::incident_)::iterator edgeIncIt{};
        decltype(FacePool::cellFaces_)::iterator faceCfIt{};
        bool edgeIncPlaced = false, faceCfPlaced = false;
    };
    struct TetRemoveRow {
        CellHandle cell{};
        decltype(EdgePool::incident_)::iterator edgeIncIt{};
        decltype(FacePool::cellFaces_)::iterator faceCfIt{};
        bool hasEdgeInc = false, hasFaceCf = false;
    };
    struct TetPlan {
        std::vector<TetRemoveRow> removes;
        std::vector<TetAddRow> adds;
        std::vector<TetEdgeRow> edges;
        std::vector<TetFaceRow> faces;
        std::size_t targetFaceRow = std::numeric_limits<std::size_t>::max();
        std::size_t targetEdgeRow = std::numeric_limits<std::size_t>::max();
        mesh::Classification targetClassification{};
        bool edgeTargetMode = false;
        bool prepared = false;
    };

    struct TetCommitBuffers {
        std::vector<CellHandle>* resultCells = nullptr;
        FaceHandle* targetFace = nullptr;
        EdgeHandle* targetEdge = nullptr;
    };

    // T2-B target-edge metadata mode.  Existing-edge marking is one atomic
    // two-field publication on the already-authorized Tet transaction; it
    // requires no table, allocator or incidence preparation.
    struct TetEdgeMetadataPlan {
        EdgeHandle edge{};
        mesh::Classification classification{};
        bool required = true;
        bool prepared = false;
    };

    static TetEdgeMetadataPlan planTetEdgeMetadata(
        const EdgePool& edges, EdgeHandle edge,
        mesh::Classification classification, bool required) {
        if (!edges.isCurrent(edge) || classification.tag < 0 ||
            (classification.dim != 1 && classification.dim != 2) ||
            !required)
            throw std::invalid_argument("invalid Tet target-edge metadata");
        return {edge,classification,required,false};
    }

    static void prepareTetEdgeMetadata(TetEdgeMetadataPlan& plan) {
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        plan.prepared = true;
    }

    static void commitTetEdgeMetadata(
        EdgePool& edges, TetEdgeMetadataPlan& plan) noexcept {
        auto& data = edges.rec_[static_cast<std::size_t>(
            raw::slotUnchecked(plan.edge))].data;
        data.cls = plan.classification;
        data.requiredConstraint = plan.required;
        plan.prepared = false;
    }

    static TetEdgeRankKey tetEdgeKey(const NodePool& nodes, NodeHandle a,
                                     NodeHandle b) {
        auto ra = nodes.seqOf(a), rb = nodes.seqOf(b);
        if (rb < ra) {
            std::swap(a, b);
            std::swap(ra, rb);
        }
        return {ra, rb};
    }
    static TetFaceRankKey tetFaceKey(
        const NodePool& nodes,
        const std::array<NodeHandle, 3>& directed) {
        TetFaceRankKey key{{nodes.seqOf(directed[0]),
                            nodes.seqOf(directed[1]),
                            nodes.seqOf(directed[2])}};
        std::sort(key.ranks.begin(), key.ranks.end());
        return key;
    }
    static FacePool::Staged tetStagedFace(
        const NodePool& nodes,
        const std::array<NodeHandle, 3>& directed) {
        FacePool::Staged staged;
        staged.shape = FaceShape::Tri;
        staged.n = 3;
        for (std::size_t i = 0; i < 3; ++i)
            staged.outward[i] = directed[i];
        staged.key.n = 3;
        const auto key = tetFaceKey(nodes, directed);
        for (std::size_t i = 0; i < 3; ++i)
            staged.key.s[i] = key.ranks[i];
        return staged;
    }
    static std::int8_t tetFaceSenseBetween(
        const NodePool& nodes, const FacePool::Staged& stored,
        const FacePool::Staged& use) {
        std::array<std::uint64_t, 3> a{}, b{};
        for (std::size_t i = 0; i < 3; ++i) {
            a[i] = nodes.seqOf(stored.outward[i]);
            b[i] = nodes.seqOf(use.outward[i]);
        }
        for (std::size_t shift = 0; shift < 3; ++shift)
            if (b[0] == a[shift] && b[1] == a[(shift + 1) % 3] &&
                b[2] == a[(shift + 2) % 3])
                return +1;
        for (std::size_t shift = 0; shift < 3; ++shift)
            if (b[0] == a[shift] && b[1] == a[(shift + 2) % 3] &&
                b[2] == a[(shift + 1) % 3])
                return -1;
        return 0;
    }

    static TetPlan planTetCavityImpl(
        const NodePool& nodes, CellPool& cells, EdgePool& edges,
        FacePool& faces, std::span<const CellHandle> removed,
        std::span<const std::array<NodeHandle, 4>> added,
        std::span<const mesh::Classification> addedClasses,
        const std::array<NodeHandle, 3>& target,
        const mesh::Classification& targetClassification,
        const mesh::Classification& region, bool edgeTargetMode,
        const NodeHandle* stagedNode = nullptr,
        std::uint64_t stagedSequence = 0) {
        const bool refinementMode = stagedNode != nullptr;
        if (removed.empty() || added.empty() ||
            added.size() != addedClasses.size())
            throw std::invalid_argument("Tet cavity requires nonempty, paired input");
        if (((!refinementMode || edgeTargetMode) &&
             ((edgeTargetMode
                   ? (targetClassification.dim != 1 &&
                      targetClassification.dim != 2)
                   : targetClassification.dim != 2) ||
              targetClassification.tag < 0)) ||
            region.dim != 3 || region.tag < 0)
            throw std::invalid_argument("Tet cavity classification domain");
        auto nodeRank = [&](NodeHandle h) -> std::uint64_t {
            return refinementMode && h == *stagedNode
                       ? stagedSequence
                       : nodes.seqOf(h);
        };
        auto nodeCurrent = [&](NodeHandle h) {
            return (refinementMode && h == *stagedNode) || nodes.isCurrent(h);
        };
        auto edgeRankKey = [&](NodeHandle a, NodeHandle b) {
            auto ra = nodeRank(a), rb = nodeRank(b);
            if (rb < ra) std::swap(ra, rb);
            return TetEdgeRankKey{ra, rb};
        };
        auto faceRankKey = [&](const std::array<NodeHandle, 3>& directed) {
            TetFaceRankKey key{{nodeRank(directed[0]), nodeRank(directed[1]),
                                nodeRank(directed[2])}};
            std::sort(key.ranks.begin(), key.ranks.end());
            return key;
        };
        auto stagedFace = [&](const std::array<NodeHandle, 3>& directed) {
            FacePool::Staged staged;
            staged.shape = FaceShape::Tri;
            staged.n = 3;
            for (std::size_t i = 0; i < 3; ++i)
                staged.outward[i] = directed[i];
            staged.key.n = 3;
            const auto key = faceRankKey(directed);
            for (std::size_t i = 0; i < 3; ++i)
                staged.key.s[i] = key.ranks[i];
            return staged;
        };
        auto stagedFaceSense = [&](const FacePool::Staged& stored,
                                   const FacePool::Staged& use) {
            std::array<std::uint64_t, 3> a{}, b{};
            for (std::size_t i = 0; i < 3; ++i) {
                a[i] = nodeRank(stored.outward[i]);
                b[i] = nodeRank(use.outward[i]);
            }
            for (std::size_t shift = 0; shift < 3; ++shift)
                if (b[0] == a[shift] && b[1] == a[(shift + 1) % 3] &&
                    b[2] == a[(shift + 2) % 3])
                    return std::int8_t{+1};
            for (std::size_t shift = 0; shift < 3; ++shift)
                if (b[0] == a[shift] && b[1] == a[(shift + 2) % 3] &&
                    b[2] == a[(shift + 1) % 3])
                    return std::int8_t{-1};
            return std::int8_t{0};
        };
        TetPlan plan;
        plan.targetClassification = targetClassification;
        plan.edgeTargetMode = edgeTargetMode;
        std::set<CellHandle> removedSet;
        std::set<std::uint64_t> removedNodeRanks;
        std::map<std::uint64_t,NodeHandle> removedNodeByRank;
        for (const auto cell : removed) {
            if (!cells.isCurrent(cell) || !removedSet.insert(cell).second)
                throw std::invalid_argument("Tet cavity removed identity");
            const auto snap = cells.get(cell);
            if (snap.type != ElementType::Tet4 || snap.cls.dim != region.dim ||
                snap.cls.tag != region.tag ||
                !edges.hasCellEdges(cells, cell) ||
                !faces.hasCellFaces(cells, cell))
                throw std::invalid_argument("Tet cavity removed cell contract");
            std::array<std::uint64_t, 4> key{};
            for (std::size_t i = 0; i < 4; ++i)
                key[i] = nodeRank(snap.nodes[i]);
            std::sort(key.begin(), key.end());
            for (std::size_t i = 0; i < 4; ++i)
                removedNodeByRank.emplace(nodeRank(snap.nodes[i]),
                                          snap.nodes[i]);
            for (const auto& [rank,node] : removedNodeByRank) {
                (void)node;
                removedNodeRanks.insert(rank);
            }
            plan.removes.push_back(TetRemoveRow{cell});
        }
        std::set<std::array<std::uint64_t, 4>> survivingKeys;
        std::set<std::uint64_t> finalNodeRanks;
        if (!refinementMode || edgeTargetMode) {
            for (const auto cell : cells.activeBySequence()) {
                if (removedSet.count(cell)) continue;
                const auto snap = cells.get(cell);
                if (snap.type != ElementType::Tet4) continue;
                std::array<std::uint64_t, 4> key{};
                for (std::size_t i = 0; i < 4; ++i)
                    key[i] = nodeRank(snap.nodes[i]);
                std::sort(key.begin(), key.end());
                survivingKeys.insert(key);
                for (std::size_t i = 0; i < 4; ++i)
                    finalNodeRanks.insert(nodeRank(snap.nodes[i]));
            }
        }

        struct EdgeAgg {
            NodeHandle a{}, b{};
            EdgeHandle existing{};
            std::uint32_t initial = 0, removed = 0, added = 0;
        };
        struct FaceAgg {
            FacePool::Staged staged{};
            FaceHandle existing{};
            std::vector<FacePool::Incidence> surviving;
            std::uint32_t initial = 0, removed = 0, added = 0;
            std::vector<std::pair<std::size_t, std::int8_t>> addedUses;
        };
        std::map<TetEdgeRankKey, EdgeAgg> edgeAgg;
        std::map<TetFaceRankKey, FaceAgg> faceAgg;

        for (const auto cell : removed) {
            for (const auto edge : edges.cellEdges(cells, cell)) {
                const auto data = edges.get(edge);
                const auto key = edgeRankKey(data.a, data.b);
                auto& agg = edgeAgg[key];
                agg.a = data.a;
                agg.b = data.b;
                if (nodeRank(agg.b) < nodeRank(agg.a))
                    std::swap(agg.a, agg.b);
                agg.existing = edge;
                agg.initial = data.useCount;
                ++agg.removed;
            }
            for (const auto face : faces.cellFaces(cells, cell)) {
                const auto data = faces.get(face);
                std::array<NodeHandle, 3> directed{
                    data.corners[0], data.corners[1], data.corners[2]};
                const auto key = faceRankKey(directed);
                auto& agg = faceAgg[key];
                agg.existing = face;
                agg.initial = static_cast<std::uint32_t>(data.incident.size());
                agg.staged = stagedFace(directed);
                ++agg.removed;
            }
        }

        static constexpr std::array<std::array<int, 2>, 6> localEdges{{
            {{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};
        static constexpr std::array<std::array<int, 3>, 4> localFaces{{
            {{1, 2, 3}}, {{0, 3, 2}}, {{0, 1, 3}}, {{0, 2, 1}}}};
        std::vector<std::array<TetEdgeRankKey, 6>> addEdgeKeys;
        std::vector<std::array<TetFaceRankKey, 4>> addFaceKeys;
        addEdgeKeys.resize(added.size());
        addFaceKeys.resize(added.size());
        std::set<std::array<std::uint64_t, 4>> addedKeys;
        for (std::size_t ai = 0; ai < added.size(); ++ai) {
            auto conn = added[ai];
            std::array<std::uint64_t, 4> cellKey{};
            for (std::size_t i = 0; i < 4; ++i) {
                if (!nodeCurrent(conn[i]))
                    throw std::invalid_argument("Tet cavity added node identity");
                cellKey[i] = nodeRank(conn[i]);
            }
            auto sorted = cellKey;
            std::sort(sorted.begin(), sorted.end());
            if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end() ||
                !addedKeys.insert(sorted).second || survivingKeys.count(sorted))
                throw std::invalid_argument("Tet cavity added duplicate");
            if (addedClasses[ai].dim != region.dim ||
                addedClasses[ai].tag != region.tag)
                throw std::invalid_argument("Tet cavity added region mismatch");
            for (const auto node : conn)
                finalNodeRanks.insert(nodeRank(node));
            TetAddRow add;
            add.nodes = conn;
            add.cls = addedClasses[ai];
            plan.adds.push_back(add);
            for (std::size_t e = 0; e < localEdges.size(); ++e) {
                NodeHandle a = conn[static_cast<std::size_t>(localEdges[e][0])];
                NodeHandle b = conn[static_cast<std::size_t>(localEdges[e][1])];
                const auto key = edgeRankKey(a, b);
                addEdgeKeys[ai][e] = key;
                auto& agg = edgeAgg[key];
                if (!agg.existing.isSpecified()) {
                    const bool hasStaged = refinementMode &&
                        (a == *stagedNode || b == *stagedNode);
                    const auto found = hasStaged
                        ? std::optional<EdgeHandle>{}
                        : edges.find(nodes, a, b);
                    if (found) {
                        agg.existing = *found;
                        const auto data = edges.get(*found);
                        agg.initial = data.useCount;
                        agg.a = data.a;
                        agg.b = data.b;
                    } else {
                        if (nodeRank(b) < nodeRank(a)) std::swap(a, b);
                        agg.a = a;
                        agg.b = b;
                    }
                }
                ++agg.added;
            }
            for (std::size_t f = 0; f < localFaces.size(); ++f) {
                std::array<NodeHandle, 3> directed{
                    conn[static_cast<std::size_t>(localFaces[f][0])],
                    conn[static_cast<std::size_t>(localFaces[f][1])],
                    conn[static_cast<std::size_t>(localFaces[f][2])]};
                const auto key = faceRankKey(directed);
                addFaceKeys[ai][f] = key;
                auto staged = stagedFace(directed);
                auto& agg = faceAgg[key];
                if (!agg.existing.isSpecified()) {
                    const bool hasStaged = refinementMode &&
                        (directed[0] == *stagedNode ||
                         directed[1] == *stagedNode ||
                         directed[2] == *stagedNode);
                    const auto found = hasStaged
                        ? std::optional<FaceHandle>{}
                        : faces.find(nodes, directed);
                    if (found) {
                        agg.existing = *found;
                        const auto data = faces.get(*found);
                        agg.initial =
                            static_cast<std::uint32_t>(data.incident.size());
                        std::array<NodeHandle, 3> stored{
                            data.corners[0], data.corners[1], data.corners[2]};
                        agg.staged = stagedFace(stored);
                    } else if (agg.added == 0) {
                        agg.staged = staged;
                    }
                }
                std::int8_t sense = agg.existing.isSpecified()
                    ? faceSense(faces, nodes, agg.existing, staged)
                    : stagedFaceSense(agg.staged, staged);
                if (sense == 0)
                    throw std::invalid_argument("Tet cavity face winding");
                agg.addedUses.push_back({ai, sense});
                ++agg.added;
            }
        }

        const auto targetEdgeKey = edgeTargetMode
            ? edgeRankKey(target[0], target[1])
            : TetEdgeRankKey{};
        std::map<TetEdgeRankKey, std::size_t> edgeIndex;
        for (const auto& [key, agg] : edgeAgg) {
            TetEdgeRow row;
            row.key = EdgeKey{key.lo, key.hi};
            row.a = agg.a;
            row.b = agg.b;
            row.existing = agg.existing;
            row.initialUse = agg.initial;
            row.removedUse = agg.removed;
            row.addedUse = agg.added;
            if (row.removedUse > row.initialUse)
                throw std::invalid_argument("Tet cavity edge delta underflow");
            row.finalUse = row.initialUse - row.removedUse + row.addedUse;
            row.create = !row.existing.isSpecified() && row.finalUse > 0;
            row.retire = row.existing.isSpecified() && row.finalUse == 0;
            row.target = edgeTargetMode &&
                         !(key < targetEdgeKey) &&
                         !(targetEdgeKey < key);
            if (row.target && row.existing.isSpecified()) {
                const auto data = edges.get(row.existing);
                const bool exact =
                    data.cls.dim == targetClassification.dim &&
                    data.cls.tag == targetClassification.tag;
                const bool temporaryDefault =
                    !data.requiredConstraint && data.cls.dim == -1 &&
                    data.cls.tag == -1;
                if ((!exact && !temporaryDefault) ||
                    (data.requiredConstraint && !exact))
                    throw std::invalid_argument(
                        "Tet cavity target edge metadata conflict");
            }
            if (row.retire && edges.get(row.existing).requiredConstraint)
                throw std::invalid_argument("Tet cavity protected edge conflict");
            edgeIndex.emplace(key, plan.edges.size());
            if (row.target) plan.targetEdgeRow = plan.edges.size();
            plan.edges.push_back(row);
        }

        const auto targetKey = (edgeTargetMode || refinementMode)
                                   ? TetFaceRankKey{}
                                   : faceRankKey(target);
        std::map<TetFaceRankKey, std::size_t> faceIndex;
        for (auto& [key, agg] : faceAgg) {
            TetFaceRow row;
            row.key.n = 3;
            for (std::size_t i = 0; i < 3; ++i) row.key.s[i] = key.ranks[i];
            row.staged = agg.staged;
            row.existing = agg.existing;
            row.initialUse = agg.initial;
            row.removedUse = agg.removed;
            row.addedUse = agg.added;
            if (row.removedUse > row.initialUse)
                throw std::invalid_argument("Tet cavity face delta underflow");
            row.finalUse = row.initialUse - row.removedUse + row.addedUse;
            row.create = !row.existing.isSpecified() && row.finalUse > 0;
            row.retire = row.existing.isSpecified() && row.finalUse == 0;
            row.target = !refinementMode && !edgeTargetMode &&
                         !(key < targetKey) &&
                         !(targetKey < key);
            if (row.existing.isSpecified()) {
                const auto data = faces.get(row.existing);
                for (const auto& incidence : data.incident) {
                    const CellHandle owner{incidence.cellSlot,
                                           incidence.cellGen,
                                           cells.poolId()};
                    if (!removedSet.count(owner)) row.surviving.push_back(incidence);
                }
                const bool protectedFace = data.cls.dim == 2 && data.cls.tag >= 0;
                if (protectedFace && row.removedUse > 0 && row.finalUse == 0)
                    throw std::invalid_argument("Tet cavity protected face conflict");
                if (row.initialUse == 1 && row.removedUse == 1 &&
                    row.finalUse != 1)
                    throw std::invalid_argument("Tet cavity boundary face mismatch");
            }
            std::vector<std::int8_t> finalSenses;
            for (const auto& incidence : row.surviving)
                finalSenses.push_back(incidence.sense);
            for (const auto& use : agg.addedUses)
                finalSenses.push_back(use.second);
            if (finalSenses.size() != row.finalUse || row.finalUse > 2 ||
                (row.finalUse == 2 &&
                 finalSenses[0] * finalSenses[1] != -1))
                throw std::invalid_argument("Tet cavity final face incidence");
            faceIndex.emplace(key, plan.faces.size());
            if (row.target) plan.targetFaceRow = plan.faces.size();
            plan.faces.push_back(std::move(row));
        }
        if (!refinementMode && !edgeTargetMode &&
            (plan.targetFaceRow == std::numeric_limits<std::size_t>::max() ||
             plan.faces[plan.targetFaceRow].finalUse != 2))
            throw std::invalid_argument("Tet cavity target face incidence");
        if (edgeTargetMode &&
            (plan.targetEdgeRow == std::numeric_limits<std::size_t>::max() ||
             plan.edges[plan.targetEdgeRow].finalUse == 0))
            throw std::invalid_argument("Tet cavity target edge incidence");
        if (refinementMode) {
            // Ordinary refinement uses the staged node in every cell. Edge
            // recovery can retain existing-node cells, checked against the
            // surviving keys above. Preserve a removed node used elsewhere.
            for (const auto rank : removedNodeRanks) {
                if (finalNodeRanks.count(rank)) continue;
                const NodeHandle node=removedNodeByRank.at(rank);
                std::set<CellHandle> pending,visited;
                for (const auto cell : removed) {
                    const auto snap=cells.get(cell);
                    if (std::find(snap.nodes.begin(),
                                  snap.nodes.begin()+snap.count,node)!=
                            snap.nodes.begin()+snap.count)
                        pending.insert(cell);
                }
                bool survives=false;
                while (!pending.empty()&&!survives) {
                    const auto first=pending.begin();
                    const CellHandle cell=*first;
                    pending.erase(first);
                    if (!visited.insert(cell).second) continue;
                    for (const auto face : faces.cellFaces(cells,cell)) {
                        const auto data=faces.get(face);
                        bool contains=false;
                        for (std::size_t i=0;i<data.cornerCount;++i)
                            contains=contains||data.corners[i]==node;
                        if (!contains) continue;
                        for (const auto& incidence:data.incident) {
                            const CellHandle owner{
                                incidence.cellSlot,incidence.cellGen,
                                cells.poolId()};
                            if (!removedSet.count(owner)) {
                                survives=true;
                                break;
                            }
                            if (!visited.count(owner)) pending.insert(owner);
                        }
                        if (survives) break;
                    }
                }
                if (survives) finalNodeRanks.insert(rank);
            }
        }
        for (const auto rank : removedNodeRanks)
            if (!finalNodeRanks.count(rank))
                throw std::invalid_argument("Tet cavity would orphan a node");

        for (std::size_t ai = 0; ai < plan.adds.size(); ++ai) {
            for (std::size_t e = 0; e < 6; ++e)
                plan.adds[ai].edges[e] = edgeIndex.at(addEdgeKeys[ai][e]);
            for (std::size_t f = 0; f < 4; ++f) {
                const auto rowIndex = faceIndex.at(addFaceKeys[ai][f]);
                plan.adds[ai].faces[f] = rowIndex;
                const auto& uses = faceAgg.at(addFaceKeys[ai][f]).addedUses;
                auto found = std::find_if(uses.begin(), uses.end(),
                    [&](const auto& use) { return use.first == ai; });
                plan.adds[ai].faceSenses[f] = found->second;
            }
        }
        if (remainingSequences(cells.alloc_) < plan.adds.size())
            throw std::length_error("Tet cavity cell sequence exhausted");
        std::size_t createdEdgeCount = 0;
        for (const auto& row : plan.edges)
            createdEdgeCount += row.create ? 1U : 0U;
        if (remainingSequences(edges.alloc_) < createdEdgeCount)
            throw std::length_error("Tet cavity edge sequence exhausted");
        std::size_t createdFaceCount = 0;
        for (const auto& row : plan.faces)
            createdFaceCount += row.create ? 1U : 0U;
        if (remainingSequences(faces.alloc_) < createdFaceCount)
            throw std::length_error("Tet cavity face sequence exhausted");
        return plan;
    }

    static TetPlan planTetCavity(
        const NodePool& nodes, CellPool& cells, EdgePool& edges,
        FacePool& faces, std::span<const CellHandle> removed,
        std::span<const std::array<NodeHandle, 4>> added,
        std::span<const mesh::Classification> addedClasses,
        const std::array<NodeHandle, 3>& target,
        const mesh::Classification& targetClassification,
        const mesh::Classification& region) {
        return planTetCavityImpl(nodes, cells, edges, faces, removed, added,
                                 addedClasses, target, targetClassification,
                                 region, false);
    }

    static TetPlan planTetEdgeCavity(
        const NodePool& nodes, CellPool& cells, EdgePool& edges,
        FacePool& faces, std::span<const CellHandle> removed,
        std::span<const std::array<NodeHandle, 4>> added,
        std::span<const mesh::Classification> addedClasses,
        const std::array<NodeHandle, 2>& target,
        const mesh::Classification& targetClassification,
        const mesh::Classification& region) {
        const std::array<NodeHandle, 3> expanded{target[0], target[1],
                                                 target[0]};
        return planTetCavityImpl(nodes, cells, edges, faces, removed, added,
                                 addedClasses, expanded,
                                 targetClassification, region, true);
    }

    // T3-A: the one predicted staged node is treated as current only for
    // canonical rank/key planning. No pool record is materialized until the
    // composite noexcept commit, and no target edge/face metadata is changed.
    static TetPlan planTetRefinement(
        const NodePool& nodes, CellPool& cells, EdgePool& edges,
        FacePool& faces, std::span<const CellHandle> removed,
        std::span<const std::array<NodeHandle, 4>> added,
        std::span<const mesh::Classification> addedClasses,
        NodeHandle stagedNode, std::uint64_t stagedSequence,
        const mesh::Classification& region) {
        const std::array<NodeHandle, 3> noTarget{};
        return planTetCavityImpl(nodes, cells, edges, faces, removed, added,
                                 addedClasses, noTarget,
                                 mesh::Classification{}, region, false,
                                 &stagedNode, stagedSequence);
    }

    // V0.10 bounded non-star recovery: the same one-node composite commit as
    // refinement, with one prescribed edge marked/classified in that exact
    // prepared net delta.  No node or topology is visible before commit.
    static TetPlan planTetEdgeRefinement(
        const NodePool& nodes, CellPool& cells, EdgePool& edges,
        FacePool& faces, std::span<const CellHandle> removed,
        std::span<const std::array<NodeHandle, 4>> added,
        std::span<const mesh::Classification> addedClasses,
        NodeHandle stagedNode, std::uint64_t stagedSequence,
        const std::array<NodeHandle, 2>& target,
        const mesh::Classification& targetClassification,
        const mesh::Classification& region) {
        const std::array<NodeHandle, 3> expanded{
            target[0], target[1], target[0]};
        return planTetCavityImpl(
            nodes, cells, edges, faces, removed, added, addedClasses,
            expanded, targetClassification, region, true, &stagedNode,
            stagedSequence);
    }

    static void prepareTetCavity(CellPool& cells, EdgePool& edges,
                                 FacePool& faces, TetPlan& plan) {
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        const auto cellSlots = predictSlots(cells.alloc_, plan.adds.size());
        cells.alloc_.reserveMaterializedCapacity(
            account(cells.alloc_, plan.adds.size(), plan.removes.size())
                .slotCapacity);
        cells.alloc_.reserveFreeCapacity(plan.removes.size());
        std::size_t need = cells.rec_.size();
        for (const auto slot : cellSlots)
            need = std::max(need, static_cast<std::size_t>(slot) + 1);
        if (need > cells.rec_.size()) {
            detail::reserveFor(cells.rec_, need);
            cells.rec_.resize(need);
        }
        for (std::size_t i = 0; i < plan.adds.size(); ++i) {
            plan.adds[i].predictedSlot = cellSlots[i];
            plan.adds[i].predictedGen = slotGen(cells.alloc_, cellSlots[i]);
        }

        std::size_t newEdges = 0, retiredEdges = 0;
        for (const auto& row : plan.edges) {
            newEdges += row.create ? 1 : 0;
            retiredEdges += row.retire ? 1 : 0;
        }
        const auto edgeSlots = predictSlots(edges.alloc_, newEdges);
        edges.alloc_.reserveMaterializedCapacity(
            account(edges.alloc_, newEdges, retiredEdges).slotCapacity);
        edges.alloc_.reserveFreeCapacity(retiredEdges);
        need = edges.rec_.size();
        for (const auto slot : edgeSlots)
            need = std::max(need, static_cast<std::size_t>(slot) + 1);
        if (need > edges.rec_.size()) {
            detail::reserveFor(edges.rec_, need);
            edges.rec_.resize(need);
        }
        edges.table_.reserve(edges.table_.size() + newEdges);
        std::size_t ei = 0;
        for (auto& row : plan.edges)
            if (row.create) {
                row.predictedSlot = edgeSlots[ei++];
                row.predictedGen = slotGen(edges.alloc_, row.predictedSlot);
                auto inserted = edges.table_.emplace(row.key, EdgeHandle{});
                if (!inserted.second)
                    throw std::logic_error("Tet cavity edge placeholder conflict");
                row.tabIt = inserted.first;
                row.tabPlaced = true;
#if defined(PAACS_TOPO_TESTING)
                maybeFault();
#endif
            }

        edges.incident_.reserve(edges.incident_.size() + plan.adds.size());
        for (auto& add : plan.adds) {
            auto inserted = edges.incident_.emplace(
                add.predictedSlot, EdgePool::IncidenceEntry{});
            if (!inserted.second)
                throw std::logic_error("Tet cavity edge incidence conflict");
            // Record ownership of the placeholder before its embedded vector
            // reserve: that reserve is the next throwing operation.
            add.edgeIncIt = inserted.first;
            add.edgeIncPlaced = true;
            inserted.first->second.gen = add.predictedGen;
            inserted.first->second.edges.reserve(6);
#if defined(PAACS_TOPO_TESTING)
            maybeFault();
#endif
        }
        for (auto& remove : plan.removes) {
            remove.edgeIncIt = edges.incident_.find(remove.cell.slot);
            remove.hasEdgeInc = remove.edgeIncIt != edges.incident_.end();
            if (!remove.hasEdgeInc)
                throw std::logic_error("Tet cavity missing removed edge incidence");
        }

        std::size_t newFaces = 0, retiredFaces = 0;
        for (const auto& row : plan.faces) {
            newFaces += row.create ? 1 : 0;
            retiredFaces += row.retire ? 1 : 0;
        }
        const auto faceSlots = predictSlots(faces.alloc_, newFaces);
        faces.alloc_.reserveMaterializedCapacity(
            account(faces.alloc_, newFaces, retiredFaces).slotCapacity);
        faces.alloc_.reserveFreeCapacity(retiredFaces);
        need = faces.rec_.size();
        for (const auto slot : faceSlots)
            need = std::max(need, static_cast<std::size_t>(slot) + 1);
        if (need > faces.rec_.size()) {
            detail::reserveFor(faces.rec_, need);
            faces.rec_.resize(need);
        }
        faces.table_.reserve(faces.table_.size() + newFaces);
        std::size_t fi = 0;
        for (auto& row : plan.faces) {
            if (row.create) {
                row.predictedSlot = faceSlots[fi++];
                row.predictedGen = slotGen(faces.alloc_, row.predictedSlot);
                faces.rec_[static_cast<std::size_t>(row.predictedSlot)]
                    .data.incident.reserve(row.finalUse);
                auto inserted = faces.table_.emplace(row.key, FaceHandle{});
                if (!inserted.second)
                    throw std::logic_error("Tet cavity face placeholder conflict");
                row.tabIt = inserted.first;
                row.tabPlaced = true;
            } else if (row.existing.isSpecified()) {
                faces.rec_[static_cast<std::size_t>(row.existing.slot)]
                    .data.incident.reserve(row.finalUse);
            }
#if defined(PAACS_TOPO_TESTING)
            maybeFault();
#endif
        }
        faces.cellFaces_.reserve(faces.cellFaces_.size() + plan.adds.size());
        for (auto& add : plan.adds) {
            auto inserted = faces.cellFaces_.emplace(
                add.predictedSlot, FacePool::CellFaceEntry{});
            if (!inserted.second)
                throw std::logic_error("Tet cavity face incidence conflict");
            // Same immediate-tracking rule as the edge-incidence table.
            add.faceCfIt = inserted.first;
            add.faceCfPlaced = true;
            inserted.first->second.gen = add.predictedGen;
            inserted.first->second.faces.reserve(4);
#if defined(PAACS_TOPO_TESTING)
            maybeFault();
#endif
        }
        for (auto& remove : plan.removes) {
            remove.faceCfIt = faces.cellFaces_.find(remove.cell.slot);
            remove.hasFaceCf = remove.faceCfIt != faces.cellFaces_.end();
            if (!remove.hasFaceCf)
                throw std::logic_error("Tet cavity missing removed face incidence");
        }
        // Retiring table iterators are captured only after the final reserve.
        for (auto& row : plan.edges)
            if (row.retire) row.tabIt = edges.table_.find(row.key);
        for (auto& row : plan.faces)
            if (row.retire) row.tabIt = faces.table_.find(row.key);
        plan.prepared = true;
    }

    static void abandonTetCavity(EdgePool& edges, FacePool& faces,
                                 TetPlan& plan) noexcept {
        for (auto& row : plan.edges)
            if (row.create && row.tabPlaced) {
                edges.table_.erase(row.tabIt);
                row.tabPlaced = false;
            }
        for (auto& row : plan.faces)
            if (row.create && row.tabPlaced) {
                faces.table_.erase(row.tabIt);
                row.tabPlaced = false;
            }
        for (auto& add : plan.adds) {
            if (add.edgeIncPlaced) {
                edges.incident_.erase(add.edgeIncIt);
                add.edgeIncPlaced = false;
            }
            if (add.faceCfPlaced) {
                faces.cellFaces_.erase(add.faceCfIt);
                add.faceCfPlaced = false;
            }
        }
        plan.prepared = false;
    }

    static FaceHandle predictedTetTargetFace(const FacePool& faces,
                                             const TetPlan& plan) noexcept {
        const auto& row = plan.faces[plan.targetFaceRow];
        return row.create
            ? FaceHandle{row.predictedSlot, row.predictedGen,
                         faces.alloc_.poolId()}
            : row.existing;
    }

    static EdgeHandle predictedTetTargetEdge(const EdgePool& edges,
                                             const TetPlan& plan) noexcept {
        const auto& row = plan.edges[plan.targetEdgeRow];
        return row.create
            ? EdgeHandle{row.predictedSlot, row.predictedGen,
                         edges.alloc_.poolId()}
            : row.existing;
    }

    static void commitTetCavity(NodePool& nodes, CellPool& cells,
                                EdgePool& edges, FacePool& faces,
                                TetPlan& plan,
                                const TetCommitBuffers& buffers) noexcept {
        // Materialize all additions before retiring any source slot: retired
        // slots are deliberately not reused within this commit.
        for (auto& add : plan.adds) {
            const auto cell = allocPrepared(cells.alloc_);
#if defined(PAACS_TOPO_TESTING)
            if (cell.slot != add.predictedSlot || cell.gen != add.predictedGen)
                std::abort();
#endif
            add.materialized = cell;
            auto& record = cells.rec_[static_cast<std::size_t>(cell.slot)];
            record.type = ElementType::Tet4;
            record.count = 4;
            for (std::size_t i = 0; i < 4; ++i)
                record.conn[i] = CellPool::LocalNodeRef{
                    add.nodes[i].slot, add.nodes[i].gen};
            record.cls = add.cls;
            record.seq = cells.alloc_.slots_[static_cast<std::size_t>(cell.slot)].seq;
            record.gen = cell.gen;
            cells.linkTail(cell.slot);
            buffers.resultCells->push_back(cell);
        }
        cells.bumpPeaks();

        for (auto& remove : plan.removes) {
            edges.incident_.erase(remove.edgeIncIt);
            faces.cellFaces_.erase(remove.faceCfIt);
        }

        for (auto& row : plan.edges) {
            if (row.create) {
                const auto edge = allocPrepared(edges.alloc_);
#if defined(PAACS_TOPO_TESTING)
                if (edge.slot != row.predictedSlot || edge.gen != row.predictedGen)
                    std::abort();
#endif
                row.finalHandle = edge;
                auto& record = edges.rec_[static_cast<std::size_t>(edge.slot)];
                record.data.a = row.a;
                record.data.b = row.b;
                record.data.cls = row.target ? plan.targetClassification
                                             : mesh::Classification{};
                record.data.requiredConstraint = row.target;
                record.data.useCount = row.finalUse;
                record.key = row.key;
                record.seq = edges.alloc_.slots_[static_cast<std::size_t>(edge.slot)].seq;
                record.gen = edge.gen;
                edges.linkTail(edge.slot);
                row.tabIt->second = edge;
                row.tabPlaced = false;
            } else if (!row.retire) {
                row.finalHandle = row.existing;
                auto& data = edges.rec_[static_cast<std::size_t>(
                    row.existing.slot)].data;
                data.useCount = row.finalUse;
                if (row.target) {
                    data.cls = plan.targetClassification;
                    data.requiredConstraint = true;
                }
            }
        }
        for (auto& add : plan.adds) {
            auto& incidence = add.edgeIncIt->second.edges;
            for (const auto rowIndex : add.edges)
                incidence.push_back(plan.edges[rowIndex].finalHandle);
        }
        edges.bumpPeaks();

        for (auto& row : plan.faces) {
            FaceHandle face;
            if (row.create) {
                face = allocPrepared(faces.alloc_);
#if defined(PAACS_TOPO_TESTING)
                if (face.slot != row.predictedSlot || face.gen != row.predictedGen)
                    std::abort();
#endif
                auto& record = faces.rec_[static_cast<std::size_t>(face.slot)];
                record.data.incident.clear();
                record.data.shape = FaceShape::Tri;
                record.data.cornerCount = 3;
                for (std::size_t i = 0; i < 3; ++i)
                    record.data.corners[i] = row.staged.outward[i];
                record.data.cls = row.target ? plan.targetClassification
                                             : mesh::Classification{};
                record.key = row.key;
                record.seq = faces.alloc_.slots_[static_cast<std::size_t>(face.slot)].seq;
                record.gen = face.gen;
                faces.linkTail(face.slot);
                row.tabIt->second = face;
                row.tabPlaced = false;
            } else if (!row.retire) {
                face = row.existing;
                auto& record = faces.rec_[static_cast<std::size_t>(face.slot)];
                record.data.incident.clear();
                for (const auto& incidence : row.surviving)
                    record.data.incident.push_back(incidence);
                if (row.target) record.data.cls = plan.targetClassification;
            }
            row.finalHandle = face;
        }
        for (auto& add : plan.adds) {
            auto& cellFaces = add.faceCfIt->second.faces;
            for (std::size_t f = 0; f < 4; ++f) {
                auto& row = plan.faces[add.faces[f]];
                auto& record = faces.rec_[static_cast<std::size_t>(
                    row.finalHandle.slot)];
                record.data.incident.push_back(FacePool::Incidence{
                    add.materialized.slot, add.materialized.gen,
                    add.faceSenses[f]});
                cellFaces.push_back(row.finalHandle);
            }
        }
        faces.bumpPeaks();

        for (auto& row : plan.faces)
            if (row.retire) {
                faces.table_.erase(row.tabIt);
                faces.unlink(row.existing.slot);
                freePrepared(faces.alloc_, row.existing);
            }
        for (auto& row : plan.edges)
            if (row.retire) {
                edges.table_.erase(row.tabIt);
                edges.unlink(row.existing.slot);
                freePrepared(edges.alloc_, row.existing);
            }
        for (auto& remove : plan.removes) {
            cells.unlink(remove.cell.slot);
            freePrepared(cells.alloc_, remove.cell);
        }
        if (buffers.targetFace)
            *buffers.targetFace = plan.faces[plan.targetFaceRow].finalHandle;
        if (buffers.targetEdge)
            *buffers.targetEdge = plan.edges[plan.targetEdgeRow].finalHandle;
        plan.prepared = false;
        (void)nodes;
    }

    // T3-A composite preparation/commit. Node capacity is prepared first;
    // Tet placeholders are then prepared by the accepted net-delta machinery.
    // The source node allocator is not consumed until the composite commit.
    static void prepareTetRefinement(
        NodePool& nodes, CellPool& cells, EdgePool& edges, FacePool& faces,
        const std::vector<StagedNodePredict>& nodePrediction,
        TetPlan& tetPlan) {
        reserveNodeBatch(nodes, nodePrediction);
        prepareTetCavity(cells, edges, faces, tetPlan);
    }

    static void abandonTetRefinementPrepared(
        EdgePool& edges, FacePool& faces, TetPlan& tetPlan) noexcept {
        abandonTetCavity(edges, faces, tetPlan);
    }

    static void commitTetRefinementPrepared(
        NodePool& nodes, CellPool& cells, EdgePool& edges, FacePool& faces,
        const std::vector<StagedNodeSpec>& nodeSpecs,
        const std::vector<StagedNodePredict>& nodePrediction,
        std::vector<NodeHandle>& createdNodes, TetPlan& tetPlan,
        const TetCommitBuffers& buffers) noexcept {
        commitNodeBatch(nodes, nodeSpecs, nodePrediction, createdNodes);
        commitTetCavity(nodes, cells, edges, faces, tetPlan, buffers);
    }
};

} // namespace paacs::topo::detail
