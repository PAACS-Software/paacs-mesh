// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// paacs::topo - M4-B.3: explicit reusable edge topology.
//
// Edges are first-class entities derived from cell connectivity, for
// constraint recovery, refinement encroachment, and hex-layer logic.
// Composed ON the frozen M4-B.1 HandleAllocator<EdgeTag> and the frozen
// M4-B.2 NodePool / CellPool. Consumes ONLY storage that passes
// storageEdgeReady(nodes, cells).
//
// IDENTITY (directive 2.1: topological identity is node identity, never
// geometry): an edge's canonical key is the unordered pair of its two
// endpoint nodes, ordered (lo, hi) by NODE CREATION SEQUENCE - a stable,
// reuse-independent identity that survives far translation and does not
// depend on coordinates. Runtime pool ids / slots / generations never
// enter the key, serialization, or hashes.
//
// INVARIANTS (constitution M4-B.3):
//   * key uses node identities only;
//   * an edge exists iff it has >= 1 incident-cell use;
//   * use counts stay consistent with incident cells;
//   * no duplicate edge record for one identity pair;
//   * self-edge (a, a) is rejected.
//
// SCOPE (reviewer): canonical undirected edge identity; EdgeHandle
// ownership; unique reusable edge records; sequence-ordered iteration;
// cell<->edge incidence; traits-driven mixed-element edge extraction;
// structured rejection of stale/foreign/duplicate endpoints; three-state
// accounting + bounded churn; pool tokens excluded from persistence.
// DEFERRED: explicit faces, transactions, Delaunay migration, public Mesh
// replacement, curved-edge caches, parallel edge builds.
//
#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "paacs/core/Types.hpp"
#include "paacs/diag/Diagnostics.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/topo/ElementTraits.hpp"
#include "paacs/topo/Handles.hpp"
#include "paacs/topo/Storage.hpp"

namespace paacs::topo {

// Canonical undirected edge key: endpoint node sequences (lo <= hi). Node
// SEQUENCE, not slot/gen/pool - a stable identity independent of storage
// reuse. Self-edges are rejected before a key is ever built, so lo < hi
// always holds for a valid key.
struct EdgeKey {
    std::uint64_t lo = 0;
    std::uint64_t hi = 0;
    bool operator==(const EdgeKey& o) const {
        return lo == o.lo && hi == o.hi;
    }
};

struct EdgeKeyHash {
    std::size_t operator()(const EdgeKey& k) const noexcept {
        // 64-bit mix of the two sequences (FNV-1a style over 16 bytes).
        std::uint64_t h = 1469598103934665603ULL;
        for (std::uint64_t v : {k.lo, k.hi})
            for (int b = 0; b < 8; ++b) {
                h ^= (v >> (8 * b)) & 0xFF;
                h *= 1099511628211ULL;
            }
        return static_cast<std::size_t>(h);
    }
};

struct EdgePoolCompleteBuildMetrics {
    std::uint64_t validationCellVisits = 0;
    std::uint64_t buildCellVisits = 0;
    std::uint64_t edgeUses = 0;
    std::uint64_t uniqueEdges = 0;
    std::uint64_t finalEdgeAudits = 0;
    std::uint64_t perCellGlobalAudits = 0;
    bool runComplete = false;
};

static_assert(std::is_standard_layout_v<EdgePoolCompleteBuildMetrics>);
static_assert(std::is_trivially_copyable_v<EdgePoolCompleteBuildMetrics>);

class EdgePool;
inline diag::Report auditEdgeTopology(
    const NodePool&,const CellPool&,const EdgePool&);

// ---------------------------------------------------------------------------
// EdgePool - reusable edge records with canonical identity and cell<->edge
// incidence. Bound to one NodePool and one CellPool.
// ---------------------------------------------------------------------------
class EdgePool {
    friend class detail::TransactionAccess;   // M4-B.6 prepared commit
public:
    struct EdgeData {
        NodeHandle a;                     // endpoints, ordered lo->hi by seq
        NodeHandle b;
        mesh::Classification cls{};
        bool requiredConstraint = false;
        std::uint32_t useCount = 0;       // # incident cell-edges
    };

    EdgePool() = default;
    EdgePool(const EdgePool&) = delete;
    EdgePool& operator=(const EdgePool&) = delete;
    EdgePool(EdgePool&& o)
        : EdgePool(TokenlessAdoptionTag{}) {
#if defined(PAACS_TOPO_TESTING)
        if (adoptionBeginHook_ != nullptr) adoptionBeginHook_();
#endif
        adoptFrom(o);
    }
    EdgePool& operator=(EdgePool&&) = delete;

    static EdgePool buildComplete(
        const NodePool& nodes,const CellPool& cells,
        EdgePoolCompleteBuildMetrics* metrics = nullptr) {
        return EdgePool(CompleteBuildTag{},nodes,cells,metrics);
    }

    std::uint64_t poolId() const { return alloc_.poolId(); }
    std::uint64_t boundNodePool() const { return boundNodePool_; }
    std::uint64_t boundCellPool() const { return boundCellPool_; }

    // Cell-edge incidence queries (M4-B.4 gate requirement): does THIS
    // specific cell have edge incidence recorded here (by slot AND
    // generation)? Distinguishes "this cell's edges were inserted" from
    // "some geometrically identical cell created the same edge records".
    // Pure const reads of existing incidence state - additive, no behavior
    // change.
    bool hasCellEdges(const CellPool& cells, const CellHandle& cell) const {
        if (boundCellPool_ != 0 && cells.poolId() != boundCellPool_)
            return false;
        const auto it = incident_.find(cell.slot);
        return it != incident_.end() && it->second.gen == cell.gen;
    }
    std::span<const EdgeHandle> cellEdges(const CellPool& cells,
                                          const CellHandle& cell) const {
        if (boundCellPool_ != 0 && cells.poolId() != boundCellPool_)
            return {};
        const auto it = incident_.find(cell.slot);
        if (it == incident_.end() || it->second.gen != cell.gen) return {};
        return std::span<const EdgeHandle>(it->second.edges.data(),
                                           it->second.edges.size());
    }
    bool isCurrent(const EdgeHandle& h) const { return alloc_.isCurrent(h); }

    // Look up an existing edge by endpoints; nullopt if absent. Endpoints
    // may be given in any order.
    std::optional<EdgeHandle> find(const NodePool& nodes,
                                   const NodeHandle& x,
                                   const NodeHandle& y) const {
        // The edge key deliberately excludes pool ids, so runtime lookup
        // MUST gate on pool binding separately (review blocker 2: a foreign
        // NodePool with matching node sequences otherwise returned an edge
        // owned by the bound pool).
        if (boundNodePool_ != 0 && nodes.poolId() != boundNodePool_)
            return std::nullopt;
        EdgeKey k;
        if (!makeKey(nodes, x, y, k)) return std::nullopt;
        const auto it = table_.find(k);
        if (it == table_.end()) return std::nullopt;
        return it->second;
    }

    EdgeData get(const EdgeHandle& h) const {
        const Rec& r = rec_[static_cast<std::size_t>(alloc_.deref(h))];
        return r.data;
    }
    std::uint64_t seqOf(const EdgeHandle& h) const { return alloc_.seqOf(h); }

    // Insert every edge of `cell` (traits-driven), creating edge records on
    // first use and incrementing use counts / incidence on each. STRONGLY
    // EXCEPTION-SAFE: on any throw the pool is left exactly as before the
    // call (newly created records rolled back, use-count increments undone,
    // no incidence recorded). Rejects a stale cell, a reused cell slot that
    // still carries old (different-generation) incidence, and re-insertion
    // of the same cell.
    std::vector<EdgeHandle> insertCellEdges(const NodePool& nodes,
                                            const CellPool& cells,
                                            const CellHandle& cell) {
        bind(nodes, cells);
        if (!cells.isCurrent(cell))
            throw std::invalid_argument(
                "EdgePool::insertCellEdges: stale cell handle");
        // Slot-indexed incidence with generation validation: reject if this
        // slot still holds ANY incidence (same cell => already inserted;
        // different generation => stale/unremoved incidence from a
        // since-reused slot). The pool is implicit (boundCellPool_).
        {
            const auto it = incident_.find(cell.slot);
            if (it != incident_.end()) {
                if (it->second.gen == cell.gen)
                    throw std::invalid_argument(
                        "EdgePool::insertCellEdges: cell edges already "
                        "inserted");
                throw std::invalid_argument(
                    "EdgePool::insertCellEdges: stale edge incidence for a "
                    "reused cell slot; remove the old cell's edges first");
            }
        }
        const auto snap = cells.get(cell);
        const ElementTraits& tr = checkedTraits(snap.type);
        // Stage + validate ALL endpoint pairs before any mutation.
        struct Pair { NodeHandle x, y; EdgeKey key; };
        std::vector<Pair> pairs;
        pairs.reserve(tr.edgeCount);
        for (int e = 0; e < tr.edgeCount; ++e) {
            const NodeHandle x = snap.nodes[tr.edges[e].a];
            const NodeHandle y = snap.nodes[tr.edges[e].b];
            if (!nodes.isCurrent(x) || !nodes.isCurrent(y))
                throw std::invalid_argument(
                    "EdgePool::insertCellEdges: stale endpoint");
            EdgeKey k;
            if (!makeKey(nodes, x, y, k))
                throw std::invalid_argument(
                    "EdgePool::insertCellEdges: degenerate self-edge");
            pairs.push_back({x, y, k});
        }
        // Reserve every container that may allocate, up front, so the
        // commit phase below is effectively no-throw for container growth.
        std::vector<EdgeHandle> cellEdges;
        cellEdges.reserve(pairs.size());
        std::vector<EdgeHandle> out;
        out.reserve(pairs.size());
        table_.reserve(table_.size() + pairs.size());
        rec_.reserve(rec_.size() + pairs.size());
        // Commit with rollback. `created` is populated INSIDE createNew, in
        // its noexcept tail, so a throw from any table/allocator step can
        // always be rolled back (review: the old code registered the handle
        // only after createOrGet returned, so a throw from table_.emplace
        // inside it leaked the edge). `bumped` records existing edges whose
        // use count we incremented.
        std::vector<EdgeHandle> created;      // retired on rollback
        std::vector<EdgeHandle> bumped;       // use-count decremented on rb
        created.reserve(pairs.size());
        bumped.reserve(pairs.size());
        // Guarantee rollback's retireEdge -> alloc_.free() cannot allocate:
        // pre-reserve free-list capacity for every edge we might create.
        alloc_.reserveFreeCapacity(pairs.size());
        // Pre-size the incidence edge list so its growth can't throw during
        // the commit tail.
        cellEdges.reserve(pairs.size());
        auto rollback = [&]() noexcept {
            for (auto it = bumped.rbegin(); it != bumped.rend(); ++it) {
                Rec& r = rec_[static_cast<std::size_t>(
                    raw::slotUnchecked(*it))];
                if (r.data.useCount > 0) --r.data.useCount;
            }
            for (auto it = created.rbegin(); it != created.rend(); ++it)
                retireEdge(*it);              // noexcept: see retireEdge
        };
        try {
#if defined(PAACS_TOPO_TESTING)
            maybeFault();                     // fault before any mutation
#endif
            for (const Pair& p : pairs) {
                EdgeHandle eh = lookup(p.key);
                if (!alloc_.isCurrent(eh)) {          // absent -> create
                    eh = createNew(nodes, p.x, p.y, p.key, created);
                }
                Rec& r =
                    rec_[static_cast<std::size_t>(raw::slotUnchecked(eh))];
                ++r.data.useCount;
                bumped.push_back(eh);         // reserved: noexcept
                cellEdges.push_back(eh);      // reserved: noexcept
                out.push_back(eh);            // reserved: noexcept
#if defined(PAACS_TOPO_TESTING)
                maybeFault();                 // fault mid-commit
#endif
            }
            // Incidence is the last thing recorded, only after all edge
            // uses are committed.
            incident_.emplace(cell.slot,
                              IncidenceEntry{cell.gen, std::move(cellEdges)});
        } catch (...) {
            rollback();
            throw;
        }
        return out;
    }

    // Remove every edge use contributed by `cell`; edges whose use count
    // reaches zero are retired. Validates the cell against its CellPool so
    // the caller cannot silently desync the pools. The cell must still be
    // current AND owned by the bound CellPool.
    void removeCellEdges(const CellPool& cells, const CellHandle& cell) {
        if (boundCellPool_ != 0 && cells.poolId() != boundCellPool_)
            throw std::invalid_argument(
                "EdgePool::removeCellEdges: foreign CellPool");
        // The cell must still be current: edge incidence is defined only for
        // current cells, so lifecycle order is enforced (edges first, THEN
        // cell) rather than allowing repair after deletion (review blocker
        // 1). A separate explicit recovery path can be added later if a use
        // case needs it.
        if (!cells.isCurrent(cell))
            throw std::invalid_argument(
                "EdgePool::removeCellEdges: stale cell handle (remove edges "
                "before removing the cell)");
        const auto it = incident_.find(cell.slot);
        if (it == incident_.end() || it->second.gen != cell.gen)
            throw std::invalid_argument(
                "EdgePool::removeCellEdges: cell edges were not inserted");
        // Pre-reserve free-list capacity for every edge that could be
        // retired here, BEFORE mutating anything (review blocker: ordinary
        // removal shares retireEdge -> alloc_.free() -> freeList_.push_back,
        // which can allocate; retireEdge is noexcept, so an allocation
        // failure there would terminate). This throwing step runs before
        // any state change, preserving logical state on failure; after it
        // succeeds the retireEdge calls below cannot allocate.
        alloc_.reserveFreeCapacity(it->second.edges.size());
        for (const EdgeHandle& eh : it->second.edges) {
            if (!alloc_.isCurrent(eh)) continue;   // already retired
            Rec& r = rec_[static_cast<std::size_t>(raw::slotUnchecked(eh))];
            if (r.data.useCount > 0) --r.data.useCount;
            if (r.data.useCount == 0) retireEdge(eh);
        }
        incident_.erase(it);
    }

    void setClassification(const EdgeHandle& h,
                           const mesh::Classification& c) {
        rec_[static_cast<std::size_t>(alloc_.deref(h))].data.cls = c;
    }
    void setRequiredConstraint(const EdgeHandle& h, bool required) {
        rec_[static_cast<std::size_t>(alloc_.deref(h))]
            .data.requiredConstraint = required;
    }

    std::vector<EdgeHandle> activeBySequence() const {
        std::vector<EdgeHandle> out;
        out.reserve(liveLinks_);
        for (std::int32_t s = head_; s != detail::kNil; s = rec_[s].next)
            out.push_back(EdgeHandle{s, rec_[s].gen, alloc_.poolId()});
        return out;
    }

    std::size_t liveCount() const { return alloc_.liveCount(); }

    PoolStats stats() const {
        PoolStats st;
        st.materializedSlots = alloc_.materializedSlots();
        st.live = alloc_.liveCount();
        st.reusableFree = alloc_.reusableFreeSlots();
        st.retired = alloc_.retiredSlots();
        st.totalCreated = alloc_.totalCreated();
        st.peakLive = peakLive_;
        st.peakMaterialized = peakMaterialized_;
        // Complete accounting (review item 6): record array + hash table
        // (buckets + nodes) + incidence map (buckets + nodes) + each
        // incidence entry's edge vector.
        std::size_t tableBytes =
            table_.bucket_count() * sizeof(void*) +
            table_.size() * (sizeof(EdgeKey) + sizeof(EdgeHandle) +
                             sizeof(void*));
        std::size_t incidenceBytes =
            incident_.bucket_count() * sizeof(void*) +
            incident_.size() * (sizeof(std::int32_t) +
                                sizeof(IncidenceEntry) + sizeof(void*));
        for (const auto& [slot, entry] : incident_)
            incidenceBytes += entry.edges.capacity() * sizeof(EdgeHandle);
        st.recordBytes = rec_.capacity() * sizeof(Rec);
        st.estimatedOtherBytes = tableBytes + incidenceBytes;
        st.payloadBytes = st.recordBytes + st.estimatedOtherBytes;
        st.allocatorBytes = alloc_.approxBytes();
        st.totalBytes = st.payloadBytes + st.allocatorBytes;
        return st;
    }

    // Bounds-safe invariant check. The nodes-only form checks edge-record
    // structure and table agreement; the (nodes, cells) form is the REAL
    // topology audit that also proves incidence is backed by current cells
    // and use counts equal incidence counts (review: the nodes-only form
    // passed on ghost edges from an uncoordinated cell removal).
    bool checkInvariants(const NodePool& nodes) const {
        return checkListAndTable(nodes);
    }
    bool checkInvariants(const NodePool& nodes, const CellPool& cells) const {
        return checkListAndTable(nodes) && checkIncidence(nodes, cells);
    }
    bool checkListAndTable(const NodePool& nodes) const {
        if (alloc_.materializedSlots() !=
            alloc_.liveCount() + alloc_.reusableFreeSlots() +
                alloc_.retiredSlots())
            return false;
        if (liveLinks_ != alloc_.liveCount()) return false;
        if (rec_.size() < alloc_.materializedSlots()) return false;
        const auto validSlot = [&](std::int32_t s) {
            return s >= 0 && static_cast<std::size_t>(s) < rec_.size();
        };
        if (liveLinks_ == 0) {
            if (!(head_ == detail::kNil && tail_ == detail::kNil))
                return false;
            return table_.empty();
        }
        if (!validSlot(head_) || !validSlot(tail_)) return false;
        if (rec_[head_].prev != detail::kNil) return false;
        if (rec_[tail_].next != detail::kNil) return false;
        std::vector<std::uint8_t> seen(rec_.size(), 0);
        std::int32_t prevSlot = detail::kNil;
        std::uint64_t prevSeq = 0;
        std::size_t walked = 0;
        for (std::int32_t s = head_; s != detail::kNil; ) {
            if (!validSlot(s)) return false;
            if (seen[static_cast<std::size_t>(s)]) return false;
            seen[static_cast<std::size_t>(s)] = 1;
            const Rec& r = rec_[static_cast<std::size_t>(s)];
            if (r.prev != detail::kNil && !validSlot(r.prev)) return false;
            if (r.prev != prevSlot) return false;
            if (r.next != detail::kNil && !validSlot(r.next)) return false;
            const EdgeHandle h{s, r.gen, alloc_.poolId()};
            if (!alloc_.isCurrent(h)) return false;
            if (alloc_.seqOf(h) != r.seq) return false;
            if (r.data.useCount == 0) return false;   // exists iff used
            if (!nodes.isCurrent(r.data.a) || !nodes.isCurrent(r.data.b))
                return false;
            if (!(nodes.seqOf(r.data.a) < nodes.seqOf(r.data.b)))
                return false;                          // ordered, non-self
            if (walked > 0 && !(r.seq > prevSeq)) return false;
            prevSeq = r.seq;
            prevSlot = s;
            if (++walked > liveLinks_) return false;
            s = r.next;
        }
        if (prevSlot != tail_) return false;
        if (walked != liveLinks_) return false;
        // Table must have exactly one entry per live edge...
        if (table_.size() != liveLinks_) return false;
        // ...and every live edge's stored key must map back to THAT handle
        // in the table (no stale, mis-keyed, or dangling table entries).
        for (std::int32_t s = head_; s != detail::kNil; s = rec_[s].next) {
            const Rec& r = rec_[static_cast<std::size_t>(s)];
            const auto it = table_.find(r.key);
            if (it == table_.end()) return false;
            if (raw::slotUnchecked(it->second) != s) return false;
            if (!alloc_.isCurrent(it->second)) return false;
            // Stored key must match the endpoints' current sequences.
            const std::uint64_t lo =
                std::min(nodes.seqOf(r.data.a), nodes.seqOf(r.data.b));
            const std::uint64_t hi =
                std::max(nodes.seqOf(r.data.a), nodes.seqOf(r.data.b));
            if (r.key.lo != lo || r.key.hi != hi) return false;
        }
        return true;
    }

    // Incidence audit (review blocker): every incidence entry points to a
    // CURRENT cell of the bound CellPool, its recorded edges are current
    // and match the cell's traits-extracted edges, and each live edge's
    // use count equals the number of incident cell-edge uses actually
    // reachable from current cells. This is what catches ghost edges.
    bool checkIncidence(const NodePool& nodes, const CellPool& cells) const {
        if (boundCellPool_ != 0 && cells.poolId() != boundCellPool_)
            return false;
        if (boundNodePool_ != 0 && nodes.poolId() != boundNodePool_)
            return false;
        // Recompute use counts from incidence and compare to stored counts.
        std::unordered_map<std::int32_t, std::uint32_t> counted;
        for (const auto& [slot, entry] : incident_) {
            const CellHandle ch{slot, entry.gen, cells.poolId()};
            if (!cells.isCurrent(ch)) return false;       // ghost incidence
            const auto snap = cells.get(ch);
            if (!isValidElementType(snap.type)) return false;
            const ElementTraits& tr = checkedTraits(snap.type);
            if (static_cast<int>(entry.edges.size()) != tr.edgeCount)
                return false;
            // Each recorded edge must be current and correspond to the
            // cell's e-th trait edge endpoints.
            for (int e = 0; e < tr.edgeCount; ++e) {
                const EdgeHandle eh = entry.edges[static_cast<std::size_t>(e)];
                if (!alloc_.isCurrent(eh)) return false;
                const Rec& r =
                    rec_[static_cast<std::size_t>(raw::slotUnchecked(eh))];
                const NodeHandle x = snap.nodes[tr.edges[e].a];
                const NodeHandle y = snap.nodes[tr.edges[e].b];
                const std::uint64_t lo =
                    std::min(nodes.seqOf(x), nodes.seqOf(y));
                const std::uint64_t hi =
                    std::max(nodes.seqOf(x), nodes.seqOf(y));
                if (r.key.lo != lo || r.key.hi != hi) return false;
                ++counted[raw::slotUnchecked(eh)];
            }
        }
        // Every live edge's stored use count must equal its incidence count,
        // and there must be no live edge with zero reachable uses.
        std::size_t liveSeen = 0;
        for (std::int32_t s = head_; s != detail::kNil; s = rec_[s].next) {
            const Rec& r = rec_[static_cast<std::size_t>(s)];
            const auto it = counted.find(s);
            const std::uint32_t inc =
                (it == counted.end()) ? 0u : it->second;
            if (r.data.useCount != inc) return false;     // count mismatch
            if (inc == 0) return false;                   // ghost edge
            ++liveSeen;
        }
        // No incidence may reference an edge slot that is not live.
        return liveSeen == liveLinks_;
    }

#if defined(PAACS_TOPO_TESTING)
    using AdoptionBeginHook = void (*)() noexcept;
    static void testSetAdoptionBeginHook(AdoptionBeginHook hook) noexcept {
        adoptionBeginHook_ = hook;
    }
    std::uint32_t testUseCount(const EdgeHandle& h) const {
        return rec_[static_cast<std::size_t>(alloc_.deref(h))].data.useCount;
    }
    // Corruption seams for the edge-level invariant probes (mirrors the
    // M4-B.2 pool seams): let tests inject malformed list endpoints, links,
    // table entries and use counts and confirm checkInvariants rejects them
    // without OOB access.
    void testSetHead(std::int32_t s) { head_ = s; }
    void testSetTail(std::int32_t s) { tail_ = s; }
    void testSetPrev(std::int32_t slot, std::int32_t prev) {
        rec_[static_cast<std::size_t>(slot)].prev = prev;
    }
    void testSetNext(std::int32_t slot, std::int32_t next) {
        rec_[static_cast<std::size_t>(slot)].next = next;
    }
    std::int32_t testHead() const { return head_; }
    std::int32_t testMaterialized() const {
        return static_cast<std::int32_t>(alloc_.materializedSlots());
    }
    // Point a live edge's table entry at the wrong handle.
    void testCorruptTableEntry(const EdgeHandle& liveEdge,
                               const EdgeHandle& wrong) {
        const Rec& r =
            rec_[static_cast<std::size_t>(raw::slotUnchecked(liveEdge))];
        table_[r.key] = wrong;
    }
    // Erase a live edge's table entry (missing-entry corruption).
    void testDropTableEntry(const EdgeHandle& liveEdge) {
        const Rec& r =
            rec_[static_cast<std::size_t>(raw::slotUnchecked(liveEdge))];
        table_.erase(r.key);
    }
    // Swap two cells' COMPLETE incidence edge lists (generations stay
    // with their own slots, so hasCellEdges remains true for both).
    // Testing-only addition to the frozen file, disclosed: lets the hash
    // instrument prove sensitivity to cell->edge ASSIGNMENT with nodes,
    // cells, edge records and use counts all held constant.
    void testSwapCellIncidence(std::int32_t slotA, std::int32_t slotB) {
        auto a = incident_.find(slotA), b = incident_.find(slotB);
        if (a != incident_.end() && b != incident_.end())
            a->second.edges.swap(b->second.edges);
    }
    // Force a use count that disagrees with incidence.
    void testSetUseCount(const EdgeHandle& h, std::uint32_t v) {
        rec_[static_cast<std::size_t>(raw::slotUnchecked(h))]
            .data.useCount = v;
    }
#endif

private:
#if defined(PAACS_TOPO_TESTING)
    inline static AdoptionBeginHook adoptionBeginHook_ = nullptr;
#endif
    struct TokenlessAdoptionTag {};
    struct CompleteBuildTag {};

    struct Rec {
        EdgeData data{};
        EdgeKey key{};                    // canonical identity, for O(1)
                                          // table erase at retirement
        std::uint64_t seq = 0;
        std::uint32_t gen = 0;
        std::int32_t prev = detail::kNil;
        std::int32_t next = detail::kNil;
    };
    struct IncidenceEntry {
        std::uint32_t gen = 0;            // cell generation when inserted
        std::vector<EdgeHandle> edges;
    };

    // Stage 1 may throw while ordinary vector/map members establish their
    // empty target state. The source is not referenced until delegation has
    // completed. The allocator alone uses its private tokenless constructor,
    // so target preparation consumes no pool identity.
    explicit EdgePool(TokenlessAdoptionTag)
        : alloc_(HandleAllocator<EdgeTag>::TokenlessAdoptionTag{}) {}

    EdgePool(CompleteBuildTag,const NodePool& nodes,const CellPool& cells,
             EdgePoolCompleteBuildMetrics* metrics) {
        populateComplete(nodes,cells,metrics);
    }

    struct CompleteKeyLess {
        bool operator()(const EdgeKey& a,const EdgeKey& b) const noexcept {
            return a.lo < b.lo || (a.lo == b.lo && a.hi < b.hi);
        }
    };

    static void completeMetricAdd(std::uint64_t& value,
                                  std::uint64_t amount) {
        if (amount > (std::numeric_limits<std::uint64_t>::max)() - value)
            throw std::length_error(
                "EdgePool::buildComplete: metric count overflow");
        value += amount;
    }

    void populateComplete(const NodePool& nodes,const CellPool& cells,
                          EdgePoolCompleteBuildMetrics* publishedMetrics) {
        EdgePoolCompleteBuildMetrics metrics;
        const auto active = cells.activeBySequence();
        std::map<EdgeKey,std::uint64_t,CompleteKeyLess> multiplicity;

        // Pass 1: validate every traits edge and compute exact unique/use
        // counts. Map iteration never determines handle creation order.
        for (const CellHandle cell : active) {
            ++metrics.validationCellVisits;
            if (!cells.isCurrent(cell))
                throw std::invalid_argument(
                    "EdgePool::buildComplete: stale cell handle");
            const auto snap = cells.get(cell);
            const ElementTraits& traits = checkedTraits(snap.type);
            for (int edge = 0; edge < traits.edgeCount; ++edge) {
                const NodeHandle a = snap.nodes[traits.edges[edge].a];
                const NodeHandle b = snap.nodes[traits.edges[edge].b];
                if (!nodes.isCurrent(a) || !nodes.isCurrent(b))
                    throw std::invalid_argument(
                        "EdgePool::buildComplete: stale endpoint");
                EdgeKey key;
                if (!makeKey(nodes,a,b,key))
                    throw std::invalid_argument(
                        "EdgePool::buildComplete: degenerate self-edge");
                auto [it,inserted] = multiplicity.emplace(key,0);
                (void)inserted;
                if (it->second ==
                    (std::numeric_limits<std::uint32_t>::max)())
                    throw std::length_error(
                        "EdgePool::buildComplete: edge incidence overflow");
                ++it->second;
                completeMetricAdd(metrics.edgeUses,1);
            }
        }
        metrics.uniqueEdges = multiplicity.size();
        if (multiplicity.size() > static_cast<std::size_t>(
                (std::numeric_limits<std::int32_t>::max)()))
            throw std::length_error(
                "EdgePool::buildComplete: edge slot limit exceeded");

#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        // Prepare every destination container once. Incidence vectors are
        // then sized exactly once per cell during pass 2.
        alloc_.reserveMaterializedCapacity(multiplicity.size());
        rec_.resize(multiplicity.size());
        table_.reserve(multiplicity.size());
        incident_.reserve(active.size());
        if (!active.empty()) {
            boundNodePool_ = nodes.poolId();
            boundCellPool_ = cells.poolId();
        }

        // Pass 2: cell sequence then local traits order is the sole handle
        // creation order. No hash-table iteration participates.
        for (const CellHandle cell : active) {
            ++metrics.buildCellVisits;
            const auto snap = cells.get(cell);
            const ElementTraits& traits = checkedTraits(snap.type);
            IncidenceEntry entry;
            entry.gen = cell.gen;
            entry.edges.resize(static_cast<std::size_t>(traits.edgeCount));
            for (int edge = 0; edge < traits.edgeCount; ++edge) {
                const NodeHandle a = snap.nodes[traits.edges[edge].a];
                const NodeHandle b = snap.nodes[traits.edges[edge].b];
                EdgeKey key;
                if (!makeKey(nodes,a,b,key))
                    throw std::logic_error(
                        "EdgePool::buildComplete: pass mismatch");
                auto found = table_.find(key);
                EdgeHandle handle;
                if (found == table_.end()) {
                    auto [tableIt,inserted] =
                        table_.emplace(key,EdgeHandle{});
                    if (!inserted)
                        throw std::logic_error(
                            "EdgePool::buildComplete: duplicate table row");
                    handle = alloc_.alloc();
                    const auto slot = static_cast<std::size_t>(
                        raw::slotUnchecked(handle));
                    if (slot >= rec_.size())
                        throw std::logic_error(
                            "EdgePool::buildComplete: unprepared slot");
                    Rec& record = rec_[slot];
                    record.data = EdgeData{};
                    if (nodes.seqOf(a) < nodes.seqOf(b)) {
                        record.data.a = a; record.data.b = b;
                    } else {
                        record.data.a = b; record.data.b = a;
                    }
                    record.key = key;
                    record.seq = alloc_.seqOf(handle);
                    record.gen = handle.gen;
                    linkTail(static_cast<std::int32_t>(slot));
                    tableIt->second = handle;
                    bumpPeaks();
                } else {
                    handle = found->second;
                }
                Rec& record = rec_[static_cast<std::size_t>(
                    raw::slotUnchecked(handle))];
                ++record.data.useCount;
                entry.edges[static_cast<std::size_t>(edge)] = handle;
            }
            const auto [it,inserted] =
                incident_.emplace(cell.slot,std::move(entry));
            (void)it;
            if (!inserted)
                throw std::logic_error(
                    "EdgePool::buildComplete: duplicate cell incidence");
        }
        for (const auto& [key,count] : multiplicity) {
            const auto found = table_.find(key);
            if (found == table_.end() ||
                rec_[static_cast<std::size_t>(raw::slotUnchecked(
                    found->second))].data.useCount != count)
                throw std::logic_error(
                    "EdgePool::buildComplete: multiplicity mismatch");
        }
        ++metrics.finalEdgeAudits;
        if (!auditEdgeTopology(nodes,cells,*this)
                 .acceptable(diag::Severity::Warning))
            throw std::invalid_argument(
                "EdgePool::buildComplete: final audit failed");
        metrics.runComplete = true;
        if (publishedMetrics != nullptr) *publishedMetrics = metrics;
    }

    // Stage 2 is allocation-free and noexcept. Every ownership field moves
    // by exact container/scalar swap; no hashing, reinsertion or element copy
    // occurs after the first source mutation.
    void adoptFrom(EdgePool& o) noexcept {
        alloc_.adoptFrom(o.alloc_);
        static_assert(noexcept(rec_.swap(o.rec_)));
        static_assert(noexcept(table_.swap(o.table_)));
        static_assert(noexcept(incident_.swap(o.incident_)));
        rec_.swap(o.rec_);
        table_.swap(o.table_);
        incident_.swap(o.incident_);
        using std::swap;
        swap(head_, o.head_);
        swap(tail_, o.tail_);
        swap(liveLinks_, o.liveLinks_);
        swap(boundNodePool_, o.boundNodePool_);
        swap(boundCellPool_, o.boundCellPool_);
        swap(peakLive_, o.peakLive_);
        swap(peakMaterialized_, o.peakMaterialized_);
    }

    void bind(const NodePool& nodes, const CellPool& cells) {
        if (boundNodePool_ == 0) {
            boundNodePool_ = nodes.poolId();
            boundCellPool_ = cells.poolId();
        } else if (boundNodePool_ != nodes.poolId() ||
                   boundCellPool_ != cells.poolId()) {
            throw std::invalid_argument(
                "EdgePool: bound to different node/cell pools");
        }
    }
    // Build the canonical (lo,hi)-by-sequence key; false on self-edge.
    bool makeKey(const NodePool& nodes, const NodeHandle& x,
                 const NodeHandle& y, EdgeKey& out) const {
        const std::uint64_t sx = nodes.seqOf(x);
        const std::uint64_t sy = nodes.seqOf(y);
        if (sx == sy) return false;               // self-edge
        out.lo = std::min(sx, sy);
        out.hi = std::max(sx, sy);
        return true;
    }
    // Create a NEW edge record for `key` (caller guarantees key is absent).
    // Exception-safety contract: the ONLY throwing steps happen BEFORE
    // alloc_.alloc(); everything after allocation is noexcept, and the new
    // handle is pushed to `createdOut` in the same noexcept tail so the
    // caller can always roll it back. table_ node is materialized by an
    // emplace of a placeholder BEFORE alloc, then overwritten with the real
    // handle (no post-alloc map allocation).
    EdgeHandle createNew(const NodePool& nodes, const NodeHandle& x,
                         const NodeHandle& y, const EdgeKey& key,
                         std::vector<EdgeHandle>& createdOut) {
        // Slot-space exhaustion precheck (before any mutation) when no free
        // slot is available.
        if (alloc_.reusableFreeSlots() == 0 &&
            HandleAllocator<EdgeTag>::wouldExhaustSlots(
                alloc_.materializedSlots()))
            throw std::length_error("EdgePool: slot space exhausted");
        // Grow payload for the worst-case slot BEFORE alloc (throwing step).
        const std::size_t worst = alloc_.materializedSlots() + 1;
        if (worst > rec_.size()) {
            detail::reserveFor(rec_, worst);
            rec_.resize(worst);                       // may throw: pre-alloc
        }
        // Materialize the table node NOW (throwing step) with a placeholder;
        // we overwrite the mapped value after alloc with no further
        // allocation. If a later step throws we erase this key (noexcept).
        auto [tabIt, inserted] = table_.emplace(key, EdgeHandle{});
        (void)inserted;                               // guaranteed new key
        try {
            const EdgeHandle h = alloc_.alloc();      // may throw: table key
                                                      // erased in catch
            // ---- noexcept tail ----
            const auto s = static_cast<std::size_t>(raw::slotUnchecked(h));
            if (s >= rec_.size()) rec_.resize(s + 1);  // reused slot < size
            Rec& r = rec_[s];
            r.data = EdgeData{};
            if (nodes.seqOf(x) < nodes.seqOf(y)) { r.data.a = x; r.data.b = y; }
            else { r.data.a = y; r.data.b = x; }
            r.data.useCount = 0;
            r.key = key;
            r.seq = alloc_.seqOf(h);
            r.gen = h.gen;
            linkTail(static_cast<std::int32_t>(s));
            tabIt->second = h;                        // overwrite placeholder
            createdOut.push_back(h);                  // reserved: noexcept
            bumpPeaks();
            return h;
        } catch (...) {
            table_.erase(tabIt);                      // noexcept: undo node
            throw;
        }
    }
    EdgeHandle lookup(const EdgeKey& key) const {
        const auto it = table_.find(key);
        return (it == table_.end()) ? EdgeHandle{} : it->second;
    }
    // Retire an edge. Used both on use-count->0 and during insertion
    // rollback. The table entry and list link (both noexcept) are removed
    // first, then the slot is freed. During insertion rollback the caller
    // has pre-reserved free-list capacity (alloc_.reserveFreeCapacity), so
    // free() cannot allocate and this is genuinely noexcept.
    void retireEdge(const EdgeHandle& h) noexcept {
        const auto s = static_cast<std::int32_t>(raw::slotUnchecked(h));
        table_.erase(rec_[static_cast<std::size_t>(s)].key);  // noexcept
        unlink(s);                                            // noexcept
        alloc_.free(h);          // noexcept here: free list pre-reserved
    }

    void linkTail(std::int32_t s) {
        rec_[s].prev = tail_;
        rec_[s].next = detail::kNil;
        if (tail_ != detail::kNil) rec_[tail_].next = s;
        else head_ = s;
        tail_ = s;
        ++liveLinks_;
    }
    void unlink(std::int32_t s) {
        Rec& r = rec_[static_cast<std::size_t>(s)];
        if (r.prev != detail::kNil) rec_[r.prev].next = r.next;
        else head_ = r.next;
        if (r.next != detail::kNil) rec_[r.next].prev = r.prev;
        else tail_ = r.prev;
        r.prev = r.next = detail::kNil;
        --liveLinks_;
    }
    void bumpPeaks() {
        peakLive_ = std::max(peakLive_, alloc_.liveCount());
        peakMaterialized_ = std::max(peakMaterialized_,
                                     alloc_.materializedSlots());
    }

    HandleAllocator<EdgeTag> alloc_;
    std::vector<Rec> rec_;
    std::unordered_map<EdgeKey, EdgeHandle, EdgeKeyHash> table_;
    // Incidence is a SLOT-INDEXED table with generation validation (NOT a
    // full (pool,slot,gen) key): the value stores the cell generation, and
    // every lookup checks it, so a reused slot with stale incidence is
    // rejected rather than silently overwritten. The pool is implicit
    // (boundCellPool_). This is sufficient because the EdgePool binds to
    // exactly one CellPool; a full composite key would add no safety here.
    std::unordered_map<std::int32_t, IncidenceEntry> incident_;
    std::int32_t head_ = detail::kNil;
    std::int32_t tail_ = detail::kNil;
    std::size_t liveLinks_ = 0;
    std::uint64_t boundNodePool_ = 0;
    std::uint64_t boundCellPool_ = 0;
    std::size_t peakLive_ = 0;
    std::size_t peakMaterialized_ = 0;

    static_assert(std::is_nothrow_swappable_v<decltype(rec_)>);
    static_assert(std::is_nothrow_swappable_v<decltype(table_)>);
    static_assert(std::is_nothrow_swappable_v<decltype(incident_)>);
};

static_assert(!std::is_nothrow_move_constructible_v<EdgePool>);
static_assert(std::is_nothrow_move_constructible_v<
              HandleAllocator<EdgeTag>>);

// ---------------------------------------------------------------------------
// Structured edge-topology audit (the M4-B.4 prerequisite gate). Short-
// circuits: if the storage pair is not edge-ready, or the edge records /
// table are inconsistent, it does NOT proceed to incidence traversal.
// ---------------------------------------------------------------------------
inline diag::Report auditEdgeTopology(const NodePool& nodes,
                                      const CellPool& cells,
                                      const EdgePool& edges) {
    diag::Report r;
    auto add = [&](diag::Code code, diag::Severity sev, std::string msg) {
        diag::Diagnostic d;
        d.code = code;
        d.severity = sev;
        d.stage = "edge-audit";
        d.message = std::move(msg);
        r.add(std::move(d));
    };
    // Pool binding must agree.
    if ((edges.boundNodePool() != 0 &&
         edges.boundNodePool() != nodes.poolId()) ||
        (edges.boundCellPool() != 0 &&
         edges.boundCellPool() != cells.poolId())) {
        add(diag::Code::EdgePoolBindingMismatch,
            diag::Severity::InvalidOutput,
            "EdgePool is bound to different node/cell pools than audited");
        return r;
    }
    // Underlying storage must be sound first.
    if (!storageEdgeReady(nodes, cells)) {
        add(diag::Code::StorageInvariantViolation,
            diag::Severity::InvalidOutput,
            "underlying node/cell storage is not edge-ready");
        return r;
    }
    // Edge record + list + table structure.
    if (!edges.checkInvariants(nodes)) {
        add(diag::Code::EdgeInvariantViolation,
            diag::Severity::InvalidOutput,
            "edge records / list / table invariants violated");
        return r;
    }
    // Incidence: ghost edges, stale incident cells, use-count mismatch.
    if (!edges.checkInvariants(nodes, cells))
        add(diag::Code::GhostEdge,
            diag::Severity::InvalidOutput,
            "edge incidence is inconsistent: a live edge lacks current-cell "
            "backing, an incident cell is stale, or a use count disagrees "
            "with incidence");
    return r;
}

inline bool edgeTopologyReady(const NodePool& nodes, const CellPool& cells,
                              const EdgePool& edges) {
    return auditEdgeTopology(nodes, cells, edges)
        .acceptable(diag::Severity::Warning);
}

} // namespace paacs::topo
