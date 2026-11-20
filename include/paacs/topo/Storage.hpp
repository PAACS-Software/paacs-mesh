// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// paacs::topo - M4-B.2: generic node and mixed-cell storage composed ON
// the frozen M4-B.1 HandleAllocator. Adds PAYLOAD only; the allocator is
// used as-is and never modified.
//
// DATA STRUCTURE:
//   * Each pool is a slot-indexed record array with INTRUSIVE prev/next
//     links (by slot) forming a doubly-linked list in CREATION-SEQUENCE
//     order. add() amortized O(1), remove() O(1), iterate O(live),
//     audit O(materialized).
//   * Cell connectivity is stored INLINE in each slot record as a fixed
//     std::array<LocalNodeRef, kMaxCellNodes>, where LocalNodeRef is
//     {slot, gen} (8 bytes; the pool token is implicit, the CellPool being
//     bound to one NodePool). On slot reuse the inline block is simply
//     overwritten, so cell memory is bounded by PEAK live+free+retired
//     slots - there is NO connectivity arena and NO historical retention.
//     A Tet4 record carries 10 * 8 = 80 B of connectivity capacity; only
//     `count` entries are meaningful. (An earlier draft used an
//     append-only arena that retained connectivity for every historical
//     insertion - removed.)
//
// API SAFETY (round-5/7 review): ordinary accessors return payload BY VALUE
// (NodeData / CellSnapshot), and CellSnapshot::connectivity() is
// lvalue-only (rvalue overload deleted) so no span escapes a temporary.
// forEachActive intentionally provides a callback-scoped BORROW of the
// payload. Precise contract: no ordinary accessor returns an escaping
// reference or span; callback views (forEachActive) are valid ONLY for the
// callback's duration and must not mutate the pool. Ordinary insertion may
// reallocate, but because no ordinary accessor hands out an interior view,
// it can never dangle a caller.
//
// Correctness contracts (unchanged): materialized == live + reusableFree +
// retired; cells own node identity and bind to one NodePool;
// dangling/foreign connectivity detectable; auditStorage() is the M4-B.3
// gate; strong-safety insertion (throwing staging precedes alloc; logical
// state preserved on failure, capacity may grow); runtime pool
// ids/slots/generations never enter hashes; allocator+payload move as one.
//
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <functional>
#include <new>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paacs/core/Types.hpp"
#include "paacs/diag/Diagnostics.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/topo/ElementTraits.hpp"
#include "paacs/topo/Handles.hpp"

namespace paacs::topo {

inline constexpr int kMaxCellNodes = 10;

#if defined(PAACS_TOPO_TESTING)
inline int& faultCountdown() { static thread_local int n = 0; return n; }
inline void maybeFault() {
    int& n = faultCountdown();
    if (n > 0 && --n == 0) throw std::bad_alloc();
}
#endif

struct PoolStats {
    std::size_t materializedSlots = 0;
    std::size_t live = 0;
    std::size_t reusableFree = 0;
    std::size_t retired = 0;
    std::uint64_t totalCreated = 0;
    std::size_t peakLive = 0;
    std::size_t peakMaterialized = 0;
    // payloadBytes keeps its ORIGINAL meaning (records + any
    // table/incidence estimates) for compatibility. M4-B.7 C1 review
    // added the narrower classes below, because payloadBytes conflates
    // known record capacity with hash/incidence ESTIMATES and made a
    // subtraction-based split identically zero:
    //   payloadBytes == recordBytes + estimatedOtherBytes
    //   totalBytes   == recordBytes + allocatorBytes + estimatedOtherBytes
    std::size_t payloadBytes = 0;
    std::size_t recordBytes = 0;          // record-array capacity ONLY
    std::size_t estimatedOtherBytes = 0;  // hash/incidence ESTIMATES
    std::size_t allocatorBytes = 0;
    std::size_t totalBytes = 0;
    bool consistent() const {
        return materializedSlots == live + reusableFree + retired;
    }
    bool bytesConsistent() const {
        return payloadBytes == recordBytes + estimatedOtherBytes &&
               totalBytes ==
                   recordBytes + allocatorBytes + estimatedOtherBytes;
    }
};

namespace detail {
constexpr std::int32_t kNil = -1;
// Reserve geometrically so appends are amortized O(1).
template <class V>
inline void reserveFor(V& v, std::size_t need) {
    if (need > v.capacity()) {
        std::size_t cap = v.capacity();
        cap = cap < 8 ? 8 : cap + cap / 2;
        v.reserve(std::max(cap, need));
    }
}
} // namespace detail

// ---------------------------------------------------------------------------
// NodePool
// ---------------------------------------------------------------------------
class NodePool {
    friend class detail::TransactionAccess;   // M4-B.6 prepared commit
public:
    struct NodeData {
        Vec3 coord{0, 0, 0};
        mesh::Classification cls{};
    };

    NodePool() = default;
    NodePool(const NodePool&) = delete;
    NodePool& operator=(const NodePool&) = delete;
    NodePool(NodePool&& o) noexcept
        : alloc_(std::move(o.alloc_)), rec_(std::move(o.rec_)),
          head_(o.head_), tail_(o.tail_), liveLinks_(o.liveLinks_),
          peakLive_(o.peakLive_), peakMaterialized_(o.peakMaterialized_) {
        o.rec_.clear();
        o.head_ = o.tail_ = detail::kNil;
        o.liveLinks_ = 0;
        o.peakLive_ = o.peakMaterialized_ = 0;
    }
    NodePool& operator=(NodePool&&) = delete;

    NodeHandle add(const Vec3& p, const mesh::Classification& c) {
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        // Prepare only the rows actually required (review defect 2): if a
        // reusable free slot exists the allocator will reuse it (no new
        // materialized slot), so no payload growth is needed. If NO free
        // slot exists, precheck the allocator's documented slot-space
        // exhaustion BEFORE touching payload capacity, so a bad_alloc from
        // payload growth can't pre-empt the structured length_error the
        // M4-B.1 contract promises.
        if (alloc_.reusableFreeSlots() == 0) {
            if (HandleAllocator<NodeTag>::wouldExhaustSlots(
                    alloc_.materializedSlots()))
                throw std::length_error(
                    "NodePool::add: slot space exhausted");
            const std::size_t required = alloc_.materializedSlots() + 1;
            if (required > rec_.size()) {
                detail::reserveFor(rec_, required);
                rec_.resize(required);
            }
        }
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        const NodeHandle h = alloc_.alloc();
        const auto s = static_cast<std::size_t>(raw::slotUnchecked(h));
        if (s >= rec_.size()) rec_.resize(s + 1);   // reuse-safety net
        Rec& r = rec_[s];
        r.data = NodeData{p, c};
        r.seq = alloc_.seqOf(h);
        r.gen = h.gen;
        linkTail(static_cast<std::int32_t>(s));
        bumpPeaks();
        return h;
    }

    bool isCurrent(const NodeHandle& h) const { return alloc_.isCurrent(h); }
    std::uint64_t poolId() const { return alloc_.poolId(); }

    // BY VALUE (review: returning a reference into rec_ dangled after a
    // reallocating add()).
    NodeData get(const NodeHandle& h) const {
        return rec_[static_cast<std::size_t>(alloc_.deref(h))].data;
    }
    void setCoord(const NodeHandle& h, const Vec3& p) {
        rec_[static_cast<std::size_t>(alloc_.deref(h))].data.coord = p;
    }
    void setClassification(const NodeHandle& h,
                           const mesh::Classification& c) {
        rec_[static_cast<std::size_t>(alloc_.deref(h))].data.cls = c;
    }

    void remove(const NodeHandle& h) {
        const std::int32_t s = static_cast<std::int32_t>(alloc_.deref(h));
        alloc_.free(h);
        unlink(s);
    }

    std::uint64_t seqOf(const NodeHandle& h) const { return alloc_.seqOf(h); }

    std::vector<NodeHandle> activeBySequence() const {
        std::vector<NodeHandle> out;
        out.reserve(liveLinks_);
        for (std::int32_t s = head_; s != detail::kNil; s = rec_[s].next)
            out.push_back(NodeHandle{s, rec_[s].gen, alloc_.poolId()});
        return out;
    }

    // Borrow-scoped iteration: the NodeData is valid ONLY during the
    // callback; the pool must not be mutated inside it.
    template <class Fn>
    void forEachActive(Fn&& fn) const {
        for (std::int32_t s = head_; s != detail::kNil; s = rec_[s].next)
            fn(NodeHandle{s, rec_[s].gen, alloc_.poolId()}, rec_[s].data);
    }

    bool checkInvariants() const { return checkList(); }

    PoolStats stats() const {
        PoolStats st;
        st.materializedSlots = alloc_.materializedSlots();
        st.live = alloc_.liveCount();
        st.reusableFree = alloc_.reusableFreeSlots();
        st.retired = alloc_.retiredSlots();
        st.totalCreated = alloc_.totalCreated();
        st.peakLive = peakLive_;
        st.peakMaterialized = peakMaterialized_;
        st.recordBytes = rec_.capacity() * sizeof(Rec);
        st.estimatedOtherBytes = 0;      // no tables/incidence here
        st.payloadBytes = st.recordBytes + st.estimatedOtherBytes;
        st.allocatorBytes = alloc_.approxBytes();
        st.totalBytes = st.payloadBytes + st.allocatorBytes;
        return st;
    }

#if defined(PAACS_TOPO_TESTING)
    void testForceSeqMax() { alloc_.testForceSequenceMax(); }
    void testRetire(const NodeHandle& h) {
        const std::int32_t s = static_cast<std::int32_t>(alloc_.deref(h));
        alloc_.testForceGeneration(h, std::numeric_limits<std::uint32_t>::max());
        rec_[static_cast<std::size_t>(s)].gen =
            std::numeric_limits<std::uint32_t>::max();
        alloc_.free(NodeHandle{s, std::numeric_limits<std::uint32_t>::max(),
                               alloc_.poolId()});
        unlink(s);
    }
    // Corruption probes for the invariant checker (see checkList).
    void testSetTail(std::int32_t s) { tail_ = s; }
    void testSetHead(std::int32_t s) { head_ = s; }
    void testSetPrev(std::int32_t slot, std::int32_t prev) {
        rec_[static_cast<std::size_t>(slot)].prev = prev;
    }
    void testSetNext(std::int32_t slot, std::int32_t next) {
        rec_[static_cast<std::size_t>(slot)].next = next;
    }
    void testSetSeq(const NodeHandle& h, std::uint64_t seq) {
        rec_[static_cast<std::size_t>(h.slot)].seq = seq;
    }
    std::int32_t testHead() const { return head_; }
    std::int32_t testMaterialized() const {
        return static_cast<std::int32_t>(alloc_.materializedSlots());
    }
#endif

private:
    struct Rec {
        NodeData data{};
        std::uint64_t seq = 0;
        std::uint32_t gen = 0;
        std::int32_t prev = detail::kNil;
        std::int32_t next = detail::kNil;
    };

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
    // COMPLETE doubly-linked-list verification (review item 1): empty
    // consistency, head.prev==nil, each prev matches the previous visited
    // record, last visited == tail_, tail.next==nil, visited==live, stored
    // seq == allocator seq, and every current handle valid, all O(n) via a
    // slot bitmap.
    // Bounds-SAFE list auditor (review blocker): a corruption auditor must
    // NEVER dereference an unvalidated index. head_, tail_, and every prev
    // and next link is range-checked BEFORE it indexes rec_. A malformed
    // endpoint or link returns false; it can never segfault.
    bool checkList() const {
        if (alloc_.materializedSlots() !=
            alloc_.liveCount() + alloc_.reusableFreeSlots() +
                alloc_.retiredSlots())
            return false;
        if (liveLinks_ != alloc_.liveCount()) return false;
        // rec_ must be at least as large as the materialized slot count, or
        // a valid slot index could still be out of the payload array.
        if (rec_.size() < alloc_.materializedSlots()) return false;
        const auto validSlot = [&](std::int32_t s) {
            return s >= 0 && static_cast<std::size_t>(s) < rec_.size();
        };
        if (liveLinks_ == 0)
            return head_ == detail::kNil && tail_ == detail::kNil;
        // Nonempty: BOTH endpoints must be valid record indices (not just
        // "not nil") before any dereference.
        if (!validSlot(head_) || !validSlot(tail_)) return false;
        if (rec_[head_].prev != detail::kNil) return false;
        if (rec_[tail_].next != detail::kNil) return false;
        std::vector<std::uint8_t> seen(rec_.size(), 0);
        std::int32_t prevSlot = detail::kNil;
        std::uint64_t prevSeq = 0;
        std::size_t walked = 0;
        std::int32_t s = head_;
        while (s != detail::kNil) {
            if (!validSlot(s)) return false;                 // guard current
            if (seen[static_cast<std::size_t>(s)]) return false;  // cycle
            seen[static_cast<std::size_t>(s)] = 1;
            const Rec& r = rec_[static_cast<std::size_t>(s)];
            // prev must be nil or a valid slot, and must equal the last
            // visited slot.
            if (r.prev != detail::kNil && !validSlot(r.prev)) return false;
            if (r.prev != prevSlot) return false;
            // next must be nil or a valid slot BEFORE we follow it.
            if (r.next != detail::kNil && !validSlot(r.next)) return false;
            const NodeHandle h{s, r.gen, alloc_.poolId()};
            if (!alloc_.isCurrent(h)) return false;
            if (alloc_.seqOf(h) != r.seq) return false;
            if (walked > 0 && !(r.seq > prevSeq)) return false;
            prevSeq = r.seq;
            prevSlot = s;
            if (++walked > liveLinks_) return false;         // overrun
            s = r.next;
        }
        if (prevSlot != tail_) return false;                 // last == tail
        return walked == liveLinks_;
    }
    void bumpPeaks() {
        peakLive_ = std::max(peakLive_, alloc_.liveCount());
        peakMaterialized_ = std::max(peakMaterialized_,
                                     alloc_.materializedSlots());
    }

    HandleAllocator<NodeTag> alloc_;
    std::vector<Rec> rec_;
    std::int32_t head_ = detail::kNil;
    std::int32_t tail_ = detail::kNil;
    std::size_t liveLinks_ = 0;
    std::size_t peakLive_ = 0;
    std::size_t peakMaterialized_ = 0;
};

// ---------------------------------------------------------------------------
// CellPool - compact connectivity (inline std::array<LocalNodeRef,10> per
// record; NO arena), intrusive sequence list, value-returning API.
// ---------------------------------------------------------------------------
class CellPool {
    friend class detail::TransactionAccess;   // M4-B.6 prepared commit
public:
    // Compact per-cell node reference: slot+generation only, 8 bytes. The
    // pool token is implicit (boundNodePool_), so it is not repeated per
    // connectivity entry (review item 3).
    struct LocalNodeRef {
        std::int32_t slot = -1;
        std::uint32_t gen = 0;
    };

    // Isolated snapshot returned BY VALUE (review item 2). Holds up to
    // kMaxCellNodes reconstructed NodeHandles plus a count; no borrow into
    // pool memory.
    struct CellSnapshot {
        ElementType type{};
        mesh::Classification cls{};
        std::array<NodeHandle, kMaxCellNodes> nodes{};
        std::uint8_t count = 0;
        // Ref-qualified so a span is NEVER returned from a temporary
        // snapshot (review defect 2: cells.get(h).connectivity() left a
        // dangling span into the dead temporary - ASan stack-use-after-
        // scope). Callable only on a named lvalue snapshot; the rvalue
        // overload is deleted. Callers on a temporary must use
        // get(h).nodes with get(h).count instead.
        std::span<const NodeHandle> connectivity() const& noexcept {
            return std::span<const NodeHandle>(nodes.data(), count);
        }
        std::span<const NodeHandle> connectivity() const&& = delete;
    };

    CellPool() = default;
    CellPool(const CellPool&) = delete;
    CellPool& operator=(const CellPool&) = delete;
    CellPool(CellPool&& o) noexcept
        : alloc_(std::move(o.alloc_)), rec_(std::move(o.rec_)),
          head_(o.head_), tail_(o.tail_),
          liveLinks_(o.liveLinks_), boundNodePool_(o.boundNodePool_),
          peakLive_(o.peakLive_), peakMaterialized_(o.peakMaterialized_) {
        o.rec_.clear();
        o.head_ = o.tail_ = detail::kNil;
        o.liveLinks_ = 0;
        o.boundNodePool_ = 0;
        o.peakLive_ = o.peakMaterialized_ = 0;
    }
    CellPool& operator=(CellPool&&) = delete;

    // span overload (review carry-item): callers need not heap-allocate a
    // vector per cell.
    CellHandle add(ElementType type, std::span<const NodeHandle> nodes,
                   const mesh::Classification& c, const NodePool& owner) {
        (void)checkedTraits(type);
        if (boundNodePool_ != 0 && boundNodePool_ != owner.poolId())
            throw std::invalid_argument(
                "CellPool::add: cell references a different NodePool than "
                "this pool is bound to");
        if (nodes.size() > static_cast<std::size_t>(kMaxCellNodes))
            throw std::invalid_argument(
                "CellPool::add: connectivity exceeds kMaxCellNodes");
        const bool ok = validateHandleConnectivity<NodeHandle>(
            type, nodes, [&](const NodeHandle& nh) {
                return nh.pool == owner.poolId() && owner.isCurrent(nh);
            });
        if (!ok)
            throw std::invalid_argument(
                "CellPool::add: connectivity fails traits validation "
                "(arity / stale or foreign node / repeated node)");
        const std::size_t k = nodes.size();
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        if (alloc_.reusableFreeSlots() == 0) {
            if (HandleAllocator<CellTag>::wouldExhaustSlots(
                    alloc_.materializedSlots()))
                throw std::length_error(
                    "CellPool::add: slot space exhausted");
            const std::size_t required = alloc_.materializedSlots() + 1;
            if (required > rec_.size()) {
                detail::reserveFor(rec_, required);
                rec_.resize(required);
            }
        }
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        const CellHandle h = alloc_.alloc();
        const auto s = static_cast<std::size_t>(raw::slotUnchecked(h));
        if (s >= rec_.size()) rec_.resize(s + 1);   // reuse-safety net
        Rec& r = rec_[s];
        r.type = type;
        r.count = static_cast<std::uint8_t>(k);
        for (std::size_t i = 0; i < k; ++i)
            r.conn[i] = LocalNodeRef{nodes[i].slot, nodes[i].gen};  // INLINE
        r.cls = c;
        r.seq = alloc_.seqOf(h);
        r.gen = h.gen;
        linkTail(static_cast<std::int32_t>(s));
        boundNodePool_ = owner.poolId();
        bumpPeaks();
        return h;
    }
    // Convenience for vector callers.
    CellHandle add(ElementType type, const std::vector<NodeHandle>& nodes,
                   const mesh::Classification& c, const NodePool& owner) {
        return add(type, std::span<const NodeHandle>(nodes.data(),
                                                     nodes.size()),
                   c, owner);
    }

    bool isCurrent(const CellHandle& h) const { return alloc_.isCurrent(h); }
    std::uint64_t poolId() const { return alloc_.poolId(); }
    std::uint64_t boundNodePool() const { return boundNodePool_; }
    bool isBound() const { return boundNodePool_ != 0; }
    // BINDING SEMANTICS (review carry-item): a CellPool binds to its
    // NodePool on the first insert and STAYS bound for its lifetime, even
    // if every cell is later removed (removal does not unbind). This keeps
    // the bound-pool invariant stable for audits and edge topology. There
    // is deliberately no unbind/reset in M4-B.2 - re-targeting a populated
    // or previously-populated pool is a mutation the transaction layer
    // (later) should own, not a silent storage side effect. A never-used
    // CellPool (boundNodePool_ == 0) accepts any NodePool on first insert.

    // BY VALUE (review item 2): a self-contained snapshot with full
    // NodeHandles reconstructed from the bound pool id.
    CellSnapshot get(const CellHandle& h) const {
        const Rec& r = rec_[static_cast<std::size_t>(alloc_.deref(h))];
        CellSnapshot snap;
        snap.type = r.type;
        snap.cls = r.cls;
        snap.count = r.count;
        for (std::uint8_t i = 0; i < r.count; ++i)
            snap.nodes[i] = NodeHandle{r.conn[i].slot, r.conn[i].gen,
                                       boundNodePool_};
        return snap;
    }

    void setClassification(const CellHandle& h,
                           const mesh::Classification& c) {
        rec_[static_cast<std::size_t>(alloc_.deref(h))].cls = c;
    }
    void remove(const CellHandle& h) {
        const std::int32_t s = static_cast<std::int32_t>(alloc_.deref(h));
        alloc_.free(h);
        unlink(s);
        // Connectivity is inline in the slot record; on reuse it is simply
        // overwritten. No arena, no tombstones, no historical retention:
        // cell memory is bounded by PEAK live+free+retired slots, never by
        // total insertions (review found the old append-only arena retained
        // ~40 MB after 1M create/remove of a single cell).
    }
    std::uint64_t seqOf(const CellHandle& h) const { return alloc_.seqOf(h); }

    bool connectivityIntact(const NodePool& owner) const {
        if (boundNodePool_ != 0 && boundNodePool_ != owner.poolId())
            return false;
        for (std::int32_t s = head_; s != detail::kNil; s = rec_[s].next) {
            const Rec& r = rec_[static_cast<std::size_t>(s)];
            for (std::uint8_t i = 0; i < r.count; ++i)
                if (!owner.isCurrent(NodeHandle{r.conn[i].slot,
                                                r.conn[i].gen,
                                                boundNodePool_}))
                    return false;
        }
        return true;
    }

    std::vector<CellHandle> activeBySequence() const {
        std::vector<CellHandle> out;
        out.reserve(liveLinks_);
        for (std::int32_t s = head_; s != detail::kNil; s = rec_[s].next)
            out.push_back(CellHandle{s, rec_[s].gen, alloc_.poolId()});
        return out;
    }

    // Borrow-scoped iteration: the span is valid ONLY during the callback;
    // the pool must not be mutated inside it. Uses a small stack buffer, no
    // heap-borrow leaks.
    template <class Fn>
    void forEachActive(Fn&& fn) const {
        std::array<NodeHandle, kMaxCellNodes> buf{};
        for (std::int32_t s = head_; s != detail::kNil; s = rec_[s].next) {
            const Rec& r = rec_[static_cast<std::size_t>(s)];
            for (std::uint8_t i = 0; i < r.count; ++i)
                buf[i] = NodeHandle{r.conn[i].slot, r.conn[i].gen,
                                    boundNodePool_};
            fn(CellHandle{s, r.gen, alloc_.poolId()}, r.type,
               std::span<const NodeHandle>(buf.data(), r.count), r.cls);
        }
    }

    bool checkInvariants() const { return checkList(); }

    PoolStats stats() const {
        PoolStats st;
        st.materializedSlots = alloc_.materializedSlots();
        st.live = alloc_.liveCount();
        st.reusableFree = alloc_.reusableFreeSlots();
        st.retired = alloc_.retiredSlots();
        st.totalCreated = alloc_.totalCreated();
        st.peakLive = peakLive_;
        st.peakMaterialized = peakMaterialized_;
        st.recordBytes = rec_.capacity() * sizeof(Rec);
        st.estimatedOtherBytes = 0;      // no tables/incidence here
        st.payloadBytes = st.recordBytes + st.estimatedOtherBytes;
        st.allocatorBytes = alloc_.approxBytes();
        st.totalBytes = st.payloadBytes + st.allocatorBytes;
        return st;
    }

#if defined(PAACS_TOPO_TESTING)
    void testRetire(const CellHandle& h) {
        const std::int32_t s = static_cast<std::int32_t>(alloc_.deref(h));
        alloc_.testForceGeneration(h, std::numeric_limits<std::uint32_t>::max());
        rec_[static_cast<std::size_t>(s)].gen =
            std::numeric_limits<std::uint32_t>::max();
        alloc_.free(CellHandle{s, std::numeric_limits<std::uint32_t>::max(),
                               alloc_.poolId()});
        unlink(s);
    }
    void testSetTail(std::int32_t s) { tail_ = s; }
    void testSetHead(std::int32_t s) { head_ = s; }
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
#endif

private:
    struct Rec {
        std::array<LocalNodeRef, kMaxCellNodes> conn{};  // INLINE, no arena
        mesh::Classification cls{};
        std::uint64_t seq = 0;
        std::uint32_t gen = 0;
        ElementType type{};
        std::uint8_t count = 0;
        std::int32_t prev = detail::kNil;
        std::int32_t next = detail::kNil;
    };

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
    // Bounds-SAFE list auditor (review blocker): identical safety contract
    // to NodePool::checkList - every endpoint and link is range-checked
    // before it indexes rec_, so corrupt state returns false, never
    // segfaults.
    bool checkList() const {
        if (alloc_.materializedSlots() !=
            alloc_.liveCount() + alloc_.reusableFreeSlots() +
                alloc_.retiredSlots())
            return false;
        if (liveLinks_ != alloc_.liveCount()) return false;
        if (rec_.size() < alloc_.materializedSlots()) return false;
        const auto validSlot = [&](std::int32_t s) {
            return s >= 0 && static_cast<std::size_t>(s) < rec_.size();
        };
        if (liveLinks_ == 0)
            return head_ == detail::kNil && tail_ == detail::kNil;
        if (!validSlot(head_) || !validSlot(tail_)) return false;
        if (rec_[head_].prev != detail::kNil) return false;
        if (rec_[tail_].next != detail::kNil) return false;
        std::vector<std::uint8_t> seen(rec_.size(), 0);
        std::int32_t prevSlot = detail::kNil;
        std::uint64_t prevSeq = 0;
        std::size_t walked = 0;
        std::int32_t s = head_;
        while (s != detail::kNil) {
            if (!validSlot(s)) return false;
            if (seen[static_cast<std::size_t>(s)]) return false;
            seen[static_cast<std::size_t>(s)] = 1;
            const Rec& r = rec_[static_cast<std::size_t>(s)];
            if (r.prev != detail::kNil && !validSlot(r.prev)) return false;
            if (r.prev != prevSlot) return false;
            if (r.next != detail::kNil && !validSlot(r.next)) return false;
            const CellHandle h{s, r.gen, alloc_.poolId()};
            if (!alloc_.isCurrent(h)) return false;
            if (alloc_.seqOf(h) != r.seq) return false;
            if (!isValidElementType(r.type)) return false;
            if (r.count != traits(r.type).nodeCount) return false;
            if (boundNodePool_ == 0) return false;
            for (std::uint8_t i = 0; i < r.count; ++i)
                for (std::uint8_t j = i + 1; j < r.count; ++j)
                    if (r.conn[i].slot == r.conn[j].slot &&
                        r.conn[i].gen == r.conn[j].gen)
                        return false;
            if (walked > 0 && !(r.seq > prevSeq)) return false;
            prevSeq = r.seq;
            prevSlot = s;
            if (++walked > liveLinks_) return false;
            s = r.next;
        }
        if (prevSlot != tail_) return false;
        return walked == liveLinks_;
    }
    void bumpPeaks() {
        peakLive_ = std::max(peakLive_, alloc_.liveCount());
        peakMaterialized_ = std::max(peakMaterialized_,
                                     alloc_.materializedSlots());
    }

    HandleAllocator<CellTag> alloc_;
    std::vector<Rec> rec_;
    std::int32_t head_ = detail::kNil;
    std::int32_t tail_ = detail::kNil;
    std::size_t liveLinks_ = 0;
    std::uint64_t boundNodePool_ = 0;
    std::size_t peakLive_ = 0;
    std::size_t peakMaterialized_ = 0;
};

// ---------------------------------------------------------------------------
// Storage validity gate (short-circuits before connectivity traversal).
// ---------------------------------------------------------------------------
inline diag::Report auditStorage(const NodePool& nodes,
                                 const CellPool& cells) {
    diag::Report r;
    auto add = [&](diag::Code code, diag::Severity sev, std::string msg) {
        diag::Diagnostic d;
        d.code = code;
        d.severity = sev;
        d.stage = "storage-audit";
        d.message = std::move(msg);
        r.add(std::move(d));
    };
    const bool nodesOk = nodes.checkInvariants();
    const bool cellsOk = cells.checkInvariants();
    if (!nodesOk)
        add(diag::Code::StorageInvariantViolation,
            diag::Severity::InvalidOutput, "NodePool storage invariants violated");
    if (!cellsOk)
        add(diag::Code::StorageInvariantViolation,
            diag::Severity::InvalidOutput, "CellPool storage invariants violated");
    if (!nodesOk || !cellsOk) return r;
    if (cells.boundNodePool() != 0 &&
        cells.boundNodePool() != nodes.poolId())
        add(diag::Code::CellPoolNodePoolMismatch,
            diag::Severity::InvalidOutput,
            "CellPool is bound to a different NodePool than audited");
    else if (!cells.connectivityIntact(nodes))
        add(diag::Code::DanglingCellConnectivity,
            diag::Severity::InvalidOutput,
            "a live cell references a node no longer current in the NodePool");
    return r;
}

inline bool storageEdgeReady(const NodePool& nodes, const CellPool& cells) {
    return auditStorage(nodes, cells).acceptable(diag::Severity::Warning);
}

// ---------------------------------------------------------------------------
// Ordered-storage content hash (domain+version marked). Requires edge-ready
// storage; excludes pool ids/slots/generations.
// ---------------------------------------------------------------------------
struct StorageContentHash { std::uint64_t value = 0; };

namespace detail {
inline void hashMix(std::uint64_t& h, std::uint64_t v) {
    for (int b = 0; b < 8; ++b) { h ^= (v >> (8 * b)) & 0xFF; h *= 1099511628211ULL; }
}
inline void hashStr(std::uint64_t& h, const char* s) {
    for (const char* p = s; *p; ++p) { h ^= (unsigned char)*p; h *= 1099511628211ULL; }
}
inline void hashBits(std::uint64_t& h, double d) {
    std::uint64_t b = 0;
    static_assert(sizeof(b) == sizeof(d));
    std::memcpy(&b, &d, sizeof(b));
    hashMix(h, b);
}
} // namespace detail

inline StorageContentHash contentHash(const NodePool& nodes,
                                      const CellPool& cells) {
    if (!storageEdgeReady(nodes, cells))
        throw std::invalid_argument(
            "contentHash: storage is not edge-ready");
    std::uint64_t h = 1469598103934665603ULL;
    detail::hashStr(h, "PAACS_STORAGE_HASH_V1");
    const auto nseq = nodes.activeBySequence();
    detail::hashMix(h, nseq.size());
    std::vector<std::pair<std::uint64_t, std::size_t>> rankBySeq;
    rankBySeq.reserve(nseq.size());
    for (std::size_t i = 0; i < nseq.size(); ++i) {
        const auto nd = nodes.get(nseq[i]);
        detail::hashBits(h, nd.coord.x);
        detail::hashBits(h, nd.coord.y);
        detail::hashBits(h, nd.coord.z);
        detail::hashMix(h, static_cast<std::uint64_t>(nd.cls.dim));
        detail::hashMix(h, static_cast<std::uint64_t>(nd.cls.tag));
        rankBySeq.push_back({nodes.seqOf(nseq[i]), i});
    }
    std::sort(rankBySeq.begin(), rankBySeq.end());
    auto rankOf = [&](const NodeHandle& nh) -> std::size_t {
        const std::uint64_t s = nodes.seqOf(nh);
        const auto it = std::lower_bound(
            rankBySeq.begin(), rankBySeq.end(),
            std::pair<std::uint64_t, std::size_t>{s, 0},
            [](const auto& a, const auto& b) { return a.first < b.first; });
        return it->second;
    };
    detail::hashStr(h, "CELLS");
    const auto cseq = cells.activeBySequence();
    detail::hashMix(h, cseq.size());
    for (const CellHandle& ch : cseq) {
        const auto snap = cells.get(ch);
        detail::hashMix(h, static_cast<std::uint64_t>(snap.type));
        detail::hashMix(h, snap.count);
        for (std::uint8_t i = 0; i < snap.count; ++i)
            detail::hashMix(h, rankOf(snap.nodes[i]));
        detail::hashMix(h, static_cast<std::uint64_t>(snap.cls.dim));
        detail::hashMix(h, static_cast<std::uint64_t>(snap.cls.tag));
    }
    return {h};
}

} // namespace paacs::topo
