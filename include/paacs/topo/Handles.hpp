// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// paacs::topo - M4-B.1: generation-protected generic entity handles.
//
// SCOPE (deliberately narrow, per the issue plan): the handle types and a
// payload-free allocator that manages slot lifetimes, generations, and
// deterministic creation sequence. Generic node/edge/face/cell STORAGE
// composes these at M4-B.2; TetArena keeps its own specialized handle
// until its migration (M4-B.9) and is unchanged by this issue.
//
// Contracts:
//   * A default-constructed handle is INVALID.
//   * Handles of different entity kinds are DISTINCT TYPES with no
//     implicit conversion between each other or to raw indices - passing
//     a NodeHandle where a CellHandle is expected is a compile error
//     (statically asserted in the test suite).
//   * Freeing a slot bumps its generation: every previously issued
//     handle to that slot becomes stale; stale isCurrent() is false and
//     stale deref throws std::logic_error deterministically.
//   * GENERATION-EXHAUSTION POLICY: a slot's 32-bit generation does NOT
//     wrap. On the 2^32-th free of one slot the slot is RETIRED (removed
//     from reuse permanently) so no ancient handle can ever re-validate.
//     Retirement is counted in the three-state accounting invariant.
//   * Equality/ordering/hashing are (pool, slot, gen)-based and stable:
//     they never change for a given handle value, and handles from
//     different allocators are always distinct (pool participates).
//   * Deterministic order: allocations carry a strictly increasing
//     sequence number independent of slot reuse; public iteration in
//     later storage layers orders by sequence, never by slot.
//   * Unchecked raw-slot access lives ONLY in topo::raw and is forbidden
//     in public APIs, transactions, and audits (directive 7.5).
//
#include <atomic>
#include <cstdint>
#include <limits>
#include <functional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace paacs::topo {

class TopologyTransaction;                    // M4-B.6 coordinator (fwd)
class EdgePool;                               // reviewed adoption friends
class FacePool;
namespace detail { class TransactionAccess; } // M4-B.6 prepared commit

template <class Tag>
struct Handle {
    std::int32_t slot = -1;
    std::uint32_t gen = 0;
    // Pool identity: same-type handles from DIFFERENT allocators must not
    // validate against each other. Each allocator stamps its handles with
    // its own 64-bit id; isCurrent() rejects foreign handles. 0 =
    // default/null and moved-from (matches no live allocator; ids start
    // at 1). RUNTIME-ONLY OWNERSHIP TOKEN (review): pool ids depend on
    // construction order and MUST be excluded from serialization,
    // semantic hashes, classifications, stable external ids, and
    // deterministic public ordering (enforced at M4-B.2).
    std::uint64_t pool = 0;

    // "Has a slot" - i.e. not the default/null handle. This is NOT a
    // usability check: a freed handle is still specified. Use
    // HandleAllocator::isCurrent(h) to test whether a handle is USABLE
    // (review finding: valid() invited exactly that confusion).
    constexpr bool isSpecified() const { return slot >= 0; }

    friend constexpr bool operator==(const Handle& a, const Handle& b) {
        return a.slot == b.slot && a.gen == b.gen && a.pool == b.pool;
    }
    friend constexpr bool operator!=(const Handle& a, const Handle& b) {
        return !(a == b);
    }
    friend constexpr bool operator<(const Handle& a, const Handle& b) {
        if (a.pool != b.pool) return a.pool < b.pool;
        if (a.slot != b.slot) return a.slot < b.slot;
        return a.gen < b.gen;
    }
};

struct NodeTag {};
struct EdgeTag {};
struct FaceTag {};
struct CellTag {};

using NodeHandle = Handle<NodeTag>;
using EdgeHandle = Handle<EdgeTag>;
using FaceHandle = Handle<FaceTag>;
using CellHandle = Handle<CellTag>;

// Payload-free slot/generation/sequence manager. Storage layers (M4-B.2)
// hold their payload arrays parallel to the allocator's slots.
//
// STATE MACHINE (round-6 review: the previous version traded one broken
// invariant for several subtler ones, hidden because tests only checked
// ordinary lifecycle - not allocator state AFTER a throw). Each slot is
// exactly one of: live, reusable-free, or retired (generation-exhausted,
// never reissued). The accounting invariant is therefore:
//     materializedSlots == live + reusableFree + retired
// Every mutating operation is either fully applied or leaves this
// invariant and all handle validity unchanged (strong exception safety).
template <class Tag>
class HandleAllocator {
public:
    using H = Handle<Tag>;

    HandleAllocator() : pool_(nextPoolId()) {}
    HandleAllocator(const HandleAllocator&) = delete;
    HandleAllocator& operator=(const HandleAllocator&) = delete;
    HandleAllocator(HandleAllocator&& o) noexcept
        : HandleAllocator(TokenlessAdoptionTag{}) {
        adoptFrom(o);
    }
    // Move ASSIGNMENT is DELETED (round-6 review): assigning onto a
    // non-empty allocator would silently invalidate all its handles and,
    // once payload storage exists, detach lifetime metadata from data.
    // Allocator + payload arrays must move together as a unit (the M4-B.2
    // storage type owns both and is itself move-constructible). Move
    // construction is retained.
    HandleAllocator& operator=(HandleAllocator&&) = delete;

    H alloc() {
        // A moved-from allocator (pool 0) mints nothing (review: alloc()
        // previously produced a specified-but-never-current handle).
        if (pool_ == 0)
            throw std::logic_error(
                "HandleAllocator::alloc on a moved-from allocator");
        // Exhaustion is checked BEFORE any state mutation (review: the
        // sequence check ran after popping/appending, losing a slot on
        // throw). Both checks are pure reads.
        if (nextSeq_ == std::numeric_limits<std::uint64_t>::max())
            throw std::length_error(
                "HandleAllocator: creation sequence exhausted (UINT64_MAX)");
        if (freeList_.empty() && atSlotCapacity(slots_.size()))
            throw std::length_error(
                "HandleAllocator: slot space exhausted (INT32_MAX)");
        std::int32_t s;
        if (!freeList_.empty()) {
            s = freeList_.back();
            freeList_.pop_back();
        } else {
            s = static_cast<std::int32_t>(slots_.size());
            slots_.push_back(Slot{});   // may throw: no state changed yet
        }
        Slot& sl = slots_[static_cast<std::size_t>(s)];
        sl.alive = true;
        sl.seq = nextSeq_++;
        ++live_;
        return H{s, sl.gen, pool_};
    }

    void free(const H& h) {
        requireCurrent(h, "HandleAllocator::free");
        Slot& sl = slots_[static_cast<std::size_t>(h.slot)];
        if (sl.gen == std::numeric_limits<std::uint32_t>::max()) {
            // Generation exhausted: RETIRE in place (review: the previous
            // path pushed to a retiredSlots_ vector AFTER mutating the
            // slot, so a throw there corrupted state). Retirement is now a
            // Slot flag - noexcept, O(1), no secondary allocation.
            sl.alive = false;
            sl.retired = true;
            ++retired_;
            --live_;
            return;
        }
        // Normal path: reserve the free-list entry FIRST (only throwing
        // step), then apply noexcept slot mutations.
        freeList_.push_back(h.slot);
        sl.alive = false;
        ++sl.gen;                       // invalidate outstanding handles
        --live_;
    }

    bool isCurrent(const H& h) const {
        return h.pool == pool_ && pool_ != 0 && h.slot >= 0 &&
               h.slot < static_cast<std::int32_t>(slots_.size()) &&
               slots_[static_cast<std::size_t>(h.slot)].alive &&
               slots_[static_cast<std::size_t>(h.slot)].gen == h.gen;
    }
    std::uint64_t poolId() const { return pool_; }

    std::int32_t deref(const H& h) const {
        requireCurrent(h, "HandleAllocator::deref");
        return h.slot;
    }
    std::uint64_t seqOf(const H& h) const {
        requireCurrent(h, "HandleAllocator::seqOf");
        return slots_[static_cast<std::size_t>(h.slot)].seq;
    }

    std::size_t liveCount() const { return live_; }
    std::size_t materializedSlots() const { return slots_.size(); }
    std::size_t reusableFreeSlots() const { return freeList_.size(); }
    std::size_t retiredSlots() const { return retired_; }
    std::uint64_t totalCreated() const { return nextSeq_; }
    // Pure capacity HINT: ensure the free list can hold at least `n` more
    // entries without reallocating, so a subsequent free() is guaranteed
    // not to allocate (used by composed pools to make transactional
    // rollback noexcept). This changes NO observable behavior, sequence,
    // hash, or accounting - only the free list's capacity. Frozen contract
    // preserved.
    void reserveFreeCapacity(std::size_t n) {
        freeList_.reserve(freeList_.size() + n);
    }
    // Absolute-target slot-capacity HINT (M4-B.6 prepared commit,
    // reviewer-sanctioned additive frozen-file edit; absolute form
    // preferred over an additive delta to avoid size+delta overflow):
    // ensure the slot table can hold at least `requiredSlots`
    // materialized slots without reallocating, so a later fresh alloc()
    // cannot throw from slots_.push_back. Capacity ONLY: does not
    // materialize slots, advance sequences, change generations, change
    // live/free/retired counts, or change any hash-visible state.
    void reserveMaterializedCapacity(std::size_t requiredSlots) {
        if (requiredSlots <= slots_.capacity()) return;   // already held
        // Reviewer checklist: INT32_MAX slot limit, then max_size; the
        // absolute-target form has no size+delta addition to overflow.
        if (requiredSlots >
            static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
            throw std::length_error(
                "HandleAllocator::reserveMaterializedCapacity: beyond the "
                "INT32_MAX slot limit");
        if (requiredSlots > slots_.max_size())
            throw std::length_error(
                "HandleAllocator::reserveMaterializedCapacity: beyond "
                "vector max_size");
        slots_.reserve(requiredSlots);                    // capacity only
    }
    // Approximate heap footprint of the allocator's own arrays (slot
    // table + free list). Lets storage report TOTAL pool memory, not just
    // payload (review item 6).
    std::size_t approxBytes() const {
        return slots_.capacity() * sizeof(Slot) +
               freeList_.capacity() * sizeof(std::int32_t);
    }
    // Pure, constexpr slot-capacity predicate (public: it leaks no state
    // and lets callers and tests reason about the INT32_MAX boundary).
    static constexpr bool wouldExhaustSlots(std::size_t materialized) {
        return atSlotCapacity(materialized);
    }
    // The full accounting invariant, exposed for tests and M4-B.2 counters.
    bool accountingConsistent() const {
        return materializedSlots() ==
               liveCount() + reusableFreeSlots() + retiredSlots();
    }

#if defined(PAACS_TOPO_TESTING)
    // Test-only seam (round-6 review: this must NOT be in the production
    // API). Compiled only under PAACS_TOPO_TESTING. Forces a live slot's
    // generation so retirement can be exercised without 2^32 real frees.
    void testForceGeneration(const H& h, std::uint32_t g) {
        requireCurrent(h, "testForceGeneration");
        slots_[static_cast<std::size_t>(h.slot)].gen = g;
    }
    bool testSlotRetired(std::int32_t s) const {
        return s >= 0 && s < static_cast<std::int32_t>(slots_.size()) &&
               slots_[static_cast<std::size_t>(s)].retired;
    }
    // Force the creation sequence to its max so the next alloc() must hit
    // the pre-mutation exhaustion check.
    void testForceSequenceMax() {
        nextSeq_ = std::numeric_limits<std::uint64_t>::max();
    }
#endif

private:
    friend class EdgePool;
    friend class FacePool;

    // Reviewed Increment-5 tokenless adoption target. Unlike the public
    // constructor, this state does not consume a runtime pool identity.
    // Access is limited to allocator self-move and the authorized composed
    // EdgePool/FacePool move paths.
    struct TokenlessAdoptionTag {};
    explicit HandleAllocator(TokenlessAdoptionTag) noexcept {}

    // M4-B.6 prepared-commit access (reviewer-sanctioned narrow friend;
    // TransactionAccess's primitives are PRIVATE and callable only by
    // TopologyTransaction).
    friend class detail::TransactionAccess;

    struct Slot {
        std::uint32_t gen = 0;
        std::uint64_t seq = 0;
        bool alive = false;
        bool retired = false;
    };

    // Pure capacity decision (review finding: extracted so the INT32_MAX
    // boundary is testable without materializing billions of slots).
    static constexpr bool atSlotCapacity(std::size_t materialized) {
        return materialized >= static_cast<std::size_t>(
                                   std::numeric_limits<std::int32_t>::max());
    }

    void requireCurrent(const H& h, const char* who) const {
        if (!isCurrent(h))
            throw std::logic_error(
                std::string(who) + ": stale, freed, or invalid handle");
    }
    void adoptFrom(HandleAllocator& o) noexcept {
        // The target is private, tokenless and empty. Exact container swaps
        // preserve slot generations, creation sequence, reusable free-list
        // order and retired state. Scalar swaps transfer the complete
        // accounting and pool token while leaving `o` in that target state.
        static_assert(noexcept(slots_.swap(o.slots_)));
        static_assert(noexcept(freeList_.swap(o.freeList_)));
        slots_.swap(o.slots_);
        freeList_.swap(o.freeList_);
        using std::swap;
        swap(live_, o.live_);
        swap(retired_, o.retired_);
        swap(nextSeq_, o.nextSeq_);
        swap(pool_, o.pool_);
    }
    static std::uint64_t nextPoolId() {
        // 64-bit ownership token; a wrap would let a handle from a
        // destroyed allocator validate against a later one - the same
        // resurrection class retirement prevents per-slot. Saturating:
        // after the final id is issued the counter stays at UINT64_MAX
        // PERMANENTLY and every later request throws, rather than a plain
        // ++counter that would wrap to 0 and reissue id 1. At 1e9
        // allocators/sec this ceiling is ~585 years away. The transition
        // logic lives in the pure helper poolIdStep() so it is testable
        // without 2^64 real constructions (review finding: the saturating
        // path was previously asserted-tested only in prose).
        static std::atomic<std::uint64_t> counter{0};
        std::uint64_t cur = counter.load(std::memory_order_relaxed);
        for (;;) {
            const std::uint64_t next = poolIdStep(cur);   // throws at max
            if (counter.compare_exchange_weak(cur, next,
                                              std::memory_order_relaxed))
                return next;
        }
    }

public:
    // Pure pool-id transition (public, stateless, no allocator state
    // touched): given the current counter value, returns the next id or
    // throws length_error at saturation. ids start at 1; 0 is reserved
    // for null/moved-from. Extracted for direct exhaustion testing.
    static constexpr std::uint64_t poolIdStep(std::uint64_t cur) {
        if (cur == std::numeric_limits<std::uint64_t>::max())
            throw std::length_error(
                "HandleAllocator: pool-id space permanently exhausted");
        return cur + 1;
    }

private:

    std::vector<Slot> slots_;
    std::vector<std::int32_t> freeList_;
    std::size_t live_ = 0;
    std::size_t retired_ = 0;
    std::uint64_t nextSeq_ = 0;
    std::uint64_t pool_ = 0;

    static_assert(
        std::is_nothrow_default_constructible_v<decltype(slots_)>);
    static_assert(
        std::is_nothrow_default_constructible_v<decltype(freeList_)>);
    static_assert(std::is_nothrow_swappable_v<decltype(slots_)>);
    static_assert(std::is_nothrow_swappable_v<decltype(freeList_)>);
};

// Compile-time guards for the 64-bit ownership-token contract (review
// finding: poolId() had silently truncated to uint32). Placed after the
// template so both the member and the method signature are complete.
static_assert(
    std::is_same_v<std::uint64_t,
                   decltype(std::declval<const HandleAllocator<NodeTag>>()
                                .poolId())>,
    "HandleAllocator::poolId() must return the full 64-bit token");
static_assert(
    std::is_same_v<std::uint64_t,
                   decltype(std::declval<Handle<NodeTag>>().pool)>,
    "Handle::pool must be the full 64-bit token");

// Unchecked raw-slot access: INTERNAL construction hot paths only
// (directive 7.5). Never in public APIs, transactions, or audits.
namespace raw {
template <class Tag>
inline std::int32_t slotUnchecked(const Handle<Tag>& h) { return h.slot; }
} // namespace raw

} // namespace paacs::topo

// Stable hashing for unordered containers (explicit namespace std block
// for maximum compiler portability - review suggestion).
namespace std {
template <class Tag>
struct hash<paacs::topo::Handle<Tag>> {
    std::size_t operator()(const paacs::topo::Handle<Tag>& h) const {
        auto mix = [](std::uint64_t x) {
            x += 0x9e3779b97f4a7c15ULL;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            return x ^ (x >> 31);
        };
        const std::uint64_t a =
            (std::uint64_t(std::uint32_t(h.slot)) << 32) | h.gen;
        return static_cast<std::size_t>(mix(a) ^ mix(h.pool + 1));  // 64-bit pool
    }
};
} // namespace std
