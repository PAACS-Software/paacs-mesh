// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// paacs::topo - M4 construction topology, stage A: a reusable tetrahedron
// arena for incremental construction algorithms (Delaunay, and later
// constrained recovery and optimization).
//
// Problem it solves: the M3-A Delaunay arena was append-only with alive
// flags, so every cavity replacement leaked its dead cells - ~145-208 MB
// peak RSS at 192k output tets, with the live complex a small fraction of
// that. This arena reuses freed slots via a free list, so capacity tracks
// the LIVE cell count (plus transient cavity churn), not construction
// history.
//
// Facilities:
//   * Slot pool with free list; O(1) alloc/free; slots reused immediately.
//   * Generation counters: every free bumps the slot's generation, so a
//     stale external Handle can never silently alias a reused slot.
//     Internal hot paths use raw slots (documented; the construction
//     algorithm maintains its own consistency and is audited).
//   * Epoch-stamped visit marks: constant-time "seen this insertion?"
//     without clearing between insertions and WITHOUT abusing liveness
//     (the M3-A code repurposed the alive flag as a cavity marker, which
//     made freeing-while-marking impossible).
//   * Creation sequence numbers: deterministic output ordering that is
//     independent of slot reuse.
//   * Memory counters: capacity, live, free, peak-live, byte estimate -
//     the explicit instrumentation the benchmark's RSS column cannot
//     provide portably.
//
// Stage B (separate milestone step): explicit face objects with two-cell
// ownership, node/edge adjacency, and transactional multi-cell edits for
// the mutable public mesh. This arena deliberately stays minimal: same
// meshes as M3-A, materially lower memory.
//
#include <array>
#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "paacs/core/Types.hpp"

namespace paacs::topo {

struct TetHandle {
    // Staleness detection is generation-based. NOTE (M4-B.1): this is the
    // LEGACY construction-arena handle; its 32-bit generation still wraps
    // after 2^32 frees of one slot. The generic topo::Handle allocator
    // (Handles.hpp) instead RETIRES a slot on generation exhaustion. This
    // arena is scheduled to migrate onto the generic allocator at M4-B.9;
    // until then its wrap limitation stands and is harmless at realistic
    // construction scales.
    std::int32_t slot = -1;
    std::uint32_t gen = 0;
    bool valid() const { return slot >= 0; }
};

class TetArena {
public:
    struct Cell {
        std::array<Index, 4> v{};
        std::array<std::int32_t, 4> nb{-1, -1, -1, -1};
        std::uint64_t seq = 0;        // creation sequence (deterministic)
        std::uint32_t gen = 0;        // bumped on every free
        std::uint32_t mark = 0;       // epoch visit stamp
        bool alive = false;
    };

    struct MemoryStats {
        // Slot accounting (review-corrected semantics):
        //   materializedSlots = slots actually constructed (cells_.size());
        //   reservedSlots     = vector capacity, which additionally includes
        //                       reserved-but-never-constructed storage;
        //   freeSlots         = materialized slots currently reusable.
        // materializedSlots == liveCells + freeSlots always holds;
        // reservedSlots - liveCells does NOT equal freeSlots and the two
        // must not be conflated.
        std::size_t materializedSlots = 0;
        std::size_t reservedSlots = 0;
        std::size_t liveCells = 0;
        std::size_t freeSlots = 0;
        std::size_t peakLiveCells = 0;
        std::size_t totalCellsCreated = 0; // == slots an append-only arena
                                           //    (M3-A) would have retained
        std::size_t bytes = 0;             // storage estimate (reserved)
    };

    // ---- allocation ------------------------------------------------------
    std::int32_t alloc(const std::array<Index, 4>& v) {
        std::int32_t s;
        if (!freeList_.empty()) {
            s = freeList_.back();
            freeList_.pop_back();
        } else {
            s = static_cast<std::int32_t>(cells_.size());
            cells_.emplace_back();
        }
        Cell& c = cells_[static_cast<size_t>(s)];
        c.v = v;
        c.nb = {-1, -1, -1, -1};
        c.seq = nextSeq_++;
        c.mark = 0;
        c.alive = true;
        ++live_;
        if (live_ > peakLive_) peakLive_ = live_;
        return s;
    }

    void free(std::int32_t s) {
        assert(s >= 0 && s < slotCount());
        Cell& c = cells_[static_cast<size_t>(s)];
        if (!c.alive)
            throw std::logic_error("TetArena::free: slot already free");
        c.alive = false;
        ++c.gen;                      // invalidate outstanding handles
        --live_;
        freeList_.push_back(s);
    }

    // ---- access (raw slots) -------------------------------------------------
    // CONTRACT: raw-slot access is unchecked in Release and reserved for
    // audited internal construction loops that maintain their own
    // consistency (Delaunay3D). Out-of-range slots are undefined behavior
    // in Release; Debug builds assert. All external consumers must use the
    // generation-protected Handle API below.
    Cell& cell(std::int32_t s) {
        assert(s >= 0 && s < slotCount());
        return cells_[static_cast<size_t>(s)];
    }
    const Cell& cell(std::int32_t s) const {
        assert(s >= 0 && s < slotCount());
        return cells_[static_cast<size_t>(s)];
    }
    bool alive(std::int32_t s) const {
        return s >= 0 && cells_[static_cast<size_t>(s)].alive;
    }
    std::int32_t slotCount() const {
        return static_cast<std::int32_t>(cells_.size());
    }

    // ---- generation-protected external handles ----------------------------
    TetHandle handleOf(std::int32_t s) const {
        assert(s >= 0 && s < slotCount());
        return TetHandle{s, cells_[static_cast<size_t>(s)].gen};
    }
    bool isCurrent(const TetHandle& h) const {
        return h.slot >= 0 && h.slot < slotCount() &&
               cells_[static_cast<size_t>(h.slot)].alive &&
               cells_[static_cast<size_t>(h.slot)].gen == h.gen;
    }
    const Cell& deref(const TetHandle& h) const {
        if (!isCurrent(h))
            throw std::logic_error("TetArena::deref: stale or invalid handle");
        return cells_[static_cast<size_t>(h.slot)];
    }

    // ---- epoch visit marks -------------------------------------------------
    // beginEpoch() invalidates all previous marks in O(1).
    void beginEpoch() {
        if (++epoch_ == 0) {           // wrap: hard-clear once per 2^32
            for (Cell& c : cells_) c.mark = 0;
            epoch_ = 1;
        }
    }
    void markVisited(std::int32_t s) {
        cells_[static_cast<size_t>(s)].mark = epoch_;
    }
    void unmark(std::int32_t s) { cells_[static_cast<size_t>(s)].mark = 0; }
    bool visited(std::int32_t s) const {
        return cells_[static_cast<size_t>(s)].mark == epoch_;
    }

    // ---- instrumentation ----------------------------------------------------
    MemoryStats memory() const {
        MemoryStats m;
        m.materializedSlots = cells_.size();
        m.reservedSlots = cells_.capacity();
        m.liveCells = live_;
        m.freeSlots = freeList_.size();
        m.peakLiveCells = peakLive_;
        m.totalCellsCreated = static_cast<std::size_t>(nextSeq_);
        m.bytes = cells_.capacity() * sizeof(Cell) +
                  freeList_.capacity() * sizeof(std::int32_t);
        return m;
    }

private:
    std::vector<Cell> cells_;
    std::vector<std::int32_t> freeList_;
    std::size_t live_ = 0, peakLive_ = 0;
    std::uint64_t nextSeq_ = 0;
    std::uint32_t epoch_ = 0;
};

} // namespace paacs::topo
