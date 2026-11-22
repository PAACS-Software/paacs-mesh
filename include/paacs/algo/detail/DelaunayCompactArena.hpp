// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/core/Types.hpp"
#include "paacs/algo/detail/TetConstructionBudget.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#if defined(_MSC_VER)
#  define PAACS_DELAUNAY_ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#  define PAACS_DELAUNAY_ALWAYS_INLINE inline __attribute__((always_inline))
#else
#  define PAACS_DELAUNAY_ALWAYS_INLINE inline
#endif

namespace paacs::algo::detail {

// Internal exception used to restart safely through the 64-bit legacy arena.
// It never escapes the public Delaunay API.
class CompactEnvelopeExceeded final : public std::runtime_error {
public:
    explicit CompactEnvelopeExceeded(const char* message)
        : std::runtime_error(message) {}
};

// Delaunay-private construction arena.  The 32-byte hot record contains only
// the data touched by point-location, cavity traversal and rewiring.  Liveness,
// epoch marks and deterministic live-order links are held in cold side arrays.
//
// A packed neighbour stores (slot << 2) | reciprocalFace.  The all-ones value
// is reserved for the hull/no-neighbour sentinel, so slot 0x3fffffff is outside
// the compact envelope.  The caller catches CompactEnvelopeExceeded and
// restarts with the unchanged public-Index legacy implementation.
class CompactTetArena {
public:
    using Vertex = std::uint32_t;
    using Slot = std::uint32_t;
    static constexpr Slot kInvalidSlot = std::numeric_limits<Slot>::max();
    static constexpr std::uint32_t kNoNeighbour =
        std::numeric_limits<std::uint32_t>::max();
    static constexpr Slot kMaximumPackedSlot = 0x3ffffffeU;

    struct HotTet {
        std::array<Vertex, 4> vertices{};
        std::array<std::uint32_t, 4> neighbours{
            kNoNeighbour, kNoNeighbour, kNoNeighbour, kNoNeighbour};
    };
    static_assert(sizeof(HotTet) == 32,
                  "compact Delaunay hot tetrahedron must remain 32 bytes");

    struct MemoryStats {
        std::uint64_t materializedSlots = 0;
        std::uint64_t reservedSlots = 0;
        std::uint64_t liveCells = 0;
        std::uint64_t freeSlots = 0;
        std::uint64_t peakLiveCells = 0;
        std::uint64_t totalCellsCreated = 0;
        std::uint64_t hotBytes = 0;
        std::uint64_t coldBytes = 0;
        std::uint64_t bytes = 0;
    };

    CompactTetArena() = default;
    CompactTetArena(const CompactTetArena&) = delete;
    CompactTetArena& operator=(const CompactTetArena&) = delete;
    CompactTetArena(CompactTetArena&&) noexcept = default;
    CompactTetArena& operator=(CompactTetArena&&) noexcept = default;

    void setConstructionBudget(TetConstructionBudget* budget, Vertex infinity) noexcept {
        budget_=budget; infinity_=infinity;
        hot_.setConstructionBudget(budget); state_.setConstructionBudget(budget);
        links_.setConstructionBudget(budget); freeList_.setConstructionBudget(budget);
        versions_.setConstructionBudget(budget);
    }

    void reserveSlots(std::size_t slots) {
        if (slots > static_cast<std::size_t>(kMaximumPackedSlot) + 1U)
            throw CompactEnvelopeExceeded(
                "compact Delaunay reserve envelope exceeded");
        hot_.reserve(slots);
        state_.reserve(slots);
        links_.reserve(slots);
        if (versionsEnabled_) versions_.reserve(slots);
    }

    void reset() noexcept {
        // Retain materialized storage.  The next job overwrites slots from
        // zero upward, avoiding both allocation and value-initialization of
        // every HotTet.
        activeSlots_ = 0;
        finiteLive_ = 0;
        freeList_.clear();
        liveHead_ = kInvalidSlot;
        liveTail_ = kInvalidSlot;
        liveCount_ = 0;
        peakLive_ = 0;
        totalCreated_ = 0;
        epoch_ = 0;
        if (versionsEnabled_) {
            if (versions_.size() < hot_.size()) versions_.resize(hot_.size(), 0U);
            std::fill(versions_.begin(), versions_.end(), 0U);
        }
    }

    void release() noexcept {
        CompactTetArena empty;
        *this = std::move(empty);
    }

    PAACS_DELAUNAY_ALWAYS_INLINE Slot allocate(
        const std::array<Vertex, 4>& vertices) {
        if (budget_) {
            budget_->step();
            if (vertices[3] != infinity_) budget_->cells(finiteLive_+1);
        }
        Slot slot = kInvalidSlot;
        if (!freeList_.empty()) {
            slot = freeList_.back();
            freeList_.pop_back();
        } else if (activeSlots_ < hot_.size()) {
            slot = activeSlots_++;
        } else {
            if (hot_.size() > static_cast<std::size_t>(kMaximumPackedSlot))
                throw CompactEnvelopeExceeded(
                    "compact Delaunay tetrahedron-slot envelope exceeded");
            slot = static_cast<Slot>(hot_.size());
            hot_.emplace_back();
            state_.push_back(0);
            links_.push_back({kInvalidSlot, kInvalidSlot});
            if (versionsEnabled_) versions_.push_back(0U);
            ++activeSlots_;
        }

        HotTet& tet = hot_[slot];
        tet.vertices = vertices;
        tet.neighbours = {kNoNeighbour, kNoNeighbour,
                          kNoNeighbour, kNoNeighbour};
        state_[slot] = kAliveBit;
        links_[slot].previous = liveTail_;
        links_[slot].next = kInvalidSlot;
        if (liveTail_ != kInvalidSlot) links_[liveTail_].next = slot;
        else liveHead_ = slot;
        liveTail_ = slot;

        bumpVersion(slot);
        ++liveCount_;
        if (budget_ && vertices[3] != infinity_) ++finiteLive_;
        if (liveCount_ > peakLive_) peakLive_ = liveCount_;
        ++totalCreated_;
        return slot;
    }

    PAACS_DELAUNAY_ALWAYS_INLINE void free(Slot slot) {
        if (slot >= hot_.size() || !alive(slot))
            throw std::logic_error(
                "CompactTetArena::free: invalid or already-free slot");

        if (budget_) {
            budget_->step();
            if (hot_[slot].vertices[3] != infinity_) --finiteLive_;
        }
        const Slot previous = links_[slot].previous;
        const Slot next = links_[slot].next;
        if (previous != kInvalidSlot) links_[previous].next = next;
        else liveHead_ = next;
        if (next != kInvalidSlot) links_[next].previous = previous;
        else liveTail_ = previous;

        links_[slot] = {kInvalidSlot, kInvalidSlot};
        state_[slot] = 0;
        bumpVersion(slot);
        --liveCount_;
        freeList_.push_back(slot);
    }

    [[nodiscard]] HotTet& hot(Slot slot) noexcept { return hot_[slot]; }
    [[nodiscard]] const HotTet& hot(Slot slot) const noexcept {
        return hot_[slot];
    }

    [[nodiscard]] bool alive(Slot slot) const noexcept {
        return slot < state_.size() && (state_[slot] & kAliveBit) != 0;
    }
    [[nodiscard]] Slot slotCount() const noexcept { return activeSlots_; }
    [[nodiscard]] Slot liveHead() const noexcept { return liveHead_; }
    [[nodiscard]] Slot liveTail() const noexcept { return liveTail_; }
    [[nodiscard]] Slot nextLive(Slot slot) const noexcept {
        return links_[slot].next;
    }

    static constexpr std::uint32_t packNeighbour(Slot slot,
                                                  std::uint32_t face) {
        return (slot << 2U) | (face & 3U);
    }
    static constexpr Slot neighbourSlot(std::uint32_t packed) noexcept {
        return packed == kNoNeighbour ? kInvalidSlot : packed >> 2U;
    }
    static constexpr std::uint32_t neighbourFace(
        std::uint32_t packed) noexcept {
        return packed & 3U;
    }

    [[nodiscard]] std::uint32_t neighbourRef(Slot slot,
                                              std::uint32_t face) const noexcept {
        return hot_[slot].neighbours[face];
    }
    [[nodiscard]] Slot neighbour(Slot slot,
                                 std::uint32_t face) const noexcept {
        return neighbourSlot(neighbourRef(slot, face));
    }
    [[nodiscard]] std::uint32_t reciprocalFace(
        Slot slot, std::uint32_t face) const noexcept {
        return neighbourFace(neighbourRef(slot, face));
    }

    void setBoundary(Slot slot, std::uint32_t face) noexcept {
        if (hot_[slot].neighbours[face] != kNoNeighbour) {
            hot_[slot].neighbours[face] = kNoNeighbour;
            bumpVersion(slot);
        }
    }
    void setNeighbour(Slot slot, std::uint32_t face, Slot neighbour,
                      std::uint32_t reciprocalFace) {
        if (slot > kMaximumPackedSlot || neighbour > kMaximumPackedSlot ||
            face > 3U || reciprocalFace > 3U)
            throw CompactEnvelopeExceeded(
                "compact Delaunay packed-neighbour envelope exceeded");
        const std::uint32_t packed = packNeighbour(neighbour, reciprocalFace);
        if (hot_[slot].neighbours[face] != packed) {
            hot_[slot].neighbours[face] = packed;
            bumpVersion(slot);
        }
    }
    void connect(Slot a, std::uint32_t faceA, Slot b,
                 std::uint32_t faceB) {
        setNeighbour(a, faceA, b, faceB);
        setNeighbour(b, faceB, a, faceA);
    }

    void beginEpoch() {
        ++epoch_;
        if (epoch_ == 0 || epoch_ > kEpochMask) {
            // Rare epoch wrap clears retained storage as one admitted chunk.
            if (budget_) budget_->step(state_.size());
            for (std::uint32_t& state : state_) state &= kAliveBit;
            epoch_ = 1;
        }
    }
    void markVisited(Slot slot) noexcept {
        state_[slot] = kAliveBit | epoch_;
    }
    void unmark(Slot slot) noexcept { state_[slot] = kAliveBit; }
    [[nodiscard]] bool visited(Slot slot) const noexcept {
        return (state_[slot] & kEpochMask) == epoch_;
    }


    // V0.9.3 speculative-discovery support. Version storage is enabled only
    // for deterministic-parallel jobs, so the frozen serial compact path
    // retains its V0.9.2 memory footprint. Versions cover all topology words
    // that can affect point location, conflict traversal or shell ownership.
    void enableVersions() {
        versionsEnabled_ = true;
        versions_.assign(hot_.size(), 0U);
        versions_.reserve(hot_.capacity());
    }
    void disableVersions() noexcept {
        versionsEnabled_ = false;
        versions_.releaseStorage();
    }
    [[nodiscard]] bool versionsEnabled() const noexcept { return versionsEnabled_; }
    [[nodiscard]] std::uint32_t version(Slot slot) const noexcept {
        return versionsEnabled_ && slot < versions_.size() ? versions_[slot] : 0U;
    }

    [[nodiscard]] MemoryStats memory() const noexcept {
        MemoryStats result;
        result.materializedSlots = activeSlots_;
        result.reservedSlots = hot_.capacity();
        result.liveCells = liveCount_;
        result.freeSlots = freeList_.size();
        result.peakLiveCells = peakLive_;
        result.totalCellsCreated = totalCreated_;
        result.hotBytes = hot_.capacity() * sizeof(HotTet);
        result.coldBytes = state_.capacity() * sizeof(std::uint32_t) +
                           links_.capacity() * sizeof(LiveLinks) +
                           freeList_.capacity() * sizeof(Slot) +
                           versions_.capacity() * sizeof(std::uint32_t);
        result.bytes = result.hotBytes + result.coldBytes;
        return result;
    }

private:

    void bumpVersion(Slot slot) noexcept {
        if (!versionsEnabled_) return;
        if (slot >= versions_.size()) versions_.resize(static_cast<std::size_t>(slot) + 1U, 0U);
        std::uint32_t next = versions_[slot] + 1U;
        if (next == 0U) next = 1U;
        versions_[slot] = next;
    }

    static constexpr std::uint32_t kAliveBit = 0x80000000U;
    static constexpr std::uint32_t kEpochMask = 0x7fffffffU;

    ConstructionVector<HotTet> hot_;
    struct LiveLinks {
        Slot previous = kInvalidSlot;
        Slot next = kInvalidSlot;
    };
    static_assert(sizeof(LiveLinks) == 8);

    ConstructionVector<std::uint32_t> state_;
    ConstructionVector<LiveLinks> links_;
    ConstructionVector<Slot> freeList_;
    ConstructionVector<std::uint32_t> versions_;
    TetConstructionBudget* budget_ = nullptr;
    Vertex infinity_ = 0;
    std::uint64_t finiteLive_ = 0;
    bool versionsEnabled_ = false;
    Slot activeSlots_ = 0;
    Slot liveHead_ = kInvalidSlot;
    Slot liveTail_ = kInvalidSlot;
    std::uint64_t liveCount_ = 0;
    std::uint64_t peakLive_ = 0;
    std::uint64_t totalCreated_ = 0;
    std::uint32_t epoch_ = 0;
};

} // namespace paacs::algo::detail

#undef PAACS_DELAUNAY_ALWAYS_INLINE
