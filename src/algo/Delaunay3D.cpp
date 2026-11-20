// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Delaunay3D.hpp"
#include "paacs/algo/detail/DelaunayCompactArena.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "TetCandidateDetail.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <barrier>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef PAACS_DELAUNAY_USE_RADIX
#define PAACS_DELAUNAY_USE_RADIX 1
#endif
#ifndef PAACS_DELAUNAY_DIRECT_COMMIT
#define PAACS_DELAUNAY_DIRECT_COMMIT 1
#endif
#ifndef PAACS_DELAUNAY_LINEAR_FINALIZATION
#define PAACS_DELAUNAY_LINEAR_FINALIZATION 1
#endif

namespace paacs::predicates::detail {
void resetDelaunayZeroSeen() noexcept;
bool delaunayZeroSeen() noexcept;
}

namespace paacs::algo {
namespace detail {
std::vector<Tet> delaunay3DLegacy64(const std::vector<Vec3>&,
                                    Delaunay3DStats*);
}
namespace {

using Clock = std::chrono::steady_clock;
using predicates::Sign;
using Arena = detail::CompactTetArena;
using Slot = Arena::Slot;
using Vertex = Arena::Vertex;

constexpr int kFace[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};

inline double secondsBetween(const Clock::time_point begin,
                             const Clock::time_point end) noexcept {
    return std::chrono::duration<double>(end - begin).count();
}

template <bool Enabled>
inline void detailCount(std::uint64_t& value,
                        std::uint64_t amount = 1) noexcept {
    if constexpr (Enabled) value += amount;
}

template <bool Enabled>
inline void detailMaximum(std::uint64_t& value,
                          std::uint64_t candidate) noexcept {
    if constexpr (Enabled) value = std::max(value, candidate);
}

inline double orientHeur(const Vec3& a, const Vec3& b, const Vec3& c,
                         const Vec3& d) noexcept {
    return (b - a).cross(c - a).dot(d - a);
}

std::uint64_t morton3(std::uint32_t x, std::uint32_t y, std::uint32_t z) {
    auto part = [](std::uint64_t value) {
        value &= 0x1FFFFFU;
        value = (value | value << 32U) & 0x1F00000000FFFFULL;
        value = (value | value << 16U) & 0x1F0000FF0000FFULL;
        value = (value | value << 8U) & 0x100F00F00F00F00FULL;
        value = (value | value << 4U) & 0x10C30C30C30C30C3ULL;
        value = (value | value << 2U) & 0x1249249249249249ULL;
        return value;
    };
    return part(x) | (part(y) << 1U) | (part(z) << 2U);
}

inline void sort2(Vertex& a, Vertex& b) noexcept {
    if (b < a) std::swap(a, b);
}

[[maybe_unused]] inline std::array<Vertex, 4> canonical4(
    const std::array<Vertex, 4>& input) noexcept {
    auto key = input;
    sort2(key[0], key[1]);
    sort2(key[2], key[3]);
    sort2(key[0], key[2]);
    sort2(key[1], key[3]);
    sort2(key[1], key[2]);
    return key;
}


inline bool sameSet3(const std::array<Vertex, 3>& a,
                     const std::array<Vertex, 3>& b) noexcept {
    return (a[0] == b[0] || a[0] == b[1] || a[0] == b[2]) &&
           (a[1] == b[0] || a[1] == b[1] || a[1] == b[2]) &&
           (a[2] == b[0] || a[2] == b[1] || a[2] == b[2]);
}

inline bool sameSet4(const std::array<Vertex, 4>& a,
                     const std::array<Vertex, 4>& b) noexcept {
    for (Vertex value : a) {
        if (value != b[0] && value != b[1] && value != b[2] &&
            value != b[3])
            return false;
    }
    return true;
}

struct PendingSide {
    Slot newTet = Arena::kInvalidSlot;
    std::uint8_t side = 0;
    Vertex e0 = 0;
    Vertex e1 = 0;
};

class ShellEdgePairTable {
public:
    enum class Result : std::uint8_t { First, Second, Third };

    void begin(std::size_t expected) {
        std::size_t capacity = 8;
        while (capacity < expected * 2U) capacity <<= 1U;
        if (slots_.size() < capacity) slots_.resize(capacity);
        mask_ = slots_.size() - 1U;
        unmatched_ = 0;
        if (++epoch_ == 0 || epoch_ > kEpochMask) {
            for (auto& slot : slots_) slot.stamp = 0;
            epoch_ = 1;
        }
    }

    template <bool CollectDetailed>
    [[nodiscard]] Result insert(const PendingSide& value,
                                PendingSide& first,
                                std::uint64_t& probes) noexcept {
        const std::uint64_t edgeKey =
            (static_cast<std::uint64_t>(value.e0) << 32U) |
            static_cast<std::uint64_t>(value.e1);
        std::uint64_t hashed = edgeKey;
        hashed ^= hashed >> 33U;
        hashed *= 0xff51afd7ed558ccdULL;
        hashed ^= hashed >> 33U;
        std::size_t position = static_cast<std::size_t>(hashed) & mask_;
        for (;;) {
            detailCount<CollectDetailed>(probes);
            SlotRecord& slot = slots_[position];
            const std::uint32_t slotEpoch = slot.stamp & kEpochMask;
            if (slotEpoch != epoch_) {
                slot.edgeKey = edgeKey;
                slot.packedTetSide =
                    Arena::packNeighbour(value.newTet, value.side);
                slot.stamp = epoch_;
                ++unmatched_;
                return Result::First;
            }
            if (slot.edgeKey == edgeKey) {
                first.newTet = Arena::neighbourSlot(slot.packedTetSide);
                first.side = static_cast<std::uint8_t>(
                    Arena::neighbourFace(slot.packedTetSide));
                first.e0 = value.e0;
                first.e1 = value.e1;
                if ((slot.stamp & kMatchedBit) != 0U)
                    return Result::Third;
                slot.stamp |= kMatchedBit;
                --unmatched_;
                return Result::Second;
            }
            position = (position + 1U) & mask_;
        }
    }

    [[nodiscard]] std::size_t unmatched() const noexcept { return unmatched_; }
    [[nodiscard]] std::size_t capacity() const noexcept {
        return slots_.capacity();
    }
    [[nodiscard]] std::uint64_t retainedBytes() const noexcept {
        return slots_.capacity() * sizeof(SlotRecord);
    }

private:
    static constexpr std::uint32_t kMatchedBit = 0x80000000U;
    static constexpr std::uint32_t kEpochMask = 0x7fffffffU;
    struct SlotRecord {
        std::uint64_t edgeKey = 0;
        std::uint32_t packedTetSide = Arena::kNoNeighbour;
        std::uint32_t stamp = 0;
    };
    static_assert(sizeof(SlotRecord) == 16);

    std::vector<SlotRecord> slots_;
    std::uint32_t epoch_ = 0;
    std::size_t mask_ = 0;
    std::size_t unmatched_ = 0;
};

struct ShellFace {
    Vertex f0 = 0;
    Vertex f1 = 0;
    Vertex f2 = 0;
    std::uint32_t outerRef = Arena::kNoNeighbour;
    bool directFiniteWinding = false;
};

class UniqueTetTable {
public:
    void begin(std::size_t expected) {
        std::size_t capacity = 8;
        const std::size_t target = expected + expected / 2U + 1U;
        while (capacity < target) capacity <<= 1U;
        if (slots_.size() < capacity) slots_.resize(capacity);
        mask_ = slots_.size() - 1U;
        if (++epoch_ == 0) {
            for (auto& slot : slots_) slot.epoch = 0;
            epoch_ = 1;
        }
    }

    bool insert(const std::array<Vertex, 4>& vertices,
                const std::vector<Tet>& output) noexcept {
        const std::uint64_t fullHash = unorderedHash(vertices);
        const std::uint16_t fingerprint = static_cast<std::uint16_t>(
            (fullHash >> 32U) ^ fullHash);
        std::size_t position = static_cast<std::size_t>(fullHash) & mask_;
        for (;;) {
            SlotRecord& slot = slots_[position];
            if (slot.epoch != epoch_) {
                slot.epoch = epoch_;
                slot.fingerprint = fingerprint;
                slot.outputIndex = static_cast<std::uint32_t>(output.size());
                return true;
            }
            if (slot.fingerprint == fingerprint) {
                const Tet& existing = output[slot.outputIndex];
                const std::array<Vertex, 4> existingVertices{
                    static_cast<Vertex>(existing.v[0]),
                    static_cast<Vertex>(existing.v[1]),
                    static_cast<Vertex>(existing.v[2]),
                    static_cast<Vertex>(existing.v[3])};
                if (sameSet4(existingVertices, vertices)) return false;
            }
            position = (position + 1U) & mask_;
        }
    }

    [[nodiscard]] std::uint64_t retainedBytes() const noexcept {
        return slots_.capacity() * sizeof(SlotRecord);
    }

private:
    struct SlotRecord {
        std::uint32_t outputIndex = 0;
        std::uint16_t epoch = 0;
        std::uint16_t fingerprint = 0;
    };
    static_assert(sizeof(SlotRecord) == 8);

    static std::uint64_t mix(std::uint64_t value) noexcept {
        value ^= value >> 30U;
        value *= 0xbf58476d1ce4e5b9ULL;
        value ^= value >> 27U;
        value *= 0x94d049bb133111ebULL;
        return value ^ (value >> 31U);
    }

    static std::uint64_t unorderedHash(
        const std::array<Vertex, 4>& vertices) noexcept {
        std::uint64_t sum = 0;
        std::uint64_t xored = 0;
        for (Vertex value : vertices) {
            const std::uint64_t hashed = mix(
                static_cast<std::uint64_t>(value) +
                0x9e3779b97f4a7c15ULL);
            sum += hashed;
            xored ^= hashed;
        }
        return sum ^ ((xored << 29U) | (xored >> 35U));
    }

    std::vector<SlotRecord> slots_;
    std::uint16_t epoch_ = 0;
    std::size_t mask_ = 0;
};

void stableRadixMorton(const std::vector<std::uint64_t>& keys,
                       std::vector<Vertex>& order,
                       std::vector<Vertex>& scratch) {
    const std::size_t count = order.size();
    scratch.resize(count);
    std::vector<Vertex>* input = &order;
    std::vector<Vertex>* output = &scratch;
    for (unsigned shift = 0; shift < 64U; shift += 8U) {
        std::array<std::size_t, 256> offsets{};
        for (Vertex index : *input)
            ++offsets[(keys[index] >> shift) & 0xffU];
        std::size_t running = 0;
        for (std::size_t& value : offsets) {
            const std::size_t size = value;
            value = running;
            running += size;
        }
        for (Vertex index : *input) {
            const std::size_t bucket = (keys[index] >> shift) & 0xffU;
            (*output)[offsets[bucket]++] = index;
        }
        std::swap(input, output);
    }
    if (input != &order) order.swap(*input);
}

} // namespace


enum class SpecCandidateStatus : std::uint8_t {
    Ready,
    Duplicate,
    Skipped,
    Seeded
};

struct VersionedRead {
    Slot slot = Arena::kInvalidSlot;
    std::uint32_t version = 0;
};

struct SpecCandidate {
    Vertex pointIndex = 0;
    SpecCandidateStatus status = SpecCandidateStatus::Skipped;
    Slot located = Arena::kInvalidSlot;
    std::vector<Slot> cavity;
    std::vector<ShellFace> shell;
    // Squared radius of a conservative ball centered at shell.f0 that
    // contains the complete replacement Tet circumsphere. Infinity means
    // no safe broad-phase certificate was available.
    std::vector<long double> replacementInfluenceRadius2;
    std::vector<VersionedRead> reads;
    std::vector<Slot> writes;
    bool touchesInfinity = false;
    std::array<long double, 3> influenceMin{};
    std::array<long double, 3> influenceMax{};
    bool hasInfluenceBounds = false;
    std::uint64_t predicateCalls = 0;
    std::uint64_t cavityVisits = 0;

    void reset(Vertex point) {
        pointIndex = point;
        status = SpecCandidateStatus::Skipped;
        located = Arena::kInvalidSlot;
        cavity.clear();
        shell.clear();
        replacementInfluenceRadius2.clear();
        reads.clear();
        writes.clear();
        touchesInfinity = false;
        influenceMin = {0.0L, 0.0L, 0.0L};
        influenceMax = {0.0L, 0.0L, 0.0L};
        hasInfluenceBounds = false;
        predicateCalls = 0;
        cavityVisits = 0;
    }
};

struct SpecScratch {
    std::vector<std::uint32_t> visited;
    std::vector<std::uint32_t> readMarked;
    std::vector<std::uint32_t> writeMarked;
    std::uint32_t visitEpoch = 0;
    std::uint32_t readEpoch = 0;
    std::uint32_t writeEpoch = 0;
    std::vector<Slot> stack;
    std::vector<Slot> work;
    std::vector<Slot> filtered;

    void begin(std::size_t slots) {
        if (visited.size() < slots) visited.resize(slots, 0U);
        if (readMarked.size() < slots) readMarked.resize(slots, 0U);
        if (writeMarked.size() < slots) writeMarked.resize(slots, 0U);
        if (++visitEpoch == 0U) {
            std::fill(visited.begin(), visited.end(), 0U);
            visitEpoch = 1U;
        }
        if (++readEpoch == 0U) {
            std::fill(readMarked.begin(), readMarked.end(), 0U);
            readEpoch = 1U;
        }
        if (++writeEpoch == 0U) {
            std::fill(writeMarked.begin(), writeMarked.end(), 0U);
            writeEpoch = 1U;
        }
        stack.clear();
        work.clear();
        filtered.clear();
    }

    [[nodiscard]] bool isVisited(Slot slot) const noexcept {
        return slot < visited.size() && visited[slot] == visitEpoch;
    }
    void markVisited(Slot slot) noexcept { visited[slot] = visitEpoch; }
    void unmark(Slot slot) noexcept { visited[slot] = 0U; }

    void resetReads(SpecCandidate& candidate) {
        if (++readEpoch == 0U) {
            std::fill(readMarked.begin(), readMarked.end(), 0U);
            readEpoch = 1U;
        }
        candidate.reads.clear();
    }

    void recordRead(Slot slot, const Arena& arena, SpecCandidate& candidate) {
        if (slot == Arena::kInvalidSlot) return;
        if (slot >= readMarked.size()) readMarked.resize(static_cast<std::size_t>(slot) + 1U, 0U);
        if (readMarked[slot] == readEpoch) return;
        readMarked[slot] = readEpoch;
        candidate.reads.push_back({slot, arena.version(slot)});
    }
    void recordWrite(Slot slot, SpecCandidate& candidate) {
        if (slot == Arena::kInvalidSlot) return;
        if (slot >= writeMarked.size())
            writeMarked.resize(static_cast<std::size_t>(slot) + 1U, 0U);
        if (writeMarked[slot] == writeEpoch) return;
        writeMarked[slot] = writeEpoch;
        candidate.writes.push_back(slot);
    }
};



class DelaunayRoundTeam {
public:
    explicit DelaunayRoundTeam(int requestedThreads)
        : threadCount_(resolve(requestedThreads)), barrier_(threadCount_) {
        for (int worker = 1; worker < threadCount_; ++worker) {
            workers_.emplace_back([this, worker] { workerLoop(worker); });
        }
    }

    DelaunayRoundTeam(const DelaunayRoundTeam&) = delete;
    DelaunayRoundTeam& operator=(const DelaunayRoundTeam&) = delete;

    ~DelaunayRoundTeam() {
        stopping_ = true;
        barrier_.arrive_and_wait();
        barrier_.arrive_and_wait();
        for (std::thread& worker : workers_) worker.join();
    }

    [[nodiscard]] int threadCount() const noexcept { return threadCount_; }

    void run(std::size_t count, const std::function<void(std::size_t)>& fn) {
        if (threadCount_ == 1) {
            for (std::size_t i = 0; i < count; ++i) fn(i);
            return;
        }
        count_ = count;
        job_ = &fn;
        {
            std::lock_guard<std::mutex> lock(errorMutex_);
            error_ = nullptr;
        }
        barrier_.arrive_and_wait();
        execute(0);
        barrier_.arrive_and_wait();
        job_ = nullptr;
        std::exception_ptr error;
        {
            std::lock_guard<std::mutex> lock(errorMutex_);
            error = error_;
        }
        if (error) std::rethrow_exception(error);
    }

private:
    static int resolve(int requested) {
        if (requested > 0) return std::max(1, requested);
        const unsigned hw = std::thread::hardware_concurrency();
        return std::max(1, static_cast<int>(hw == 0 ? 1U : hw));
    }

    void recordError() noexcept {
        std::lock_guard<std::mutex> lock(errorMutex_);
        if (!error_) error_ = std::current_exception();
    }

    void execute(int worker) noexcept {
        try {
            const auto* job = job_;
            if (!job) return;
            for (std::size_t i = static_cast<std::size_t>(worker);
                 i < count_; i += static_cast<std::size_t>(threadCount_)) {
                (*job)(i);
            }
        } catch (...) {
            recordError();
        }
    }

    void workerLoop(int worker) {
        for (;;) {
            barrier_.arrive_and_wait();
            if (stopping_) {
                barrier_.arrive_and_wait();
                return;
            }
            execute(worker);
            barrier_.arrive_and_wait();
        }
    }

    int threadCount_ = 1;
    std::barrier<> barrier_;
    std::vector<std::thread> workers_;
    const std::function<void(std::size_t)>* job_ = nullptr;
    std::size_t count_ = 0;
    bool stopping_ = false;
    std::mutex errorMutex_;
    std::exception_ptr error_;
};

struct DelaunayWorkspace::Impl {
    std::vector<Vec3> normalized;
    std::vector<std::uint64_t> mortonKeys;
    std::vector<Vertex> order;
    std::vector<Vertex> orderScratch;
    Arena arena;
    std::vector<Slot> bad;
    std::vector<Slot> stack;
    std::vector<Slot> work;
    std::vector<Slot> stillBad;
    std::vector<ShellFace> shell;
    ShellEdgePairTable pairTable;
    UniqueTetTable uniqueTable;
    std::vector<std::uint8_t> usedPoints;
    std::vector<SpecCandidate> parallelCandidates;
    bool completedRun = false;

    [[nodiscard]] std::uint64_t retainedBytes() const noexcept {
        return normalized.capacity() * sizeof(Vec3) +
               mortonKeys.capacity() * sizeof(std::uint64_t) +
               order.capacity() * sizeof(Vertex) +
               orderScratch.capacity() * sizeof(Vertex) +
               bad.capacity() * sizeof(Slot) +
               stack.capacity() * sizeof(Slot) +
               work.capacity() * sizeof(Slot) +
               stillBad.capacity() * sizeof(Slot) +
               shell.capacity() * sizeof(ShellFace) +
               pairTable.retainedBytes() + uniqueTable.retainedBytes() +
               usedPoints.capacity() * sizeof(std::uint8_t) +
               arena.memory().bytes;
    }

    void release() noexcept {
        Impl empty;
        *this = std::move(empty);
    }
};

DelaunayWorkspace::DelaunayWorkspace() : impl_(std::make_unique<Impl>()) {}
DelaunayWorkspace::~DelaunayWorkspace() = default;
DelaunayWorkspace::DelaunayWorkspace(DelaunayWorkspace&&) noexcept = default;
DelaunayWorkspace& DelaunayWorkspace::operator=(DelaunayWorkspace&&) noexcept =
    default;

std::uint64_t DelaunayWorkspace::retainedBytes() const noexcept {
    return impl_ ? impl_->retainedBytes() : 0;
}

void DelaunayWorkspace::release() noexcept {
    if (impl_) impl_->release();
}

namespace {

template <bool CollectDetailed, bool ExactInputCoordinates = false>
std::vector<Tet> compactDelaunay(const std::vector<Vec3>& pts,
                                 DelaunayWorkspace::Impl& ws,
                                 const Delaunay3DOptions& options,
                                 Delaunay3DStats& st) {
    const Clock::time_point totalBegin = Clock::now();
    st = Delaunay3DStats{};
    st.inputPoints = static_cast<std::uint64_t>(pts.size());
    st.usedCompactArena = true;
    st.usedStableRadixOrder = PAACS_DELAUNAY_USE_RADIX != 0;
    st.usedLinearFinalization = PAACS_DELAUNAY_LINEAR_FINALIZATION != 0;
    st.workspaceReused = ws.completedRun;
    st.collectedDetailedOperationCounts = CollectDetailed;
    st.requestedThreads = options.requestedThreads;
    st.packetSize = options.packetSize;
    const DelaunayParallelPolicy policy =
        options.parallelPolicy != DelaunayParallelPolicy::SerialExact
            ? options.parallelPolicy
            : (options.deterministicParallel
                   ? DelaunayParallelPolicy::UniqueDelaunayParallel
                   : DelaunayParallelPolicy::SerialExact);
    const bool useParallel = policy != DelaunayParallelPolicy::SerialExact;
    const bool strictSerialEquivalent =
        policy == DelaunayParallelPolicy::StrictSerialEquivalent;
    st.usedParallelDiscovery = useParallel;

    std::vector<Tet> output;
    std::uint64_t nonFinite = 0;
    for (const Vec3& point : pts) {
        if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
            !std::isfinite(point.z))
            ++nonFinite;
    }
    if (nonFinite != 0)
        throw std::invalid_argument(
            "delaunay3D: " + std::to_string(nonFinite) +
            " input point(s) have non-finite coordinates");
    if (pts.size() < 4U) {
        st.pointsSkipped = static_cast<std::uint64_t>(pts.size());
        st.runComplete = true;
        st.totalSeconds = secondsBetween(totalBegin, Clock::now());
        ws.completedRun = true;
        st.workspaceRetainedBytes = ws.retainedBytes();
        return output;
    }
    if (pts.size() >= static_cast<std::size_t>(
                           std::numeric_limits<Vertex>::max()))
        throw detail::CompactEnvelopeExceeded(
            "compact Delaunay point-index envelope exceeded");

    const Vertex n = static_cast<Vertex>(pts.size());
    const Vertex infinity = n;

    const Clock::time_point orderingBegin = Clock::now();
    Vec3 low{1e300, 1e300, 1e300};
    Vec3 high{-1e300, -1e300, -1e300};
    for (const Vec3& point : pts) {
        low.x = std::min(low.x, point.x);
        low.y = std::min(low.y, point.y);
        low.z = std::min(low.z, point.z);
        high.x = std::max(high.x, point.x);
        high.y = std::max(high.y, point.y);
        high.z = std::max(high.z, point.z);
    }
    const Vec3 middle = (low + high) * 0.5;
    const double span = std::max(
        {high.x - low.x, high.y - low.y, high.z - low.z, 1e-300});
    const double inverseSpan = 1.0 / span;
    int exactScaleExponent = 0;
    if constexpr (ExactInputCoordinates) {
        const auto largest = std::max({std::abs(low.x),std::abs(low.y),std::abs(low.z),
                                       std::abs(high.x),std::abs(high.y),std::abs(high.z)});
        std::frexp(largest,&exactScaleExponent);
    }

    ws.normalized.resize(static_cast<std::size_t>(n) + 1U);
    ws.mortonKeys.resize(n);
    ws.order.resize(n);
    for (Vertex index = 0; index < n; ++index) {
        const Vec3 normalized = (pts[index] - middle) * inverseSpan;
        Vec3 ordering = normalized;
        if constexpr (ExactInputCoordinates) {
            const Vec3 scaled{std::scalbn(pts[index].x,-exactScaleExponent),
                              std::scalbn(pts[index].y,-exactScaleExponent),
                              std::scalbn(pts[index].z,-exactScaleExponent)};
            // Refuse exponent ranges that would round away a coordinate.
            // Exact power-of-two scaling preserves all determinant signs.
            if (std::scalbn(scaled.x,exactScaleExponent)!=pts[index].x ||
                std::scalbn(scaled.y,exactScaleExponent)!=pts[index].y ||
                std::scalbn(scaled.z,exactScaleExponent)!=pts[index].z)
                throw std::invalid_argument("delaunay3D: exact input scaling range exceeded");
            ws.normalized[index] = scaled;
            // Morton ordering is not geometric authority. Extreme finite
            // inputs can overflow the old centering arithmetic; use the
            // bounded scaled values only for those ordering keys.
            if (!std::isfinite(ordering.x) || !std::isfinite(ordering.y) ||
                !std::isfinite(ordering.z))
                ordering = scaled;
        } else {
            if (!std::isfinite(ordering.x) || !std::isfinite(ordering.y) ||
                !std::isfinite(ordering.z))
                throw std::invalid_argument("delaunay3D: normalization range exceeded");
            ws.normalized[index] = normalized;
        }
        auto quantize = [](double value) {
            return static_cast<std::uint32_t>(std::clamp(
                (value + 0.5) * 2097151.0, 0.0, 2097151.0));
        };
        ws.mortonKeys[index] = morton3(quantize(ordering.x),
                                          quantize(ordering.y),
                                          quantize(ordering.z));
        ws.order[index] = index;
    }
    const double quietNan = std::numeric_limits<double>::quiet_NaN();
    ws.normalized[infinity] = {quietNan, quietNan, quietNan};
#if PAACS_DELAUNAY_USE_RADIX
    stableRadixMorton(ws.mortonKeys, ws.order, ws.orderScratch);
#else
    std::stable_sort(ws.order.begin(), ws.order.end(),
                     [&](Vertex a, Vertex b) {
                         return ws.mortonKeys[a] < ws.mortonKeys[b];
                     });
#endif
    st.normalizationOrderSeconds =
        secondsBetween(orderingBegin, Clock::now());

    const Clock::time_point insertionBegin = Clock::now();
    Arena& arena = ws.arena;
    if (useParallel && strictSerialEquivalent)
        arena.enableVersions();
    else
        arena.disableVersions();
    arena.reset();
    // General random clouds stabilize close to 6.7 finite tets per input
    // point.  Reserving seven slots per point avoids the geometric-capacity
    // overshoot of five parallel vectors without imposing a hard envelope;
    // harder inputs can still grow normally.
    if (static_cast<std::uint64_t>(n) * 7ULL <=
        static_cast<std::uint64_t>(Arena::kMaximumPackedSlot) + 1ULL)
        arena.reserveSlots(static_cast<std::size_t>(n) * 7U + 32U);
    ws.bad.clear();
    ws.stack.clear();
    ws.work.clear();
    ws.stillBad.clear();
    ws.shell.clear();
    const auto& points = ws.normalized;

    auto isInfinite = [&](Slot slot) noexcept {
        return arena.hot(slot).vertices[3] == infinity;
    };

    auto makeGenericTet = [&](Vertex a, Vertex b, Vertex c,
                              Vertex d) -> Slot {
        std::array<Vertex, 4> vertices{a, b, c, d};
        for (int local = 0; local < 3; ++local) {
            if (vertices[local] == infinity)
                std::swap(vertices[local], vertices[3]);
        }
        if (vertices[3] != infinity) {
            const Sign orientation = predicates::orient3dFiltered(
                points[vertices[0]], points[vertices[1]], points[vertices[2]],
                points[vertices[3]]);
            if (orientation == Sign::Negative)
                std::swap(vertices[2], vertices[3]);
            if (orientation == Sign::Zero)
                throw std::logic_error(
                    "delaunay3D: degenerate tet requested");
        }
        return arena.allocate(vertices);
    };

    auto makeDirectFiniteShellTet = [&](const ShellFace& face,
                                        Vertex inserted) -> Slot {
#if PAACS_DELAUNAY_DIRECT_COMMIT
        // shellOk() proved orient3d(f0,f1,f2,inserted) is Negative for a
        // finite cavity owner.  The legacy makeTet path therefore always
        // swaps positions 2/3 and emits exactly {f0,f1,inserted,f2}.
        detailCount<CollectDetailed>(st.directShellTets);
        return arena.allocate(
            {face.f0, face.f1, inserted, face.f2});
#else
        detailCount<CollectDetailed>(st.genericShellTets);
        return makeGenericTet(face.f0, face.f1, face.f2, inserted);
#endif
    };

    Vertex bootstrapA = ws.order[0];
    Vertex bootstrapB = std::numeric_limits<Vertex>::max();
    Vertex bootstrapC = std::numeric_limits<Vertex>::max();
    Vertex bootstrapD = std::numeric_limits<Vertex>::max();
    std::size_t bootstrapBIndex = 0;
    for (std::size_t index = 1;
         index < ws.order.size() &&
         bootstrapB == std::numeric_limits<Vertex>::max();
         ++index) {
        detailCount<CollectDetailed>(st.bootstrapDistinctTests);
        const Vec3& candidate = points[ws.order[index]];
        const Vec3& first = points[bootstrapA];
        if (!(candidate.x == first.x && candidate.y == first.y &&
              candidate.z == first.z)) {
            bootstrapB = ws.order[index];
            bootstrapBIndex = index;
        }
    }
    for (std::size_t i = bootstrapBIndex + 1U;
         i < ws.order.size() &&
         bootstrapD == std::numeric_limits<Vertex>::max();
         ++i) {
        for (std::size_t j = i + 1U;
             j < ws.order.size() &&
             bootstrapD == std::numeric_limits<Vertex>::max();
             ++j) {
            detailCount<CollectDetailed>(st.bootstrapOrientationTests);
            if (predicates::orient3dFiltered(points[bootstrapA], points[bootstrapB],
                                     points[ws.order[i]],
                                     points[ws.order[j]]) != Sign::Zero) {
                bootstrapC = ws.order[i];
                bootstrapD = ws.order[j];
            }
        }
    }
    if (bootstrapD == std::numeric_limits<Vertex>::max()) {
        st.pointsSkipped = static_cast<std::uint64_t>(n);
        st.runComplete = true;
        st.insertionSeconds = secondsBetween(insertionBegin, Clock::now());
        st.totalSeconds = secondsBetween(totalBegin, Clock::now());
        ws.completedRun = true;
        st.workspaceRetainedBytes = ws.retainedBytes();
        return output;
    }

    const Slot bootstrap =
        makeGenericTet(bootstrapA, bootstrapB, bootstrapC, bootstrapD);
    std::array<Slot, 4> infiniteCells{};
    for (int face = 0; face < 4; ++face) {
        const auto& vertices = arena.hot(bootstrap).vertices;
        infiniteCells[face] = makeGenericTet(
            vertices[kFace[face][0]], vertices[kFace[face][1]],
            vertices[kFace[face][2]], infinity);
        arena.connect(infiniteCells[face], 3U, bootstrap,
                      static_cast<std::uint32_t>(face));
    }
    auto sideFace = [&](Slot tet, Vertex edge0, Vertex edge1) {
        const auto& vertices = arena.hot(tet).vertices;
        for (int face = 0; face < 3; ++face) {
            const Vertex x = vertices[kFace[face][0]];
            const Vertex y = vertices[kFace[face][1]];
            const Vertex z = vertices[kFace[face][2]];
            const bool has0 = x == edge0 || y == edge0 || z == edge0;
            const bool has1 = x == edge1 || y == edge1 || z == edge1;
            const bool hasInfinity =
                x == infinity || y == infinity || z == infinity;
            if (has0 && has1 && hasInfinity) return face;
        }
        return -1;
    };
    for (int first = 0; first < 4; ++first) {
        for (int second = first + 1; second < 4; ++second) {
            std::array<Vertex, 2> shared{};
            int count = 0;
            const auto& bootstrapVertices = arena.hot(bootstrap).vertices;
            for (Vertex vertex : bootstrapVertices) {
                if (vertex != bootstrapVertices[first] &&
                    vertex != bootstrapVertices[second])
                    shared[count++] = vertex;
            }
            const int faceFirst =
                sideFace(infiniteCells[first], shared[0], shared[1]);
            const int faceSecond =
                sideFace(infiniteCells[second], shared[0], shared[1]);
            if (faceFirst < 0 || faceSecond < 0)
                throw std::logic_error(
                    "delaunay3D: bootstrap side-face lookup failed");
            arena.connect(infiniteCells[first],
                          static_cast<std::uint32_t>(faceFirst),
                          infiniteCells[second],
                          static_cast<std::uint32_t>(faceSecond));
        }
    }
    st.pointsInserted += 4;
    const std::array<Vertex, 4> seeded = arena.hot(bootstrap).vertices;

    auto conflicts = [&](Slot slot, const Vec3& point) {
        detailCount<CollectDetailed>(st.conflictQueries);
        const auto& tet = arena.hot(slot);
        if (tet.vertices[3] != infinity) {
            detailCount<CollectDetailed>(st.conflictInsphereTests);
            return predicates::insphereFiltered(
                       points[tet.vertices[0]], points[tet.vertices[1]],
                       points[tet.vertices[2]], points[tet.vertices[3]],
                       point) != Sign::Negative;
        }
        const std::uint32_t finiteReference = arena.neighbourRef(slot, 3U);
        const Slot finite = Arena::neighbourSlot(finiteReference);
        const std::uint32_t witnessFace =
            Arena::neighbourFace(finiteReference);
        const auto& finiteTet = arena.hot(finite);
        const Vertex witness = finiteTet.vertices[witnessFace];
        detailCount<CollectDetailed>(st.reciprocalFaceScansAvoided);
        detailCount<CollectDetailed>(st.conflictOrientationTests, 2);
        const Sign witnessOrientation = predicates::orient3dFiltered(
            points[tet.vertices[0]], points[tet.vertices[1]],
            points[tet.vertices[2]], points[witness]);
        const Sign pointOrientation = predicates::orient3dFiltered(
            points[tet.vertices[0]], points[tet.vertices[1]],
            points[tet.vertices[2]], point);
        if (pointOrientation == witnessOrientation) return false;
        if (pointOrientation != Sign::Zero) return true;
        detailCount<CollectDetailed>(st.conflictInsphereTests);
        if (witnessOrientation == Sign::Positive) {
            return predicates::insphereFiltered(
                       points[tet.vertices[0]], points[tet.vertices[1]],
                       points[tet.vertices[2]], points[witness], point) !=
                   Sign::Negative;
        }
        return predicates::insphereFiltered(
                   points[tet.vertices[0]], points[tet.vertices[2]],
                   points[tet.vertices[1]], points[witness], point) !=
               Sign::Negative;
    };

    auto commitPreparedCavity = [&](Vertex pointIndex,
                                    const std::vector<Slot>& badCells,
                                    const std::vector<ShellFace>& shellFaces)
                                    -> Slot {
        for (Slot slot : badCells) arena.free(slot);

        ws.pairTable.begin(shellFaces.size() * 3U);
        detailMaximum<CollectDetailed>(
            st.peakShellPairTableSlots,
            static_cast<std::uint64_t>(ws.pairTable.capacity()));
        Slot firstNew = Arena::kInvalidSlot;
        for (const ShellFace& face : shellFaces) {
            Slot newTet = Arena::kInvalidSlot;
            int pointFace = -1;
            if (face.directFiniteWinding) {
                newTet = makeDirectFiniteShellTet(face, pointIndex);
                pointFace = 2;
            } else {
                detailCount<CollectDetailed>(st.genericShellTets);
                newTet = makeGenericTet(face.f0, face.f1, face.f2, pointIndex);
                const auto& vertices = arena.hot(newTet).vertices;
                for (int local = 0; local < 4; ++local) {
                    if (vertices[local] == pointIndex) {
                        pointFace = local;
                        break;
                    }
                }
            }
            if (pointFace < 0)
                throw std::logic_error(
                    "delaunay3D: inserted-point face not found");
            if (firstNew == Arena::kInvalidSlot) firstNew = newTet;

            if (face.outerRef == Arena::kNoNeighbour) {
                arena.setBoundary(newTet,
                                  static_cast<std::uint32_t>(pointFace));
            } else {
                const Slot outer = Arena::neighbourSlot(face.outerRef);
                const std::uint32_t outerFace =
                    Arena::neighbourFace(face.outerRef);
                arena.connect(newTet, static_cast<std::uint32_t>(pointFace),
                              outer, outerFace);
            }

#if PAACS_DELAUNAY_DIRECT_COMMIT
            if (face.directFiniteWinding) {
                const std::array<PendingSide, 3> sides{
                    PendingSide{newTet, 0U,
                                std::min(face.f1, face.f2),
                                std::max(face.f1, face.f2)},
                    PendingSide{newTet, 1U,
                                std::min(face.f0, face.f2),
                                std::max(face.f0, face.f2)},
                    PendingSide{newTet, 3U,
                                std::min(face.f0, face.f1),
                                std::max(face.f0, face.f1)}};
                for (const PendingSide& side : sides) {
                    PendingSide first{};
                    detailCount<CollectDetailed>(st.shellEdgePairInsertions);
                    const auto result = ws.pairTable.insert<CollectDetailed>(
                        side, first, st.shellEdgePairProbes);
                    if (result == ShellEdgePairTable::Result::Third)
                        throw std::logic_error(
                            "delaunay3D: cavity shell edge used more than twice");
                    if (result == ShellEdgePairTable::Result::Second)
                        arena.connect(first.newTet, first.side,
                                      side.newTet, side.side);
                }
                continue;
            }
#endif
            const auto& vertices = arena.hot(newTet).vertices;
            for (int sideFaceIndex = 0; sideFaceIndex < 4; ++sideFaceIndex) {
                if (sideFaceIndex == pointFace) continue;
                std::array<Vertex, 2> edge{};
                int edgeCount = 0;
                for (int local = 0; local < 3; ++local) {
                    const Vertex vertex =
                        vertices[kFace[sideFaceIndex][local]];
                    if (vertex != pointIndex) edge[edgeCount++] = vertex;
                }
                sort2(edge[0], edge[1]);
                PendingSide side{newTet,
                                 static_cast<std::uint8_t>(sideFaceIndex),
                                 edge[0], edge[1]};
                PendingSide first{};
                detailCount<CollectDetailed>(st.shellEdgePairInsertions);
                const auto result = ws.pairTable.insert<CollectDetailed>(
                        side, first, st.shellEdgePairProbes);
                if (result == ShellEdgePairTable::Result::Third)
                    throw std::logic_error(
                        "delaunay3D: cavity shell edge used more than twice");
                if (result == ShellEdgePairTable::Result::Second)
                    arena.connect(first.newTet, first.side,
                                  side.newTet, side.side);
            }
        }
        if (ws.pairTable.unmatched() != 0)
            throw std::logic_error(
                "delaunay3D: cavity shell is not closed");
        return firstNew;
    };


    auto replacementInfluenceRadius2 = [&](Vertex a, Vertex b, Vertex c,
                                           Vertex d) -> long double {
        // Tight conservative circumradius bound using Cramer's rule.  Inputs
        // are binary64 coordinates represented exactly in long double.  Each
        // 3x3 determinant is accompanied by a deliberately conservative
        // roundoff envelope proportional to the sum of its six triple-product
        // magnitudes.  If the denominator interval contains zero we disable
        // the broad phase and keep the exact insphere test authoritative.
        const Vec3& p0 = points[a];
        const Vec3& p1 = points[b];
        const Vec3& p2 = points[c];
        const Vec3& p3 = points[d];
        const long double x1 = static_cast<long double>(p1.x) - p0.x;
        const long double y1 = static_cast<long double>(p1.y) - p0.y;
        const long double z1 = static_cast<long double>(p1.z) - p0.z;
        const long double x2 = static_cast<long double>(p2.x) - p0.x;
        const long double y2 = static_cast<long double>(p2.y) - p0.y;
        const long double z2 = static_cast<long double>(p2.z) - p0.z;
        const long double x3 = static_cast<long double>(p3.x) - p0.x;
        const long double y3 = static_cast<long double>(p3.y) - p0.y;
        const long double z3 = static_cast<long double>(p3.z) - p0.z;
        const long double h1 = (x1*x1 + y1*y1 + z1*z1) * 0.5L;
        const long double h2 = (x2*x2 + y2*y2 + z2*z2) * 0.5L;
        const long double h3 = (x3*x3 + y3*y3 + z3*z3) * 0.5L;
        const long double eps = std::numeric_limits<long double>::epsilon();
        auto detBound = [&](long double a00, long double a01, long double a02,
                            long double a10, long double a11, long double a12,
                            long double a20, long double a21, long double a22) {
            const std::array<long double, 6> terms{
                a00*a11*a22, a01*a12*a20, a02*a10*a21,
                -a02*a11*a20, -a01*a10*a22, -a00*a12*a21};
            long double value = 0.0L, magnitude = 0.0L;
            for (long double term : terms) {
                value += term;
                magnitude += std::abs(term);
            }
            // Covers three input-square sums, triple products and the five
            // determinant additions with ample margin on IEEE extended/quad.
            const long double error = 128.0L * eps * magnitude +
                                      16.0L * std::numeric_limits<long double>::min();
            return std::pair<long double,long double>{value, error};
        };
        const auto [det, detErr] =
            detBound(x1,y1,z1, x2,y2,z2, x3,y3,z3);
        const long double detLower = std::abs(det) - detErr;
        if (!(detLower > 0.0L))
            return std::numeric_limits<long double>::infinity();
        const auto [nx, ex] =
            detBound(h1,y1,z1, h2,y2,z2, h3,y3,z3);
        const auto [ny, ey] =
            detBound(x1,h1,z1, x2,h2,z2, x3,h3,z3);
        const auto [nz, ez] =
            detBound(x1,y1,h1, x2,y2,h2, x3,y3,h3);
        const long double cx = (std::abs(nx) + ex) / detLower;
        const long double cy = (std::abs(ny) + ey) / detLower;
        const long double cz = (std::abs(nz) + ez) / detLower;
        const long double rUpper = std::sqrt(cx*cx + cy*cy + cz*cz) *
                                   (1.0L + 256.0L*eps);
        if (!std::isfinite(rUpper))
            return std::numeric_limits<long double>::infinity();
        const long double diameter = 2.0L * rUpper;
        return diameter * diameter;
    };

    auto discoverSpecCandidate = [&](Vertex pointIndex, Slot startSeed,
                                     SpecCandidate& candidate) {
        static thread_local SpecScratch scratch;
        scratch.begin(arena.slotCount());
        candidate.reset(pointIndex);
        const Vec3 point = points[pointIndex];

        const bool trackSerialReads = strictSerialEquivalent;
        auto readSlot = [&](Slot slot) {
            if (trackSerialReads) scratch.recordRead(slot, arena, candidate);
        };
        auto localConflicts = [&](Slot slot) {
            readSlot(slot);
            ++candidate.predicateCalls;
            const auto& tet = arena.hot(slot);
            if (tet.vertices[3] != infinity) {
                return predicates::insphereFiltered(
                           points[tet.vertices[0]], points[tet.vertices[1]],
                           points[tet.vertices[2]], points[tet.vertices[3]],
                           point) != Sign::Negative;
            }
            const std::uint32_t finiteReference = arena.neighbourRef(slot, 3U);
            const Slot finite = Arena::neighbourSlot(finiteReference);
            if (finite == Arena::kInvalidSlot) return false;
            readSlot(finite);
            const std::uint32_t witnessFace = Arena::neighbourFace(finiteReference);
            const auto& finiteTet = arena.hot(finite);
            const Vertex witness = finiteTet.vertices[witnessFace];
            candidate.predicateCalls += 2;
            const Sign witnessOrientation = predicates::orient3dFiltered(
                points[tet.vertices[0]], points[tet.vertices[1]],
                points[tet.vertices[2]], points[witness]);
            const Sign pointOrientation = predicates::orient3dFiltered(
                points[tet.vertices[0]], points[tet.vertices[1]],
                points[tet.vertices[2]], point);
            if (pointOrientation == witnessOrientation) return false;
            if (pointOrientation != Sign::Zero) return true;
            ++candidate.predicateCalls;
            if (witnessOrientation == Sign::Positive) {
                return predicates::insphereFiltered(
                           points[tet.vertices[0]], points[tet.vertices[1]],
                           points[tet.vertices[2]], points[witness], point) !=
                       Sign::Negative;
            }
            return predicates::insphereFiltered(
                       points[tet.vertices[0]], points[tet.vertices[2]],
                       points[tet.vertices[1]], points[witness], point) !=
                   Sign::Negative;
        };

        Slot current = startSeed;
        if (current == Arena::kInvalidSlot || !arena.alive(current))
            current = arena.liveTail();
        Slot located = Arena::kInvalidSlot;
        const std::uint64_t maximumSteps =
            4ULL * static_cast<std::uint64_t>(arena.slotCount()) + 64ULL;
        for (std::uint64_t step = 0;
             current != Arena::kInvalidSlot && step < maximumSteps; ++step) {
            readSlot(current);
            const auto& tet = arena.hot(current);
            if (tet.vertices[3] == infinity) {
                if (localConflicts(current)) {
                    located = current;
                    break;
                }
                current = arena.neighbour(current, 3U);
                continue;
            }
            int bestFace = -1;
            double bestValue = -1e300;
            for (int face = 0; face < 4; ++face) {
                ++candidate.predicateCalls;
                const Vec3& f0 = points[tet.vertices[kFace[face][0]]];
                const Vec3& f1 = points[tet.vertices[kFace[face][1]]];
                const Vec3& f2 = points[tet.vertices[kFace[face][2]]];
                if (predicates::orient3dFiltered(f0, f1, f2, point) !=
                    Sign::Positive)
                    continue;
                const double heuristic = orientHeur(f0, f1, f2, point);
                if (heuristic > bestValue) {
                    bestValue = heuristic;
                    bestFace = face;
                }
            }
            if (bestFace < 0) {
                located = current;
                break;
            }
            current = arena.neighbour(current,
                                       static_cast<std::uint32_t>(bestFace));
        }
        candidate.located = located;

        if (located != Arena::kInvalidSlot) {
            readSlot(located);
            bool duplicate = false;
            for (Vertex vertex : arena.hot(located).vertices) {
                if (vertex == infinity) continue;
                const Vec3& existing = points[vertex];
                duplicate = duplicate ||
                            (existing.x == point.x && existing.y == point.y &&
                             existing.z == point.z);
            }
            if (duplicate) {
                candidate.status = SpecCandidateStatus::Duplicate;
                return;
            }
        }

        Slot seed = located;
        if (seed == Arena::kInvalidSlot) {
            for (Slot slot = 0; slot < arena.slotCount(); ++slot) {
                if (!arena.alive(slot)) continue;
                if (localConflicts(slot)) {
                    seed = slot;
                    break;
                }
            }
        }
        if (seed == Arena::kInvalidSlot) {
            candidate.status = SpecCandidateStatus::Skipped;
            return;
        }

        // The point-location walk is a search mechanism, not a dependency of
        // an already discovered cavity. Revalidation therefore starts from
        // the exact conflict seed/cavity rather than conservatively treating
        // every transient walk cell as a serial-order dependency.
        if (trackSerialReads) scratch.resetReads(candidate);
        scratch.stack.push_back(seed);
        scratch.markVisited(seed);
        while (!scratch.stack.empty()) {
            const Slot slot = scratch.stack.back();
            scratch.stack.pop_back();
            candidate.cavity.push_back(slot);
            ++candidate.cavityVisits;
            readSlot(slot);
            for (int face = 0; face < 4; ++face) {
                const Slot neighbour = arena.neighbour(
                    slot, static_cast<std::uint32_t>(face));
                if (neighbour == Arena::kInvalidSlot ||
                    scratch.isVisited(neighbour))
                    continue;
                if (!localConflicts(neighbour)) continue;
                scratch.markVisited(neighbour);
                scratch.stack.push_back(neighbour);
            }
        }

        auto shellOkayLocal = [&](Slot slot) {
            readSlot(slot);
            const auto& tet = arena.hot(slot);
            for (int face = 0; face < 4; ++face) {
                const Slot neighbour = arena.neighbour(
                    slot, static_cast<std::uint32_t>(face));
                if (neighbour != Arena::kInvalidSlot &&
                    scratch.isVisited(neighbour))
                    continue;
                const Vertex f0 = tet.vertices[kFace[face][0]];
                const Vertex f1 = tet.vertices[kFace[face][1]];
                const Vertex f2 = tet.vertices[kFace[face][2]];
                if (f0 == infinity || f1 == infinity || f2 == infinity)
                    continue;
                ++candidate.predicateCalls;
                const Sign pointOrientation = predicates::orient3dFiltered(
                    points[f0], points[f1], points[f2], point);
                if (tet.vertices[3] != infinity) {
                    if (pointOrientation != Sign::Negative) return false;
                } else {
                    if (neighbour == Arena::kInvalidSlot) return false;
                    const std::uint32_t reference = arena.neighbourRef(
                        slot, static_cast<std::uint32_t>(face));
                    const std::uint32_t witnessFace =
                        Arena::neighbourFace(reference);
                    readSlot(neighbour);
                    const auto& finiteTet = arena.hot(neighbour);
                    const Vertex witness = finiteTet.vertices[witnessFace];
                    ++candidate.predicateCalls;
                    const Sign witnessOrientation = predicates::orient3dFiltered(
                        points[f0], points[f1], points[f2], points[witness]);
                    if (pointOrientation == Sign::Zero ||
                        pointOrientation == witnessOrientation)
                        return false;
                }
            }
            return true;
        };

        candidate.touchesInfinity = false;
        for (Slot slot : candidate.cavity)
            candidate.touchesInfinity = candidate.touchesInfinity ||
                                        arena.hot(slot).vertices[3] == infinity;
        if (candidate.touchesInfinity) {
            for (Slot slot : candidate.cavity) {
                if (arena.hot(slot).vertices[3] == infinity)
                    scratch.work.push_back(slot);
            }
            bool resurrected = false;
            while (!scratch.work.empty()) {
                const Slot slot = scratch.work.back();
                scratch.work.pop_back();
                if (!scratch.isVisited(slot)) continue;
                if (shellOkayLocal(slot)) continue;
                scratch.unmark(slot);
                resurrected = true;
                for (int face = 0; face < 4; ++face) {
                    const Slot neighbour = arena.neighbour(
                        slot, static_cast<std::uint32_t>(face));
                    if (neighbour != Arena::kInvalidSlot &&
                        scratch.isVisited(neighbour))
                        scratch.work.push_back(neighbour);
                }
            }
            if (resurrected) {
                scratch.filtered.clear();
                for (Slot slot : candidate.cavity) {
                    if (scratch.isVisited(slot)) scratch.filtered.push_back(slot);
                }
                candidate.cavity.assign(scratch.filtered.begin(),
                                        scratch.filtered.end());
            }
        }
        if (candidate.cavity.empty()) {
            candidate.status = SpecCandidateStatus::Skipped;
            return;
        }

        candidate.shell.clear();
        candidate.writes.clear();
        candidate.writes.reserve(candidate.cavity.size() + 16U);
        for (Slot slot : candidate.cavity) scratch.recordWrite(slot, candidate);
        for (Slot slot : candidate.cavity) {
            readSlot(slot);
            const auto& tet = arena.hot(slot);
            const bool finiteOwner = tet.vertices[3] != infinity;
            for (int face = 0; face < 4; ++face) {
                const Slot neighbour = arena.neighbour(
                    slot, static_cast<std::uint32_t>(face));
                if (neighbour != Arena::kInvalidSlot &&
                    scratch.isVisited(neighbour))
                    continue;
                const std::uint32_t outerRef = arena.neighbourRef(
                    slot, static_cast<std::uint32_t>(face));
                candidate.shell.push_back(
                    {tet.vertices[kFace[face][0]],
                     tet.vertices[kFace[face][1]],
                     tet.vertices[kFace[face][2]], outerRef, finiteOwner});
                if (neighbour != Arena::kInvalidSlot) {
                    readSlot(neighbour);
                    scratch.recordWrite(neighbour, candidate);
                }
            }
        }
        // write footprint is already unique through SpecScratch::recordWrite.
        std::sort(candidate.writes.begin(), candidate.writes.end());
        if (trackSerialReads) {
            std::sort(candidate.reads.begin(), candidate.reads.end(),
                      [](const VersionedRead& a, const VersionedRead& b) {
                          return a.slot < b.slot;
                      });
        }
        candidate.replacementInfluenceRadius2.clear();
        candidate.hasInfluenceBounds = false;
        candidate.status = SpecCandidateStatus::Ready;
    };

    auto prepareInfluence = [&](SpecCandidate& candidate) {
        if (!candidate.replacementInfluenceRadius2.empty() ||
            candidate.shell.empty())
            return;
        candidate.replacementInfluenceRadius2.reserve(candidate.shell.size());
        candidate.hasInfluenceBounds = true;
        candidate.influenceMin = {
            std::numeric_limits<long double>::infinity(),
            std::numeric_limits<long double>::infinity(),
            std::numeric_limits<long double>::infinity()};
        candidate.influenceMax = {
            -std::numeric_limits<long double>::infinity(),
            -std::numeric_limits<long double>::infinity(),
            -std::numeric_limits<long double>::infinity()};
        for (const ShellFace& face : candidate.shell) {
            long double radius2 = std::numeric_limits<long double>::infinity();
            if (face.directFiniteWinding) {
                radius2 = replacementInfluenceRadius2(
                    face.f0, face.f1, candidate.pointIndex, face.f2);
            }
            candidate.replacementInfluenceRadius2.push_back(radius2);
            if (!std::isfinite(radius2)) {
                candidate.hasInfluenceBounds = false;
                continue;
            }
            const long double radius = std::sqrt(radius2);
            const Vec3& origin = points[face.f0];
            const std::array<long double,3> q{
                static_cast<long double>(origin.x),
                static_cast<long double>(origin.y),
                static_cast<long double>(origin.z)};
            for (int axis = 0; axis < 3; ++axis) {
                candidate.influenceMin[axis] =
                    std::min(candidate.influenceMin[axis], q[axis] - radius);
                candidate.influenceMax[axis] =
                    std::max(candidate.influenceMax[axis], q[axis] + radius);
            }
        }
    };

    auto versionsUnchanged = [&](const SpecCandidate& candidate) {
        for (const VersionedRead& read : candidate.reads) {
            if (arena.version(read.slot) != read.version) return false;
        }
        return true;
    };

    auto writeReadIntersect = [](const std::vector<Slot>& writes,
                                     const std::vector<VersionedRead>& reads) {
        std::size_t wi = 0, ri = 0;
        while (wi < writes.size() && ri < reads.size()) {
            if (writes[wi] < reads[ri].slot) ++wi;
            else if (reads[ri].slot < writes[wi]) ++ri;
            else return true;
        }
        return false;
    };

    auto replacementMayAffect = [&](const SpecCandidate& earlier,
                                    Vertex laterPoint) {
        if (earlier.status != SpecCandidateStatus::Ready) return false;
        if (earlier.touchesInfinity) return true;
        const Vec3& p = points[laterPoint];
        if (earlier.hasInfluenceBounds &&
            (static_cast<long double>(p.x) < earlier.influenceMin[0] ||
             static_cast<long double>(p.x) > earlier.influenceMax[0] ||
             static_cast<long double>(p.y) < earlier.influenceMin[1] ||
             static_cast<long double>(p.y) > earlier.influenceMax[1] ||
             static_cast<long double>(p.z) < earlier.influenceMin[2] ||
             static_cast<long double>(p.z) > earlier.influenceMax[2])) {
            ++st.parallelCandidateAabbCulls;
            return false;
        }
        for (std::size_t i = 0; i < earlier.shell.size(); ++i) {
            const ShellFace& face = earlier.shell[i];
            // A finite shell replacement is exactly
            // {f0, f1, insertedPoint, f2} in the compact direct path.
            if (!face.directFiniteWinding) return true;
            const long double dx = static_cast<long double>(p.x) -
                                   points[face.f0].x;
            const long double dy = static_cast<long double>(p.y) -
                                   points[face.f0].y;
            const long double dz = static_cast<long double>(p.z) -
                                   points[face.f0].z;
            const long double distance2 = dx*dx + dy*dy + dz*dz;
            if (distance2 > earlier.replacementInfluenceRadius2[i]) {
                ++st.parallelHazardBroadPhaseCulls;
                continue;
            }
            ++st.parallelHazardPredicateTests;
            if (predicates::insphereFiltered(
                    points[face.f0], points[face.f1],
                    points[earlier.pointIndex], points[face.f2], p) !=
                Sign::Negative)
                return true;
        }
        return false;
    };

    Slot lastAlive = bootstrap;
    if (useParallel) {
        DelaunayRoundTeam roundTeam(options.requestedThreads);
        if (strictSerialEquivalent) {
            const std::size_t packetCapacity = std::max<std::size_t>(
                1U, std::min<std::size_t>(options.packetSize, 512U));
            ws.parallelCandidates.resize(packetCapacity);
            std::vector<std::size_t> packetRanks(packetCapacity, ws.order.size());
            std::size_t nextRank = 0U;
            auto isBootstrapVertex = [&](Vertex pointIndex) {
                return std::find(seeded.begin(), seeded.end(), pointIndex) !=
                       seeded.end();
            };
            while (nextRank < ws.order.size()) {
                while (nextRank < ws.order.size() &&
                       isBootstrapVertex(ws.order[nextRank]))
                    ++nextRank;
                if (nextRank >= ws.order.size()) break;
                ++st.parallelRounds;
                std::size_t packetCount = 0U;
                std::size_t scanRank = nextRank;
                while (scanRank < ws.order.size() &&
                       packetCount < packetCapacity) {
                    if (!isBootstrapVertex(ws.order[scanRank]))
                        packetRanks[packetCount++] = scanRank;
                    ++scanRank;
                }
                st.speculativeDiscoveries += packetCount;
                const Clock::time_point discoveryBegin = Clock::now();
                const Slot snapshotSeed = lastAlive;
                roundTeam.run(packetCount, [&](std::size_t i) {
                    discoverSpecCandidate(ws.order[packetRanks[i]],
                                          snapshotSeed,
                                          ws.parallelCandidates[i]);
                });
                st.parallelDiscoverySeconds +=
                    secondsBetween(discoveryBegin, Clock::now());

                std::vector<std::size_t> accepted;
                accepted.reserve(packetCount);
                bool blocked = false;
                for (std::size_t i = 0; i < packetCount; ++i) {
                    SpecCandidate& candidate = ws.parallelCandidates[i];
                    if (candidate.status == SpecCandidateStatus::Ready)
                        ++st.parallelReadyObservations;
                    else if (candidate.status == SpecCandidateStatus::Skipped)
                        ++st.parallelSkippedObservations;
                    if (candidate.status == SpecCandidateStatus::Duplicate) {
                        accepted.push_back(i);
                        continue;
                    }
                    if (candidate.status != SpecCandidateStatus::Ready) {
                        blocked = true;
                        break;
                    }
                    bool equivalent = true;
                    for (std::size_t earlierIndex : accepted) {
                        ++st.parallelConflictPairChecks;
                        const SpecCandidate& earlier =
                            ws.parallelCandidates[earlierIndex];
                        if (earlier.status != SpecCandidateStatus::Ready)
                            continue;
                        if (writeReadIntersect(earlier.writes, candidate.reads) ||
                            writeReadIntersect(candidate.writes, earlier.reads) ||
                            replacementMayAffect(earlier,
                                                 candidate.pointIndex)) {
                            equivalent = false;
                            break;
                        }
                    }
                    if (!equivalent) {
                        blocked = true;
                        break;
                    }
                    prepareInfluence(candidate);
                    accepted.push_back(i);
                }
                if (accepted.empty()) {
                    // Exact serial progress at the first unresolved Morton rank.
                    SpecCandidate& candidate = ws.parallelCandidates[0];
                    discoverSpecCandidate(ws.order[packetRanks[0]], lastAlive,
                                          candidate);
                    ++st.speculativeCandidatesRecomputed;
                    accepted.push_back(0U);
                }

                const Clock::time_point commitBegin = Clock::now();
                std::size_t committed = 0U;
                for (std::size_t i : accepted) {
                    SpecCandidate& candidate = ws.parallelCandidates[i];
                    ++st.speculativeCandidates;
                    st.speculativePredicateCalls += candidate.predicateCalls;
                    st.speculativeCavityVisits += candidate.cavityVisits;
                    if (candidate.status == SpecCandidateStatus::Duplicate) {
                        ++st.pointsDuplicate;
                    } else {
                        if (!versionsUnchanged(candidate)) {
                            // The directional certificate was conservative;
                            // if a read changed, stop the prefix and rediscover
                            // this Morton point on the next round.
                            ++st.speculativeCandidatesRecomputed;
                            break;
                        }
                        const Slot firstNew = commitPreparedCavity(
                            candidate.pointIndex, candidate.cavity,
                            candidate.shell);
                        if (firstNew != Arena::kInvalidSlot) lastAlive = firstNew;
                        ++st.pointsInserted;
                    }
                    ++st.speculativeCandidatesReused;
                    st.serialEquivalentPredicateCalls += candidate.predicateCalls;
                    st.serialEquivalentCavityVisits += candidate.cavityVisits;
                    ++committed;
                }
                st.canonicalCommitSeconds +=
                    secondsBetween(commitBegin, Clock::now());
                if (committed == 0U)
                    throw std::logic_error(
                        "delaunay3D: strict packet made no progress");
                // Advance exactly over the committed Morton-prefix entries,
                // skipping bootstrap vertices but never crossing an unresolved
                // candidate.
                for (std::size_t k = 0; k < committed; ++k) {
                    do { ++nextRank; }
                    while (nextRank < ws.order.size() &&
                           isBootstrapVertex(ws.order[nextRank]));
                }
                (void)blocked;
            }
        } else {
        // V0.9.3 uses deterministic virtual Morton partitions. Each lane
        // advances in its own canonical Morton order; lane heads are evaluated
        // concurrently from one immutable snapshot. Only pairwise independent
        // finite-cavity candidates commit in a round. Exact degeneracies are
        // handled by the production serial fallback (see qualification probe).
        // Seed the whole normalized domain, not a Morton prefix. A Morton
        // prefix occupies only one spatial region and leaves almost every
        // distant lane head coupled to the symbolic hull. The deterministic
        // sparse skeleton establishes a global hull before concurrent waves.
        // General-position parallel path: seed the global hull with a small,
        // deterministic set of directional extreme points.  The seed set is
        // independent of worker count and packet size.  Exact degeneracies are
        // handled by the serial fallback gate below, so reordering in this path
        // does not change the unique Delaunay complex on its supported domain.
        std::vector<std::size_t> pointRank(n, ws.order.size());
        for (std::size_t rank = 0; rank < ws.order.size(); ++rank)
            pointRank[ws.order[rank]] = rank;
        std::vector<Vertex> extremeSeeds;
        extremeSeeds.reserve(192U);
        constexpr int kDirectionRadius = 2;
        for (int dx = -kDirectionRadius; dx <= kDirectionRadius; ++dx) {
            for (int dy = -kDirectionRadius; dy <= kDirectionRadius; ++dy) {
                for (int dz = -kDirectionRadius; dz <= kDirectionRadius; ++dz) {
                    if (dx == 0 && dy == 0 && dz == 0) continue;
                    const int gcdXY = std::gcd(std::abs(dx), std::abs(dy));
                    const int gcdXYZ = std::gcd(gcdXY, std::abs(dz));
                    if (gcdXYZ > 1) continue;
                    double best = -std::numeric_limits<double>::infinity();
                    Vertex bestVertex = 0;
                    for (Vertex vertex = 0; vertex < n; ++vertex) {
                        const Vec3& q = points[vertex];
                        const double projection =
                            static_cast<double>(dx) * q.x +
                            static_cast<double>(dy) * q.y +
                            static_cast<double>(dz) * q.z;
                        if (projection > best ||
                            (projection == best && vertex < bestVertex)) {
                            best = projection;
                            bestVertex = vertex;
                        }
                    }
                    extremeSeeds.push_back(bestVertex);
                }
            }
        }
        std::sort(extremeSeeds.begin(), extremeSeeds.end(),
                  [&](Vertex a, Vertex b) {
                      return pointRank[a] < pointRank[b];
                  });
        extremeSeeds.erase(
            std::unique(extremeSeeds.begin(), extremeSeeds.end()),
            extremeSeeds.end());
        std::vector<std::uint8_t> preinsertedRank(ws.order.size(), 0U);
        for (Vertex pointIndex : extremeSeeds) {
            const std::size_t rank = pointRank[pointIndex];
            if (std::find(seeded.begin(), seeded.end(), pointIndex) != seeded.end()) {
                preinsertedRank[rank] = 1U;
                continue;
            }
            SpecCandidate candidate;
            discoverSpecCandidate(pointIndex, lastAlive, candidate);
            if (candidate.status == SpecCandidateStatus::Duplicate) {
                ++st.pointsDuplicate;
                preinsertedRank[rank] = 1U;
                continue;
            }
            if (candidate.status != SpecCandidateStatus::Ready) continue;
            const Slot firstNew = commitPreparedCavity(
                pointIndex, candidate.cavity, candidate.shell);
            if (firstNew != Arena::kInvalidSlot) lastAlive = firstNew;
            ++st.pointsInserted;
            preinsertedRank[rank] = 1U;
        }

        // Exact half-space hull closure.  Infinite tetrahedra encode the
        // current convex-hull triangles.  Points classified on the witness side
        // of every hull triangle are inside the current hull and remain inside
        // after subsequent hull expansion.  Only points outside (or exactly on
        // a hull plane, conservatively) are serially inserted here.
        struct HullPlane {
            Vertex f0 = 0, f1 = 0, f2 = 0;
            Sign witnessOrientation = Sign::Zero;
        };
        std::vector<HullPlane> hullPlanes;
        hullPlanes.reserve(512U);
        for (Slot slot = arena.liveHead(); slot != Arena::kInvalidSlot;
             slot = arena.nextLive(slot)) {
            if (!isInfinite(slot)) continue;
            const auto& tet = arena.hot(slot);
            const std::uint32_t finiteRef = arena.neighbourRef(slot, 3U);
            const Slot finite = Arena::neighbourSlot(finiteRef);
            if (finite == Arena::kInvalidSlot) continue;
            const std::uint32_t witnessFace = Arena::neighbourFace(finiteRef);
            const Vertex witness = arena.hot(finite).vertices[witnessFace];
            const Sign witnessOrientation = predicates::orient3dFiltered(
                points[tet.vertices[0]], points[tet.vertices[1]],
                points[tet.vertices[2]], points[witness]);
            if (witnessOrientation == Sign::Zero)
                throw std::logic_error("delaunay3D: degenerate hull witness");
            hullPlanes.push_back({tet.vertices[0], tet.vertices[1],
                                  tet.vertices[2], witnessOrientation});
        }
        std::vector<std::uint8_t> outsideHull(ws.order.size(), 0U);
        const Clock::time_point hullBegin = Clock::now();
        roundTeam.run(ws.order.size(), [&](std::size_t rank) {
            if (preinsertedRank[rank] != 0U ||
                std::find(seeded.begin(), seeded.end(), ws.order[rank]) !=
                    seeded.end())
                return;
            const Vertex pointIndex = ws.order[rank];
            for (const HullPlane& face : hullPlanes) {
                const Sign side = predicates::orient3dFiltered(
                    points[face.f0], points[face.f1], points[face.f2],
                    points[pointIndex]);
                if (side != face.witnessOrientation) {
                    outsideHull[rank] = 1U;
                    return;
                }
            }
        });
        st.parallelHullPrepassSeconds =
            secondsBetween(hullBegin, Clock::now());
        for (std::size_t rank = 0; rank < ws.order.size(); ++rank) {
            if (preinsertedRank[rank] != 0U ||
                std::find(seeded.begin(), seeded.end(), ws.order[rank]) !=
                    seeded.end())
                continue;
            ++st.parallelHullPrepassCandidates;
            if (outsideHull[rank] == 0U) continue;
            ++st.parallelHullPrepassOutside;
            SpecCandidate candidate;
            discoverSpecCandidate(ws.order[rank], lastAlive, candidate);
            if (candidate.status == SpecCandidateStatus::Duplicate) {
                ++st.pointsDuplicate;
                preinsertedRank[rank] = 1U;
                continue;
            }
            if (candidate.status != SpecCandidateStatus::Ready) {
                ++st.pointsSkipped;
                preinsertedRank[rank] = 1U;
                continue;
            }
            const Slot firstNew = commitPreparedCavity(
                candidate.pointIndex, candidate.cavity, candidate.shell);
            if (firstNew != Arena::kInvalidSlot) lastAlive = firstNew;
            ++st.pointsInserted;
            preinsertedRank[rank] = 1U;
        }

        // Spatially separated deterministic packet schedule.  The current hull
        // is already closed, so all remaining points are interior.  On the
        // general-position path the exact Delaunay complex is unique; packet
        // size and physical worker count are therefore execution parameters,
        // not mathematical inputs.  Exact-zero topology events are handled by
        // the serial fallback gate before this mode is promoted.
        constexpr std::size_t kGridDim = 16U;
        constexpr std::size_t kGridCells = kGridDim * kGridDim * kGridDim;
        std::array<std::vector<std::uint16_t>, 8> colorCells;
        std::vector<std::vector<std::size_t>> binRanks(kGridCells);
        std::vector<std::size_t> binCursor(kGridCells, 0U);
        std::vector<Slot> binSeed(kGridCells, lastAlive);
        auto gridCoord = [](double value) -> std::size_t {
            const double scaled = std::clamp(value + 0.5, 0.0,
                                             std::nextafter(1.0, 0.0));
            return std::min<std::size_t>(kGridDim - 1U,
                static_cast<std::size_t>(scaled * static_cast<double>(kGridDim)));
        };
        auto cellIndex = [&](Vertex pointIndex) -> std::size_t {
            const Vec3& q = points[pointIndex];
            const std::size_t ix = gridCoord(q.x);
            const std::size_t iy = gridCoord(q.y);
            const std::size_t iz = gridCoord(q.z);
            return (iz * kGridDim + iy) * kGridDim + ix;
        };
        for (std::size_t iz = 0; iz < kGridDim; ++iz) {
            for (std::size_t iy = 0; iy < kGridDim; ++iy) {
                for (std::size_t ix = 0; ix < kGridDim; ++ix) {
                    const std::size_t cell =
                        (iz * kGridDim + iy) * kGridDim + ix;
                    const std::size_t color =
                        (ix & 1U) | ((iy & 1U) << 1U) | ((iz & 1U) << 2U);
                    colorCells[color].push_back(static_cast<std::uint16_t>(cell));
                }
            }
        }
        std::vector<std::uint8_t> doneRank(ws.order.size(), 0U);
        for (std::size_t rank = 0; rank < ws.order.size(); ++rank) {
            if (preinsertedRank[rank] != 0U ||
                std::find(seeded.begin(), seeded.end(), ws.order[rank]) !=
                    seeded.end()) {
                doneRank[rank] = 1U;
                continue;
            }
            binRanks[cellIndex(ws.order[rank])].push_back(rank);
        }
        const std::size_t packetCapacity = std::max<std::size_t>(
            1U, std::min<std::size_t>(options.packetSize, 512U));
        ws.parallelCandidates.resize(packetCapacity);
        std::vector<std::size_t> packetRanks(packetCapacity, ws.order.size());
        std::vector<std::uint16_t> packetCells(packetCapacity, 0U);
        std::vector<std::uint32_t> selectedWriteMarks(arena.slotCount(), 0U);
        std::uint32_t selectedWriteEpoch = 0U;

        for (std::size_t color = 0; color < colorCells.size(); ++color) {
            std::size_t pendingInColor = 0U;
            for (std::uint16_t cell : colorCells[color])
                pendingInColor += binRanks[cell].size();
            std::size_t scanStart = 0U;
            while (pendingInColor != 0U) {
                ++st.parallelRounds;
                std::size_t packetCount = 0U;
                const auto& cells = colorCells[color];
                for (std::size_t scan = 0;
                     scan < cells.size() && packetCount < packetCapacity; ++scan) {
                    const std::size_t ci = (scanStart + scan) % cells.size();
                    const std::uint16_t cell = cells[ci];
                    if (binCursor[cell] >= binRanks[cell].size()) continue;
                    packetRanks[packetCount] = binRanks[cell][binCursor[cell]];
                    packetCells[packetCount] = cell;
                    ++packetCount;
                }
                if (packetCount == 0U)
                    throw std::logic_error(
                        "delaunay3D: spatial packet made no progress");
                scanStart = (scanStart + packetCount) % cells.size();
                st.speculativeDiscoveries += packetCount;

                const Clock::time_point discoveryBegin = Clock::now();
                roundTeam.run(packetCount, [&](std::size_t i) {
                    SpecCandidate& candidate = ws.parallelCandidates[i];
                    const std::size_t rank = packetRanks[i];
                    discoverSpecCandidate(ws.order[rank],
                                          binSeed[packetCells[i]], candidate);
                });
                st.parallelDiscoverySeconds +=
                    secondsBetween(discoveryBegin, Clock::now());

                if (selectedWriteMarks.size() < arena.slotCount())
                    selectedWriteMarks.resize(arena.slotCount(), 0U);
                if (++selectedWriteEpoch == 0U) {
                    std::fill(selectedWriteMarks.begin(),
                              selectedWriteMarks.end(), 0U);
                    selectedWriteEpoch = 1U;
                }
                std::vector<std::size_t> accepted;
                accepted.reserve(packetCount);
                for (std::size_t i = 0; i < packetCount; ++i) {
                    SpecCandidate& candidate = ws.parallelCandidates[i];
                    if (candidate.status == SpecCandidateStatus::Ready) {
                        ++st.parallelReadyObservations;
                        if (candidate.touchesInfinity)
                            ++st.parallelHullObservations;
                    } else if (candidate.status == SpecCandidateStatus::Skipped) {
                        ++st.parallelSkippedObservations;
                    }
                    if (candidate.status == SpecCandidateStatus::Duplicate) {
                        accepted.push_back(i);
                        continue;
                    }
                    if (candidate.status != SpecCandidateStatus::Ready ||
                        candidate.touchesInfinity)
                        continue;
                    bool independent = true;
                    for (Slot slot : candidate.writes) {
                        if (slot < selectedWriteMarks.size() &&
                            selectedWriteMarks[slot] == selectedWriteEpoch) {
                            independent = false;
                            ++st.parallelWriteMarkConflicts;
                            break;
                        }
                    }
                    if (!independent) continue;
                    for (std::size_t earlierIndex : accepted) {
                        ++st.parallelConflictPairChecks;
                        const SpecCandidate& earlier =
                            ws.parallelCandidates[earlierIndex];
                        if (earlier.status != SpecCandidateStatus::Ready)
                            continue;
                        if (replacementMayAffect(earlier,
                                                 candidate.pointIndex)) {
                            independent = false;
                            break;
                        }
                    }
                    if (!independent) continue;
                    prepareInfluence(candidate);
                    accepted.push_back(i);
                    for (Slot slot : candidate.writes) {
                        if (slot < selectedWriteMarks.size())
                            selectedWriteMarks[slot] = selectedWriteEpoch;
                    }

                }

                if (accepted.empty()) {
                    // Deterministic progress fallback for a dense local
                    // interaction: execute the first packet candidate exactly.
                    accepted.push_back(0U);
                    SpecCandidate& candidate = ws.parallelCandidates[0];
                    discoverSpecCandidate(candidate.pointIndex, lastAlive,
                                          candidate);
                    ++st.speculativeCandidatesRecomputed;
                }

                const Clock::time_point commitBegin = Clock::now();
                for (std::size_t i : accepted) {
                    const std::size_t rank = packetRanks[i];
                    if (rank >= ws.order.size() || doneRank[rank] != 0U)
                        continue;
                    const std::uint16_t cell = packetCells[i];
                    SpecCandidate& candidate = ws.parallelCandidates[i];
                    ++st.speculativeCandidates;
                    st.speculativePredicateCalls += candidate.predicateCalls;
                    st.speculativeCavityVisits += candidate.cavityVisits;
                    if (candidate.status == SpecCandidateStatus::Duplicate) {
                        ++st.pointsDuplicate;
                    } else if (candidate.status == SpecCandidateStatus::Ready) {
                        const Slot firstNew = commitPreparedCavity(
                            candidate.pointIndex, candidate.cavity,
                            candidate.shell);
                        if (firstNew != Arena::kInvalidSlot) {
                            lastAlive = firstNew;
                            binSeed[cell] = firstNew;
                        }
                        ++st.pointsInserted;
                    } else {
                        continue;
                    }
                    ++st.speculativeCandidatesReused;
                    doneRank[rank] = 1U;
                    ++binCursor[cell];
                    --pendingInColor;
                    st.serialEquivalentPredicateCalls += candidate.predicateCalls;
                    st.serialEquivalentCavityVisits += candidate.cavityVisits;
                }
                st.canonicalCommitSeconds +=
                    secondsBetween(commitBegin, Clock::now());
            }
        }
        }
    } else {
    for (Vertex orderIndex = 0; orderIndex < n; ++orderIndex) {
        const Vertex pointIndex = ws.order[orderIndex];
        if (std::find(seeded.begin(), seeded.end(), pointIndex) != seeded.end())
            continue;
        const Vec3 point = points[pointIndex];

        Slot current = lastAlive;
        if (!arena.alive(current)) current = arena.liveTail();
        Slot located = Arena::kInvalidSlot;
        const std::uint64_t maximumSteps =
            4ULL * static_cast<std::uint64_t>(arena.slotCount()) + 64ULL;
        for (std::uint64_t step = 0;
             current != Arena::kInvalidSlot && step < maximumSteps; ++step) {
            detailCount<CollectDetailed>(st.walkCellVisits);
            const auto& tet = arena.hot(current);
            if (isInfinite(current)) {
                if (conflicts(current, point)) {
                    located = current;
                    break;
                }
                current = arena.neighbour(current, 3U);
                continue;
            }
            int bestFace = -1;
            double bestValue = -1e300;
            for (int face = 0; face < 4; ++face) {
                detailCount<CollectDetailed>(st.walkOrientationTests);
                const Vec3& f0 = points[tet.vertices[kFace[face][0]]];
                const Vec3& f1 = points[tet.vertices[kFace[face][1]]];
                const Vec3& f2 = points[tet.vertices[kFace[face][2]]];
                if (predicates::orient3dFiltered(f0, f1, f2, point) !=
                    Sign::Positive)
                    continue;
                const double heuristic = orientHeur(f0, f1, f2, point);
                if (heuristic > bestValue) {
                    bestValue = heuristic;
                    bestFace = face;
                }
            }
            if (bestFace < 0) {
                located = current;
                break;
            }
            current = arena.neighbour(current,
                                       static_cast<std::uint32_t>(bestFace));
        }

        if (located != Arena::kInvalidSlot) {
            bool duplicate = false;
            for (Vertex vertex : arena.hot(located).vertices) {
                if (vertex == infinity) continue;
                const Vec3& existing = points[vertex];
                duplicate = duplicate ||
                            (existing.x == point.x && existing.y == point.y &&
                             existing.z == point.z);
            }
            if (duplicate) {
                ++st.pointsDuplicate;
                continue;
            }
        }

        // A successful walk either stopped in an infinite cell whose
        // conflict was already proven, or in the closure of a positive finite
        // tetrahedron.  The latter lies wholly inside its circumsphere, so the
        // located cell is a certified conflict seed without another insphere.
        Slot seed = located;
        if (seed == Arena::kInvalidSlot) {
            ++st.walkFallbackScans;
            for (Slot slot = 0; slot < arena.slotCount(); ++slot) {
                if (!arena.alive(slot)) continue;
                detailCount<CollectDetailed>(st.fallbackCellTests);
                if (conflicts(slot, point)) {
                    seed = slot;
                    break;
                }
            }
        }
        if (seed == Arena::kInvalidSlot) {
            ++st.pointsSkipped;
            continue;
        }

        ws.bad.clear();
        ws.stack.clear();
        ws.stack.push_back(seed);
        arena.beginEpoch();
        arena.markVisited(seed);
        while (!ws.stack.empty()) {
            const Slot slot = ws.stack.back();
            ws.stack.pop_back();
            ws.bad.push_back(slot);
            detailCount<CollectDetailed>(st.cavityCellVisits);
            for (int face = 0; face < 4; ++face) {
                detailCount<CollectDetailed>(st.cavityNeighbourVisits);
                const Slot neighbour = arena.neighbour(
                    slot, static_cast<std::uint32_t>(face));
                if (neighbour == Arena::kInvalidSlot ||
                    arena.visited(neighbour))
                    continue;
                if (!conflicts(neighbour, point)) continue;
                arena.markVisited(neighbour);
                ws.stack.push_back(neighbour);
            }
        }

        auto shellOkay = [&](Slot slot) {
            detailCount<CollectDetailed>(st.trimCellVisits);
            const auto& tet = arena.hot(slot);
            for (int face = 0; face < 4; ++face) {
                const Slot neighbour = arena.neighbour(
                    slot, static_cast<std::uint32_t>(face));
                if (neighbour != Arena::kInvalidSlot &&
                    arena.visited(neighbour))
                    continue;
                const Vertex f0 = tet.vertices[kFace[face][0]];
                const Vertex f1 = tet.vertices[kFace[face][1]];
                const Vertex f2 = tet.vertices[kFace[face][2]];
                if (f0 == infinity || f1 == infinity || f2 == infinity)
                    continue;
                detailCount<CollectDetailed>(st.trimOrientationTests);
                const Sign pointOrientation = predicates::orient3dFiltered(
                    points[f0], points[f1], points[f2], point);
                if (tet.vertices[3] != infinity) {
                    if (pointOrientation != Sign::Negative) return false;
                } else {
                    if (neighbour == Arena::kInvalidSlot) return false;
                    const std::uint32_t reference = arena.neighbourRef(
                        slot, static_cast<std::uint32_t>(face));
                    const std::uint32_t witnessFace =
                        Arena::neighbourFace(reference);
                    const auto& finiteTet = arena.hot(neighbour);
                    const Vertex witness = finiteTet.vertices[witnessFace];
                    detailCount<CollectDetailed>(st.reciprocalFaceScansAvoided);
                    detailCount<CollectDetailed>(st.trimOrientationTests);
                    const Sign witnessOrientation = predicates::orient3dFiltered(
                        points[f0], points[f1], points[f2], points[witness]);
                    if (pointOrientation == Sign::Zero ||
                        pointOrientation == witnessOrientation)
                        return false;
                }
            }
            return true;
        };

        bool touchesInfinity = false;
        for (Slot slot : ws.bad) {
            touchesInfinity = touchesInfinity || isInfinite(slot);
        }
        if (!touchesInfinity) {
            // Bowyer-Watson theorem: for an interior insertion into a valid
            // finite Delaunay complex, the connected conflict region is
            // star-shaped with respect to the inserted point.  Every finite
            // boundary face is therefore already a legal replacement face.
            // The legacy visibility trim is retained unchanged whenever the
            // symbolic hull participates, where coplanar extension rules can
            // require resurrection.
            detailCount<CollectDetailed>(st.finiteCavityCertificates);
            detailCount<CollectDetailed>(
                st.finiteCavityTrimTestsAvoided,
                static_cast<std::uint64_t>(ws.bad.size()));
        } else {
            ws.work.clear();
            ws.work.reserve(ws.bad.size());
            for (Slot slot : ws.bad) {
                if (isInfinite(slot)) ws.work.push_back(slot);
            }
            bool resurrected = false;
            while (!ws.work.empty()) {
                const Slot slot = ws.work.back();
                ws.work.pop_back();
                if (!arena.visited(slot)) continue;
                if (shellOkay(slot)) continue;
                arena.unmark(slot);
                resurrected = true;
                detailCount<CollectDetailed>(st.cavityTrimResurrections);
                for (int face = 0; face < 4; ++face) {
                    const Slot neighbour = arena.neighbour(
                        slot, static_cast<std::uint32_t>(face));
                    if (neighbour != Arena::kInvalidSlot &&
                        arena.visited(neighbour))
                        ws.work.push_back(neighbour);
                }
            }
            if (resurrected) {
                ws.stillBad.clear();
                ws.stillBad.reserve(ws.bad.size());
                for (Slot slot : ws.bad) {
                    if (arena.visited(slot)) ws.stillBad.push_back(slot);
                }
                ws.bad.swap(ws.stillBad);
            }
        }
        detailMaximum<CollectDetailed>(
            st.peakCavityCells,
            static_cast<std::uint64_t>(ws.bad.size()));
        if (ws.bad.empty()) {
            ++st.pointsSkipped;
            continue;
        }

        ws.shell.clear();
        for (Slot slot : ws.bad) {
            const auto& tet = arena.hot(slot);
            const bool finiteOwner = tet.vertices[3] != infinity;
            for (int face = 0; face < 4; ++face) {
                const Slot neighbour = arena.neighbour(
                    slot, static_cast<std::uint32_t>(face));
                if (neighbour != Arena::kInvalidSlot &&
                    arena.visited(neighbour))
                    continue;
                ws.shell.push_back(
                    {tet.vertices[kFace[face][0]],
                     tet.vertices[kFace[face][1]],
                     tet.vertices[kFace[face][2]],
                     arena.neighbourRef(slot, static_cast<std::uint32_t>(face)),
                     finiteOwner});
                detailCount<CollectDetailed>(st.shellFacesEmitted);
                if (neighbour != Arena::kInvalidSlot)
                    detailCount<CollectDetailed>(st.reciprocalFaceScansAvoided);
            }
        }
        detailMaximum<CollectDetailed>(
            st.peakShellFaces,
            static_cast<std::uint64_t>(ws.shell.size()));
        const Slot firstNew =
            commitPreparedCavity(pointIndex, ws.bad, ws.shell);
        if (firstNew != Arena::kInvalidSlot) lastAlive = firstNew;
        ++st.pointsInserted;
    }
    }
    st.insertionSeconds = secondsBetween(insertionBegin, Clock::now());

    const Clock::time_point auditBegin = Clock::now();
    bool adjacencyOkay = true;
    for (Slot slot = 0; slot < arena.slotCount() && adjacencyOkay; ++slot) {
        if (!arena.alive(slot)) continue;
        const auto& tet = arena.hot(slot);
        for (int face = 0; face < 4; ++face) {
            detailCount<CollectDetailed>(st.adjacencyAuditFaceVisits);
            const std::uint32_t reference =
                tet.neighbours[static_cast<std::size_t>(face)];
            if (reference == Arena::kNoNeighbour) continue;
            const Slot neighbour = Arena::neighbourSlot(reference);
            const std::uint32_t reciprocal = Arena::neighbourFace(reference);
            if (!arena.alive(neighbour) || reciprocal > 3U ||
                arena.neighbourRef(neighbour, reciprocal) !=
                    Arena::packNeighbour(slot,
                                         static_cast<std::uint32_t>(face))) {
                adjacencyOkay = false;
                break;
            }
            const std::array<Vertex, 3> mine{
                tet.vertices[kFace[face][0]],
                tet.vertices[kFace[face][1]],
                tet.vertices[kFace[face][2]]};
            const auto& other = arena.hot(neighbour).vertices;
            const std::array<Vertex, 3> theirs{
                other[kFace[reciprocal][0]],
                other[kFace[reciprocal][1]],
                other[kFace[reciprocal][2]]};
            if (!sameSet3(mine, theirs)) {
                adjacencyOkay = false;
                break;
            }
        }
    }
    st.auditPassed = adjacencyOkay;
    if (!adjacencyOkay)
        throw std::logic_error(
            "delaunay3D: adjacency audit failed (shared-face identity)");
    st.adjacencyAuditSeconds = secondsBetween(auditBegin, Clock::now());

    const Clock::time_point finalizationBegin = Clock::now();
#if PAACS_DELAUNAY_LINEAR_FINALIZATION
    std::size_t finiteCount = 0;
    for (Slot slot = arena.liveHead(); slot != Arena::kInvalidSlot;
         slot = arena.nextLive(slot)) {
        if (arena.hot(slot).vertices[3] != infinity) ++finiteCount;
    }
    output.reserve(finiteCount);
    ws.uniqueTable.begin(finiteCount);
    for (Slot slot = arena.liveHead(); slot != Arena::kInvalidSlot;
         slot = arena.nextLive(slot)) {
        const auto& compactVertices = arena.hot(slot).vertices;
        if (compactVertices[3] == infinity) continue;
        if (!ws.uniqueTable.insert(compactVertices, output)) {
            st.duplicateTetDetected = true;
            continue;
        }
        output.push_back(Tet{{static_cast<Index>(compactVertices[0]),
                              static_cast<Index>(compactVertices[1]),
                              static_cast<Index>(compactVertices[2]),
                              static_cast<Index>(compactVertices[3])}});
    }
#else
    std::vector<std::pair<std::array<Vertex, 4>, std::uint32_t>> keyed;
    for (Slot slot = arena.liveHead(); slot != Arena::kInvalidSlot;
         slot = arena.nextLive(slot)) {
        const auto& vertices = arena.hot(slot).vertices;
        if (vertices[3] == infinity) continue;
        keyed.push_back({canonical4(vertices),
                         static_cast<std::uint32_t>(output.size())});
        output.push_back(Tet{{static_cast<Index>(vertices[0]),
                              static_cast<Index>(vertices[1]),
                              static_cast<Index>(vertices[2]),
                              static_cast<Index>(vertices[3])}});
    }
    std::sort(keyed.begin(), keyed.end());
    std::vector<Tet> deduplicated;
    deduplicated.reserve(output.size());
    std::vector<std::uint8_t> keep(output.size(), 1U);
    for (std::size_t index = 1; index < keyed.size(); ++index) {
        if (keyed[index].first == keyed[index - 1].first) {
            keep[keyed[index].second] = 0U;
            st.duplicateTetDetected = true;
        }
    }
    for (std::size_t index = 0; index < output.size(); ++index) {
        if (keep[index] != 0U) deduplicated.push_back(output[index]);
    }
    output.swap(deduplicated);
#endif
    if (policy == DelaunayParallelPolicy::UniqueDelaunayParallel) {
        std::sort(output.begin(), output.end(), [](const Tet& a, const Tet& b) {
            std::array<Vertex,4> ak{static_cast<Vertex>(a.v[0]),
                                    static_cast<Vertex>(a.v[1]),
                                    static_cast<Vertex>(a.v[2]),
                                    static_cast<Vertex>(a.v[3])};
            std::array<Vertex,4> bk{static_cast<Vertex>(b.v[0]),
                                    static_cast<Vertex>(b.v[1]),
                                    static_cast<Vertex>(b.v[2]),
                                    static_cast<Vertex>(b.v[3])};
            return canonical4(ak) < canonical4(bk);
        });
    }
    st.outputCells = static_cast<std::uint64_t>(output.size());

    if (!output.empty()) {
        ws.usedPoints.assign(n, 0U);
        for (const Tet& tet : output) {
            for (Index vertex : tet.v)
                ws.usedPoints[static_cast<std::size_t>(vertex)] = 1U;
        }
        std::uint64_t usedCount = 0;
        for (std::uint8_t used : ws.usedPoints)
            usedCount += static_cast<std::uint64_t>(used != 0U);
        if (usedCount < st.pointsInserted) {
            st.auditPassed = false;
            throw std::logic_error(
                "delaunay3D: output audit failed (" +
                std::to_string(st.pointsInserted - usedCount) +
                " inserted point(s) not used by any finite tetrahedron)");
        }
    }
    st.finalizationSeconds =
        secondsBetween(finalizationBegin, Clock::now());

    const auto memory = arena.memory();
    st.arenaPeakLiveCells = memory.peakLiveCells;
    st.arenaMaterializedSlots = memory.materializedSlots;
    st.arenaBytes = memory.bytes;
    st.arenaHotBytes = memory.hotBytes;
    st.arenaColdBytes = memory.coldBytes;
    st.arenaTotalCellsCreated = memory.totalCellsCreated;
    st.outputBytes = output.capacity() * sizeof(Tet);
    st.scratchBytes =
        ws.normalized.capacity() * sizeof(Vec3) +
        ws.mortonKeys.capacity() * sizeof(std::uint64_t) +
        ws.order.capacity() * sizeof(Vertex) +
        ws.orderScratch.capacity() * sizeof(Vertex) +
        ws.bad.capacity() * sizeof(Slot) +
        ws.stack.capacity() * sizeof(Slot) +
        ws.work.capacity() * sizeof(Slot) +
        ws.stillBad.capacity() * sizeof(Slot) +
        ws.shell.capacity() * sizeof(ShellFace) +
        ws.pairTable.retainedBytes() + ws.uniqueTable.retainedBytes() +
        ws.usedPoints.capacity() * sizeof(std::uint8_t);
    st.runComplete = true;
    st.totalSeconds = secondsBetween(totalBegin, Clock::now());
    ws.completedRun = true;
    st.workspaceRetainedBytes = ws.retainedBytes();
    return output;
}

} // namespace

namespace tet_candidate_detail {
std::uint64_t measuredDelaunayConstructionWork(const Delaunay3DStats& s) noexcept {
    std::uint64_t total=0;
    const auto add=[&](std::uint64_t value){total=value>UINT64_MAX-total?UINT64_MAX:total+value;};
    add(s.inputPoints);
    add(s.pointsInserted);
    add(s.pointsDuplicate);
    add(s.pointsSkipped);
    add(s.walkFallbackScans);
    add(s.cavityTrimResurrections);
    add(s.bootstrapDistinctTests);
    add(s.bootstrapOrientationTests);
    add(s.walkCellVisits);
    add(s.walkOrientationTests);
    add(s.conflictQueries);
    add(s.conflictOrientationTests);
    add(s.conflictInsphereTests);
    add(s.fallbackCellTests);
    add(s.cavityCellVisits);
    add(s.cavityNeighbourVisits);
    add(s.trimCellVisits);
    add(s.trimOrientationTests);
    add(s.shellFacesEmitted);
    add(s.shellEdgePairInsertions);
    add(s.shellEdgePairProbes);
    add(s.adjacencyAuditFaceVisits);
    add(s.outputCells);
    add(s.speculativePredicateCalls);
    add(s.speculativeCavityVisits);
    add(s.parallelHazardPredicateTests);
    add(s.parallelConflictPairChecks);
    add(s.parallelWriteMarkConflicts);
    add(s.originalOrientationAuditTests);
    add(s.discardedConstructionWork);
    return total;
}
void retainDiscardedDelaunayAttempt(Delaunay3DStats& accepted,const Delaunay3DStats& discarded) noexcept {
    const auto work=measuredDelaunayConstructionWork(discarded);
    accepted.discardedConstructionWork=work>UINT64_MAX-accepted.discardedConstructionWork?
        UINT64_MAX:accepted.discardedConstructionWork+work;
    accepted.originalCoordinateRetryUsed=true;
    // Attempts are sequential; these fields describe a conservative high-water
    // estimate, not the sum of native/allocator RSS across attempts.
    accepted.arenaPeakLiveCells=std::max(accepted.arenaPeakLiveCells,discarded.arenaPeakLiveCells);
    accepted.arenaMaterializedSlots=std::max(accepted.arenaMaterializedSlots,discarded.arenaMaterializedSlots);
    accepted.arenaBytes=std::max(accepted.arenaBytes,discarded.arenaBytes);
    accepted.arenaHotBytes=std::max(accepted.arenaHotBytes,discarded.arenaHotBytes);
    accepted.arenaColdBytes=std::max(accepted.arenaColdBytes,discarded.arenaColdBytes);
    accepted.outputBytes=std::max(accepted.outputBytes,discarded.outputBytes);
    accepted.scratchBytes=std::max(accepted.scratchBytes,discarded.scratchBytes);
    accepted.workspaceRetainedBytes=std::max(accepted.workspaceRetainedBytes,discarded.workspaceRetainedBytes);
    accepted.normalizationOrderSeconds+=discarded.normalizationOrderSeconds;
    accepted.insertionSeconds+=discarded.insertionSeconds;
    accepted.adjacencyAuditSeconds+=discarded.adjacencyAuditSeconds;
    accepted.finalizationSeconds+=discarded.finalizationSeconds;
    accepted.totalSeconds+=discarded.totalSeconds;
}
std::vector<Tet> delaunay3DExactInputCoordinates(const std::vector<Vec3>& points,Delaunay3DStats* stats) {
    DelaunayWorkspace::Impl workspace;
    Delaunay3DOptions serial;serial.collectDetailedOperationCounts=true;
    Delaunay3DStats local;auto& out=stats?*stats:local;
    auto result=compactDelaunay<true,true>(points,workspace,serial,out);
    for(const auto& tet:result){
        ++out.originalOrientationAuditTests;
        if(predicates::orient3d(points[tet.v[0]],points[tet.v[1]],points[tet.v[2]],points[tet.v[3]])!=Sign::Positive)
            throw std::logic_error("delaunay3D: original-coordinate retry orientation audit failed");
    }
    out.originalCoordinateRetryUsed=true;
    return result;
}
} // namespace tet_candidate_detail

std::vector<Tet> delaunay3D(const std::vector<Vec3>& points,
                            DelaunayWorkspace& workspace,
                            const Delaunay3DOptions& options,
                            Delaunay3DStats* stats) {
    Delaunay3DStats localStats;
    Delaunay3DStats& resultStats = stats ? *stats : localStats;
    if (!workspace.impl_)
        throw std::logic_error("delaunay3D: moved-from workspace");
    try {
        const DelaunayParallelPolicy requestedPolicy =
            options.parallelPolicy != DelaunayParallelPolicy::SerialExact
                ? options.parallelPolicy
                : (options.deterministicParallel
                       ? DelaunayParallelPolicy::UniqueDelaunayParallel
                       : DelaunayParallelPolicy::SerialExact);
        const bool uniqueParallel =
            requestedPolicy == DelaunayParallelPolicy::UniqueDelaunayParallel;
        if (uniqueParallel) predicates::detail::resetDelaunayZeroSeen();
        std::vector<Tet> result;
        if (options.collectDetailedOperationCounts)
            result = compactDelaunay<true>(points, *workspace.impl_, options, resultStats);
        else
            result = compactDelaunay<false>(points, *workspace.impl_, options, resultStats);
        if (uniqueParallel && predicates::detail::delaunayZeroSeen()) {
            Delaunay3DOptions serial = options;
            serial.deterministicParallel = false;
            serial.parallelPolicy = DelaunayParallelPolicy::SerialExact;
            Delaunay3DStats parallelAttempt = resultStats;
            if (serial.collectDetailedOperationCounts)
                result = compactDelaunay<true>(points, *workspace.impl_, serial, resultStats);
            else
                result = compactDelaunay<false>(points, *workspace.impl_, serial, resultStats);
            resultStats.parallelDegeneracyFallback = true;
            resultStats.usedParallelDiscovery = true;
            resultStats.requestedThreads = options.requestedThreads;
            resultStats.packetSize = options.packetSize;
            // Retain enough discarded-work accounting to make the fallback
            // visible without contaminating the accepted serial result.
            resultStats.speculativeDiscoveries = parallelAttempt.speculativeDiscoveries;
            resultStats.parallelRounds = parallelAttempt.parallelRounds;
        }
        const auto auditBegin=Clock::now();
        bool originalPositive=true;
        for(const auto& tet:result){
            ++resultStats.originalOrientationAuditTests;
            if(predicates::orient3d(points[tet.v[0]],points[tet.v[1]],points[tet.v[2]],points[tet.v[3]])!=Sign::Positive){originalPositive=false;break;}
        }
        const auto auditSeconds=secondsBetween(auditBegin,Clock::now());
        resultStats.finalizationSeconds+=auditSeconds;resultStats.totalSeconds+=auditSeconds;
        if(!originalPositive){
            const auto discarded=resultStats;
            // Release rejected physical output before the exact-coordinate
            // attempt; reuse the caller's retained arena and ordering storage.
            std::vector<Tet>{}.swap(result);
            Delaunay3DOptions serial=options;serial.deterministicParallel=false;
            serial.parallelPolicy=DelaunayParallelPolicy::SerialExact;
            try {
                if(serial.collectDetailedOperationCounts)
                    result=compactDelaunay<true,true>(points,*workspace.impl_,serial,resultStats);
                else result=compactDelaunay<false,true>(points,*workspace.impl_,serial,resultStats);
                const auto retryAuditBegin=Clock::now();
                for(const auto& tet:result){
                    ++resultStats.originalOrientationAuditTests;
                    if(predicates::orient3d(points[tet.v[0]],points[tet.v[1]],points[tet.v[2]],points[tet.v[3]])!=Sign::Positive)
                        throw std::logic_error("delaunay3D: original-coordinate retry orientation audit failed");
                }
                const auto retryAuditSeconds=secondsBetween(retryAuditBegin,Clock::now());
                resultStats.finalizationSeconds+=retryAuditSeconds;
                resultStats.totalSeconds+=retryAuditSeconds;
            } catch (...) {
                // Retain the rejected attempt even when the exact retry runs
                // out of memory or refuses its reversible scaling range.
                tet_candidate_detail::retainDiscardedDelaunayAttempt(resultStats,discarded);
                throw;
            }
            tet_candidate_detail::retainDiscardedDelaunayAttempt(resultStats,discarded);
        }
        return result;
    } catch (const detail::CompactEnvelopeExceeded&) {
        std::vector<Tet> result = detail::delaunay3DLegacy64(points, &resultStats);
        resultStats.compactEnvelopeFallback = true;
        resultStats.workspaceRetainedBytes = workspace.retainedBytes();
        // The legacy 64-bit envelope is not the compact exact-coordinate
        // retry. It must still refuse physically nonpositive output.
        for(const auto& tet:result){
            ++resultStats.originalOrientationAuditTests;
            if(predicates::orient3d(points[tet.v[0]],points[tet.v[1]],points[tet.v[2]],points[tet.v[3]])!=Sign::Positive)
                throw std::logic_error("delaunay3D: legacy original-coordinate orientation audit failed");
        }
        return result;
    }
}

std::vector<Tet> delaunay3D(const std::vector<Vec3>& points,
                            DelaunayWorkspace& workspace,
                            Delaunay3DStats* stats) {
    Delaunay3DOptions options;
    options.collectDetailedOperationCounts = stats != nullptr;
    return delaunay3D(points, workspace, options, stats);
}

std::vector<Tet> delaunay3D(const std::vector<Vec3>& points,
                            const Delaunay3DOptions& options,
                            Delaunay3DStats* stats) {
    DelaunayWorkspace workspace;
    return delaunay3D(points, workspace, options, stats);
}

std::vector<Tet> delaunay3D(const std::vector<Vec3>& points,
                            Delaunay3DStats* stats) {
    DelaunayWorkspace workspace;
    return delaunay3D(points, workspace, stats);
}

} // namespace paacs::algo
