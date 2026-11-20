// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/TetRefinementMesher.hpp"
#include "paacs/algo/TetQualityMesher.hpp"

#include "paacs/bnd/TetRecoveryIndex.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/predicates/Predicates.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace paacs::algo {
namespace {

using NodePair = std::array<topo::NodeHandle,2>;

#if defined(PAACS_TOPO_TESTING)
struct Phase3TailGuard {
    Phase3TailGuard() noexcept { topo::detail::phase3Active()=true; }
    ~Phase3TailGuard() { stop(); }
    void stop() noexcept {
        if (active) {
            topo::detail::phase3Active()=false;
            active=false;
        }
    }
    bool active=true;
};
class TestingT4AllocationScope {
public:
    TestingT4AllocationScope() noexcept
        : previous_(t4a_detail::allocationScopeForTesting()) {
        t4a_detail::allocationScopeForTesting()=true;
    }
    ~TestingT4AllocationScope() {
        t4a_detail::allocationScopeForTesting()=previous_;
    }
    void suspend() noexcept {
        t4a_detail::allocationScopeForTesting()=false;
    }
    void resume() noexcept {
        t4a_detail::allocationScopeForTesting()=true;
    }
private:
    bool previous_=false;
};
class TestingWallScope {
public:
    explicit TestingWallScope(t3a_detail::T3AWallPhase phase) noexcept
        : phase_(phase),enabled_(
              t3a_detail::wallTimingTraceForTesting().enabled) {
        if (enabled_) start_=std::chrono::steady_clock::now();
    }
    ~TestingWallScope() {
        if (!enabled_) return;
        const auto elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now()-start_).count();
        auto& trace=t3a_detail::wallTimingTraceForTesting();
        const auto index=static_cast<std::size_t>(phase_);
        trace.nanoseconds[index]+=static_cast<std::uint64_t>(elapsed);
        ++trace.calls[index];
    }
private:
    t3a_detail::T3AWallPhase phase_;
    bool enabled_=false;
    std::chrono::steady_clock::time_point start_{};
};
#endif

bool sameClass(const mesh::Classification& a,
               const mesh::Classification& b) noexcept {
    return a.dim == b.dim && a.tag == b.tag;
}

bool finite(const Vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

bool validPolicy(const sizing::MeshSpec& spec,
                 const TetRefinementPolicy& p) noexcept {
    const auto finitePositive = [](double v) {
        return std::isfinite(v) && v > 0.0;
    };
    const auto finiteNonnegative = [](double v) {
        return std::isfinite(v) && v >= 0.0;
    };
    if (spec.order != 1 || !finitePositive(p.maxEdgeToTarget) ||
        !finitePositive(p.maxRadiusEdge) ||
        !finiteNonnegative(p.minimumCreatedTetVolume) ||
        !finiteNonnegative(p.sliverRadiusRatio) ||
        !finiteNonnegative(p.minimumSliverImprovement) ||
        p.maxCandidateAttemptsPerViolation == 0 || p.maxWorkUnits == 0 ||
        p.maxOutputNodes == 0 || p.maxOutputTets == 0 ||
        p.maxReportRows == 0 || p.maxReportTopologyHandles == 0)
        return false;
    if (p.sliverPolicy == TetSliverPolicy::BoundedInsertion &&
        (p.sliverRadiusRatio <= 0.0 || p.minimumSliverImprovement <= 0.0 ||
         p.maxSliverCandidateAttempts == 0))
        return false;
    return true;
}

bool checkedAdd(std::uint64_t a,std::uint64_t b,std::uint64_t& out) noexcept {
    if (a > (std::numeric_limits<std::uint64_t>::max)() - b) return false;
    out = a + b;
    return true;
}
bool checkedMul(std::uint64_t a,std::uint64_t b,std::uint64_t& out) noexcept {
    if (a != 0 && b > (std::numeric_limits<std::uint64_t>::max)() / a)
        return false;
    out = a * b;
    return true;
}

struct WorkBudgetStop final {};

struct WorkBudgetContext {
    std::uint64_t limit=(std::numeric_limits<std::uint64_t>::max)();
    bool counting=true;
};

thread_local WorkBudgetContext* activeWorkBudget=nullptr;

class WorkBudgetScope {
public:
    explicit WorkBudgetScope(std::uint64_t limit) noexcept
        : previous_(activeWorkBudget) {
        context_.limit=limit;
        activeWorkBudget=&context_;
    }
    ~WorkBudgetScope() { activeWorkBudget=previous_; }
    WorkBudgetScope(const WorkBudgetScope&)=delete;
    WorkBudgetScope& operator=(const WorkBudgetScope&)=delete;
    void stopCounting() noexcept { context_.counting=false; }
private:
    WorkBudgetContext context_{};
    WorkBudgetContext* previous_=nullptr;
};

void recordConstructionWork(TetRefinementMetrics& metrics,
                            std::uint64_t amount=1) {
    if (activeWorkBudget&&!activeWorkBudget->counting) return;
    if (activeWorkBudget&&
        (metrics.workUnitsConsumed>activeWorkBudget->limit||
         amount>activeWorkBudget->limit-metrics.workUnitsConsumed))
        throw WorkBudgetStop{};
    const auto maximum=(std::numeric_limits<std::uint64_t>::max)();
    metrics.workUnitsConsumed=amount>maximum-metrics.workUnitsConsumed
        ? maximum:metrics.workUnitsConsumed+amount;
}
bool fitsSize(std::uint64_t value) noexcept {
    return value <= static_cast<std::uint64_t>(
        (std::numeric_limits<std::size_t>::max)());
}

std::uint64_t workspaceAccountedBytes(
    const topo::TopologyWorkspace& workspace) noexcept {
    return static_cast<std::uint64_t>(workspace.nodes().stats().totalBytes)+
           static_cast<std::uint64_t>(workspace.cells().stats().totalBytes)+
           static_cast<std::uint64_t>(workspace.edges().stats().totalBytes)+
           static_cast<std::uint64_t>(workspace.faces().stats().totalBytes);
}

std::uint64_t flatMeshAccountedBytes(const mesh::Mesh& mesh) noexcept {
    std::uint64_t connectivity=0;
    for (Index i=0;i<mesh.cellCount();++i)
        connectivity+=static_cast<std::uint64_t>(
            mesh::nodesPerCell(mesh.cell(i).type));
    return static_cast<std::uint64_t>(mesh.nodeCount())*
               (sizeof(Vec3)+sizeof(mesh::Classification))+
           static_cast<std::uint64_t>(mesh.cellCount())*sizeof(mesh::Cell)+
           connectivity*sizeof(Index);
}

std::uint64_t orderedNodeAccountedBytes(std::size_t rows,
                                        std::size_t valueBytes) noexcept {
    // Accounted convention, not a live-heap claim: payload plus four pointer
    // words for the ordered-tree node/link/color implementation envelope.
    return static_cast<std::uint64_t>(rows)*
           static_cast<std::uint64_t>(valueBytes+4*sizeof(void*));
}

void addAccountedBytes(std::uint64_t& total,std::uint64_t amount) {
    std::uint64_t next=0;
    if (!checkedAdd(total,amount,next))
        throw std::length_error("refinement memory accounting overflow");
    total=next;
}

template<class T>
std::uint64_t vectorAccountedBytes(const std::vector<T>& rows) noexcept {
    return static_cast<std::uint64_t>(rows.capacity())*sizeof(T);
}

std::uint64_t reportAccountedBytes(const T3AReport& report) {
    std::uint64_t result=sizeof(T3AReport);
    addAccountedBytes(result,vectorAccountedBytes(report.attemptStorage));
    addAccountedBytes(result,
        vectorAccountedBytes(report.protectedWitnessStorage));
    addAccountedBytes(result,vectorAccountedBytes(report.blockedStorage));
    addAccountedBytes(result,vectorAccountedBytes(report.unattemptedStorage));
    addAccountedBytes(result,
        vectorAccountedBytes(report.committedPrefixStorage));
    addAccountedBytes(result,
        vectorAccountedBytes(report.topologyHandleStorage));
    return result;
}

std::uint64_t reportAccountedBytes(const T4AReport& report) {
    std::uint64_t result=sizeof(T4AReport);
    addAccountedBytes(result,vectorAccountedBytes(report.attemptStorage));
    addAccountedBytes(result,vectorAccountedBytes(report.blockedStorage));
    addAccountedBytes(result,vectorAccountedBytes(report.commitStorage));
    addAccountedBytes(result,vectorAccountedBytes(report.dependencyStorage));
    addAccountedBytes(result,vectorAccountedBytes(report.witnessStorage));
    addAccountedBytes(result,vectorAccountedBytes(report.topologyStorage));
    addAccountedBytes(result,vectorAccountedBytes(report.unattemptedStorage));
    addAccountedBytes(result,vectorAccountedBytes(report.tailStorage));
    return result;
}

std::array<std::uint64_t,3> coordinateBits(Vec3 p) noexcept {
    const auto bits = [](double v) {
        if (v == 0.0) v = 0.0; // canonicalize signed zero
        return std::bit_cast<std::uint64_t>(v);
    };
    return {bits(p.x),bits(p.y),bits(p.z)};
}

std::uint64_t severityBits(double value) noexcept {
    if (value == 0.0) value = 0.0;
    return std::bit_cast<std::uint64_t>(value);
}

bool vertexLess(const RefinementVertexKey& a,
                const RefinementVertexKey& b) noexcept { return a < b; }

bool violationKeyLess(const TetViolationKey& a,
                      const TetViolationKey& b) noexcept {
    if (a.kind != b.kind) return a.kind < b.kind;
    if (a.tet < b.tet) return true;
    if (b.tet < a.tet) return false;
    if (a.hasSizeEdge != b.hasSizeEdge) return a.hasSizeEdge > b.hasSizeEdge;
    return a.sizeEdge < b.sizeEdge;
}

bool queueLess(const TetViolationRow& a,const TetViolationRow& b) noexcept {
    if (a.key.kind != b.key.kind) return a.key.kind < b.key.kind;
    if (a.severityBits != b.severityBits)
        return a.severityBits > b.severityBits;
    if (violationKeyLess(a.key,b.key)) return true;
    if (violationKeyLess(b.key,a.key)) return false;
    return a.generation < b.generation;
}

struct QueuedViolation {
    TetViolationRow row{};
    topo::CellHandle cell{};
    std::array<topo::NodeHandle,4> nodes{};
    std::uint8_t nextDescriptor=0;
};

bool queueEntryLess(const QueuedViolation& a,
                    const QueuedViolation& b) noexcept {
    if (queueLess(a.row,b.row)) return true;
    if (queueLess(b.row,a.row)) return false;
    return a.cell<b.cell;
}

struct QueuedViolationLess {
    bool operator()(const QueuedViolation& a,
                    const QueuedViolation& b) const noexcept {
        return queueEntryLess(a,b);
    }
};

class ViolationHeap {
public:
    void reserve(std::size_t rowCapacity,std::size_t slotCapacity) {
        rows_.reserve(rowCapacity);
        position_.assign(slotCapacity,kMissing);
    }
    bool empty() const noexcept { return rows_.empty(); }
    std::size_t size() const noexcept { return rows_.size(); }
    std::size_t capacity() const noexcept { return rows_.capacity(); }
    std::size_t indexedSize() const noexcept { return position_.capacity(); }
    std::span<const QueuedViolation> rows() const noexcept { return rows_; }
    const QueuedViolation& top() const noexcept {
        if (rows_.empty()) std::abort();
        return rows_.front();
    }
    const QueuedViolation* find(topo::CellHandle cell) const noexcept {
        const auto slot=slotOf(cell);
        if (slot>=position_.size()||position_[slot]==kMissing) return nullptr;
        const auto index=position_[slot];
        return index<rows_.size()&&rows_[index].cell==cell
            ? &rows_[index]:nullptr;
    }

    void push(QueuedViolation value) {
        const std::size_t slot=slotOf(value.cell);
        if (slot>=position_.size()) {
            throw std::length_error("prepared queue index exhausted");
        }
        if (position_[slot]!=kMissing) {
            throw std::logic_error("duplicate queued Tet violation");
        }
        if (rows_.size()==rows_.capacity()) {
            throw std::length_error("prepared queue capacity exhausted");
        }
        rows_.push_back(std::move(value));
        position_[slot]=rows_.size()-1;
        siftUp(rows_.size()-1);
    }

    QueuedViolation pop() noexcept {
        if (rows_.empty()) std::abort();
        QueuedViolation result=rows_.front();
        position_[slotOf(result.cell)]=kMissing;
        if (rows_.size()==1) {
            rows_.pop_back();
            return result;
        }
        rows_.front()=rows_.back();
        rows_.pop_back();
        position_[slotOf(rows_.front().cell)]=0;
        siftDown(0);
        return result;
    }

    void erase(topo::CellHandle cell) noexcept {
        const std::size_t slot=slotOf(cell);
        if (slot>=position_.size()||position_[slot]==kMissing) return;
        const std::size_t index=position_[slot];
        if (rows_[index].cell!=cell) return;
        position_[slot]=kMissing;
        if (index==rows_.size()-1) {
            rows_.pop_back();
            return;
        }
        rows_[index]=rows_.back();
        rows_.pop_back();
        position_[slotOf(rows_[index].cell)]=index;
        if (index>0&&higherPriority(rows_[index],rows_[(index-1)/2]))
            siftUp(index);
        else
            siftDown(index);
    }

    bool canApply(std::span<const topo::CellHandle> removals,
                  std::span<const QueuedViolation> additions) const noexcept {
        std::size_t removedRows=0;
        for (const auto cell:removals) {
            const auto slot=slotOf(cell);
            if (slot<position_.size()&&position_[slot]!=kMissing&&
                rows_[position_[slot]].cell==cell) ++removedRows;
        }
        if (additions.size()>capacity()-(size()-removedRows)) return false;
        for (std::size_t i=0;i<additions.size();++i) {
            const auto slot=slotOf(additions[i].cell);
            if (slot>=position_.size()) return false;
            for (std::size_t j=0;j<i;++j)
                if (slotOf(additions[j].cell)==slot) return false;
            if (position_[slot]==kMissing) continue;
            const auto occupied=rows_[position_[slot]].cell;
            if (occupied==additions[i].cell) return false;
            bool removed=false;
            for (const auto cell:removals)
                removed=removed||(cell==occupied);
            if (!removed) return false;
        }
        return true;
    }

    void applyPrepared(std::span<const topo::CellHandle> removals,
                       std::span<const QueuedViolation> additions) noexcept {
        static_assert(std::is_nothrow_copy_constructible_v<QueuedViolation>);
        static_assert(std::is_nothrow_copy_assignable_v<QueuedViolation>);
        for (const auto cell:removals) erase(cell);
        for (const auto& value:additions) {
            const auto slot=slotOf(value.cell);
            if (slot>=position_.size()||position_[slot]!=kMissing||
                rows_.size()==rows_.capacity()) std::abort();
            rows_.push_back(value); // capacity and element move are prepared
            position_[slot]=rows_.size()-1;
            siftUp(rows_.size()-1);
        }
    }

private:
    static constexpr std::size_t kMissing=
        (std::numeric_limits<std::size_t>::max)();
    static std::size_t slotOf(topo::CellHandle cell) noexcept {
        if (cell.slot<0) std::abort();
        return static_cast<std::size_t>(cell.slot);
    }
    static bool higherPriority(const QueuedViolation& a,
                               const QueuedViolation& b) noexcept {
        return queueEntryLess(a,b);
    }
    void swapRows(std::size_t a,std::size_t b) noexcept {
        std::swap(rows_[a],rows_[b]);
        position_[slotOf(rows_[a].cell)]=a;
        position_[slotOf(rows_[b].cell)]=b;
    }
    void siftUp(std::size_t index) noexcept {
        while (index>0) {
            const std::size_t parent=(index-1)/2;
            if (!higherPriority(rows_[index],rows_[parent])) break;
            swapRows(index,parent);
            index=parent;
        }
    }
    void siftDown(std::size_t index) noexcept {
        for (;;) {
            const std::size_t left=index*2+1;
            if (left>=rows_.size()) break;
            const std::size_t right=left+1;
            std::size_t best=left;
            if (right<rows_.size()&&higherPriority(rows_[right],rows_[left]))
                best=right;
            if (!higherPriority(rows_[best],rows_[index])) break;
            swapRows(index,best);
            index=best;
        }
    }

    std::vector<QueuedViolation> rows_;
    std::vector<std::size_t> position_;
};

struct CurrentBlockedEntry {
    topo::CellHandle cell{};
    TetViolationRow row{};
    std::uint64_t blockedRowIndex=0;
};

class CurrentBlockedHeap {
public:
    struct Position {
        topo::CellHandle cell{};
        std::size_t heapIndex=kMissing;
        bool occupied=false;
    };

    void reserve(std::size_t rowCapacity,std::size_t slotCapacity) {
        rows_.reserve(rowCapacity);
        position_.assign(slotCapacity,Position{});
    }
    bool empty() const noexcept { return rows_.empty(); }
    std::size_t size() const noexcept { return rows_.size(); }
    std::size_t capacity() const noexcept { return rows_.capacity(); }
    std::size_t indexedSize() const noexcept { return position_.capacity(); }
    const CurrentBlockedEntry& top() const noexcept {
        if (rows_.empty()) std::abort();
        return rows_.front();
    }
    const CurrentBlockedEntry* find(topo::CellHandle cell) const noexcept {
        const auto slot=slotOf(cell);
        if (slot>=position_.size()) return nullptr;
        const auto& position=position_[slot];
        if (!position.occupied||position.cell!=cell||
            position.heapIndex>=rows_.size()||
            rows_[position.heapIndex].cell!=cell) return nullptr;
        return &rows_[position.heapIndex];
    }
    bool canInsert(const CurrentBlockedEntry& value) const noexcept {
        const auto slot=slotOf(value.cell);
        return slot<position_.size()&&!position_[slot].occupied&&
               rows_.size()<rows_.capacity();
    }
    void insertPrepared(const CurrentBlockedEntry& value) noexcept {
        static_assert(
            std::is_nothrow_copy_constructible_v<CurrentBlockedEntry>);
        if (!canInsert(value)) std::abort();
        rows_.push_back(value);
        setPosition(rows_.size()-1);
        siftUp(rows_.size()-1);
    }
    std::uint64_t erasePrepared(topo::CellHandle cell,
                                std::uint64_t expectedRow) noexcept {
        const auto slot=slotOf(cell);
        if (slot>=position_.size()) std::abort();
        const auto position=position_[slot];
        if (!position.occupied||position.cell!=cell||
            position.heapIndex>=rows_.size()) std::abort();
        const std::size_t index=position.heapIndex;
        if (rows_[index].cell!=cell||
            rows_[index].blockedRowIndex!=expectedRow) std::abort();
        const auto result=rows_[index].blockedRowIndex;
        position_[slot]=Position{};
        if (index==rows_.size()-1) {
            rows_.pop_back();
            return result;
        }
        rows_[index]=rows_.back();
        rows_.pop_back();
        setPosition(index);
        if (index>0&&higherPriority(rows_[index],rows_[(index-1)/2]))
            siftUp(index);
        else
            siftDown(index);
        return result;
    }
    std::uint64_t accountedBytes() const noexcept {
        return static_cast<std::uint64_t>(rows_.capacity())*
                   sizeof(CurrentBlockedEntry)+
               static_cast<std::uint64_t>(position_.capacity())*
                   sizeof(Position);
    }
#if defined(PAACS_TOPO_TESTING)
    std::span<const CurrentBlockedEntry> rowsForTesting() const noexcept {
        return rows_;
    }
    std::span<const Position> positionsForTesting() const noexcept {
        return position_;
    }
#endif

private:
    static constexpr std::size_t kMissing=
        (std::numeric_limits<std::size_t>::max)();
    static std::size_t slotOf(topo::CellHandle cell) noexcept {
        if (cell.slot<0) std::abort();
        return static_cast<std::size_t>(cell.slot);
    }
    static bool higherPriority(const CurrentBlockedEntry& a,
                               const CurrentBlockedEntry& b) noexcept {
        return queueLess(a.row,b.row);
    }
    void setPosition(std::size_t index) noexcept {
        const auto slot=slotOf(rows_[index].cell);
        if (slot>=position_.size()) std::abort();
        position_[slot]={rows_[index].cell,index,true};
    }
    void swapRows(std::size_t a,std::size_t b) noexcept {
        std::swap(rows_[a],rows_[b]);
        setPosition(a);
        setPosition(b);
    }
    void siftUp(std::size_t index) noexcept {
        while (index>0) {
            const auto parent=(index-1)/2;
            if (!higherPriority(rows_[index],rows_[parent])) break;
            swapRows(index,parent);
            index=parent;
        }
    }
    void siftDown(std::size_t index) noexcept {
        for (;;) {
            const auto left=index*2+1;
            if (left>=rows_.size()) break;
            const auto right=left+1;
            auto best=left;
            if (right<rows_.size()&&
                higherPriority(rows_[right],rows_[left])) best=right;
            if (!higherPriority(rows_[best],rows_[index])) break;
            swapRows(index,best);
            index=best;
        }
    }

    std::vector<CurrentBlockedEntry> rows_;
    std::vector<Position> position_;
};

class T4BlockedDependencyIndex {
public:
    struct Entry {
        T4ADependencyKey key{};
        topo::CellHandle selected{};
        std::uint64_t blockedRowIndex=0;
    };

    void prepare(std::size_t capacity) { storage_.resize(capacity); }
    std::uint64_t accountedBytes() const noexcept {
        return static_cast<std::uint64_t>(storage_.capacity())*sizeof(Entry);
    }
    bool canInsert(std::span<const T4ADependencyKey> keys,
                   topo::CellHandle selected,
                   std::uint64_t rowIndex) const noexcept {
        if (keys.size()>storage_.size()-used_) return false;
        for (const auto& key:keys) {
            const Entry sought{key,selected,rowIndex};
            const auto position=lowerBound(sought);
            if (position<used_&&same(storage_[position],sought)) return false;
        }
        return true;
    }
    void insertPrepared(std::span<const T4ADependencyKey> keys,
                        topo::CellHandle selected,
                        std::uint64_t rowIndex) noexcept {
        if (!canInsert(keys,selected,rowIndex)) std::abort();
        for (const auto& key:keys) {
            const Entry value{key,selected,rowIndex};
            const auto position=lowerBound(value);
            for (std::size_t i=used_;i>position;--i)
                storage_[i]=storage_[i-1];
            storage_[position]=value;
            ++used_;
        }
    }
    template<class Function>
    void forEach(const T4ADependencyKey& key,Function&& function) const {
        std::size_t position=lowerBoundKey(key);
        while (position<used_&&storage_[position].key==key) {
            function(storage_[position]);
            ++position;
        }
    }
    void erasePrepared(std::span<const T4ADependencyKey> keys,
                       topo::CellHandle selected,
                       std::uint64_t rowIndex) noexcept {
        for (const auto& key:keys) {
            const Entry sought{key,selected,rowIndex};
            const auto position=lowerBound(sought);
            if (position>=used_||!same(storage_[position],sought))
                std::abort();
            for (std::size_t i=position+1;i<used_;++i)
                storage_[i-1]=storage_[i];
            --used_;
            storage_[used_]=Entry{};
        }
    }
#if defined(PAACS_TOPO_TESTING)
    std::span<const Entry> rowsForTesting() const noexcept {
        return {storage_.data(),used_};
    }
#endif
private:
    static bool less(const Entry& a,const Entry& b) noexcept {
        if (a.key<b.key) return true;
        if (b.key<a.key) return false;
        if (a.selected<b.selected) return true;
        if (b.selected<a.selected) return false;
        return a.blockedRowIndex<b.blockedRowIndex;
    }
    static bool same(const Entry& a,const Entry& b) noexcept {
        return a.key==b.key&&a.selected==b.selected&&
               a.blockedRowIndex==b.blockedRowIndex;
    }
    std::size_t lowerBound(const Entry& value) const noexcept {
        std::size_t first=0,count=used_;
        while (count) {
            const auto step=count/2,middle=first+step;
            if (less(storage_[middle],value)) {
                first=middle+1;
                count-=step+1;
            } else count=step;
        }
        return first;
    }
    std::size_t lowerBoundKey(const T4ADependencyKey& key) const noexcept {
        std::size_t first=0,count=used_;
        while (count) {
            const auto step=count/2,middle=first+step;
            if (storage_[middle].key<key) {
                first=middle+1;
                count-=step+1;
            } else count=step;
        }
        return first;
    }
    std::vector<Entry> storage_;
    std::size_t used_=0;
};

void digestWord(std::uint64_t& digest,std::uint64_t value) noexcept {
    digest^=value+0x9e3779b97f4a7c15ull+(digest<<6)+(digest>>2);
}
#if defined(PAACS_TOPO_TESTING)
template<class Tag>
void digestHandle(std::uint64_t& digest,topo::Handle<Tag> value) noexcept {
    digestWord(digest,value.pool);
    digestWord(digest,static_cast<std::uint32_t>(value.slot));
    digestWord(digest,value.gen);
}
template<class Tag>
void digestHandleLogical(std::uint64_t& digest,
                         topo::Handle<Tag> value) noexcept {
    digestWord(digest,static_cast<std::uint32_t>(value.slot));
    digestWord(digest,value.gen);
}
void digestVertex(std::uint64_t& digest,RefinementVertexKey value) noexcept {
    digestWord(digest,static_cast<std::uint64_t>(value.kind));
    digestWord(digest,value.value);
}
void digestViolationKey(std::uint64_t& digest,
                        const TetViolationKey& value) noexcept {
    digestWord(digest,static_cast<std::uint64_t>(value.kind));
    for (const auto vertex:value.tet.vertices) digestVertex(digest,vertex);
    digestWord(digest,value.hasSizeEdge?1:0);
    for (const auto vertex:value.sizeEdge.vertices) digestVertex(digest,vertex);
}
void digestQueued(std::uint64_t& digest,
                  const QueuedViolation& value,
                  bool withPoolIdentity=true) noexcept {
    digestViolationKey(digest,value.row.key);
    digestWord(digest,value.row.severityBits);
    digestWord(digest,value.row.generation);
    if (withPoolIdentity) {
        digestHandle(digest,value.cell);
        for (const auto node:value.nodes) digestHandle(digest,node);
    } else {
        digestHandleLogical(digest,value.cell);
        for (const auto node:value.nodes) digestHandleLogical(digest,node);
    }
    digestWord(digest,value.nextDescriptor);
}
void digestBlocked(std::uint64_t& digest,
                   const CurrentBlockedEntry& value,
                   bool withPoolIdentity=true) noexcept {
    if (withPoolIdentity) digestHandle(digest,value.cell);
    else digestHandleLogical(digest,value.cell);
    digestViolationKey(digest,value.row.key);
    digestWord(digest,value.row.severityBits);
    digestWord(digest,value.row.generation);
    digestWord(digest,value.blockedRowIndex);
}
void digestBlockedHistory(std::uint64_t& digest,
                          const T3ABlockedViolationRow& value) noexcept {
    digestViolationKey(digest,value.violation);
    digestWord(digest,value.violationGeneration);
    digestWord(digest,value.attemptOrdinal);
    digestWord(digest,value.witnessOffset);
    digestWord(digest,value.witnessCount);
    digestWord(digest,value.invalidatedByInsertionOrdinal);
}
void recordBlockerStateDigest(
    t3a_detail::T3ABlockerDigestCheckpoint checkpoint,
    const topo::TopologyWorkspace& workspace,const ViolationHeap& queue,
    const CurrentBlockedHeap& blocked,const T3AReport& report,
    const TetRefinementMetrics& metrics) noexcept {
    auto& trace=t3a_detail::blockerDigestTraceForTesting();
    const auto checkpointBit=std::uint32_t{1}<<
        static_cast<std::uint32_t>(checkpoint);
    if (!trace.enabled||(trace.checkpointMask&checkpointBit)==0||
        trace.used>=trace.rows.size()) return;
    try {
        t3a_detail::T3ABlockerStateDigest row;
        row.checkpoint=checkpoint;
        row.topologySemantic=workspace.semanticTopologyHash();
        row.topologyStrict=workspace.strictWorkspaceHash();
        row.topologyRevision=workspace.topologyRevision();
        row.poolState=0xcbf29ce484222325ull;
        for (const auto& stats:{workspace.nodes().stats(),
                               workspace.cells().stats(),
                               workspace.edges().stats(),
                               workspace.faces().stats()}) {
            digestWord(row.poolState,stats.live);
            digestWord(row.poolState,stats.reusableFree);
            digestWord(row.poolState,stats.retired);
            digestWord(row.poolState,stats.totalCreated);
            digestWord(row.poolState,stats.materializedSlots);
        }
        row.queueStructural=0xcbf29ce484222325ull;
        for (const auto& value:queue.rows())
            digestQueued(row.queueStructural,value);
        std::vector<QueuedViolation> logicalQueue(
            queue.rows().begin(),queue.rows().end());
        std::sort(logicalQueue.begin(),logicalQueue.end(),queueEntryLess);
        row.queueLogical=0xcbf29ce484222325ull;
        for (const auto& value:logicalQueue)
            digestQueued(row.queueLogical,value,false);

        row.heapStructural=0xcbf29ce484222325ull;
        for (const auto& value:blocked.rowsForTesting())
            digestBlocked(row.heapStructural,value);
        std::vector<CurrentBlockedEntry> logicalBlocked(
            blocked.rowsForTesting().begin(),blocked.rowsForTesting().end());
        std::sort(logicalBlocked.begin(),logicalBlocked.end(),
            [](const auto& a,const auto& b) {
                if (queueLess(a.row,b.row)) return true;
                if (queueLess(b.row,a.row)) return false;
                return a.cell<b.cell;
            });
        row.heapLogical=0xcbf29ce484222325ull;
        for (const auto& value:logicalBlocked)
            digestBlocked(row.heapLogical,value,false);
        row.positionStructural=0xcbf29ce484222325ull;
        for (std::size_t i=0;i<blocked.positionsForTesting().size();++i) {
            const auto& value=blocked.positionsForTesting()[i];
            if (!value.occupied) continue;
            digestWord(row.positionStructural,i);
            digestHandle(row.positionStructural,value.cell);
            digestWord(row.positionStructural,value.heapIndex);
        }
        row.history=0xcbf29ce484222325ull;
        for (std::uint64_t i=0;i<report.blockedUsed;++i)
            digestBlockedHistory(row.history,report.blockedStorage[
                static_cast<std::size_t>(i)]);
        row.attemptsUsed=report.attemptsUsed;
        row.protectedWitnessesUsed=report.protectedWitnessesUsed;
        row.blockedUsed=report.blockedUsed;
        row.unattemptedUsed=report.unattemptedUsed;
        row.committedPrefixUsed=report.committedPrefixUsed;
        row.currentBlocked=blocked.size();
        row.reportState=0xcbf29ce484222325ull;
        digestWord(row.reportState,row.attemptsUsed);
        digestWord(row.reportState,row.protectedWitnessesUsed);
        digestWord(row.reportState,row.blockedUsed);
        digestWord(row.reportState,row.unattemptedUsed);
        digestWord(row.reportState,row.committedPrefixUsed);
        digestWord(row.reportState,report.hasTerminalViolation?1:0);
        if (report.hasTerminalViolation)
            digestViolationKey(row.reportState,report.terminalViolation);
        row.blockerMetrics=0xcbf29ce484222325ull;
        digestWord(row.blockerMetrics,metrics.blockedRowsRecorded);
        digestWord(row.blockerMetrics,metrics.blockedRowsInvalidated);
        digestWord(row.blockerMetrics,metrics.blockedRowsCurrent);
        digestWord(row.blockerMetrics,metrics.blockedRowsCurrentHighWater);
        row.structural=0xcbf29ce484222325ull;
        for (const auto value:{row.topologySemantic,row.topologyStrict,
             row.topologyRevision,row.poolState,row.queueStructural,
             row.heapStructural,row.positionStructural,row.history,
             row.reportState,row.blockerMetrics})
            digestWord(row.structural,value);
        row.logical=0xcbf29ce484222325ull;
        for (const auto value:{row.topologySemantic,row.topologyRevision,
             row.poolState,row.queueLogical,row.heapLogical,row.history,
             row.reportState,row.blockerMetrics})
            digestWord(row.logical,value);
        if (checkpoint==
            t3a_detail::T3ABlockerDigestCheckpoint::Final) {
            trace.finalLiveUsed=0;
            trace.finalLiveComplete=
                queue.rows().size()<=trace.finalLive.size();
            for (const auto& value:queue.rows()) {
                if (trace.finalLiveUsed>=trace.finalLive.size()) break;
                trace.finalLive[trace.finalLiveUsed++]={
                    value.row.key,value.row.generation};
            }
            std::sort(trace.finalLive.begin(),
                      trace.finalLive.begin()+trace.finalLiveUsed,
                [](const auto& a,const auto& b) {
                    if (violationKeyLess(a.key,b.key)) return true;
                    if (violationKeyLess(b.key,a.key)) return false;
                    return a.generation<b.generation;
                });
        }
        trace.rows[trace.used++]=row;
    } catch (...) {
        trace.enabled=false;
    }
}
#endif

class SizingLedger {
public:
    struct Entry {
        std::array<std::uint64_t,3> key{};
        Vec3 point{};
        double value=0.0;
        bool occupied=false;
    };

    void prepare(std::size_t maximumSamples) {
        if (maximumSamples>((std::numeric_limits<std::size_t>::max)()-1)/2)
            throw std::length_error("sizing ledger capacity overflow");
        slots_.assign(maximumSamples*2+1,Entry{});
        occupiedSlots_.reserve(maximumSamples);
    }
    const Entry* find(const std::array<std::uint64_t,3>& key) const noexcept {
        if (slots_.empty()) return nullptr;
        std::size_t slot=hash(key)%slots_.size();
        for (std::size_t probe=0;probe<slots_.size();++probe) {
            const auto& row=slots_[slot];
            if (!row.occupied) return nullptr;
            if (row.key==key) return &row;
            slot=(slot+1)%slots_.size();
        }
        return nullptr;
    }
    void insert(const std::array<std::uint64_t,3>& key,Vec3 point,
                double value) {
        if (slots_.empty()) prepare(64);
        std::size_t slot=hash(key)%slots_.size();
        for (std::size_t probe=0;probe<slots_.size();++probe) {
            auto& row=slots_[slot];
            if (!row.occupied) {
                row={key,point,value,true};
                occupiedSlots_.push_back(slot);
                return;
            }
            if (row.key==key) throw std::logic_error("duplicate sizing sample");
            slot=(slot+1)%slots_.size();
        }
        throw std::length_error("prepared sizing ledger exhausted");
    }
    std::span<const std::size_t> occupiedSlots() const noexcept {
        return occupiedSlots_;
    }
    const Entry& at(std::size_t slot) const noexcept { return slots_[slot]; }
    std::size_t accountedBytes() const noexcept {
        return slots_.capacity()*sizeof(Entry)+
               occupiedSlots_.capacity()*sizeof(std::size_t);
    }
private:
    static std::size_t hash(
        const std::array<std::uint64_t,3>& key) noexcept {
        std::uint64_t value=0x9e3779b97f4a7c15ull;
        for (const auto word:key) {
            value^=word+0x9e3779b97f4a7c15ull+(value<<6)+(value>>2);
            value^=value>>30;
            value*=0xbf58476d1ce4e5b9ull;
            value^=value>>27;
            value*=0x94d049bb133111ebull;
            value^=value>>31;
        }
        return static_cast<std::size_t>(value);
    }
    std::vector<Entry> slots_;
    std::vector<std::size_t> occupiedSlots_;
};

struct SizingFailure : std::exception {
    T3AReason reason = T3AReason::InvalidSizingEvaluation;

    explicit SizingFailure(T3AReason value) noexcept : reason(value) {}
};

double checkedSize(const sizing::MeshSpec& spec,Vec3 p,SizingLedger& ledger,
                   TetRefinementMetrics& metrics,bool retainSample=true) {
    double value = 0.0;
    try {
        recordConstructionWork(metrics);
        ++metrics.sizingEvaluations;
        value = spec.h(p);
    } catch (const WorkBudgetStop&) {
        throw;
    } catch (const std::bad_alloc&) {
        throw;
    } catch (...) {
        throw SizingFailure{T3AReason::SizingFieldThrew};
    }
    if (!std::isfinite(value) || value <= 0.0)
        throw SizingFailure{T3AReason::InvalidSizingEvaluation};
    const auto key = coordinateBits(p);
    const auto* found=ledger.find(key);
    if (!found) {
        if (retainSample) ledger.insert(key,p,value);
    } else if (std::bit_cast<std::uint64_t>(found->value) !=
               std::bit_cast<std::uint64_t>(value))
        throw SizingFailure{T3AReason::SizingFieldChanged};
    return value;
}

bool verifySizingLedger(const sizing::MeshSpec& spec,const SizingLedger& ledger,
                        TetRefinementMetrics& metrics,T3AReason& reason) {
    try {
        for (const auto slot:ledger.occupiedSlots()) {
            const auto& sample=ledger.at(slot);
            ++metrics.sizingEvaluations;
            const double value = spec.h(sample.point);
            if (!std::isfinite(value) || value <= 0.0) {
                reason = T3AReason::InvalidSizingEvaluation;
                return false;
            }
            if (std::bit_cast<std::uint64_t>(value) !=
                std::bit_cast<std::uint64_t>(sample.value)) {
                reason = T3AReason::SizingFieldChanged;
                return false;
            }
        }
        return true;
    } catch (const std::bad_alloc&) {
        throw;
    } catch (...) {
        reason = T3AReason::SizingFieldThrew;
        return false;
    }
}

bool verifyPreparationSizing(
    const sizing::MeshSpec& spec,const TetViolationRow& selected,
    const std::optional<T3ASizingSample>& proposal,
    TetRefinementMetrics& metrics,T3AReason& reason) {
    const auto verify=[&](const T3ASizingSample& sample) {
        recordConstructionWork(metrics);
        ++metrics.sizingEvaluations;
        const double value=spec.h(sample.point);
        if (!std::isfinite(value)||value<=0.0) {
            reason=T3AReason::InvalidSizingEvaluation;
            return false;
        }
        if (std::bit_cast<std::uint64_t>(value)!=
                std::bit_cast<std::uint64_t>(sample.value)) {
            reason=T3AReason::SizingFieldChanged;
            return false;
        }
        return true;
    };
    try {
        if (selected.hasSizeSamples)
            for (const auto& sample:selected.sizeSamples)
                if (!verify(sample)) return false;
        if (proposal&&!verify(*proposal)) return false;
        return true;
    } catch (const WorkBudgetStop&) {
        throw;
    } catch (const std::bad_alloc&) {
        throw;
    } catch (...) {
        reason=T3AReason::SizingFieldThrew;
        return false;
    }
}

double tetRadiusEdge(const std::array<Vec3,4>& p) noexcept {
    const Vec3 u=p[1]-p[0],v=p[2]-p[0],w=p[3]-p[0];
    const double triple=u.cross(v).dot(w);
    if (!std::isfinite(triple) || triple == 0.0)
        return (std::numeric_limits<double>::infinity)();
    const Vec3 numerator = u.norm2()*v.cross(w) +
                           v.norm2()*w.cross(u) +
                           w.norm2()*u.cross(v);
    const double radius = numerator.norm()/(2.0*std::abs(triple));
    double shortest=(std::numeric_limits<double>::infinity)();
    static constexpr std::array<std::array<int,2>,6> edges{{
        {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
    for (const auto& edge:edges)
        shortest=(std::min)(shortest,
            (p[static_cast<std::size_t>(edge[0])]-
             p[static_cast<std::size_t>(edge[1])]).norm());
    return radius/shortest;
}

std::optional<Vec3> tetCircumcenter(const std::array<Vec3,4>& p) noexcept {
    const Vec3 u=p[1]-p[0],v=p[2]-p[0],w=p[3]-p[0];
    const double denominator=2.0*u.cross(v).dot(w);
    const Vec3 numerator = u.norm2()*v.cross(w) +
                           v.norm2()*w.cross(u) +
                           w.norm2()*u.cross(v);
    if (!std::isfinite(denominator) || denominator == 0.0 ||
        !finite(numerator)) return std::nullopt;
    const Vec3 result=p[0]+numerator/denominator;
    if (!finite(result)) return std::nullopt;
    return result;
}

struct PendingInsertionIdentity {
    topo::NodeHandle node{};
    std::uint64_t ordinal=0;
};

RefinementVertexKey vertexKey(
    topo::NodeHandle node,const std::vector<Index>& existingGlobal,
    const std::map<topo::NodeHandle,std::uint64_t>& inserted,
    const PendingInsertionIdentity* pending=nullptr) {
    const std::size_t slot=static_cast<std::size_t>(node.slot);
    if (slot<existingGlobal.size() && existingGlobal[slot]!=kInvalidIndex)
        return {RefinementVertexKind::ExistingGlobal,
                static_cast<std::uint64_t>(existingGlobal[slot])};
    const auto found=inserted.find(node);
    if (found!=inserted.end())
        return {RefinementVertexKind::InteriorInsertion,found->second};
    if (pending&&pending->node==node)
        return {RefinementVertexKind::InteriorInsertion,pending->ordinal};
    throw std::logic_error("refinement vertex identity is unmapped");
}

TetViolationKey makeTetKey(
    const std::array<topo::NodeHandle,4>& nodes,
    const std::vector<Index>& existingGlobal,
    const std::map<topo::NodeHandle,std::uint64_t>& inserted,
    const PendingInsertionIdentity* pending=nullptr) {
    TetViolationKey key;
    for (std::size_t i=0;i<4;++i)
        key.tet.vertices[i]=vertexKey(nodes[i],existingGlobal,inserted,pending);
    std::sort(key.tet.vertices.begin(),key.tet.vertices.end(),vertexLess);
    return key;
}

std::optional<TetViolationRow> evaluateViolation(
    const std::array<topo::NodeHandle,4>& nodes,
    const std::array<Vec3,4>& p,int regionTag,
    const std::vector<Index>& existingGlobal,
    const std::map<topo::NodeHandle,std::uint64_t>& inserted,
    const sizing::MeshSpec& spec,const TetRefinementPolicy& policy,
    SizingLedger& ledger,TetRefinementMetrics& metrics,
    std::uint64_t generation,
    const PendingInsertionIdentity* pending=nullptr) {
    recordConstructionWork(metrics);
    ++metrics.violationEvaluations;
    auto base=makeTetKey(nodes,existingGlobal,inserted,pending);
    static constexpr std::array<std::array<int,2>,6> localEdges{{
        {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
    bool hasSize=false;
    double worst=0.0;
    T3AEdgeKey worstKey{};
    std::array<T3ASizingSample,3> worstSamples{};
    for (const auto& edge:localEdges) {
        const std::size_t ia=static_cast<std::size_t>(edge[0]);
        const std::size_t ib=static_cast<std::size_t>(edge[1]);
        const Vec3 mid{std::midpoint(p[ia].x,p[ib].x),
                       std::midpoint(p[ia].y,p[ib].y),
                       std::midpoint(p[ia].z,p[ib].z)};
        const double ha=checkedSize(spec,p[ia],ledger,metrics);
        const double hb=checkedSize(spec,p[ib],ledger,metrics);
        const double hm=checkedSize(spec,mid,ledger,metrics);
        const double h=(std::min)({ha,hb,hm});
        const double limit=policy.maxEdgeToTarget*h;
        const double length2=(p[ib]-p[ia]).norm2();
        const double limit2=limit*limit;
        if (!std::isfinite(length2) || !std::isfinite(limit2) || limit2<=0.0)
            throw SizingFailure{T3AReason::InvalidSizingEvaluation};
        if (length2<=limit2) continue;
        T3AEdgeKey edgeKey{{vertexKey(nodes[ia],existingGlobal,inserted,pending),
                            vertexKey(nodes[ib],existingGlobal,inserted,pending)}};
        std::sort(edgeKey.vertices.begin(),edgeKey.vertices.end(),vertexLess);
        const double severity=length2/limit2;
        if (!std::isfinite(severity))
            throw SizingFailure{T3AReason::InvalidSizingEvaluation};
        if (!hasSize || severity>worst ||
            (severity==worst && edgeKey<worstKey)) {
            hasSize=true; worst=severity; worstKey=edgeKey;
            worstSamples={T3ASizingSample{p[ia],ha},
                          T3ASizingSample{p[ib],hb},
                          T3ASizingSample{mid,hm}};
        }
    }
    TetViolationRow row;
    row.key=base;
    row.regionTag=regionTag;
    row.generation=generation;
    if (hasSize) {
        row.key.kind=TetViolationKind::Size;
        row.key.hasSizeEdge=true;
        row.key.sizeEdge=worstKey;
        row.hasSizeSamples=true;
        row.sizeSamples=worstSamples;
        row.normalizedSeverity=worst;
        row.severityBits=severityBits(worst);
        return row;
    }
    const double radiusEdge=tetRadiusEdge(p);
    if (!std::isfinite(radiusEdge))
        throw SizingFailure{T3AReason::InvalidSizingEvaluation};
    if (radiusEdge>policy.maxRadiusEdge) {
        row.key.kind=TetViolationKind::RadiusEdge;
        row.normalizedSeverity=radiusEdge/policy.maxRadiusEdge;
        row.severityBits=severityBits(row.normalizedSeverity);
        return row;
    }
    if (policy.sliverPolicy==TetSliverPolicy::BoundedInsertion) {
        const double q=mesh::tetQuality(p[0],p[1],p[2],p[3]);
        if (q<policy.sliverRadiusRatio) {
            row.key.kind=TetViolationKind::Sliver;
            row.normalizedSeverity=q==0.0
                ? (std::numeric_limits<double>::max)()
                : (std::min)((std::numeric_limits<double>::max)(),
                             policy.sliverRadiusRatio/q);
            row.severityBits=severityBits(row.normalizedSeverity);
            return row;
        }
    }
    return std::nullopt;
}

struct CandidateDescriptor {
    T3ACandidateKind kind=T3ACandidateKind::TetCircumcenter;
    bool hasFace=false;
    T3AFaceKey face{};
    std::array<std::size_t,3> localFace{};
    std::size_t opposite=0;
};

T3AFaceKey reportFaceKey(
    const std::array<topo::NodeHandle,4>& nodes,
    const std::array<std::size_t,3>& localFace,
    const std::vector<Index>& existing,
    const std::map<topo::NodeHandle,std::uint64_t>& inserted) {
    T3AFaceKey result;
    for (std::size_t i=0;i<3;++i)
        result.vertices[i]=vertexKey(nodes[localFace[i]],existing,inserted);
    std::sort(result.vertices.begin(),result.vertices.end(),vertexLess);
    return result;
}

std::array<CandidateDescriptor,5> candidateDescriptors(
    const TetViolationRow& row,const std::array<topo::NodeHandle,4>& nodes,
    const std::vector<Index>& existing,
    const std::map<topo::NodeHandle,std::uint64_t>& inserted) {
    static constexpr std::array<std::array<std::size_t,3>,4> faces{{
        {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};
    std::array<CandidateDescriptor,4> faceRows{};
    for (std::size_t i=0;i<4;++i) {
        faceRows[i].kind=T3ACandidateKind::FaceOffcenter;
        faceRows[i].hasFace=true;
        faceRows[i].localFace=faces[i];
        faceRows[i].opposite=i;
        faceRows[i].face=reportFaceKey(nodes,faces[i],existing,inserted);
    }
    std::sort(faceRows.begin(),faceRows.end(),[&](const auto& a,const auto& b) {
        if (row.key.kind==TetViolationKind::Size) {
            const auto incident=[&](const CandidateDescriptor& value) {
                bool first=false,second=false;
                for (const auto& vertex:value.face.vertices) {
                    first=first||vertex==row.key.sizeEdge.vertices[0];
                    second=second||vertex==row.key.sizeEdge.vertices[1];
                }
                return first&&second;
            };
            const bool ai=incident(a),bi=incident(b);
            if (ai!=bi) return ai;
        }
        return a.face<b.face;
    });
    std::array<CandidateDescriptor,5> result{};
    if (row.key.kind==TetViolationKind::Size) {
        std::copy(faceRows.begin(),faceRows.end(),result.begin());
        result[4]=CandidateDescriptor{};
    } else {
        result[0]=CandidateDescriptor{};
        std::copy(faceRows.begin(),faceRows.end(),result.begin()+1);
    }
    return result;
}

std::optional<Vec3> faceOffcenter(
    const CandidateDescriptor& descriptor,const std::array<Vec3,4>& p,
    const sizing::MeshSpec& spec,SizingLedger& ledger,
    TetRefinementMetrics& metrics,
    std::optional<T3ASizingSample>* consumedSample=nullptr) {
    const Vec3 a=p[descriptor.localFace[0]];
    const Vec3 b=p[descriptor.localFace[1]];
    const Vec3 c=p[descriptor.localFace[2]];
    const Vec3 u=b-a,v=c-a;
    const Vec3 n=u.cross(v);
    const double nn=n.norm2();
    if (!std::isfinite(nn)||nn==0.0) return std::nullopt;
    const Vec3 relative=(u.norm2()*v.cross(n)+v.norm2()*n.cross(u))/(2.0*nn);
    const Vec3 centre=a+relative;
    const double radius2=relative.norm2();
    if (!finite(centre)||!std::isfinite(radius2)) return std::nullopt;
    const double h=checkedSize(spec,centre,ledger,metrics);
    if (consumedSample) *consumedSample=T3ASizingSample{centre,h};
    const double height2=h*h-radius2;
    if (!std::isfinite(height2)||!(height2>0.0)) return std::nullopt;
    Vec3 normal=n/std::sqrt(nn);
    const Vec3 opposite=p[descriptor.opposite];
    if ((opposite-centre).dot(normal)<0.0) normal=normal*(-1.0);
    const Vec3 result=centre+normal*std::sqrt(height2);
    return finite(result)?std::optional<Vec3>{result}:std::nullopt;
}

NodePair nodePair(topo::NodeHandle a,topo::NodeHandle b) noexcept {
    if (b<a) std::swap(a,b);
    return {a,b};
}

using CellMap=std::map<topo::CellHandle,std::array<topo::NodeHandle,4>>;
using EdgeStars=std::map<NodePair,std::vector<topo::CellHandle>>;
using NodeStars=std::map<topo::NodeHandle,std::vector<topo::CellHandle>>;
using CoordinateMap=std::map<topo::NodeHandle,Vec3>;
using InsertedMap=std::map<topo::NodeHandle,std::uint64_t>;
using CoordinateIndex=
    std::map<std::array<std::uint64_t,3>,topo::NodeHandle>;
using CreatedCellRows=std::vector<std::pair<
    topo::CellHandle,std::array<topo::NodeHandle,4>>>;

template<class Map>
std::uint64_t orderedMapAccountedBytes(const Map& rows) noexcept {
    return orderedNodeAccountedBytes(rows.size(),sizeof(typename Map::value_type));
}

template<class Set>
std::uint64_t orderedSetAccountedBytes(const Set& rows) noexcept {
    return orderedNodeAccountedBytes(rows.size(),sizeof(typename Set::value_type));
}

template<class StarMap>
std::uint64_t starMapAccountedBytes(const StarMap& rows) {
    std::uint64_t result=orderedMapAccountedBytes(rows);
    for (const auto& [key,owners]:rows) {
        (void)key;
        addAccountedBytes(result,vectorAccountedBytes(owners));
    }
    return result;
}

template<class StarMap,class Keys>
std::uint64_t affectedStarAccountedBytes(
    const StarMap& rows,const Keys& keys) {
    std::uint64_t result=0;
    for (const auto& key:keys) {
        const auto found=rows.find(key);
        if (found==rows.end()) continue;
        addAccountedBytes(result,
            orderedNodeAccountedBytes(1,sizeof(typename StarMap::value_type)));
        addAccountedBytes(result,vectorAccountedBytes(found->second));
    }
    return result;
}

void addCellToNodeStars(NodeStars& stars,topo::CellHandle cell,
                        const std::array<topo::NodeHandle,4>& nodes) {
    for (const auto node:nodes) {
        auto& star=stars[node];
        star.push_back(cell);
        std::sort(star.begin(),star.end());
        star.erase(std::unique(star.begin(),star.end()),star.end());
    }
}

void removeCellFromNodeStars(NodeStars& stars,topo::CellHandle cell,
                             const std::array<topo::NodeHandle,4>& nodes) {
    for (const auto node:nodes) {
        const auto found=stars.find(node);
        if (found==stars.end()) throw std::logic_error("missing node star");
        auto& star=found->second;
        star.erase(std::remove(star.begin(),star.end(),cell),star.end());
        if (star.empty()) stars.erase(found);
    }
}

void addCellToStars(EdgeStars& stars,topo::CellHandle cell,
                    const std::array<topo::NodeHandle,4>& nodes) {
    static constexpr std::array<std::array<int,2>,6> edges{{
        {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
    for (const auto& edge:edges) {
        auto& star=stars[nodePair(nodes[static_cast<std::size_t>(edge[0])],
                                  nodes[static_cast<std::size_t>(edge[1])])];
        star.push_back(cell);
        std::sort(star.begin(),star.end());
        star.erase(std::unique(star.begin(),star.end()),star.end());
    }
}

void removeCellFromStars(EdgeStars& stars,topo::CellHandle cell,
                         const std::array<topo::NodeHandle,4>& nodes) {
    static constexpr std::array<std::array<int,2>,6> edges{{
        {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
    for (const auto& edge:edges) {
        const auto key=nodePair(nodes[static_cast<std::size_t>(edge[0])],
                                nodes[static_cast<std::size_t>(edge[1])]);
        auto found=stars.find(key);
        if (found==stars.end()) throw std::logic_error("missing edge star");
        auto& star=found->second;
        star.erase(std::remove(star.begin(),star.end(),cell),star.end());
        if (star.empty()) stars.erase(found);
    }
}

struct BoundaryFace {
    std::array<topo::NodeHandle,3> directed{};
    topo::NodeHandle opposite{};
};

using FaceKey=std::array<topo::NodeHandle,3>;
FaceKey faceKey(std::array<topo::NodeHandle,3> face) noexcept {
    std::sort(face.begin(),face.end());
    return face;
}

using FaceOwners=std::map<FaceKey,std::vector<topo::CellHandle>>;
static constexpr std::array<std::array<std::size_t,3>,4> kLocalTetFaces{{
    {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};

void addCellToFaces(FaceOwners& owners,topo::CellHandle cell,
                    const std::array<topo::NodeHandle,4>& nodes) {
    for (const auto& local:kLocalTetFaces) {
        const auto key=faceKey({nodes[local[0]],nodes[local[1]],nodes[local[2]]});
        auto& row=owners[key];
        row.push_back(cell);
        std::sort(row.begin(),row.end());
        row.erase(std::unique(row.begin(),row.end()),row.end());
        if (row.size()>2) throw std::logic_error("nonmanifold face owners");
    }
}
void removeCellFromFaces(FaceOwners& owners,topo::CellHandle cell,
                         const std::array<topo::NodeHandle,4>& nodes) {
    for (const auto& local:kLocalTetFaces) {
        const auto key=faceKey({nodes[local[0]],nodes[local[1]],nodes[local[2]]});
        auto found=owners.find(key);
        if (found==owners.end()) throw std::logic_error("missing face owners");
        auto& row=found->second;
        row.erase(std::remove(row.begin(),row.end(),cell),row.end());
        if (row.empty()) owners.erase(found);
    }
}

struct T4TopologyNodeKey {
    bool staged=false;
    topo::NodeHandle existing{};

    friend bool operator==(const T4TopologyNodeKey& a,
                           const T4TopologyNodeKey& b) noexcept {
        return a.staged==b.staged&&a.existing==b.existing;
    }
    friend bool operator<(const T4TopologyNodeKey& a,
                          const T4TopologyNodeKey& b) noexcept {
        if (a.staged!=b.staged) return a.staged<b.staged;
        return a.existing<b.existing;
    }
};

std::optional<std::uint64_t> t4TopologyTransitionCount(
    const topo::TetRefinementCavitySpec& cavity,const CellMap& cells) {
    using EdgeKey=std::array<T4TopologyNodeKey,2>;
    using LocalFaceKey=std::array<T4TopologyNodeKey,3>;
    std::uint64_t tetCount=0,edgeCapacity=0,faceCapacity=0;
    if (!checkedAdd(cavity.removed.size(),cavity.added.size(),tetCount)||
        !checkedMul(tetCount,6,edgeCapacity)||
        !checkedMul(tetCount,4,faceCapacity)||
        !fitsSize(edgeCapacity)||!fitsSize(faceCapacity))
        return std::nullopt;
    std::vector<EdgeKey> edges;
    std::vector<LocalFaceKey> faces;
    edges.reserve(static_cast<std::size_t>(edgeCapacity));
    faces.reserve(static_cast<std::size_t>(faceCapacity));
    static constexpr std::array<std::array<std::size_t,2>,6> localEdges{{
        {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
    const auto append=[&](std::array<T4TopologyNodeKey,4> nodes) {
        for (const auto& local:localEdges) {
            EdgeKey key{nodes[local[0]],nodes[local[1]]};
            std::sort(key.begin(),key.end());
            edges.push_back(key);
        }
        for (const auto& local:kLocalTetFaces) {
            LocalFaceKey key{
                nodes[local[0]],nodes[local[1]],nodes[local[2]]};
            std::sort(key.begin(),key.end());
            faces.push_back(key);
        }
    };
    for (const auto cell:cavity.removed) {
        const auto found=cells.find(cell);
        if (found==cells.end()) return std::nullopt;
        std::array<T4TopologyNodeKey,4> nodes{};
        for (std::size_t i=0;i<4;++i) nodes[i].existing=found->second[i];
        append(nodes);
    }
    for (const auto& cell:cavity.added) {
        std::array<T4TopologyNodeKey,4> nodes{};
        for (std::size_t i=0;i<4;++i) {
            const auto& ref=cell.vertices[i];
            nodes[i].staged=
                ref.kind==topo::TetRefinementNodeRefKind::Staged;
            nodes[i].existing=nodes[i].staged?topo::NodeHandle{}:ref.existing;
        }
        append(nodes);
    }
    const auto canonicalize=[](auto& values) {
        std::sort(values.begin(),values.end());
        values.erase(std::unique(values.begin(),values.end()),values.end());
    };
    canonicalize(edges);
    canonicalize(faces);
    std::uint64_t result=1;
    if (!checkedAdd(result,tetCount,result)||
        !checkedAdd(result,edges.size(),result)||
        !checkedAdd(result,faces.size(),result))
        return std::nullopt;
    return result;
}

bnd::TetAabb pointAabb(Vec3 point) noexcept { return {point,point}; }

bnd::TetAabb tetAabb(const std::array<topo::NodeHandle,4>& nodes,
                     const std::map<topo::NodeHandle,Vec3>& coordinates) {
    const auto first=coordinates.at(nodes[0]);
    bnd::TetAabb result{first,first};
    for (std::size_t i=1;i<4;++i) {
        const auto p=coordinates.at(nodes[i]);
        result.lo.x=(std::min)(result.lo.x,p.x);
        result.lo.y=(std::min)(result.lo.y,p.y);
        result.lo.z=(std::min)(result.lo.z,p.z);
        result.hi.x=(std::max)(result.hi.x,p.x);
        result.hi.y=(std::max)(result.hi.y,p.y);
        result.hi.z=(std::max)(result.hi.z,p.z);
    }
    return result;
}

bool overlaps(const bnd::TetAabb& a,const bnd::TetAabb& b) noexcept {
    return a.lo.x<=b.hi.x&&b.lo.x<=a.hi.x&&
           a.lo.y<=b.hi.y&&b.lo.y<=a.hi.y&&
           a.lo.z<=b.hi.z&&b.lo.z<=a.hi.z;
}

class RefinementLocationIndex {
public:
    struct CreatedRow {
        topo::CellHandle cell{};
        bnd::TetAabb bounds{};
    };

    void prepare(std::size_t maximumRows,std::size_t cellSlots) {
        if (maximumRows==0||maximumRows>
                (std::numeric_limits<std::size_t>::max)()/2)
            throw std::length_error("location tree capacity overflow");
        nodes_.resize(maximumRows*2);
        byCellSlot_.assign(cellSlots,kMissing);
        for (std::size_t i=0;i<nodes_.size();++i)
            nodes_[i].next=i+1<nodes_.size()?i+1:kMissing;
        freeHead_=nodes_.empty()?kMissing:0;
        freeCount_=nodes_.size();
        maximumLeaves_=maximumRows;
    }

    std::vector<topo::CellHandle> query(
        const bnd::DynamicTetAabbIndex& base,const bnd::TetAabb& bounds,
        const topo::TopologyWorkspace& workspace,
        std::uint64_t& indexVisits) const {
        auto result=base.query(bounds,&indexVisits);
        result.erase(std::remove_if(result.begin(),result.end(),
            [&](const auto cell){return !workspace.cells().isCurrent(cell);}),
            result.end());
        queryOverlay(root_,bounds,result,indexVisits);
        std::sort(result.begin(),result.end(),[&](const auto a,const auto b){
            const auto as=workspace.cells().seqOf(a);
            const auto bs=workspace.cells().seqOf(b);
            return as<bs||(as==bs&&a<b);
        });
        result.erase(std::unique(result.begin(),result.end()),result.end());
        return result;
    }

    bool canApply(std::span<const topo::CellHandle> removals,
                  std::span<const CreatedRow> creations) const noexcept {
        std::size_t removedLeaves=0;
        for (const auto cell:removals) {
            const auto slot=cellSlot(cell);
            if (slot>=byCellSlot_.size()) return false;
            const auto node=byCellSlot_[slot];
            if (node!=kMissing) {
                if (nodes_[node].cell!=cell) return false;
                ++removedLeaves;
            }
        }
        for (std::size_t i=0;i<creations.size();++i) {
            const auto slot=cellSlot(creations[i].cell);
            if (slot>=byCellSlot_.size()||!finite(creations[i].bounds.lo)||
                !finite(creations[i].bounds.hi)) return false;
            for (std::size_t j=0;j<i;++j)
                if (cellSlot(creations[j].cell)==slot) return false;
            const auto leaf=byCellSlot_[slot];
            if (leaf==kMissing) continue;
            bool removed=false;
            for (const auto cell:removals)
                removed=removed||(cell==nodes_[leaf].cell);
            if (!removed) return false;
        }
        if (removedLeaves>activeLeaves_) return false;
        const auto finalLeaves=activeLeaves_-removedLeaves+creations.size();
        if (finalLeaves>maximumLeaves_) return false;
        const auto requiredNodes=finalLeaves==0?0:finalLeaves*2-1;
        if (requiredNodes>nodes_.size()) return false;
        return true;
    }

    std::uint64_t preparedUpdateVisits(
        std::span<const topo::CellHandle> removals,
        std::span<const CreatedRow> creations) const noexcept {
        std::uint64_t changed=static_cast<std::uint64_t>(creations.size());
        for (const auto cell:removals) {
            const auto slot=cellSlot(cell);
            if (slot<byCellSlot_.size()&&byCellSlot_[slot]!=kMissing&&
                nodes_[byCellSlot_[slot]].cell==cell) ++changed;
        }
        const auto logarithmicBound=static_cast<std::uint64_t>(
            std::bit_width(maximumLeaves_))+1;
        // One sibling descent plus one AVL refit/balance ascent is bounded
        // by four visits per balanced-tree level (including rotations).
        return changed*4*logarithmicBound;
    }

    void applyPrepared(std::span<const topo::CellHandle> removals,
                       std::span<const CreatedRow> creations,
                       std::uint64_t& updateVisits) noexcept {
        for (const auto cell:removals) {
            const auto slot=cellSlot(cell);
            if (slot>=byCellSlot_.size()) std::abort();
            const auto leaf=byCellSlot_[slot];
            if (leaf==kMissing) continue; // cell belonged to immutable base BVH
            if (nodes_[leaf].cell!=cell) std::abort();
            removeLeaf(leaf,updateVisits);
            byCellSlot_[slot]=kMissing;
            releaseNode(leaf);
            --activeLeaves_;
        }
        for (const auto& created:creations) {
            const auto slot=cellSlot(created.cell);
            if (slot>=byCellSlot_.size()||byCellSlot_[slot]!=kMissing)
                std::abort();
            const auto leafIndex=acquireNode();
            auto& leaf=nodes_[leafIndex];
            leaf.bounds=created.bounds;
            leaf.cell=created.cell;
            leaf.height=0;
            insertLeaf(leafIndex,updateVisits);
            byCellSlot_[slot]=leafIndex;
            ++activeLeaves_;
        }
    }

    std::size_t accountedBytes() const noexcept {
        return nodes_.capacity()*sizeof(Node)+
               byCellSlot_.capacity()*sizeof(std::size_t);
    }

#if defined(PAACS_TOPO_TESTING)
    std::vector<topo::CellHandle> queryOverlayForTesting(
        const bnd::TetAabb& bounds,std::uint64_t& visits) const {
        std::vector<topo::CellHandle> result;
        queryOverlay(root_,bounds,result,visits);
        std::sort(result.begin(),result.end());
        return result;
    }
    std::uint64_t stateDigestForTesting() const noexcept {
        std::uint64_t digest=0xcbf29ce484222325ull;
        const auto addPoint=[&](Vec3 point) noexcept {
            digestWord(digest,std::bit_cast<std::uint64_t>(point.x));
            digestWord(digest,std::bit_cast<std::uint64_t>(point.y));
            digestWord(digest,std::bit_cast<std::uint64_t>(point.z));
        };
        digestWord(digest,root_);
        digestWord(digest,freeHead_);
        digestWord(digest,freeCount_);
        digestWord(digest,activeLeaves_);
        digestWord(digest,maximumLeaves_);
        for (const auto& node:nodes_) {
            addPoint(node.bounds.lo);
            addPoint(node.bounds.hi);
            digestHandle(digest,node.cell);
            digestWord(digest,node.parent);
            digestWord(digest,node.left);
            digestWord(digest,node.right);
            digestWord(digest,node.next);
            digestWord(digest,static_cast<std::uint32_t>(node.height));
        }
        for (const auto slot:byCellSlot_) digestWord(digest,slot);
        return digest;
    }
#endif

private:
    struct Node {
        bnd::TetAabb bounds{};
        topo::CellHandle cell{};
        std::size_t parent=kMissing,left=kMissing,right=kMissing,next=kMissing;
        int height=-1;
    };

    static constexpr std::size_t kMissing=
        (std::numeric_limits<std::size_t>::max)();

    static std::size_t cellSlot(topo::CellHandle cell) noexcept {
        if (cell.slot<0) std::abort();
        return static_cast<std::size_t>(cell.slot);
    }

    static bnd::TetAabb joined(const bnd::TetAabb& a,
                               const bnd::TetAabb& b) noexcept {
        return {{(std::min)(a.lo.x,b.lo.x),(std::min)(a.lo.y,b.lo.y),
                 (std::min)(a.lo.z,b.lo.z)},
                {(std::max)(a.hi.x,b.hi.x),(std::max)(a.hi.y,b.hi.y),
                 (std::max)(a.hi.z,b.hi.z)}};
    }

    static long double measure(const bnd::TetAabb& bounds) noexcept {
        const long double dx=static_cast<long double>(bounds.hi.x)-bounds.lo.x;
        const long double dy=static_cast<long double>(bounds.hi.y)-bounds.lo.y;
        const long double dz=static_cast<long double>(bounds.hi.z)-bounds.lo.z;
        const long double value=dx*dy+dy*dz+dz*dx;
        return std::isfinite(value)?value:
            (std::numeric_limits<long double>::max)();
    }

    bool isLeaf(std::size_t node) const noexcept {
        return nodes_[node].left==kMissing;
    }

    std::size_t acquireNode() noexcept {
        if (freeHead_==kMissing||freeCount_==0) std::abort();
        const auto result=freeHead_;
        freeHead_=nodes_[result].next;
        nodes_[result]=Node{};
        --freeCount_;
        return result;
    }

    void releaseNode(std::size_t node) noexcept {
        nodes_[node]=Node{};
        nodes_[node].next=freeHead_;
        freeHead_=node;
        ++freeCount_;
    }

    void refitNode(std::size_t node) noexcept {
        auto& current=nodes_[node];
        if (isLeaf(node)) {
            current.height=0;
            return;
        }
        const auto& left=nodes_[current.left];
        const auto& right=nodes_[current.right];
        current.bounds=joined(left.bounds,right.bounds);
        current.height=1+(std::max)(left.height,right.height);
    }

    std::size_t balance(std::size_t a,std::uint64_t& visits) noexcept {
        ++visits;
        auto& A=nodes_[a];
        if (isLeaf(a)||A.height<2) return a;
        const auto b=A.left,c=A.right;
        const int imbalance=nodes_[c].height-nodes_[b].height;
        if (imbalance>1) {
            const auto f=nodes_[c].left,g=nodes_[c].right;
            nodes_[c].left=a;
            nodes_[c].parent=A.parent;
            A.parent=c;
            if (nodes_[c].parent==kMissing) root_=c;
            else {
                auto& parent=nodes_[nodes_[c].parent];
                if (parent.left==a) parent.left=c;
                else if (parent.right==a) parent.right=c;
                else std::abort();
            }
            if (nodes_[f].height>nodes_[g].height) {
                nodes_[c].right=f;
                A.right=g;
                nodes_[g].parent=a;
                nodes_[f].parent=c;
            } else {
                nodes_[c].right=g;
                A.right=f;
                nodes_[f].parent=a;
                nodes_[g].parent=c;
            }
            refitNode(a);
            refitNode(c);
            return c;
        }
        if (imbalance< -1) {
            const auto d=nodes_[b].left,e=nodes_[b].right;
            nodes_[b].left=a;
            nodes_[b].parent=A.parent;
            A.parent=b;
            if (nodes_[b].parent==kMissing) root_=b;
            else {
                auto& parent=nodes_[nodes_[b].parent];
                if (parent.left==a) parent.left=b;
                else if (parent.right==a) parent.right=b;
                else std::abort();
            }
            if (nodes_[d].height>nodes_[e].height) {
                nodes_[b].right=d;
                A.left=e;
                nodes_[e].parent=a;
                nodes_[d].parent=b;
            } else {
                nodes_[b].right=e;
                A.left=d;
                nodes_[d].parent=a;
                nodes_[e].parent=b;
            }
            refitNode(a);
            refitNode(b);
            return b;
        }
        return a;
    }

    void refitUpward(std::size_t node,std::uint64_t& visits) noexcept {
        while (node!=kMissing) {
            node=balance(node,visits);
            refitNode(node);
            node=nodes_[node].parent;
        }
    }

    void insertLeaf(std::size_t leaf,std::uint64_t& visits) noexcept {
        if (root_==kMissing) {
            root_=leaf;
            nodes_[leaf].parent=kMissing;
            return;
        }
        std::size_t sibling=root_;
        while (!isLeaf(sibling)) {
            ++visits;
            const auto& current=nodes_[sibling];
            const auto combined=joined(current.bounds,nodes_[leaf].bounds);
            const auto inheritance=measure(combined)-measure(current.bounds);
            const auto childCost=[&](std::size_t child) {
                const auto expanded=joined(nodes_[child].bounds,
                                           nodes_[leaf].bounds);
                return (isLeaf(child)?measure(expanded):
                    measure(expanded)-measure(nodes_[child].bounds))+inheritance;
            };
            const auto leftCost=childCost(current.left);
            const auto rightCost=childCost(current.right);
            const auto parentCost=2*measure(combined);
            if (parentCost<leftCost&&parentCost<rightCost) break;
            sibling=leftCost<rightCost?current.left:
                    rightCost<leftCost?current.right:
                    (std::min)(current.left,current.right);
        }
        const auto oldParent=nodes_[sibling].parent;
        const auto newParent=acquireNode();
        nodes_[newParent].parent=oldParent;
        nodes_[newParent].left=sibling;
        nodes_[newParent].right=leaf;
        nodes_[newParent].height=nodes_[sibling].height+1;
        nodes_[newParent].bounds=joined(
            nodes_[sibling].bounds,nodes_[leaf].bounds);
        nodes_[sibling].parent=newParent;
        nodes_[leaf].parent=newParent;
        if (oldParent==kMissing) root_=newParent;
        else if (nodes_[oldParent].left==sibling)
            nodes_[oldParent].left=newParent;
        else if (nodes_[oldParent].right==sibling)
            nodes_[oldParent].right=newParent;
        else std::abort();
        refitUpward(newParent,visits);
    }

    void removeLeaf(std::size_t leaf,std::uint64_t& visits) noexcept {
        if (leaf==root_) {
            root_=kMissing;
            return;
        }
        const auto parent=nodes_[leaf].parent;
        const auto grand=nodes_[parent].parent;
        const auto sibling=nodes_[parent].left==leaf?
            nodes_[parent].right:nodes_[parent].left;
        if (grand==kMissing) {
            root_=sibling;
            nodes_[sibling].parent=kMissing;
        } else {
            if (nodes_[grand].left==parent) nodes_[grand].left=sibling;
            else if (nodes_[grand].right==parent) nodes_[grand].right=sibling;
            else std::abort();
            nodes_[sibling].parent=grand;
        }
        releaseNode(parent);
        nodes_[leaf].parent=kMissing;
        if (grand!=kMissing) refitUpward(grand,visits);
    }

    void queryOverlay(std::size_t node,const bnd::TetAabb& bounds,
                      std::vector<topo::CellHandle>& result,
                      std::uint64_t& visits) const {
        if (node==kMissing) return;
        ++visits;
        const auto& current=nodes_[node];
        if (!overlaps(current.bounds,bounds)) return;
        if (isLeaf(node)) {
            result.push_back(current.cell);
            return;
        }
        queryOverlay(current.left,bounds,result,visits);
        queryOverlay(current.right,bounds,result,visits);
    }

    std::vector<Node> nodes_;
    std::vector<std::size_t> byCellSlot_;
    std::size_t root_=kMissing,freeHead_=kMissing,freeCount_=0,
                activeLeaves_=0,maximumLeaves_=0;
};

std::uint64_t completeContextAccountedBytes(
    const CoordinateMap& coordinates,const InsertedMap& inserted,
    const CoordinateIndex& coordinateIndex,const CellMap& cells,
    const EdgeStars& edgeStars,const FaceOwners& faceOwners,
    const NodeStars& nodeStars,const ViolationHeap& queue,
    const CurrentBlockedHeap& blocked,
    const SizingLedger& ledger,const RefinementLocationIndex& location) {
    std::uint64_t result=0;
    addAccountedBytes(result,orderedMapAccountedBytes(coordinates));
    addAccountedBytes(result,orderedMapAccountedBytes(inserted));
    addAccountedBytes(result,orderedMapAccountedBytes(coordinateIndex));
    addAccountedBytes(result,orderedMapAccountedBytes(cells));
    addAccountedBytes(result,starMapAccountedBytes(edgeStars));
    addAccountedBytes(result,starMapAccountedBytes(faceOwners));
    addAccountedBytes(result,starMapAccountedBytes(nodeStars));
    addAccountedBytes(result,static_cast<std::uint64_t>(queue.capacity())*
        sizeof(QueuedViolation));
    addAccountedBytes(result,static_cast<std::uint64_t>(queue.indexedSize())*
        sizeof(std::size_t));
    addAccountedBytes(result,blocked.accountedBytes());
    addAccountedBytes(result,ledger.accountedBytes());
    addAccountedBytes(result,location.accountedBytes());
    return result;
}

std::uint64_t recoveryContextAccountedBytes(
    const bnd::TetRecoveryContext& context) {
    std::uint64_t result=sizeof(bnd::TetRecoveryContext);
    addAccountedBytes(result,vectorAccountedBytes(context.denseToNode));
    addAccountedBytes(result,vectorAccountedBytes(context.plainVertices));
    addAccountedBytes(result,orderedMapAccountedBytes(context.nodeToDense));
    addAccountedBytes(result,orderedMapAccountedBytes(context.cells));
    addAccountedBytes(result,orderedMapAccountedBytes(context.edges));
    addAccountedBytes(result,orderedMapAccountedBytes(context.faces));
    addAccountedBytes(result,vectorAccountedBytes(context.nodeUseCount));
    addAccountedBytes(result,vectorAccountedBytes(context.nodeStars));
    for (const auto& star:context.nodeStars)
        addAccountedBytes(result,vectorAccountedBytes(star));
    addAccountedBytes(result,vectorAccountedBytes(context.protectedEdges));
    addAccountedBytes(result,vectorAccountedBytes(context.protectedFaces));
    addAccountedBytes(result,context.index.accountedBytes());
    addAccountedBytes(result,vectorAccountedBytes(context.deltaCells));
    return result;
}

std::uint64_t t4RetainedAccountedBytes(
    const T4AMutableArtifact& artifact,
    const T4AT3SourceReceipt& receipt) {
    std::uint64_t result=sizeof(T4AMutableArtifact)+
                         sizeof(T4AT3SourceReceipt)+sizeof(T4AResult);
    addAccountedBytes(result,orderedMapAccountedBytes(
        artifact.constraintState.edgeBindings));
    addAccountedBytes(result,orderedMapAccountedBytes(
        artifact.constraintState.bindings));
    addAccountedBytes(result,vectorAccountedBytes(
        artifact.constraintState.canonicalEdgeSchedule));
    addAccountedBytes(result,vectorAccountedBytes(
        artifact.constraintState.canonicalSchedule));
    addAccountedBytes(result,vectorAccountedBytes(
        artifact.existingGlobalByAuthoritativeNodeSlot));
    addAccountedBytes(result,vectorAccountedBytes(
        artifact.sourceGlobalByAuthoritativeNodeSlot));
    addAccountedBytes(result,vectorAccountedBytes(
        artifact.sourceGlobalByAuthoritativeCellSlot));
    addAccountedBytes(result,vectorAccountedBytes(
        artifact.publishedGlobalByAuthoritativeNodeSlot));
    addAccountedBytes(result,vectorAccountedBytes(
        artifact.publishedGlobalByAuthoritativeCellSlot));
    addAccountedBytes(result,vectorAccountedBytes(artifact.immutableBoundary));
    addAccountedBytes(result,orderedMapAccountedBytes(
        artifact.vertexKeyByHandle));
    if (receipt.sourceReport)
        addAccountedBytes(result,reportAccountedBytes(*receipt.sourceReport));
    if (receipt.sourceReceipt) {
        addAccountedBytes(result,sizeof(BoundaryConformingV3Receipt));
        addAccountedBytes(result,
            receipt.sourceReceipt->memoryMetrics.finalRetainedOverlapBytes);
    }
    return result;
}

#if defined(PAACS_TOPO_TESTING)
void digestT4QualityKey(std::uint64_t& digest,
                        const T4AQualityKey& value) noexcept {
    digestWord(digest,value.radiusRatioBits);
    for (const auto vertex:value.tet.vertices) digestVertex(digest,vertex);
    digestWord(digest,value.generation);
}
void digestT4Dependency(std::uint64_t& digest,
                        const T4ADependencyKey& value) noexcept {
    digestWord(digest,static_cast<std::uint64_t>(value.kind));
    digestWord(digest,value.vertexCount);
    for (const auto vertex:value.vertices) digestVertex(digest,vertex);
}
void digestT4Report(std::uint64_t& digest,const T4AReport& report) noexcept {
    digestWord(digest,report.attemptsUsed);
    digestWord(digest,report.blockedUsed);
    digestWord(digest,report.commitsUsed);
    digestWord(digest,report.dependenciesUsed);
    digestWord(digest,report.witnessesUsed);
    digestWord(digest,report.topologyUsed);
    digestWord(digest,report.unattemptedUsed);
    digestWord(digest,report.tailUsed);
    digestWord(digest,report.hasTerminal?1:0);
    for (std::uint64_t i=0;i<report.attemptsUsed;++i) {
        const auto& row=report.attemptStorage[static_cast<std::size_t>(i)];
        digestWord(digest,row.attemptOrdinal);
        digestWord(digest,row.operationOrdinal);
        digestT4QualityKey(digest,row.selected);
        digestWord(digest,std::bit_cast<std::uint64_t>(row.proposedCoordinate.x));
        digestWord(digest,std::bit_cast<std::uint64_t>(row.proposedCoordinate.y));
        digestWord(digest,std::bit_cast<std::uint64_t>(row.proposedCoordinate.z));
        for (const auto value:{row.removedCells,row.createdCells,
             row.dependencyOffset,row.dependencyCount,row.witnessOffset,
             row.witnessCount,std::bit_cast<std::uint64_t>(row.selectedQuality),
             std::bit_cast<std::uint64_t>(row.createdMinimumQuality),
             static_cast<std::uint64_t>(row.disposition),
             static_cast<std::uint64_t>(row.reason),row.workBefore,
             row.workAfter}) digestWord(digest,value);
    }
    for (std::uint64_t i=0;i<report.blockedUsed;++i) {
        const auto& row=report.blockedStorage[static_cast<std::size_t>(i)];
        digestT4QualityKey(digest,row.selected);
        for (const auto value:{static_cast<std::uint64_t>(row.reason),
             row.attemptOrdinal,row.dependencyOffset,row.dependencyCount,
             row.witnessOffset,row.witnessCount,row.invalidatedByOperation})
            digestWord(digest,value);
    }
    for (std::uint64_t i=0;i<report.commitsUsed;++i) {
        const auto& row=report.commitStorage[static_cast<std::size_t>(i)];
        digestWord(digest,row.operationOrdinal);
        digestWord(digest,row.attemptOrdinal);
        digestT4QualityKey(digest,row.selected);
        digestHandle(digest,row.predictedNode);
        digestHandle(digest,row.actualNode);
        for (const auto value:{row.topologyOffset,row.topologyCount,
             std::bit_cast<std::uint64_t>(row.selectedQuality),
             std::bit_cast<std::uint64_t>(row.createdMinimumQuality),
             row.postCommitStamp.revision}) digestWord(digest,value);
    }
    for (std::uint64_t i=0;i<report.dependenciesUsed;++i)
        digestT4Dependency(digest,report.dependencyStorage[
            static_cast<std::size_t>(i)]);
    for (std::uint64_t i=0;i<report.witnessesUsed;++i) {
        const auto& witness=report.witnessStorage[static_cast<std::size_t>(i)];
        digestWord(digest,static_cast<std::uint64_t>(witness.kind));
        for (const auto node:witness.edge.canonicalGlobalNodes)
            digestWord(digest,node);
        digestWord(digest,witness.facet.facePersistentId);
        digestWord(digest,witness.facet.normalizedLocalCellOrdinal);
    }
    for (std::uint64_t i=0;i<report.topologyUsed;++i) {
        const auto& row=report.topologyStorage[static_cast<std::size_t>(i)];
        digestWord(digest,static_cast<std::uint64_t>(row.kind));
        digestWord(digest,static_cast<std::uint64_t>(row.transition));
        digestHandle(digest,row.predictedNode); digestHandle(digest,row.actualNode);
        digestHandle(digest,row.predictedCell); digestHandle(digest,row.actualCell);
        digestHandle(digest,row.predictedEdge); digestHandle(digest,row.actualEdge);
        digestHandle(digest,row.predictedFace); digestHandle(digest,row.actualFace);
    }
    for (std::uint64_t i=0;i<report.unattemptedUsed;++i)
        digestT4QualityKey(digest,report.unattemptedStorage[
            static_cast<std::size_t>(i)]);
    for (std::uint64_t i=0;i<report.tailUsed;++i)
        digestT4QualityKey(digest,report.tailStorage[static_cast<std::size_t>(i)]);
    if (report.hasTerminal) {
        digestT4QualityKey(digest,report.terminal.selected);
        digestWord(digest,static_cast<std::uint64_t>(report.terminal.reason));
        digestWord(digest,static_cast<std::uint64_t>(report.terminal.phase));
        digestWord(digest,report.terminal.attemptOrdinal);
    }
}
void digestT4Metrics(std::uint64_t& digest,const T4AMetrics& metrics) noexcept {
    for (const auto value:{metrics.initialCellsScanned,
         metrics.finalCellsAudited,metrics.qualityEvaluations,
         metrics.queuePushes,metrics.queuePops,metrics.staleRows,
         metrics.queueHighWater,metrics.selectedGenerations,
         metrics.candidateAttempts,metrics.formulaRejects,
         metrics.duplicateRejects,metrics.protectedIndexVisits,
         metrics.locationIndexVisits,metrics.locationTetTests,
         metrics.cavityCellsVisited,metrics.cavityFaceVisits,
         metrics.survivorIndexVisits,metrics.candidateSurvivorTests,
         metrics.localRowsPrepared,metrics.localRowsApplied,
         metrics.blockedRows,metrics.blockedInvalidations,
         metrics.preparedTransactions,metrics.committedOperations,
         metrics.perOperationGlobalTetScans,metrics.contextRebuilds,
         metrics.phase3Allocations,metrics.publicationTailAllocations,
         metrics.workUnitsConsumed,metrics.transientAccountedHighWaterBytes})
        digestWord(digest,value);
    digestWord(digest,metrics.runComplete?1U:0U);
    digestWord(digest,metrics.qualityComplete?1U:0U);
}
void recordT4StateDigest(
    t4a_detail::T4ADigestCheckpoint checkpoint,bool postFailure,
    const topo::TopologyWorkspace& workspace,const CoordinateMap& coordinates,
    const InsertedMap& inserted,const CoordinateIndex& coordinateIndex,
    const CellMap& cells,const EdgeStars& edgeStars,
    const FaceOwners& faceOwners,const NodeStars& nodeStars,
    const ViolationHeap& queue,const CurrentBlockedHeap& blocked,
    const T4BlockedDependencyIndex& dependencyIndex,
    const RefinementLocationIndex& location,const T4AMutableArtifact& artifact,
    const T4AReport& report,const T4AMetrics& metrics,
    std::uint64_t preparedPrediction=0) {
    auto& trace=t4a_detail::digestTraceForTesting();
    if (!trace.enabled||trace.used>=trace.rows.size()) return;
    if ((trace.faultArmed||trace.faultPending)&&
        checkpoint!=trace.faultCheckpoint&&
        !(trace.faultCheckpoint==
              t4a_detail::T4ADigestCheckpoint::Prepared&&
          checkpoint==t4a_detail::T4ADigestCheckpoint::PrePreparation))
        return;
    auto& row=trace.rows[trace.used++];
    row={};
    row.checkpoint=checkpoint;
    row.postFailure=postFailure;
    row.revision=workspace.topologyRevision();
    // The Prepared checkpoint deliberately observes a workspace while the
    // transaction owns live placeholders.  Public workspace hash entry
    // points reject that state, so its last pre-preparation topology digest
    // is the authoritative unchanged logical topology.  The exact predicted
    // staged/cell/edge/face handles are recorded separately below, while the
    // prior pool digest pins the allocator state that must be restored.
    if (checkpoint==t4a_detail::T4ADigestCheckpoint::Prepared&&
        !postFailure&&trace.used>1)
        row.topology=trace.rows[trace.used-2].topology;
    else {
        row.topology=workspace.strictWorkspaceHash();
        digestWord(row.topology,workspace.semanticTopologyHash());
    }
    row.prepared=preparedPrediction;
    if (checkpoint==t4a_detail::T4ADigestCheckpoint::Prepared&&
        !postFailure&&trace.used>1)
        row.pools=trace.rows[trace.used-2].pools;
    else {
        row.pools=0xcbf29ce484222325ull;
        for (const auto& stats:{workspace.nodes().stats(),
                               workspace.cells().stats(),
                               workspace.edges().stats(),
                               workspace.faces().stats()}) {
            digestWord(row.pools,stats.live);
            digestWord(row.pools,stats.reusableFree);
            digestWord(row.pools,stats.retired);
            digestWord(row.pools,stats.totalCreated);
            digestWord(row.pools,stats.materializedSlots);
        }
    }
    row.mappings=0xcbf29ce484222325ull;
    const auto digestIndices=[&](const auto& values) noexcept {
        digestWord(row.mappings,values.size());
        for (const auto value:values) digestWord(row.mappings,value);
    };
    digestIndices(artifact.existingGlobalByAuthoritativeNodeSlot);
    digestIndices(artifact.sourceGlobalByAuthoritativeNodeSlot);
    digestIndices(artifact.sourceGlobalByAuthoritativeCellSlot);
    digestIndices(artifact.publishedGlobalByAuthoritativeNodeSlot);
    digestIndices(artifact.publishedGlobalByAuthoritativeCellSlot);
    for (const auto& [handle,key]:artifact.vertexKeyByHandle) {
        digestHandle(row.mappings,handle); digestVertex(row.mappings,key);
    }
    row.resolvers=0xcbf29ce484222325ull;
    for (const auto& [handle,point]:coordinates) {
        digestHandle(row.resolvers,handle);
        digestWord(row.resolvers,std::bit_cast<std::uint64_t>(point.x));
        digestWord(row.resolvers,std::bit_cast<std::uint64_t>(point.y));
        digestWord(row.resolvers,std::bit_cast<std::uint64_t>(point.z));
    }
    for (const auto& [handle,ordinal]:inserted) {
        digestHandle(row.resolvers,handle); digestWord(row.resolvers,ordinal);
    }
    for (const auto& [bits,handle]:coordinateIndex) {
        for (const auto value:bits) digestWord(row.resolvers,value);
        digestHandle(row.resolvers,handle);
    }
    for (const auto& [cell,nodes]:cells) {
        digestHandle(row.resolvers,cell);
        for (const auto node:nodes) digestHandle(row.resolvers,node);
    }
    const auto digestStars=[&](const auto& stars) noexcept {
        for (const auto& [key,owners]:stars) {
            if constexpr (requires { key.pool; }) digestHandle(row.resolvers,key);
            else for (const auto node:key) digestHandle(row.resolvers,node);
            for (const auto owner:owners) digestHandle(row.resolvers,owner);
        }
    };
    digestStars(edgeStars); digestStars(faceOwners); digestStars(nodeStars);
    row.location=location.stateDigestForTesting();
    row.queue=0xcbf29ce484222325ull;
    for (const auto& value:queue.rows()) digestQueued(row.queue,value);
    row.blockers=0xcbf29ce484222325ull;
    for (const auto& value:blocked.rowsForTesting())
        digestBlocked(row.blockers,value);
    row.dependencyIndex=0xcbf29ce484222325ull;
    for (const auto& value:dependencyIndex.rowsForTesting()) {
        digestT4Dependency(row.dependencyIndex,value.key);
        digestHandle(row.dependencyIndex,value.selected);
        digestWord(row.dependencyIndex,value.blockedRowIndex);
    }
    row.report=0xcbf29ce484222325ull;
    digestT4Report(row.report,report);
    row.metrics=0xcbf29ce484222325ull;
    digestT4Metrics(row.metrics,metrics);
    row.attemptsUsed=report.attemptsUsed;
    row.blockedUsed=report.blockedUsed;
    row.commitsUsed=report.commitsUsed;
    row.dependenciesUsed=report.dependenciesUsed;
    row.witnessesUsed=report.witnessesUsed;
    row.topologyUsed=report.topologyUsed;
    row.unattemptedUsed=report.unattemptedUsed;
    row.structural=0xcbf29ce484222325ull;
    for (const auto value:{row.topology,row.pools,row.prepared,row.mappings,
         row.resolvers,row.location,row.queue,row.blockers,
         row.dependencyIndex,row.revision})
        digestWord(row.structural,value);
    if (!postFailure&&trace.faultArmed&&
        trace.faultCheckpoint==checkpoint) {
        trace.faultArmed=false;
        trace.faultPending=true;
        throw std::bad_alloc{};
    }
}
#endif

struct ExactLocation {
    topo::CellHandle seed{};
    std::vector<topo::CellHandle> star;
};

std::optional<ExactLocation> locateExact(
    Vec3 point,const bnd::DynamicTetAabbIndex& base,
    const RefinementLocationIndex& index,
    const topo::TopologyWorkspace& workspace,
    const std::map<topo::NodeHandle,Vec3>& coordinates,
    TetRefinementMetrics& metrics) {
    const auto visitsBefore=metrics.locationIndexVisits;
    auto candidates=index.query(base,pointAabb(point),workspace,
                                metrics.locationIndexVisits);
    recordConstructionWork(
        metrics,metrics.locationIndexVisits-visitsBefore);
    ExactLocation result;
    std::size_t bestZeros=(std::numeric_limits<std::size_t>::max)();
    for (const auto cell:candidates) {
        ++metrics.locationTetTests;
        recordConstructionWork(metrics);
        const auto data=workspace.cells().get(cell);
        if (data.type!=topo::ElementType::Tet4) continue;
        const std::array<topo::NodeHandle,4> nodes{
            data.nodes[0],data.nodes[1],data.nodes[2],data.nodes[3]};
        bool inside=true;
        std::size_t zeros=0;
        for (std::size_t f=0;f<4;++f) {
            const auto& local=kLocalTetFaces[f];
            const auto a=coordinates.at(nodes[local[0]]);
            const auto b=coordinates.at(nodes[local[1]]);
            const auto c=coordinates.at(nodes[local[2]]);
            const auto opposite=coordinates.at(nodes[f]);
            metrics.orient3dTests+=2;
            recordConstructionWork(metrics,2);
            const auto reference=predicates::orient3d(a,b,c,opposite);
            const auto side=predicates::orient3d(a,b,c,point);
            if (reference==predicates::Sign::Zero||
                (side!=predicates::Sign::Zero&&side!=reference)) {
                inside=false;
                break;
            }
            if (side==predicates::Sign::Zero) ++zeros;
        }
        if (!inside) continue;
        result.star.push_back(cell);
        if (!result.seed.isSpecified()||zeros<bestZeros) {
            result.seed=cell;
            bestZeros=zeros;
        }
    }
    if (result.star.empty()) return std::nullopt;
    // Candidate iteration is sequence ordered.  Equal-dimensional exact
    // face/edge/vertex ties therefore retain the canonical first owner.
    return result;
}

std::optional<std::vector<topo::CellHandle>> discoverConflictCavity(
    topo::CellHandle seed,Vec3 point,const CellMap& cells,
    const FaceOwners& owners,const std::map<topo::NodeHandle,Vec3>& coordinates,
    TetRefinementMetrics& metrics) {
    if (!cells.count(seed)) return std::nullopt;
    std::set<topo::CellHandle> pending{seed};
    std::vector<topo::CellHandle> result;
    std::set<topo::CellHandle> visited;
    while (!pending.empty()) {
        const auto first=pending.begin();
        const auto cell=*first;
        pending.erase(first);
        if (!visited.insert(cell).second) continue;
        recordConstructionWork(metrics);
        ++metrics.cavityCellsVisited;
        const auto found=cells.find(cell);
        if (found==cells.end()) return std::nullopt;
        std::array<Vec3,4> p{};
        for (std::size_t i=0;i<4;++i) p[i]=coordinates.at(found->second[i]);
        recordConstructionWork(metrics);
        ++metrics.insphereTests;
        if (predicates::insphere(p[0],p[1],p[2],p[3],point)==
                predicates::Sign::Negative)
            continue;
        result.push_back(cell);
        for (const auto& local:kLocalTetFaces) {
            recordConstructionWork(metrics);
            const auto key=faceKey({found->second[local[0]],
                                    found->second[local[1]],
                                    found->second[local[2]]});
            const auto owner=owners.find(key);
            if (owner==owners.end()||owner->second.empty()||
                owner->second.size()>2) return std::nullopt;
            for (const auto neighbour:owner->second)
                if (neighbour!=cell&&!visited.count(neighbour))
                    pending.insert(neighbour);
        }
    }
    std::sort(result.begin(),result.end());
    if (result.empty()||!std::binary_search(result.begin(),result.end(),seed))
        return std::nullopt;
    return result;
}

std::optional<std::vector<BoundaryFace>> cavityBoundary(
    const std::vector<topo::CellHandle>& removed,const CellMap& cells) {
    struct Agg { std::size_t count=0; BoundaryFace face{}; };
    std::map<FaceKey,Agg> rows;
    static constexpr std::array<std::array<int,3>,4> faces{{
        {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};
    for (const auto cell:removed) {
        const auto found=cells.find(cell);
        if (found==cells.end()) return std::nullopt;
        for (std::size_t f=0;f<4;++f) {
            BoundaryFace face;
            for (std::size_t i=0;i<3;++i)
                face.directed[i]=found->second[
                    static_cast<std::size_t>(faces[f][i])];
            face.opposite=found->second[f];
            auto& row=rows[faceKey(face.directed)];
            ++row.count;
            row.face=face;
            if (row.count>2) return std::nullopt;
        }
    }
    std::vector<BoundaryFace> boundary;
    for (const auto& [key,row]:rows) {
        (void)key;
        if (row.count==1) boundary.push_back(row.face);
    }
    if (boundary.empty()) return std::nullopt;
    return boundary;
}

bool pointStrictlyInsideBoundary(const std::vector<BoundaryFace>& boundary,
                                 Vec3 point,const topo::TopologyWorkspace& ws,
                                 TetRefinementMetrics& metrics) {
    for (const auto& face:boundary) {
        const Vec3 a=ws.nodes().get(face.directed[0]).coord;
        const Vec3 b=ws.nodes().get(face.directed[1]).coord;
        const Vec3 c=ws.nodes().get(face.directed[2]).coord;
        const Vec3 opposite=ws.nodes().get(face.opposite).coord;
        metrics.orient3dTests+=2;
        recordConstructionWork(metrics,2);
        const auto reference=predicates::orient3d(a,b,c,opposite);
        const auto candidate=predicates::orient3d(a,b,c,point);
        if (reference==predicates::Sign::Zero || candidate!=reference)
            return false;
    }
    return true;
}

std::int8_t directedFaceSense(
    const std::array<topo::NodeHandle,3>& face,
    const topo::TopologyWorkspace& workspace) {
    std::array<std::pair<std::uint64_t,topo::NodeHandle>,3> ordered{};
    for (std::size_t i=0;i<3;++i)
        ordered[i]={workspace.nodes().seqOf(face[i]),face[i]};
    std::sort(ordered.begin(),ordered.end(),[](const auto& a,const auto& b){
        return a.first<b.first||(a.first==b.first&&a.second<b.second);
    });
    std::array<std::size_t,3> permutation{};
    for (std::size_t i=0;i<3;++i)
        for (std::size_t j=0;j<3;++j)
            if (face[i]==ordered[j].second) permutation[i]=j;
    std::size_t inversions=0;
    for (std::size_t i=0;i<3;++i)
        for (std::size_t j=i+1;j<3;++j)
            inversions+=permutation[i]>permutation[j]?1U:0U;
    return inversions%2==0?std::int8_t{1}:std::int8_t{-1};
}

struct CandidateTet {
    std::array<topo::NodeHandle,3> boundary{};
};

bool pairHasIllegalContact(
    const CandidateTet& candidate,Vec3 stagedPoint,
    const std::array<topo::NodeHandle,4>& survivor,
    const std::map<topo::NodeHandle,Vec3>& coordinates,
    const topo::TopologyWorkspace& workspace,
    mesh::Classification region) {
    std::vector<topo::NodeHandle> handles(candidate.boundary.begin(),
                                          candidate.boundary.end());
    handles.insert(handles.end(),survivor.begin(),survivor.end());
    std::sort(handles.begin(),handles.end(),[&](const auto a,const auto b){
        const auto as=workspace.nodes().seqOf(a);
        const auto bs=workspace.nodes().seqOf(b);
        return as<bs||(as==bs&&a<b);
    });
    handles.erase(std::unique(handles.begin(),handles.end()),handles.end());
    if (handles.size()>=(std::numeric_limits<bnd::TetFacetVertexId>::max)())
        return true;
    std::vector<bnd::TetFacetVertex> vertices;
    vertices.reserve(handles.size()+1);
    std::map<topo::NodeHandle,bnd::TetFacetVertexId> ids;
    for (const auto handle:handles) {
        const auto id=static_cast<bnd::TetFacetVertexId>(vertices.size());
        ids.emplace(handle,id);
        vertices.push_back({coordinates.at(handle),
                            workspace.nodes().seqOf(handle)});
    }
    const auto staged=static_cast<bnd::TetFacetVertexId>(vertices.size());
    vertices.push_back({stagedPoint,
        (std::numeric_limits<std::uint64_t>::max)()});
    bnd::TetFacetCell a{{ids.at(candidate.boundary[0]),
                         ids.at(candidate.boundary[1]),
                         ids.at(candidate.boundary[2]),staged},0,region};
    bnd::TetFacetCell b{{ids.at(survivor[0]),ids.at(survivor[1]),
                         ids.at(survivor[2]),ids.at(survivor[3])},1,region};
    const auto orient=[&](bnd::TetFacetCell& cell) {
        const auto sign=predicates::orient3d(
            vertices[cell.vertices[0]].position,
            vertices[cell.vertices[1]].position,
            vertices[cell.vertices[2]].position,
            vertices[cell.vertices[3]].position);
        if (sign==predicates::Sign::Negative)
            std::swap(cell.vertices[0],cell.vertices[1]);
        return sign!=predicates::Sign::Zero;
    };
    if (!orient(a)||!orient(b)) return true;
    const std::array<bnd::TetFacetCell,2> cells{a,b};
    const bnd::TetFacetComplexView view{
        vertices,cells,{},std::span<const bnd::TetFacetEdgeKey>{},
        std::span<const bnd::TetFacetProtectedFace>{},region};
    return bnd::tet_facet_detail::tetImproperlyIntersects(view,a,b);
}

bool validateRefinementCandidate(
    const std::vector<topo::CellHandle>& removed,
    const std::vector<BoundaryFace>& boundary,Vec3 stagedPoint,
    const CellMap& cells,const EdgeStars& edgeStars,
    const NodeStars& nodeStars,
    const std::map<topo::NodeHandle,Vec3>& coordinates,
    const bnd::DynamicTetAabbIndex& baseIndex,
    const RefinementLocationIndex& locationIndex,
    const topo::TopologyWorkspace& workspace,mesh::Classification region,
    TetRefinementMetrics& metrics) {
    if (removed.empty()||boundary.empty()) return false;
    std::set<topo::CellHandle> removedSet(removed.begin(),removed.end());
    if (removedSet.size()!=removed.size()) return false;

    std::vector<topo::NodeHandle> handles;
    for (const auto cell:removed) {
        const auto found=cells.find(cell);
        if (found==cells.end()) return false;
        handles.insert(handles.end(),found->second.begin(),found->second.end());
    }
    std::sort(handles.begin(),handles.end(),[&](const auto a,const auto b){
        const auto as=workspace.nodes().seqOf(a);
        const auto bs=workspace.nodes().seqOf(b);
        return as<bs||(as==bs&&a<b);
    });
    handles.erase(std::unique(handles.begin(),handles.end()),handles.end());
    if (handles.size()+1>=
        (std::numeric_limits<bnd::TetFacetVertexId>::max)()) return false;
    std::map<topo::NodeHandle,bnd::TetFacetVertexId> ids;
    std::vector<bnd::TetFacetVertex> vertices;
    vertices.reserve(handles.size()+1);
    for (const auto handle:handles) {
        const auto id=static_cast<bnd::TetFacetVertexId>(vertices.size());
        ids.emplace(handle,id);
        vertices.push_back({coordinates.at(handle),
                            workspace.nodes().seqOf(handle)});
    }
    const auto staged=static_cast<bnd::TetFacetVertexId>(vertices.size());
    vertices.push_back({stagedPoint,
        (std::numeric_limits<std::uint64_t>::max)()});

    std::vector<bnd::TetFacetCell> removedPlain;
    removedPlain.reserve(removed.size());
    for (const auto cell:removed) {
        const auto& conn=cells.at(cell);
        bnd::TetFacetCell row{{ids.at(conn[0]),ids.at(conn[1]),
                              ids.at(conn[2]),ids.at(conn[3])},
                             workspace.cells().seqOf(cell),region};
        removedPlain.push_back(row);
    }
    const bnd::TetFacetComplexView removedView{
        vertices,removedPlain,{},std::span<const bnd::TetFacetEdgeKey>{},
        std::span<const bnd::TetFacetProtectedFace>{},region};
    std::vector<std::size_t> removedIds(removedPlain.size());
    std::iota(removedIds.begin(),removedIds.end(),0);
    if (!bnd::tet_facet_detail::auditBallImpl(
            removedView,removedIds).valid) return false;

    std::vector<bnd::TetFacetCell> candidatePlain;
    std::vector<CandidateTet> candidateTets;
    candidatePlain.reserve(boundary.size());
    candidateTets.reserve(boundary.size());
    for (const auto& face:boundary) {
        bnd::TetFacetCell row{{ids.at(face.directed[0]),
                              ids.at(face.directed[1]),
                              ids.at(face.directed[2]),staged},
                             candidatePlain.size(),region};
        const auto sign=predicates::orient3d(
            vertices[row.vertices[0]].position,
            vertices[row.vertices[1]].position,
            vertices[row.vertices[2]].position,stagedPoint);
        metrics.orient3dTests+=1;
        recordConstructionWork(metrics);
        if (sign==predicates::Sign::Zero) return false;
        if (sign==predicates::Sign::Negative)
            std::swap(row.vertices[0],row.vertices[1]);
        candidatePlain.push_back(row);
        candidateTets.push_back({face.directed});
    }
    const bnd::TetFacetComplexView candidateView{
        vertices,candidatePlain,{},std::span<const bnd::TetFacetEdgeKey>{},
        std::span<const bnd::TetFacetProtectedFace>{},region};
    std::vector<std::size_t> candidateIds(candidatePlain.size());
    std::iota(candidateIds.begin(),candidateIds.end(),0);
    if (!bnd::tet_facet_detail::auditBallImpl(
            candidateView,candidateIds).valid) return false;

    struct FaceUse {
        std::size_t count=0;
        std::int8_t sense=0;
        bool operator==(const FaceUse&) const = default;
    };
    std::map<FaceKey,FaceUse> expectedBoundary,candidateFaces;
    for (const auto& face:boundary) {
        auto& use=expectedBoundary[faceKey(face.directed)];
        if (++use.count!=1) return false;
        use.sense=directedFaceSense(face.directed,workspace);
    }
    for (const auto& cell:candidatePlain)
        for (std::size_t f=0;f<4;++f) {
            std::array<topo::NodeHandle,3> directed{};
            bool hasStaged=false;
            for (std::size_t i=0;i<3;++i) {
                const auto id=cell.vertices[kLocalTetFaces[f][i]];
                if (id==staged) { hasStaged=true; break; }
                directed[i]=handles[id];
            }
            if (hasStaged) continue;
            auto& use=candidateFaces[faceKey(directed)];
            ++use.count;
            use.sense=directedFaceSense(directed,workspace);
        }
    if (candidateFaces!=expectedBoundary) return false;

    std::set<NodePair> candidateExistingEdges;
    for (const auto& face:boundary) {
        candidateExistingEdges.insert(nodePair(face.directed[0],face.directed[1]));
        candidateExistingEdges.insert(nodePair(face.directed[1],face.directed[2]));
        candidateExistingEdges.insert(nodePair(face.directed[2],face.directed[0]));
    }
    for (const auto removedCell:removed) {
        if (!workspace.edges().hasCellEdges(workspace.cells(),removedCell)||
            !workspace.faces().hasCellFaces(workspace.cells(),removedCell))
            return false;
        for (const auto edge:workspace.edges().cellEdges(
                 workspace.cells(),removedCell)) {
            const auto data=workspace.edges().get(edge);
            if (!data.requiredConstraint) continue;
            const auto key=nodePair(data.a,data.b);
            if (!candidateExistingEdges.count(key)||!edgeStars.count(key))
                return false;
        }
        for (const auto face:workspace.faces().cellFaces(
                 workspace.cells(),removedCell)) {
            const auto data=workspace.faces().get(face);
            if (data.cls.dim!=2||data.cls.tag<0) continue;
            if (data.cornerCount!=3) return false;
            const FaceKey key=faceKey(
                {data.corners[0],data.corners[1],data.corners[2]});
            if (!candidateFaces.count(key)) return false;
        }
    }

    for (const auto node:handles) {
        bool usedByCandidate=false;
        for (const auto& face:boundary)
            usedByCandidate=usedByCandidate||
                std::find(face.directed.begin(),face.directed.end(),node)!=
                    face.directed.end();
        if (usedByCandidate) continue;
        const auto star=nodeStars.find(node);
        if (star==nodeStars.end()) return false;
        const bool survives=std::any_of(star->second.begin(),star->second.end(),
            [&](const auto cell){return !removedSet.count(cell);});
        if (!survives) return false;
    }

    for (const auto& candidate:candidateTets) {
        std::array<topo::NodeHandle,4> temporary{
            candidate.boundary[0],candidate.boundary[1],
            candidate.boundary[2],candidate.boundary[0]};
        auto bounds=tetAabb(temporary,coordinates);
        bounds.lo.x=(std::min)(bounds.lo.x,stagedPoint.x);
        bounds.lo.y=(std::min)(bounds.lo.y,stagedPoint.y);
        bounds.lo.z=(std::min)(bounds.lo.z,stagedPoint.z);
        bounds.hi.x=(std::max)(bounds.hi.x,stagedPoint.x);
        bounds.hi.y=(std::max)(bounds.hi.y,stagedPoint.y);
        bounds.hi.z=(std::max)(bounds.hi.z,stagedPoint.z);
        const auto visitsBefore=metrics.survivorIndexVisits;
        auto survivors=locationIndex.query(baseIndex,bounds,workspace,
                                            metrics.survivorIndexVisits);
        recordConstructionWork(
            metrics,metrics.survivorIndexVisits-visitsBefore);
        for (const auto survivor:survivors) {
            if (removedSet.count(survivor)) continue;
            const auto found=cells.find(survivor);
            if (found==cells.end()) return false;
            ++metrics.candidateSurvivorTests;
            recordConstructionWork(metrics);
            if (pairHasIllegalContact(candidate,stagedPoint,found->second,
                                      coordinates,workspace,region))
                return false;
        }
    }
    return true;
}

struct PreparedBlockedInvalidation {
    topo::CellHandle cell{};
    std::uint64_t blockedRowIndex=0;
    std::uint64_t insertionOrdinal=0;
};

struct PreparedTetRefinementCoordinatorDelta {
    PreparedTetRefinementCoordinatorDelta(
        const PreparedTetRefinementCoordinatorDelta&)=delete;
    PreparedTetRefinementCoordinatorDelta& operator=(
        const PreparedTetRefinementCoordinatorDelta&)=delete;
    PreparedTetRefinementCoordinatorDelta(
        PreparedTetRefinementCoordinatorDelta&&)=delete;
    PreparedTetRefinementCoordinatorDelta& operator=(
        PreparedTetRefinementCoordinatorDelta&&)=delete;

    PreparedTetRefinementCoordinatorDelta()=default;

    CoordinateMap coordinateCreates;
    InsertedMap insertedCreates;
    CoordinateIndex coordinateIndexCreates;
    CellMap cellCreates;
    std::vector<topo::CellHandle> cellRetires;
    EdgeStars edgeRows;
    std::vector<NodePair> edgeErases;
    FaceOwners faceRows;
    std::vector<FaceKey> faceErases;
    NodeStars nodeRows;
    std::vector<topo::NodeHandle> nodeErases;
    std::vector<topo::CellHandle> queueRetires;
    std::vector<QueuedViolation> queueCreates;
    std::vector<PreparedBlockedInvalidation> blockedInvalidations;
    std::vector<RefinementLocationIndex::CreatedRow> locationCreates;
    std::vector<T3ATopologyHandleRow> topologyRows;
    std::set<NodePair> affectedEdges;
    std::set<FaceKey> affectedFaces;
    std::set<topo::NodeHandle> affectedNodes;
    std::uint64_t contextBytesRemoved=0;
    std::uint64_t contextBytesAdded=0;
};

std::uint64_t preparedCoordinatorAccountedBytes(
    const PreparedTetRefinementCoordinatorDelta& prepared) {
    std::uint64_t result=0;
    addAccountedBytes(result,orderedMapAccountedBytes(
        prepared.coordinateCreates));
    addAccountedBytes(result,orderedMapAccountedBytes(
        prepared.insertedCreates));
    addAccountedBytes(result,orderedMapAccountedBytes(
        prepared.coordinateIndexCreates));
    addAccountedBytes(result,orderedMapAccountedBytes(prepared.cellCreates));
    addAccountedBytes(result,starMapAccountedBytes(prepared.edgeRows));
    addAccountedBytes(result,starMapAccountedBytes(prepared.faceRows));
    addAccountedBytes(result,starMapAccountedBytes(prepared.nodeRows));
    addAccountedBytes(result,vectorAccountedBytes(prepared.cellRetires));
    addAccountedBytes(result,vectorAccountedBytes(prepared.edgeErases));
    addAccountedBytes(result,vectorAccountedBytes(prepared.faceErases));
    addAccountedBytes(result,vectorAccountedBytes(prepared.nodeErases));
    addAccountedBytes(result,vectorAccountedBytes(prepared.queueRetires));
    addAccountedBytes(result,vectorAccountedBytes(prepared.queueCreates));
    addAccountedBytes(result,
        vectorAccountedBytes(prepared.blockedInvalidations));
    addAccountedBytes(result,vectorAccountedBytes(prepared.locationCreates));
    addAccountedBytes(result,vectorAccountedBytes(prepared.topologyRows));
    addAccountedBytes(result,orderedSetAccountedBytes(prepared.affectedEdges));
    addAccountedBytes(result,orderedSetAccountedBytes(prepared.affectedFaces));
    addAccountedBytes(result,orderedSetAccountedBytes(prepared.affectedNodes));
    return result;
}

std::uint64_t applyAccountedContextDelta(
    std::uint64_t current,
    const PreparedTetRefinementCoordinatorDelta& prepared) {
    if (current<prepared.contextBytesRemoved)
        throw std::logic_error("refinement memory accounting underflow");
    std::uint64_t result=current-prepared.contextBytesRemoved;
    addAccountedBytes(result,prepared.contextBytesAdded);
    return result;
}

Vec3 preparedCoordinate(topo::NodeHandle node,
                        topo::NodeHandle stagedNode,Vec3 stagedPoint,
                        const CoordinateMap& coordinates) {
    return node==stagedNode?stagedPoint:coordinates.at(node);
}

bnd::TetAabb preparedTetAabb(
    const std::array<topo::NodeHandle,4>& nodes,
    topo::NodeHandle stagedNode,Vec3 stagedPoint,
    const CoordinateMap& coordinates) {
    const auto first=preparedCoordinate(
        nodes[0],stagedNode,stagedPoint,coordinates);
    bnd::TetAabb result{first,first};
    for (std::size_t i=1;i<4;++i) {
        const auto p=preparedCoordinate(
            nodes[i],stagedNode,stagedPoint,coordinates);
        result.lo.x=(std::min)(result.lo.x,p.x);
        result.lo.y=(std::min)(result.lo.y,p.y);
        result.lo.z=(std::min)(result.lo.z,p.z);
        result.hi.x=(std::max)(result.hi.x,p.x);
        result.hi.y=(std::max)(result.hi.y,p.y);
        result.hi.z=(std::max)(result.hi.z,p.z);
    }
    return result;
}

template<class Map,class Keys>
void copyAffectedRows(const Map& source,const Keys& keys,Map& destination,
                      TetRefinementMetrics& metrics) {
    for (const auto& key:keys) {
        const auto found=source.find(key);
        if (found==source.end()) destination.emplace(key,typename Map::mapped_type{});
        else destination.emplace(key,found->second);
        const std::uint64_t rows=found==source.end()?0:
            static_cast<std::uint64_t>(found->second.size());
        metrics.coordinatorRowsPrepared+=1+rows;
        recordConstructionWork(metrics,1+rows);
    }
}

void prepareCoordinatorDelta(
    PreparedTetRefinementCoordinatorDelta& prepared,
    topo::NodeHandle stagedNode,Vec3 stagedPoint,
    std::uint64_t insertionOrdinal,
    std::span<const topo::CellHandle> removed,
    const CreatedCellRows& createdRows,
    std::span<const QueuedViolation> createdViolations,
    const CoordinateMap& coordinates,const CellMap& cells,
    const EdgeStars& edgeStars,const FaceOwners& faceOwners,
    const NodeStars& nodeStars,const ViolationHeap& queue,
    const CurrentBlockedHeap& blocked,const T3AReport& report,
    const RefinementLocationIndex& location,TetRefinementMetrics& metrics) {
    if (!prepared.coordinateCreates.emplace(stagedNode,stagedPoint).second||
        !prepared.insertedCreates.emplace(stagedNode,insertionOrdinal).second||
        !prepared.coordinateIndexCreates.emplace(
            coordinateBits(stagedPoint),stagedNode).second)
        throw std::logic_error("duplicate prepared coordinator key");

    prepared.cellRetires.assign(removed.begin(),removed.end());
    prepared.queueRetires=prepared.cellRetires;
    prepared.queueCreates.assign(
        createdViolations.begin(),createdViolations.end());
    prepared.locationCreates.reserve(createdRows.size());
    prepared.blockedInvalidations.reserve(removed.size());
    for (const auto cell:removed) {
        const auto* current=blocked.find(cell);
        if (!current) continue;
        if (current->blockedRowIndex>=report.blockedUsed)
            throw std::logic_error("invalid current blocker report index");
        const auto& historical=report.blockedStorage[
            static_cast<std::size_t>(current->blockedRowIndex)];
        if (!(historical.violation==current->row.key)||
            historical.violationGeneration!=current->row.generation||
            historical.invalidatedByInsertionOrdinal!=
                kCurrentBlockedInsertion)
            throw std::logic_error("stale current blocker binding");
#if defined(PAACS_TOPO_TESTING)
        auto& invalidationFaultCountdown=
            t3a_detail::blockerInvalidationFaultCountdownForTesting();
        if (invalidationFaultCountdown>0&&
            --invalidationFaultCountdown==0)
            throw std::bad_alloc{};
#endif
        prepared.blockedInvalidations.push_back(
            {cell,current->blockedRowIndex,insertionOrdinal});
    }
    std::sort(prepared.blockedInvalidations.begin(),
              prepared.blockedInvalidations.end(),
              [](const auto& a,const auto& b) { return a.cell<b.cell; });
    if (!prepared.blockedInvalidations.empty()) {
        const auto depth=static_cast<std::uint64_t>(
            std::bit_width(blocked.size()));
        std::uint64_t perRow=0,total=0;
        if (!checkedAdd(3,2*depth,perRow)||
            !checkedMul(perRow,prepared.blockedInvalidations.size(),total))
            throw std::length_error("blocked invalidation work overflow");
        metrics.coordinatorRowsPrepared+=
            prepared.blockedInvalidations.size();
        recordConstructionWork(metrics,total);
    }

    auto& affectedEdges=prepared.affectedEdges;
    auto& affectedFaces=prepared.affectedFaces;
    auto& affectedNodes=prepared.affectedNodes;
    const auto collect=[&](const std::array<topo::NodeHandle,4>& conn) {
        static constexpr std::array<std::array<std::size_t,2>,6> edges{{
            {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
        for (const auto& edge:edges)
            affectedEdges.insert(nodePair(conn[edge[0]],conn[edge[1]]));
        for (const auto& face:kLocalTetFaces)
            affectedFaces.insert(faceKey(
                {conn[face[0]],conn[face[1]],conn[face[2]]}));
        affectedNodes.insert(conn.begin(),conn.end());
    };
    for (const auto cell:removed) {
        const auto found=cells.find(cell);
        if (found==cells.end()) throw std::logic_error("missing retired cell");
        collect(found->second);
        metrics.coordinatorRowsPrepared+=1;
        recordConstructionWork(metrics);
    }
    for (const auto& [cell,conn]:createdRows) {
        if (!prepared.cellCreates.emplace(cell,conn).second)
            throw std::logic_error("duplicate predicted cell");
        collect(conn);
        prepared.locationCreates.push_back(
            {cell,preparedTetAabb(conn,stagedNode,stagedPoint,coordinates)});
        metrics.coordinatorRowsPrepared+=1;
        recordConstructionWork(metrics);
    }

    copyAffectedRows(edgeStars,affectedEdges,prepared.edgeRows,metrics);
    copyAffectedRows(faceOwners,affectedFaces,prepared.faceRows,metrics);
    copyAffectedRows(nodeStars,affectedNodes,prepared.nodeRows,metrics);
    for (const auto cell:removed) {
        const auto& conn=cells.at(cell);
        removeCellFromStars(prepared.edgeRows,cell,conn);
        removeCellFromFaces(prepared.faceRows,cell,conn);
        removeCellFromNodeStars(prepared.nodeRows,cell,conn);
    }
    for (const auto& [cell,conn]:createdRows) {
        addCellToStars(prepared.edgeRows,cell,conn);
        addCellToFaces(prepared.faceRows,cell,conn);
        addCellToNodeStars(prepared.nodeRows,cell,conn);
    }
    for (const auto& key:affectedEdges)
        if (edgeStars.count(key)&&!prepared.edgeRows.count(key))
            prepared.edgeErases.push_back(key);
    for (const auto& key:affectedFaces)
        if (faceOwners.count(key)&&!prepared.faceRows.count(key))
            prepared.faceErases.push_back(key);
    for (const auto node:affectedNodes)
        if (nodeStars.count(node)&&!prepared.nodeRows.count(node))
            prepared.nodeErases.push_back(node);

    // The context high-water is maintained by this affected-row delta.  This
    // deliberately touches only the same local keys already prepared above;
    // no insertion-time accounting pass over a global resolver is permitted.
    addAccountedBytes(prepared.contextBytesRemoved,
        orderedNodeAccountedBytes(
            prepared.cellRetires.size(),sizeof(CellMap::value_type)));
    addAccountedBytes(prepared.contextBytesAdded,
        orderedNodeAccountedBytes(
            prepared.cellCreates.size(),sizeof(CellMap::value_type)));
    addAccountedBytes(prepared.contextBytesRemoved,
        affectedStarAccountedBytes(edgeStars,affectedEdges));
    addAccountedBytes(prepared.contextBytesRemoved,
        affectedStarAccountedBytes(faceOwners,affectedFaces));
    addAccountedBytes(prepared.contextBytesRemoved,
        affectedStarAccountedBytes(nodeStars,affectedNodes));
    addAccountedBytes(prepared.contextBytesAdded,
        starMapAccountedBytes(prepared.edgeRows));
    addAccountedBytes(prepared.contextBytesAdded,
        starMapAccountedBytes(prepared.faceRows));
    addAccountedBytes(prepared.contextBytesAdded,
        starMapAccountedBytes(prepared.nodeRows));
    addAccountedBytes(prepared.contextBytesAdded,
        orderedNodeAccountedBytes(1,sizeof(CoordinateMap::value_type)));
    addAccountedBytes(prepared.contextBytesAdded,
        orderedNodeAccountedBytes(1,sizeof(InsertedMap::value_type)));
    addAccountedBytes(prepared.contextBytesAdded,
        orderedNodeAccountedBytes(1,sizeof(CoordinateIndex::value_type)));
    const auto accountingRows=static_cast<std::uint64_t>(
        affectedEdges.size()+affectedFaces.size()+affectedNodes.size());
    recordConstructionWork(metrics,accountingRows);

    metrics.coordinatorRowsPrepared+=
        prepared.queueRetires.size()+prepared.queueCreates.size()+
        prepared.locationCreates.size();
    recordConstructionWork(metrics,static_cast<std::uint64_t>(
        prepared.queueRetires.size()+prepared.queueCreates.size()+
        prepared.locationCreates.size()));
    if (!queue.canApply(prepared.queueRetires,prepared.queueCreates)||
        !location.canApply(prepared.cellRetires,prepared.locationCreates))
        throw std::logic_error("unprepared coordinator capacity");

    // Charge the allocation-free tail before topology mutation.  Ordered-map
    // node transfer/find/erase and heap repair are logarithmically bounded;
    // the spatial-tree count is exact for the prepared leaf changes.
    std::uint64_t tailWork=0;
    const auto addTail=[&](std::uint64_t rows,std::size_t population) {
        const auto depth=static_cast<std::uint64_t>(
            std::bit_width(population+1))+1;
        std::uint64_t amount=0,next=0;
        if (!checkedMul(rows,depth,amount)||
            !checkedAdd(tailWork,amount,next))
            throw std::length_error("coordinator work accounting overflow");
        tailWork=next;
    };
    addTail(3,coordinates.size()+1);
    addTail(prepared.cellRetires.size()+prepared.cellCreates.size(),
            cells.size()+prepared.cellCreates.size());
    addTail(prepared.edgeErases.size()+prepared.edgeRows.size(),
            edgeStars.size()+prepared.edgeRows.size());
    addTail(prepared.faceErases.size()+prepared.faceRows.size(),
            faceOwners.size()+prepared.faceRows.size());
    addTail(prepared.nodeErases.size()+prepared.nodeRows.size(),
            nodeStars.size()+prepared.nodeRows.size());
    addTail(prepared.queueRetires.size()+prepared.queueCreates.size(),
            queue.capacity());
    const auto locationWork=location.preparedUpdateVisits(
        prepared.cellRetires,prepared.locationCreates);
    std::uint64_t charged=0;
    if (!checkedAdd(tailWork,locationWork,charged))
        throw std::length_error("location work accounting overflow");
    metrics.coordinatorTailWorkCharged+=tailWork;
    metrics.locationIndexUpdateWorkCharged+=locationWork;
    recordConstructionWork(metrics,charged);
}

template<class Map>
void transferPreparedMapNodes(Map& prepared,Map& destination) noexcept {
    while (!prepared.empty()) {
        auto node=prepared.extract(prepared.begin());
        const auto inserted=destination.insert(std::move(node));
        if (!inserted.inserted) std::abort();
    }
}

template<class Map,class Key>
void applyPreparedAffectedRows(Map& prepared,std::span<const Key> erased,
                               Map& destination) noexcept {
    for (const auto& key:erased) {
        const auto found=destination.find(key);
        if (found==destination.end()) std::abort();
        destination.erase(found);
    }
    for (auto it=prepared.begin();it!=prepared.end();) {
        const auto found=destination.find(it->first);
        if (found!=destination.end()) {
            found->second.swap(it->second);
            ++it;
        } else {
            auto current=it++;
            auto node=prepared.extract(current);
            const auto inserted=destination.insert(std::move(node));
            if (!inserted.inserted) std::abort();
        }
    }
}

void applyCoordinatorDelta(
    PreparedTetRefinementCoordinatorDelta& prepared,
    CoordinateMap& coordinates,InsertedMap& inserted,
    CoordinateIndex& coordinateIndex,CellMap& cells,
    EdgeStars& edgeStars,FaceOwners& faceOwners,NodeStars& nodeStars,
    ViolationHeap& queue,CurrentBlockedHeap& blocked,T3AReport& report,
    RefinementLocationIndex& location,
    TetRefinementMetrics& metrics) noexcept {
    transferPreparedMapNodes(prepared.coordinateCreates,coordinates);
    transferPreparedMapNodes(prepared.insertedCreates,inserted);
    transferPreparedMapNodes(prepared.coordinateIndexCreates,coordinateIndex);
    for (const auto cell:prepared.cellRetires) {
        const auto found=cells.find(cell);
        if (found==cells.end()) std::abort();
        cells.erase(found);
    }
    transferPreparedMapNodes(prepared.cellCreates,cells);
    applyPreparedAffectedRows(
        prepared.edgeRows,std::span<const NodePair>(prepared.edgeErases),edgeStars);
    applyPreparedAffectedRows(
        prepared.faceRows,std::span<const FaceKey>(prepared.faceErases),faceOwners);
    applyPreparedAffectedRows(
        prepared.nodeRows,
        std::span<const topo::NodeHandle>(prepared.nodeErases),nodeStars);
    for (const auto& invalidation:prepared.blockedInvalidations) {
        const auto rowIndex=blocked.erasePrepared(
            invalidation.cell,invalidation.blockedRowIndex);
        if (rowIndex>=report.blockedUsed) std::abort();
        auto& row=report.blockedStorage[static_cast<std::size_t>(rowIndex)];
        if (row.invalidatedByInsertionOrdinal!=kCurrentBlockedInsertion)
            std::abort();
        row.invalidatedByInsertionOrdinal=invalidation.insertionOrdinal;
        ++metrics.blockedRowsInvalidated;
        if (metrics.blockedRowsCurrent==0) std::abort();
        --metrics.blockedRowsCurrent;
    }
    queue.applyPrepared(prepared.queueRetires,prepared.queueCreates);
    location.applyPrepared(prepared.cellRetires,prepared.locationCreates,
                           metrics.locationIndexUpdateVisits);
    if (metrics.locationIndexUpdateVisits>
        metrics.locationIndexUpdateWorkCharged) std::abort();
}

bool selectedSizeEdgeAbsent(const TetViolationRow& selected,
                            const std::vector<BoundaryFace>& boundary,
                            const std::vector<Index>& existing,
                            const std::map<topo::NodeHandle,std::uint64_t>& inserted) {
    if (selected.key.kind!=TetViolationKind::Size) return true;
    for (const auto& face:boundary) {
        std::array<RefinementVertexKey,3> keys{};
        for (std::size_t i=0;i<3;++i)
            keys[i]=vertexKey(face.directed[i],existing,inserted);
        bool a=false,b=false;
        for (const auto& key:keys) {
            a = a || key==selected.key.sizeEdge.vertices[0];
            b = b || key==selected.key.sizeEdge.vertices[1];
        }
        if (a&&b) return false;
    }
    return true;
}

struct ProtectedQuery {
    bool invalid=false;
    std::vector<T3AProtectedFeatureWitness> witnesses;
};

struct DirectedInterval {
    double lo = 0.0;
    double hi = 0.0;
};

double down(double value) noexcept {
    return std::nextafter(value,
        -(std::numeric_limits<double>::infinity)());
}
double up(double value) noexcept {
    return std::nextafter(value,
        (std::numeric_limits<double>::infinity)());
}

std::optional<DirectedInterval> interval(double value,int shift) noexcept {
    const double scaled=std::scalbn(value,shift);
    if (!std::isfinite(scaled)) return std::nullopt;
    return DirectedInterval{down(scaled),up(scaled)};
}
std::optional<DirectedInterval> add(DirectedInterval a,
                                    DirectedInterval b) noexcept {
    const double lo=a.lo+b.lo,hi=a.hi+b.hi;
    if (!std::isfinite(lo)||!std::isfinite(hi)) return std::nullopt;
    return DirectedInterval{down(lo),up(hi)};
}
std::optional<DirectedInterval> sub(DirectedInterval a,
                                    DirectedInterval b) noexcept {
    const double lo=a.lo-b.hi,hi=a.hi-b.lo;
    if (!std::isfinite(lo)||!std::isfinite(hi)) return std::nullopt;
    return DirectedInterval{down(lo),up(hi)};
}
std::optional<DirectedInterval> mul(DirectedInterval a,
                                    DirectedInterval b) noexcept {
    const std::array<double,4> products{
        a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
    if (!std::all_of(products.begin(),products.end(),
                     [](double v){return std::isfinite(v);}))
        return std::nullopt;
    const auto [mn,mx]=std::minmax_element(products.begin(),products.end());
    return DirectedInterval{down(*mn),up(*mx)};
}
std::optional<DirectedInterval> divide(DirectedInterval a,
                                       DirectedInterval b) noexcept {
    if (!(b.lo>0.0) || !std::isfinite(b.hi)) return std::nullopt;
    return mul(a,{down(1.0/b.hi),up(1.0/b.lo)});
}

using IVec3=std::array<DirectedInterval,3>;

std::optional<DirectedInterval> l1RadiusBound(
    const IVec3& relative) noexcept {
    DirectedInterval result{};
    for (const auto component:relative) {
        const double magnitude=(std::max)(std::abs(component.lo),
                                           std::abs(component.hi));
        if (!std::isfinite(magnitude)) return std::nullopt;
        const auto next=add(result,{0.0,up(magnitude)});
        if (!next) return std::nullopt;
        result=*next;
    }
    result.lo=0.0;
    return result;
}

std::optional<IVec3> ivec(Vec3 value,int shift) noexcept {
    IVec3 result{};
    const std::array<double,3> values{value.x,value.y,value.z};
    for (std::size_t i=0;i<3;++i) {
        const auto component=interval(values[i],shift);
        if (!component) return std::nullopt;
        result[i]=*component;
    }
    return result;
}
std::optional<IVec3> iadd(const IVec3& a,const IVec3& b) noexcept {
    IVec3 result{};
    for (std::size_t i=0;i<3;++i) {
        const auto component=add(a[i],b[i]);
        if (!component) return std::nullopt;
        result[i]=*component;
    }
    return result;
}
std::optional<IVec3> isub(const IVec3& a,const IVec3& b) noexcept {
    IVec3 result{};
    for (std::size_t i=0;i<3;++i) {
        const auto component=sub(a[i],b[i]);
        if (!component) return std::nullopt;
        result[i]=*component;
    }
    return result;
}
std::optional<IVec3> icross(const IVec3& a,const IVec3& b) noexcept {
    IVec3 result{};
    static constexpr std::array<std::array<std::size_t,2>,3> ij{{
        {{1,2}},{{2,0}},{{0,1}}}};
    for (std::size_t k=0;k<3;++k) {
        const auto i=ij[k][0],j=ij[k][1];
        const auto lhs=mul(a[i],b[j]);
        const auto rhs=mul(a[j],b[i]);
        if (!lhs||!rhs) return std::nullopt;
        const auto component=sub(*lhs,*rhs);
        if (!component) return std::nullopt;
        result[k]=*component;
    }
    return result;
}
std::optional<DirectedInterval> idot(const IVec3& a,
                                     const IVec3& b) noexcept {
    DirectedInterval result{};
    for (std::size_t i=0;i<3;++i) {
        const auto product=mul(a[i],b[i]);
        if (!product) return std::nullopt;
        const auto next=add(result,*product);
        if (!next) return std::nullopt;
        result=*next;
    }
    return result;
}
std::optional<IVec3> iscale(const IVec3& a,
                            DirectedInterval scale) noexcept {
    IVec3 result{};
    for (std::size_t i=0;i<3;++i) {
        const auto component=mul(a[i],scale);
        if (!component) return std::nullopt;
        result[i]=*component;
    }
    return result;
}

struct ProtectedAabb {
    Vec3 lo{},hi{};
};
struct ProtectedIndexRow {
    ProtectedAabb bounds{};
    T3AProtectedFeatureWitness witness{};
    std::array<topo::NodeHandle,3> nodes{};
    std::uint8_t nodeCount=0;
};
struct ProtectedIndexNode {
    ProtectedAabb bounds{};
    std::size_t begin=0,end=0,left=(std::numeric_limits<std::size_t>::max)(),
        right=(std::numeric_limits<std::size_t>::max)();
    bool leaf=false;
};
struct ProtectedFeatureIndex {
    std::vector<ProtectedIndexRow> rows;
    std::vector<ProtectedIndexNode> nodes;
    std::size_t root=(std::numeric_limits<std::size_t>::max)();
};
enum class ProtectedIndexBuild : std::uint8_t { Built,Invalid,Range };

int normalizationShift(std::span<const Vec3> points) noexcept {
    double magnitude=0.0;
    for (const auto p:points)
        magnitude=(std::max)({magnitude,std::abs(p.x),std::abs(p.y),
                              std::abs(p.z)});
    if (magnitude==0.0) return 0;
    int exponent=0;
    (void)std::frexp(magnitude,&exponent);
    return -exponent;
}

std::optional<ProtectedAabb> sphereBounds(
    std::span<const Vec3> points) noexcept {
    if (points.size()!=2&&points.size()!=3) return std::nullopt;
    const int shift=normalizationShift(points);
    const auto a=ivec(points[0],shift),b=ivec(points[1],shift);
    if (!a||!b) return std::nullopt;
    IVec3 center{},relative{};
    DirectedInterval radius{};
    if (points.size()==2) {
        const auto sum=iadd(*a,*b),delta=isub(*b,*a);
        if (!sum||!delta) return std::nullopt;
        const DirectedInterval half{0.5,0.5};
        const auto c=iscale(*sum,half);
        const auto length=l1RadiusBound(*delta);
        if (!c||!length) return std::nullopt;
        const auto r=mul(*length,half);
        if (!r) return std::nullopt;
        center=*c; radius=*r;
    } else {
        const auto c=ivec(points[2],shift);
        if (!c) return std::nullopt;
        const auto u=isub(*b,*a),v=isub(*c,*a);
        if (!u||!v) return std::nullopt;
        const auto n=icross(*u,*v);
        if (!n) return std::nullopt;
        const auto nn=idot(*n,*n),uu=idot(*u,*u),vv=idot(*v,*v);
        const auto vxn=icross(*v,*n),nxu=icross(*n,*u);
        if (!nn||!uu||!vv||!vxn||!nxu) return std::nullopt;
        const auto left=iscale(*vxn,*uu),right=iscale(*nxu,*vv);
        if (!left||!right) return std::nullopt;
        const auto numerator=iadd(*left,*right);
        const auto denominator=mul(*nn,{2.0,2.0});
        if (!numerator||!denominator||!(denominator->lo>0.0))
            return std::nullopt;
        for (std::size_t i=0;i<3;++i) {
            const auto component=divide((*numerator)[i],*denominator);
            if (!component) return std::nullopt;
            relative[i]=*component;
        }
        const auto absolute=iadd(*a,relative);
        const auto r=l1RadiusBound(relative);
        if (!absolute) return std::nullopt;
        if (!r) return std::nullopt;
        center=*absolute; radius=*r;
    }
    ProtectedAabb result;
    std::array<double,3> lo{},hi{};
    for (std::size_t i=0;i<3;++i) {
        const auto low=sub(center[i],{radius.hi,radius.hi});
        const auto high=add(center[i],{radius.hi,radius.hi});
        if (!low||!high) return std::nullopt;
        lo[i]=std::scalbn(low->lo,-shift);
        hi[i]=std::scalbn(high->hi,-shift);
        if (!std::isfinite(lo[i])||!std::isfinite(hi[i]))
            return std::nullopt;
        lo[i]=down(lo[i]); hi[i]=up(hi[i]);
        if (!std::isfinite(lo[i])||!std::isfinite(hi[i]))
            return std::nullopt;
    }
    result.lo={lo[0],lo[1],lo[2]};
    result.hi={hi[0],hi[1],hi[2]};
    return result;
}

ProtectedIndexBuild buildProtectedIndex(
    const bnd::TetFacetSetState& state,const topo::TopologyWorkspace& ws,
    ProtectedFeatureIndex& index,TetRefinementMetrics& metrics) {
    index.rows.reserve(state.edgeBindings.size()+state.bindings.size());
    for (const auto& [key,binding]:state.edgeBindings) {
        if (!ws.edges().isCurrent(binding.edge)) return ProtectedIndexBuild::Invalid;
        const auto edge=ws.edges().get(binding.edge);
        const std::array<Vec3,2> points{ws.nodes().get(edge.a).coord,
                                       ws.nodes().get(edge.b).coord};
        ++metrics.protectedEdgeTests;
        const auto valid=predicates::edgeDiametralPower(points[0],points[1],
                                                        points[0]);
        if (valid.status!=predicates::PowerPredicateStatus::Valid)
            return ProtectedIndexBuild::Invalid;
        const auto bounds=sphereBounds(points);
        if (!bounds) return ProtectedIndexBuild::Range;
        ProtectedIndexRow row;
        row.bounds=*bounds;
        row.witness={T3AProtectedFeatureKind::Edge,key,{}};
        row.nodes={edge.a,edge.b,{}}; row.nodeCount=2;
        index.rows.push_back(row);
    }
    for (const auto& [key,binding]:state.bindings) {
        if (!ws.faces().isCurrent(binding.recoveredFace))
            return ProtectedIndexBuild::Invalid;
        const std::array<Vec3,3> points{
            ws.nodes().get(binding.canonicalVertices[0]).coord,
            ws.nodes().get(binding.canonicalVertices[1]).coord,
            ws.nodes().get(binding.canonicalVertices[2]).coord};
        ++metrics.protectedFacetTests;
        const auto valid=predicates::triangleCircumballPower(
            points[0],points[1],points[2],points[0]);
        if (valid.status!=predicates::PowerPredicateStatus::Valid)
            return ProtectedIndexBuild::Invalid;
        const auto bounds=sphereBounds(points);
        if (!bounds) return ProtectedIndexBuild::Range;
        ProtectedIndexRow row;
        row.bounds=*bounds;
        row.witness={T3AProtectedFeatureKind::Facet,{},key};
        row.nodes=binding.canonicalVertices; row.nodeCount=3;
        index.rows.push_back(row);
    }
    std::sort(index.rows.begin(),index.rows.end(),[](const auto& a,const auto& b) {
        if (a.bounds.lo.x!=b.bounds.lo.x) return a.bounds.lo.x<b.bounds.lo.x;
        if (a.witness.kind!=b.witness.kind) return a.witness.kind<b.witness.kind;
        return a.witness.kind==T3AProtectedFeatureKind::Edge
            ? a.witness.edge<b.witness.edge:a.witness.facet<b.witness.facet;
    });
    if (!index.rows.empty()) {
        index.nodes.reserve(index.rows.size()*2);
        const auto buildNode=[&](auto&& self,std::size_t begin,
                                 std::size_t end)->std::size_t {
            ProtectedIndexNode node;
            node.begin=begin; node.end=end;
            node.bounds=index.rows[begin].bounds;
            for (std::size_t i=begin+1;i<end;++i) {
                node.bounds.lo.x=(std::min)(node.bounds.lo.x,index.rows[i].bounds.lo.x);
                node.bounds.lo.y=(std::min)(node.bounds.lo.y,index.rows[i].bounds.lo.y);
                node.bounds.lo.z=(std::min)(node.bounds.lo.z,index.rows[i].bounds.lo.z);
                node.bounds.hi.x=(std::max)(node.bounds.hi.x,index.rows[i].bounds.hi.x);
                node.bounds.hi.y=(std::max)(node.bounds.hi.y,index.rows[i].bounds.hi.y);
                node.bounds.hi.z=(std::max)(node.bounds.hi.z,index.rows[i].bounds.hi.z);
            }
            const std::size_t slot=index.nodes.size();
            index.nodes.push_back(node);
            if (end-begin<=8) {
                index.nodes[slot].leaf=true;
                return slot;
            }
            const std::array<double,3> extent{
                node.bounds.hi.x-node.bounds.lo.x,
                node.bounds.hi.y-node.bounds.lo.y,
                node.bounds.hi.z-node.bounds.lo.z};
            std::size_t axis=0;
            if (extent[1]>extent[axis]) axis=1;
            if (extent[2]>extent[axis]) axis=2;
            const auto centre=[axis](const ProtectedIndexRow& row) noexcept {
                if (axis==0) return std::midpoint(row.bounds.lo.x,row.bounds.hi.x);
                if (axis==1) return std::midpoint(row.bounds.lo.y,row.bounds.hi.y);
                return std::midpoint(row.bounds.lo.z,row.bounds.hi.z);
            };
            const auto witnessLess=[](const auto& a,const auto& b) noexcept {
                if (a.kind!=b.kind) return a.kind<b.kind;
                return a.kind==T3AProtectedFeatureKind::Edge
                    ? a.edge<b.edge:a.facet<b.facet;
            };
            std::sort(index.rows.begin()+static_cast<std::ptrdiff_t>(begin),
                      index.rows.begin()+static_cast<std::ptrdiff_t>(end),
                [&](const auto& a,const auto& b) {
                    const double ac=centre(a),bc=centre(b);
                    return ac<bc||(ac==bc&&witnessLess(a.witness,b.witness));
                });
            const std::size_t middle=begin+(end-begin)/2;
            const std::size_t left=self(self,begin,middle);
            const std::size_t right=self(self,middle,end);
            index.nodes[slot].left=left;
            index.nodes[slot].right=right;
            return slot;
        };
        index.root=buildNode(buildNode,0,index.rows.size());
    }
    return ProtectedIndexBuild::Built;
}

void protectedEncroachment(
    Vec3 point,const ProtectedFeatureIndex& index,
    const topo::TopologyWorkspace& ws,TetRefinementMetrics& metrics,
    ProtectedQuery& result,std::vector<std::size_t>& stack) {
    result.invalid=false;
    result.witnesses.clear();
    stack.clear();
    // A non-finite proposal cannot be classified by the finite AABB broad
    // phase.  Reject it before broad-phase comparisons can cull it as though
    // it were simply outside every protected feature.
    if (!finite(point)) {
        result.invalid=true;
        return;
    }
    if (index.root!=(std::numeric_limits<std::size_t>::max)())
        stack.push_back(index.root);
    while (!stack.empty()) {
        const std::size_t nodeIndex=stack.back();
        stack.pop_back();
        const auto& node=index.nodes[nodeIndex];
        recordConstructionWork(metrics);
        ++metrics.protectedIndexVisits;
        if (point.x<node.bounds.lo.x||point.x>node.bounds.hi.x||
            point.y<node.bounds.lo.y||point.y>node.bounds.hi.y||
            point.z<node.bounds.lo.z||point.z>node.bounds.hi.z) continue;
        if (!node.leaf) {
            // Push right first so the canonical left subtree is visited first.
            stack.push_back(node.right);
            stack.push_back(node.left);
            continue;
        }
        for (std::size_t i=node.begin;i<node.end;++i) {
            const auto& row=index.rows[i];
            recordConstructionWork(metrics);
            ++metrics.protectedIndexVisits;
            if (point.x<row.bounds.lo.x||point.x>row.bounds.hi.x||
                point.y<row.bounds.lo.y||point.y>row.bounds.hi.y||
                point.z<row.bounds.lo.z||point.z>row.bounds.hi.z) continue;
            predicates::PowerPredicateResult power{};
            if (row.nodeCount==2) {
                recordConstructionWork(metrics);
                ++metrics.protectedEdgeTests;
                power=predicates::edgeDiametralPower(
                    ws.nodes().get(row.nodes[0]).coord,
                    ws.nodes().get(row.nodes[1]).coord,point);
            } else {
                recordConstructionWork(metrics);
                ++metrics.protectedFacetTests;
                power=predicates::triangleCircumballPower(
                    ws.nodes().get(row.nodes[0]).coord,
                    ws.nodes().get(row.nodes[1]).coord,
                    ws.nodes().get(row.nodes[2]).coord,point);
            }
            if (power.status!=predicates::PowerPredicateStatus::Valid) {
                result.invalid=true;
                return;
            }
            if (power.sign!=predicates::Sign::Positive)
                result.witnesses.push_back(row.witness);
        }
    }
    std::sort(result.witnesses.begin(),result.witnesses.end(),
              [](const auto& a,const auto& b) {
                  if (a.kind!=b.kind) return a.kind<b.kind;
                  return a.kind==T3AProtectedFeatureKind::Edge
                             ? a.edge<b.edge : a.facet<b.facet;
              });
}

T3ATopologyTransition publicTransition(topo::TetRefinementTransition value) {
    switch (value) {
        case topo::TetRefinementTransition::Retained:
            return T3ATopologyTransition::Retained;
        case topo::TetRefinementTransition::Created:
            return T3ATopologyTransition::Created;
        case topo::TetRefinementTransition::Retired:
            return T3ATopologyTransition::Retired;
    }
    return T3ATopologyTransition::Retained;
}

bool finalProtectedAudit(const bnd::TetFacetSetState& state,
                         const topo::TopologyWorkspace& ws,
                         const std::vector<Index>& globalByNodeSlot) {
    if (!(state.workspaceStamp==topo::WorkspaceStamp::of(ws))) return false;
    for (const auto& [key,binding]:state.edgeBindings) {
        if (!(key==binding.key)) return false;
        if (!ws.edges().isCurrent(binding.edge)) return false;
        const auto edge=ws.edges().get(binding.edge);
        if (!edge.requiredConstraint || !sameClass(edge.cls,binding.classification))
            return false;
        const auto globalFor=[&](topo::NodeHandle node) {
            const std::size_t slot=static_cast<std::size_t>(node.slot);
            return slot<globalByNodeSlot.size()
                ? globalByNodeSlot[slot]:kInvalidIndex;
        };
        auto endpoints=std::array<Index,2>{globalFor(edge.a),globalFor(edge.b)};
        if (endpoints[0]==kInvalidIndex||endpoints[1]==kInvalidIndex)
            return false;
        std::sort(endpoints.begin(),endpoints.end());
        if (endpoints!=key.canonicalGlobalNodes) return false;
    }
    for (const auto& [key,binding]:state.bindings) {
        if (!(key==binding.key)) return false;
        if (!ws.faces().isCurrent(binding.recoveredFace)) return false;
        const auto face=ws.faces().get(binding.recoveredFace);
        if (!sameClass(face.cls,binding.classification) ||
            face.incident.size()!=1) return false;
        auto canonical=std::array<topo::NodeHandle,3>{
            face.corners[0],face.corners[1],face.corners[2]};
        std::sort(canonical.begin(),canonical.end(),[&](const auto a,const auto b) {
            return ws.nodes().seqOf(a)<ws.nodes().seqOf(b);
        });
        if (canonical!=binding.canonicalVertices) return false;
    }
    return true;
}

bool auditCommittedDelta(const topo::TetRefinementCommit& commit,
                         const topo::TopologyWorkspace& ws) {
    if (!(commit.stamp==topo::WorkspaceStamp::of(ws))||
        commit.actualNode!=commit.predictedNode||
        !ws.nodes().isCurrent(commit.actualNode)) return false;
    const auto resolve=[&](const topo::TetRefinementNodeRef& ref) noexcept {
        return ref.kind==topo::TetRefinementNodeRefKind::Staged
            ? commit.actualNode:ref.existing;
    };
    for (const auto& row:commit.cells) {
        if (row.actual!=row.predicted) return false;
        const bool retired=row.transition==topo::TetRefinementTransition::Retired;
        if (ws.cells().isCurrent(row.actual)==retired) return false;
        if (retired) continue;
        const auto data=ws.cells().get(row.actual);
        if (data.type!=topo::ElementType::Tet4||
            !sameClass(data.cls,row.classification)) return false;
        std::array<topo::NodeHandle,4> expected{},actual{
            data.nodes[0],data.nodes[1],data.nodes[2],data.nodes[3]};
        for (std::size_t i=0;i<4;++i) expected[i]=resolve(row.key[i]);
        std::sort(expected.begin(),expected.end());
        std::sort(actual.begin(),actual.end());
        if (expected!=actual) return false;
    }
    for (const auto& row:commit.edges) {
        if (row.actual!=row.predicted) return false;
        const bool retired=row.transition==topo::TetRefinementTransition::Retired;
        if (ws.edges().isCurrent(row.actual)==retired) return false;
        if (retired) continue;
        const auto data=ws.edges().get(row.actual);
        auto expected=std::array<topo::NodeHandle,2>{
            resolve(row.key[0]),resolve(row.key[1])};
        auto actual=std::array<topo::NodeHandle,2>{data.a,data.b};
        std::sort(expected.begin(),expected.end());
        std::sort(actual.begin(),actual.end());
        if (expected!=actual||data.useCount!=row.finalUseCount||
            data.requiredConstraint!=row.required||
            !sameClass(data.cls,row.classification)) return false;
    }
    for (const auto& row:commit.faces) {
        if (row.actual!=row.predicted) return false;
        const bool retired=row.transition==topo::TetRefinementTransition::Retired;
        if (ws.faces().isCurrent(row.actual)==retired) return false;
        if (retired) continue;
        const auto data=ws.faces().get(row.actual);
        auto expected=std::array<topo::NodeHandle,3>{
            resolve(row.key[0]),resolve(row.key[1]),resolve(row.key[2])};
        auto actual=std::array<topo::NodeHandle,3>{
            data.corners[0],data.corners[1],data.corners[2]};
        std::sort(expected.begin(),expected.end());
        std::sort(actual.begin(),actual.end());
        if (expected!=actual||data.incident.size()!=row.finalOwnerCount||
            !sameClass(data.cls,row.classification)) return false;
        for (std::size_t i=0;i<row.finalOwnerCount;++i) {
            bool found=false;
            for (const auto& incidence:data.incident)
                found=found||
                    (incidence.cellSlot==row.finalOwners[i].slot&&
                     incidence.cellGen==row.finalOwners[i].gen&&
                     incidence.sense==row.finalSenses[i]);
            if (!found) return false;
        }
    }
    return true;
}

bool buildImmutableBoundary(
    const SurfaceSetOutcome& surface,const bnd::TetFacetSetState& state,
    const mesh::Mesh& publishedSeed,
    std::vector<T3ABoundaryTriRecord>& destination) {
    if (!surface.retainedArtifact) return false;
    using SeedFaceKey=std::array<Index,3>;
    struct SeedFace {
        SeedFaceKey directed{};
        mesh::Classification classification{};
    };
    std::map<SeedFaceKey,SeedFace> seedFaces;
    for (Index ci=0;ci<publishedSeed.cellCount();++ci) {
        const auto cell=publishedSeed.cell(ci);
        if (cell.type!=mesh::CellType::Tri3) continue;
        const Index* nodes=publishedSeed.cellNodes(ci);
        SeedFace row{{nodes[0],nodes[1],nodes[2]},cell.onEntity};
        auto key=row.directed;
        std::sort(key.begin(),key.end());
        if (!seedFaces.emplace(key,row).second) return false;
    }
    destination.clear();
    destination.reserve(state.bindings.size());
    for (const auto& [key,binding]:state.bindings) {
        const SurfaceFaceResult* face=nullptr;
        for (const auto& candidate:surface.retainedArtifact->faces)
            if (candidate&&candidate->faceStamp.persistentId==
                    key.facePersistentId) {
                face=candidate.get(); break;
            }
        if (!face||key.normalizedLocalCellOrdinal>=face->localCells.size())
            return false;
        const auto& local=face->localCells[key.normalizedLocalCellOrdinal];
        T3ABoundaryTriRecord row;
        row.facetKey=key;
        row.classification=binding.classification;
        row.authoritativeFace=binding.recoveredFace;
        SeedFaceKey surfaceDirected{};
        for (std::size_t i=0;i<3;++i) {
            if (local.nodes[i]>=face->publishedGlobalByLocalNode.size())
                return false;
            surfaceDirected[i]=
                face->publishedGlobalByLocalNode[local.nodes[i]];
            if (surfaceDirected[i]==kInvalidIndex) return false;
        }
        auto canonical=surfaceDirected;
        std::sort(canonical.begin(),canonical.end());
        const auto seedFace=seedFaces.find(canonical);
        if (seedFace==seedFaces.end()||
            !sameClass(seedFace->second.classification,binding.classification))
            return false;
        row.directedGlobalNodes=seedFace->second.directed;
        destination.push_back(row);
    }
    return destination.size()==seedFaces.size();
}

bool bindConstraintStateToAuthoritative(
    bnd::TetFacetSetState& state,std::vector<T3ABoundaryTriRecord>& boundary,
    const std::vector<Index>& globalBySlot,topo::TopologyWorkspace& ws) {
    std::map<Index,topo::NodeHandle> localByGlobal;
    for (const auto node:ws.nodes().activeBySequence()) {
        const std::size_t slot=static_cast<std::size_t>(node.slot);
        if (slot>=globalBySlot.size()||globalBySlot[slot]==kInvalidIndex||
            !localByGlobal.emplace(globalBySlot[slot],node).second)
            return false;
    }
    for (auto& [key,binding]:state.edgeBindings) {
        const auto a=localByGlobal.find(key.canonicalGlobalNodes[0]);
        const auto b=localByGlobal.find(key.canonicalGlobalNodes[1]);
        if (a==localByGlobal.end()||b==localByGlobal.end()) return false;
        const auto edge=ws.edges().find(ws.nodes(),a->second,b->second);
        if (!edge) return false;
        const auto data=ws.edges().get(*edge);
        if (!data.requiredConstraint||!sameClass(data.cls,binding.classification))
            return false;
        binding.edge=*edge;
    }
    for (auto& [key,binding]:state.bindings) {
        const auto row=std::find_if(boundary.begin(),boundary.end(),
            [&](const auto& value){return value.facetKey==key;});
        if (row==boundary.end()) return false;
        std::array<topo::NodeHandle,3> nodes{};
        for (std::size_t i=0;i<3;++i) {
            const auto found=localByGlobal.find(row->directedGlobalNodes[i]);
            if (found==localByGlobal.end()) return false;
            nodes[i]=found->second;
        }
        const auto face=ws.faces().find(ws.nodes(),nodes);
        if (!face) return false;
        const auto data=ws.faces().get(*face);
        if (!sameClass(data.cls,binding.classification)||data.incident.size()!=1)
            return false;
        std::sort(nodes.begin(),nodes.end(),[&](const auto a,const auto b) {
            return ws.nodes().seqOf(a)<ws.nodes().seqOf(b);
        });
        binding.canonicalVertices=nodes;
        binding.recoveredFace=*face;
        row->authoritativeFace=*face;
    }
    state.workspaceStamp=topo::WorkspaceStamp::of(ws);
    return true;
}

void copyReceipt(BoundaryConformingV3Receipt& receipt,
                 BoundaryConformingV3Outcome& source) noexcept {
    receipt.status=source.status;
    receipt.failure=source.failure;
    receipt.destinationUnchanged=source.destinationUnchanged;
    receipt.reportComplete=source.reportComplete;
    receipt.surfaceOutcome=std::move(source.surfaceOutcome);
    receipt.report=std::move(source.report);
    if (source.retainedResult) {
        receipt.artifactPreflightMetrics=source.retainedResult->artifactPreflightMetrics;
        receipt.recoveryMetrics=source.retainedResult->recoveryMetrics;
        receipt.extractionMetrics=source.retainedResult->extractionMetrics;
        receipt.memoryMetrics=source.retainedResult->memoryMetrics;
        receipt.boundaryAudit=source.retainedResult->boundaryAudit;
    }
}

T3AStatus sourceStatus(BoundaryConformingV3Status status) noexcept {
    switch(status) {
        case BoundaryConformingV3Status::Rejected: return T3AStatus::Rejected;
        case BoundaryConformingV3Status::Unsupported: return T3AStatus::Unsupported;
        case BoundaryConformingV3Status::Refused: return T3AStatus::Incomplete;
        case BoundaryConformingV3Status::Failed: return T3AStatus::Failed;
        case BoundaryConformingV3Status::Published: return T3AStatus::Published;
    }
    return T3AStatus::Failed;
}
T3AReason sourceReason(BoundaryConformingV3Status status) noexcept {
    switch(status) {
        case BoundaryConformingV3Status::Rejected: return T3AReason::SourceRejected;
        case BoundaryConformingV3Status::Unsupported: return T3AReason::SourceRefused;
        case BoundaryConformingV3Status::Refused: return T3AReason::SourceRefused;
        case BoundaryConformingV3Status::Failed: return T3AReason::SourceFailed;
        case BoundaryConformingV3Status::Published: return T3AReason::None;
    }
    return T3AReason::SourceFailed;
}

bool fillUnattempted(T3AReport& report,const ViolationHeap& queue) noexcept {
    if (queue.size()>report.unattemptedStorage.size()) return false;
    report.unattemptedUsed=queue.size();
    std::size_t i=0;
    for (const auto& row:queue.rows())
        report.unattemptedStorage[i++]=row.row.key;
    std::sort(report.unattemptedStorage.begin(),
              report.unattemptedStorage.begin()+
                  static_cast<std::ptrdiff_t>(report.unattemptedUsed),
              violationKeyLess);
    return true;
}

bool buildPublicationMesh(
    const mesh::Mesh& seed,const topo::TopologyWorkspace& ws,
    const std::vector<Index>& existing,
    const std::map<topo::NodeHandle,std::uint64_t>& inserted,
    mesh::Mesh& out,std::vector<Index>& nodeMap,std::vector<Index>& cellMap,
    TetRefinementQualityReport& quality) {
    (void)inserted;
    const auto activeNodes=ws.nodes().activeBySequence();
    const auto activeCells=ws.cells().activeBySequence();
    std::size_t nodeSlots=0,cellSlots=0;
    for (const auto h:activeNodes) nodeSlots=(std::max)(nodeSlots,
        static_cast<std::size_t>(h.slot)+1);
    for (const auto h:activeCells) cellSlots=(std::max)(cellSlots,
        static_cast<std::size_t>(h.slot)+1);
    nodeMap.assign(nodeSlots,kInvalidIndex);
    cellMap.assign(cellSlots,kInvalidIndex);

    std::size_t nonTet=0,nonTetConn=0;
    for (Index i=0;i<seed.cellCount();++i)
        if (seed.cell(i).type!=mesh::CellType::Tet4) {
            ++nonTet; nonTetConn+=static_cast<std::size_t>(
                mesh::nodesPerCell(seed.cell(i).type));
        }
    out.reserveAppend(static_cast<std::size_t>(seed.nodeCount())+
                          activeNodes.size(),
                      nonTet+activeCells.size(),
                      nonTetConn+4*activeCells.size());
    for (Index i=0;i<seed.nodeCount();++i) out.addNode(seed.node(i),seed.nodeClass(i));
    for (Index i=0;i<seed.cellCount();++i) {
        const auto cell=seed.cell(i);
        if (cell.type==mesh::CellType::Tet4) continue;
        const int n=mesh::nodesPerCell(cell.type);
        std::vector<Index> conn(static_cast<std::size_t>(n));
        const Index* src=seed.cellNodes(i);
        for (int j=0;j<n;++j) conn[static_cast<std::size_t>(j)]=src[j];
        out.addCell(cell.type,conn,cell.onEntity);
    }
    for (const auto node:activeNodes) {
        const auto data=ws.nodes().get(node);
        const std::size_t slot=static_cast<std::size_t>(node.slot);
        if (slot<existing.size() && existing[slot]!=kInvalidIndex) {
            const Index global=existing[slot];
            if (global<0 || global>=seed.nodeCount() ||
                seed.node(global).x!=data.coord.x ||
                seed.node(global).y!=data.coord.y ||
                seed.node(global).z!=data.coord.z ||
                !sameClass(seed.nodeClass(global),data.cls)) return false;
            nodeMap[slot]=global;
        } else {
            nodeMap[slot]=out.addNode(data.coord,data.cls);
        }
    }
    struct CellRow {
        topo::CellHandle handle{};
        std::array<Index,4> nodes{};
        mesh::Classification cls{};
    };
    std::vector<CellRow> rows;
    rows.reserve(activeCells.size());
    for (const auto cell:activeCells) {
        const auto data=ws.cells().get(cell);
        if (data.type!=topo::ElementType::Tet4) continue;
        CellRow row; row.handle=cell; row.cls=data.cls;
        for (std::size_t i=0;i<4;++i) {
            const std::size_t slot=static_cast<std::size_t>(data.nodes[i].slot);
            if (slot>=nodeMap.size() || nodeMap[slot]==kInvalidIndex) return false;
            row.nodes[i]=nodeMap[slot];
        }
        auto key=row.nodes; std::sort(key.begin(),key.end());
        rows.push_back(row);
    }
    std::sort(rows.begin(),rows.end(),[](const auto& a,const auto& b) {
        auto ka=a.nodes,kb=b.nodes;
        std::sort(ka.begin(),ka.end()); std::sort(kb.begin(),kb.end());
        return ka<kb;
    });
    for (const auto& row:rows) {
        const Index global=out.addCell(mesh::CellType::Tet4,
            std::vector<Index>(row.nodes.begin(),row.nodes.end()),row.cls);
        cellMap[static_cast<std::size_t>(row.handle.slot)]=global;
    }
    const auto qr=mesh::analyzeTets(out);
    quality.nodeCount=static_cast<std::uint64_t>(out.nodeCount());
    quality.tetCount=static_cast<std::uint64_t>(qr.nTets);
    quality.inverted=static_cast<std::uint64_t>(qr.nInverted);
    quality.nonfinite=static_cast<std::uint64_t>(qr.nNonFinite);
    quality.exactPositive=quality.tetCount-quality.inverted-quality.nonfinite;
    quality.minimumSignedVolume=qr.nTets?qr.totalVolume/static_cast<double>(qr.nTets):0.0;
    quality.minimumRadiusRatio=qr.radiusRatio.n?qr.radiusRatio.min:0.0;
    for (std::size_t i=0;i<quality.radiusRatioHistogram.size();++i)
        quality.radiusRatioHistogram[i]=static_cast<std::uint64_t>(qr.radiusRatioHist[i]);
    return qr.nInverted==0 && qr.nNonFinite==0;
}

using GlobalFaceKey=std::array<Index,3>;
GlobalFaceKey globalFaceKey(GlobalFaceKey value) noexcept {
    std::sort(value.begin(),value.end());
    return value;
}

bool sameDirectedCycle(const GlobalFaceKey& a,const GlobalFaceKey& b) noexcept {
    for (std::size_t shift=0;shift<3;++shift)
        if (a[0]==b[shift]&&a[1]==b[(shift+1)%3]&&
            a[2]==b[(shift+2)%3]) return true;
    return false;
}
bool oppositeDirectedCycle(const GlobalFaceKey& a,
                           const GlobalFaceKey& b) noexcept {
    for (std::size_t shift=0;shift<3;++shift)
        if (a[0]==b[shift]&&a[1]==b[(shift+2)%3]&&
            a[2]==b[(shift+1)%3]) return true;
    return false;
}

bool auditBoundaryBijection(
    const mesh::Mesh& mesh,
    const std::vector<T3ABoundaryTriRecord>& immutable,
    const topo::TopologyWorkspace& workspace,
    const std::vector<Index>& globalByNodeSlot) {
    std::map<GlobalFaceKey,T3ABoundaryTriRecord> expected;
    for (const auto& row:immutable)
        if (!expected.emplace(globalFaceKey(row.directedGlobalNodes),row).second)
            return false;
    struct Owners {
        std::uint32_t count=0;
        std::array<GlobalFaceKey,2> directed{};
    };
    std::map<GlobalFaceKey,Owners> incidence;
    for (Index ci=0;ci<mesh.cellCount();++ci) {
        const auto cell=mesh.cell(ci);
        if (cell.type!=mesh::CellType::Tet4) continue;
        const Index* nodes=mesh.cellNodes(ci);
        for (const auto& local:kLocalTetFaces) {
            const GlobalFaceKey directed{nodes[local[0]],nodes[local[1]],
                                         nodes[local[2]]};
            auto& owners=incidence[globalFaceKey(directed)];
            if (owners.count>=2) return false;
            owners.directed[owners.count++]=directed;
        }
    }
    std::size_t boundaryCount=0;
    for (const auto& [key,owners]:incidence) {
        if (owners.count==2) {
            if (!oppositeDirectedCycle(owners.directed[0],owners.directed[1]))
                return false;
            continue;
        }
        ++boundaryCount;
        const auto found=expected.find(key);
        if (found==expected.end()||
            !workspace.faces().isCurrent(found->second.authoritativeFace))
            return false;
        const auto protectedFace=workspace.faces().get(
            found->second.authoritativeFace);
        if (protectedFace.cornerCount!=3||protectedFace.incident.size()!=1||
            protectedFace.incident[0].sense==0) return false;
        GlobalFaceKey protectedDirected{};
        for (std::size_t i=0;i<3;++i) {
            const auto node=protectedFace.corners[i];
            const std::size_t slot=static_cast<std::size_t>(node.slot);
            if (slot>=globalByNodeSlot.size()||
                globalByNodeSlot[slot]==kInvalidIndex) return false;
            protectedDirected[i]=globalByNodeSlot[slot];
        }
        if (protectedFace.incident[0].sense<0)
            std::swap(protectedDirected[1],protectedDirected[2]);
        if (!sameDirectedCycle(owners.directed[0],protectedDirected))
            return false;
    }
    if (boundaryCount!=expected.size()) return false;

    std::map<GlobalFaceKey,std::pair<GlobalFaceKey,mesh::Classification>> surface;
    for (Index ci=0;ci<mesh.cellCount();++ci) {
        const auto cell=mesh.cell(ci);
        if (cell.type!=mesh::CellType::Tri3) continue;
        const Index* nodes=mesh.cellNodes(ci);
        const GlobalFaceKey directed{nodes[0],nodes[1],nodes[2]};
        if (!surface.emplace(globalFaceKey(directed),
                             std::make_pair(directed,cell.onEntity)).second)
            return false;
    }
    if (surface.size()!=expected.size()) return false;
    for (const auto& [key,row]:expected) {
        const auto found=surface.find(key);
        if (found==surface.end()||
            !sameDirectedCycle(found->second.first,row.directedGlobalNodes)||
            !sameClass(found->second.second,row.classification))
            return false;
    }
    return true;
}

bool completeQualityReport(
    const mesh::Mesh& mesh,const sizing::MeshSpec& spec,
    const TetRefinementPolicy& policy,SizingLedger& ledger,
    TetRefinementMetrics& metrics,TetRefinementQualityReport& quality) {
    std::vector<double> radiusEdges;
    std::vector<double> radiusRatios;
    radiusEdges.reserve(static_cast<std::size_t>(mesh.cellCount()));
    radiusRatios.reserve(static_cast<std::size_t>(mesh.cellCount()));
    quality.minimumSignedVolume=(std::numeric_limits<double>::infinity)();
    quality.maximumRadiusEdge=0.0;
    quality.worstNormalizedSize=0.0;
    quality.sizeViolations=0;
    quality.sliversRemaining=0;
    quality.radiusEdgeHistogram.fill(0);
    quality.radiusRatioHistogram.fill(0);
    static constexpr std::array<std::array<int,2>,6> edges{{
        {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
    for (Index ci=0;ci<mesh.cellCount();++ci) {
        if (mesh.cell(ci).type!=mesh::CellType::Tet4) continue;
        const Index* conn=mesh.cellNodes(ci);
        std::array<Vec3,4> p{};
        for (std::size_t i=0;i<4;++i) p[i]=mesh.node(conn[i]);
        const double volume=mesh::tetVolume(p[0],p[1],p[2],p[3]);
        const double radiusEdge=tetRadiusEdge(p);
        const double radiusRatio=mesh::tetQuality(p[0],p[1],p[2],p[3]);
        if (!(volume>0.0)||!std::isfinite(volume)||
            !std::isfinite(radiusEdge)||!std::isfinite(radiusRatio)||
            radiusRatio<0.0||radiusRatio>1.0) return false;
        quality.minimumSignedVolume=(std::min)(quality.minimumSignedVolume,volume);
        quality.maximumRadiusEdge=(std::max)(quality.maximumRadiusEdge,radiusEdge);
        radiusEdges.push_back(radiusEdge);
        radiusRatios.push_back(radiusRatio);
        const std::size_t bin=(std::min<std::size_t>)(
            quality.radiusEdgeHistogram.size()-1,
            static_cast<std::size_t>((std::min)(radiusEdge,9.0)));
        ++quality.radiusEdgeHistogram[bin];
        const std::size_t ratioBin=(std::min<std::size_t>)(
            quality.radiusRatioHistogram.size()-1,
            static_cast<std::size_t>(radiusRatio*10.0));
        ++quality.radiusRatioHistogram[ratioBin];
        if (policy.sliverRadiusRatio>0.0&&
            radiusRatio<policy.sliverRadiusRatio)
            ++quality.sliversRemaining;
        for (const auto& edge:edges) {
            const std::size_t a=static_cast<std::size_t>(edge[0]);
            const std::size_t b=static_cast<std::size_t>(edge[1]);
            const Vec3 mid{std::midpoint(p[a].x,p[b].x),
                           std::midpoint(p[a].y,p[b].y),
                           std::midpoint(p[a].z,p[b].z)};
            const double h=(std::min)({
                checkedSize(spec,p[a],ledger,metrics,false),
                checkedSize(spec,p[b],ledger,metrics,false),
                checkedSize(spec,mid,ledger,metrics,false)});
            const double limit=policy.maxEdgeToTarget*h;
            const double severity=(p[b]-p[a]).norm2()/(limit*limit);
            if (!std::isfinite(severity)) return false;
            quality.worstNormalizedSize=(std::max)(
                quality.worstNormalizedSize,severity);
            if (severity>1.0) ++quality.sizeViolations;
        }
    }
    if (radiusEdges.empty()) quality.minimumSignedVolume=0.0;
    std::sort(radiusEdges.begin(),radiusEdges.end());
    const auto percentile=[&](double p) {
        if (radiusEdges.empty()) return 0.0;
        const std::size_t index=static_cast<std::size_t>(
            std::ceil(p*static_cast<double>(radiusEdges.size())))-1;
        return radiusEdges[(std::min)(index,radiusEdges.size()-1)];
    };
    quality.p95RadiusEdge=percentile(0.95);
    quality.p99RadiusEdge=percentile(0.99);
    std::sort(radiusRatios.begin(),radiusRatios.end());
    const auto lowPercentile=[&](double p) {
        if (radiusRatios.empty()) return 0.0;
        const std::size_t index=static_cast<std::size_t>(
            std::ceil(p*static_cast<double>(radiusRatios.size())))-1;
        return radiusRatios[(std::min)(index,radiusRatios.size()-1)];
    };
    quality.minimumRadiusRatio=radiusRatios.empty()?0.0:radiusRatios.front();
    quality.p1RadiusRatio=lowPercentile(0.01);
    quality.p5RadiusRatio=lowPercentile(0.05);
    return true;
}

} // namespace

#if defined(PAACS_TOPO_TESTING)
namespace t3a_detail {
bool certifiedProtectedEdgeBoundsForTesting(
    const std::array<Vec3,2>& edge,Vec3& lo,Vec3& hi) noexcept {
    const auto bounds=sphereBounds(edge);
    if (!bounds) return false;
    lo=bounds->lo;
    hi=bounds->hi;
    return true;
}
bool certifiedProtectedFacetBoundsForTesting(
    const std::array<Vec3,3>& facet,Vec3& lo,Vec3& hi) noexcept {
    const auto bounds=sphereBounds(facet);
    if (!bounds) return false;
    lo=bounds->lo;
    hi=bounds->hi;
    return true;
}
bool evaluateTetViolationForTesting(
    const std::array<Vec3,4>& points,const std::array<Index,4>& globalNodes,
    int regionTag,const sizing::MeshSpec& spec,
    const TetRefinementPolicy& policy,TetViolationRow& row) noexcept {
    try {
        std::array<topo::NodeHandle,4> nodes{};
        std::vector<Index> existing(4,kInvalidIndex);
        for (std::size_t i=0;i<4;++i) {
            nodes[i]=topo::NodeHandle{static_cast<std::int32_t>(i),1,1};
            existing[i]=globalNodes[i];
        }
        const std::map<topo::NodeHandle,std::uint64_t> inserted;
        SizingLedger ledger;
        TetRefinementMetrics metrics;
        const auto result=evaluateViolation(
            nodes,points,regionTag,existing,inserted,spec,policy,ledger,metrics,0);
        if (!result) return false;
        row=*result;
        return true;
    } catch (...) {
        return false;
    }
}
bool exactLocationForTesting(const topo::TopologyWorkspace& workspace,
                             Vec3 point,topo::CellHandle& seed) noexcept {
    try {
        const auto active=workspace.cells().activeBySequence();
        if (active.empty()) return false;
        const auto first=workspace.cells().get(active.front());
        if (first.type!=topo::ElementType::Tet4) return false;
        const auto built=bnd::buildTetRecoveryContext(
            workspace,first.cls);
        if (!built.built()) return false;
        CoordinateMap coordinates;
        for (const auto node:workspace.nodes().activeBySequence())
            coordinates.emplace(node,workspace.nodes().get(node).coord);
        RefinementLocationIndex overlay;
        overlay.prepare(active.size()*4+1,
            workspace.cells().stats().materializedSlots+active.size()*4+1);
        TetRefinementMetrics metrics;
        const auto located=locateExact(
            point,built.context->index,overlay,workspace,coordinates,metrics);
        if (!located) return false;
        seed=located->seed;
        return true;
    } catch (...) {
        return false;
    }
}
bool validateRefinementCavityForTesting(
    const topo::TopologyWorkspace& workspace,
    std::span<const topo::CellHandle> removedSpan,Vec3 stagedPoint) noexcept {
    try {
        const auto active=workspace.cells().activeBySequence();
        if (active.empty()) return false;
        const auto first=workspace.cells().get(active.front());
        if (first.type!=topo::ElementType::Tet4) return false;
        CoordinateMap coordinates;
        for (const auto node:workspace.nodes().activeBySequence())
            coordinates.emplace(node,workspace.nodes().get(node).coord);
        CellMap cells;
        EdgeStars edgeStars;
        FaceOwners faceOwners;
        NodeStars nodeStars;
        for (const auto cell:active) {
            const auto data=workspace.cells().get(cell);
            if (data.type!=topo::ElementType::Tet4||
                !sameClass(data.cls,first.cls)) return false;
            const std::array<topo::NodeHandle,4> conn{
                data.nodes[0],data.nodes[1],data.nodes[2],data.nodes[3]};
            cells.emplace(cell,conn);
            addCellToStars(edgeStars,cell,conn);
            addCellToFaces(faceOwners,cell,conn);
            addCellToNodeStars(nodeStars,cell,conn);
        }
        std::vector<topo::CellHandle> removed(
            removedSpan.begin(),removedSpan.end());
        const auto boundary=cavityBoundary(removed,cells);
        if (!boundary) return false;
        const auto built=bnd::buildTetRecoveryContext(workspace,first.cls);
        if (!built.built()) return false;
        RefinementLocationIndex overlay;
        overlay.prepare(active.size()*4+1,
            workspace.cells().stats().materializedSlots+active.size()*4+1);
        TetRefinementMetrics metrics;
        return validateRefinementCandidate(
            removed,*boundary,stagedPoint,cells,edgeStars,nodeStars,coordinates,
            built.context->index,overlay,workspace,first.cls,metrics);
    } catch (...) {
        return false;
    }
}
bool dynamicLocationScalingForTesting(
    std::size_t insertionCount,std::uint64_t& updateVisits,
    std::uint64_t& queryVisits) noexcept {
    try {
        if (insertionCount==0||insertionCount>
                static_cast<std::size_t>((std::numeric_limits<std::int32_t>::max)()))
            return false;
        struct OracleRow {
            topo::CellHandle cell{};
            bnd::TetAabb bounds{};
            bool active=true;
        };
        std::vector<OracleRow> oracle;
        oracle.reserve(insertionCount);
        const auto exhaustive=[&](const bnd::TetAabb& query) {
            std::vector<topo::CellHandle> result;
            for (const auto& row:oracle)
                if (row.active&&overlaps(row.bounds,query))
                    result.push_back(row.cell);
            std::sort(result.begin(),result.end());
            return result;
        };
        const auto matches=[&](const RefinementLocationIndex& candidate,
                               const bnd::TetAabb& query) {
            return candidate.queryOverlayForTesting(query,queryVisits)==
                   exhaustive(query);
        };
        RefinementLocationIndex index;
        index.prepare(insertionCount,insertionCount);
        updateVisits=0;
        queryVisits=0;
        for (std::size_t i=0;i<insertionCount;++i) {
            const auto slot=(i*73)%insertionCount;
            const auto cell=topo::CellHandle{
                static_cast<std::int32_t>(slot),1,1};
            const double x=static_cast<double>(i)*4.0;
            const RefinementLocationIndex::CreatedRow created{
                cell,{{x,0.0,0.0},{x+1.0,1.0,1.0}}};
            const std::array<topo::CellHandle,0> removed{};
            const std::array<RefinementLocationIndex::CreatedRow,1> added{{
                created}};
            if (!index.canApply(removed,added)) return false;
            index.applyPrepared(removed,added,updateVisits);
            oracle.push_back({cell,created.bounds,true});
        }

        // Every old leaf must remain reachable after all later insertions and
        // AVL rotations, not merely the most recently inserted leaf.
        for (const auto& row:oracle) {
            const Vec3 point{(row.bounds.lo.x+row.bounds.hi.x)*0.5,
                             0.5,0.5};
            if (!matches(index,{point,point})) return false;
        }
        const bnd::TetAabb broad{{-1.0,-1.0,-1.0},
            {static_cast<double>(insertionCount)*4.0+1.0,2.0,2.0}};
        if (!matches(index,broad)) return false;
        for (std::size_t i=0;i<insertionCount;i+=17) {
            const double x=static_cast<double>(i)*4.0;
            const bnd::TetAabb overlap{{x-2.0,-0.25,-0.25},
                                       {x+10.0,1.25,1.25}};
            if (!matches(index,overlap)) return false;
        }

        // Remove a deterministic subset, then reuse those exact slots with a
        // new generation and different bounds before repeating the oracle.
        for (std::size_t i=0;i<insertionCount;i+=3) {
            const auto slot=(i*73)%insertionCount;
            const auto cell=topo::CellHandle{
                static_cast<std::int32_t>(slot),1,1};
            const std::array<topo::CellHandle,1> removed{{cell}};
            const std::array<RefinementLocationIndex::CreatedRow,0> added{};
            if (!index.canApply(removed,added)) return false;
            index.applyPrepared(removed,added,updateVisits);
            oracle[i].active=false;
        }
        for (std::size_t i=0;i<insertionCount;i+=3) {
            const auto slot=(i*73)%insertionCount;
            const auto cell=topo::CellHandle{
                static_cast<std::int32_t>(slot),2,1};
            const double x=-static_cast<double>(i)*5.0-4.0;
            const RefinementLocationIndex::CreatedRow created{
                cell,{{x,-2.0,-1.0},{x+1.5,-0.5,0.5}}};
            const std::array<topo::CellHandle,0> removed{};
            const std::array<RefinementLocationIndex::CreatedRow,1> added{{
                created}};
            if (!index.canApply(removed,added)) return false;
            index.applyPrepared(removed,added,updateVisits);
            oracle[i]={cell,created.bounds,true};
        }
        for (const auto& row:oracle) {
            const Vec3 point{(row.bounds.lo.x+row.bounds.hi.x)*0.5,
                             (row.bounds.lo.y+row.bounds.hi.y)*0.5,
                             (row.bounds.lo.z+row.bounds.hi.z)*0.5};
            if (!matches(index,{point,point})) return false;
        }
        if (!matches(index,broad)) return false;
        const bnd::TetAabb reusedBroad{
            {-static_cast<double>(insertionCount)*5.0-5.0,-3.0,-2.0},
            {2.0,2.0,2.0}};
        if (!matches(index,reusedBroad)) return false;
        for (std::size_t i=0;i<insertionCount;i+=19) {
            const double x=static_cast<double>(i)*4.0;
            if (!matches(index,{{x-5.0,-2.5,-1.5},
                                {x+9.0,1.5,1.5}})) return false;
        }

        // The query contract is canonical even when the same final leaves are
        // inserted in reverse order and the tree takes different rotations.
        RefinementLocationIndex reversed;
        reversed.prepare(insertionCount,insertionCount);
        for (std::size_t i=oracle.size();i>0;--i) {
            const auto& row=oracle[i-1];
            const std::array<topo::CellHandle,0> removed{};
            const std::array<RefinementLocationIndex::CreatedRow,1> added{{
                {row.cell,row.bounds}}};
            if (!reversed.canApply(removed,added)) return false;
            reversed.applyPrepared(removed,added,updateVisits);
        }
        for (const auto& row:oracle) {
            const Vec3 point{(row.bounds.lo.x+row.bounds.hi.x)*0.5,
                             (row.bounds.lo.y+row.bounds.hi.y)*0.5,
                             (row.bounds.lo.z+row.bounds.hi.z)*0.5};
            const auto forward=index.queryOverlayForTesting(pointAabb(point),
                                                             queryVisits);
            const auto reverse=reversed.queryOverlayForTesting(pointAabb(point),
                                                                queryVisits);
            if (forward!=reverse||forward!=exhaustive(pointAabb(point)))
                return false;
        }
        return true;
    } catch (...) {
        return false;
    }
}
bool coordinatorDeltaScalingForTesting(
    std::size_t currentCellCount,std::uint64_t& affectedRows,
    std::uint64_t& workUnits,std::uint64_t& accountedRemovedBytes,
    std::uint64_t& accountedAddedBytes) noexcept {
    try {
        if (currentCellCount==0||currentCellCount>
                static_cast<std::size_t>(
                    (std::numeric_limits<std::int32_t>::max)()/4-4))
            return false;
        CoordinateMap coordinates;
        InsertedMap inserted;
        CoordinateIndex coordinateIndex;
        CellMap cells;
        EdgeStars edgeStars;
        FaceOwners faceOwners;
        NodeStars nodeStars;
        for (std::size_t i=0;i<currentCellCount;++i) {
            std::array<topo::NodeHandle,4> conn{};
            for (std::size_t j=0;j<4;++j) {
                conn[j]=topo::NodeHandle{
                    static_cast<std::int32_t>(i*4+j),1,1};
                coordinates.emplace(conn[j],Vec3{
                    static_cast<double>(i)*8.0+
                        (j==1?1.0:0.0),j==2?1.0:0.0,j==3?1.0:0.0});
                coordinateIndex.emplace(
                    coordinateBits(coordinates.at(conn[j])),conn[j]);
            }
            const auto cell=topo::CellHandle{
                static_cast<std::int32_t>(i),1,1};
            cells.emplace(cell,conn);
            addCellToStars(edgeStars,cell,conn);
            addCellToFaces(faceOwners,cell,conn);
            addCellToNodeStars(nodeStars,cell,conn);
        }
        const auto removedCell=topo::CellHandle{0,1,1};
        const std::array<topo::CellHandle,1> removed{{removedCell}};
        const auto original=cells.at(removedCell);
        const auto stagedNode=topo::NodeHandle{
            static_cast<std::int32_t>(currentCellCount*4),1,1};
        const Vec3 stagedPoint{0.2,0.2,0.2};
        CreatedCellRows created;
        created.reserve(4);
        for (std::size_t i=0;i<kLocalTetFaces.size();++i) {
            const auto& face=kLocalTetFaces[i];
            created.push_back({
                topo::CellHandle{
                    static_cast<std::int32_t>(currentCellCount+i),1,1},
                {original[face[0]],original[face[1]],original[face[2]],
                 stagedNode}});
        }
        ViolationHeap queue;
        queue.reserve(currentCellCount+4,currentCellCount+4);
        CurrentBlockedHeap blocked;
        blocked.reserve(currentCellCount+4,currentCellCount+4);
        T3AReport report;
        RefinementLocationIndex location;
        location.prepare(currentCellCount+4,currentCellCount+4);
        SizingLedger ledger;
        ledger.prepare(currentCellCount+4);
        const std::vector<QueuedViolation> violations;
        PreparedTetRefinementCoordinatorDelta prepared;
        TetRefinementMetrics metrics;
        const auto before=completeContextAccountedBytes(
            coordinates,inserted,coordinateIndex,cells,edgeStars,faceOwners,
            nodeStars,queue,blocked,ledger,location);
        prepareCoordinatorDelta(
            prepared,stagedNode,stagedPoint,0,removed,created,violations,
            coordinates,cells,edgeStars,faceOwners,nodeStars,queue,blocked,
            report,location,metrics);
        affectedRows=metrics.coordinatorRowsPrepared;
        workUnits=metrics.workUnitsConsumed;
        accountedRemovedBytes=prepared.contextBytesRemoved;
        accountedAddedBytes=prepared.contextBytesAdded;
        const bool localDelta=!prepared.cellCreates.empty()&&
            prepared.edgeRows.size()<edgeStars.size()&&
            accountedRemovedBytes>0&&accountedAddedBytes>0;
        const auto incrementallyUpdated=applyAccountedContextDelta(
            before,prepared);
        applyCoordinatorDelta(
            prepared,coordinates,inserted,coordinateIndex,cells,edgeStars,
            faceOwners,nodeStars,queue,blocked,report,location,metrics);
        const auto exhaustivelyRecomputed=completeContextAccountedBytes(
            coordinates,inserted,coordinateIndex,cells,edgeStars,faceOwners,
            nodeStars,queue,blocked,ledger,location);
        return localDelta&&incrementallyUpdated==exhaustivelyRecomputed;
    } catch (...) {
        return false;
    }
}

bool blockedSchedulerReuseForTesting() noexcept {
    try {
        CurrentBlockedHeap blocked;
        blocked.reserve(3,2);
        const topo::CellHandle oldCell{0,1,1};
        TetViolationRow oldRow;
        oldRow.key.tet.vertices[0]={RefinementVertexKind::ExistingGlobal,10};
        oldRow.generation=7;
        const CurrentBlockedEntry oldEntry{oldCell,oldRow,0};
        if (!blocked.canInsert(oldEntry)) return false;
        blocked.insertPrepared(oldEntry);
        if (!blocked.find(oldCell)||blocked.top().cell!=oldCell) return false;
        if (blocked.erasePrepared(oldCell,0)!=0||blocked.find(oldCell))
            return false;

        const topo::CellHandle reusedCell{0,2,1};
        TetViolationRow reusedRow;
        reusedRow.key.tet.vertices[0]=
            {RefinementVertexKind::ExistingGlobal,20};
        reusedRow.generation=8;
        const CurrentBlockedEntry reusedEntry{reusedCell,reusedRow,1};
        if (!blocked.canInsert(reusedEntry)) return false;
        blocked.insertPrepared(reusedEntry);
        if (blocked.find(oldCell)||!blocked.find(reusedCell)) return false;

        const topo::CellHandle otherCell{1,1,1};
        TetViolationRow otherRow;
        otherRow.key.tet.vertices[0]=
            {RefinementVertexKind::ExistingGlobal,5};
        otherRow.generation=9;
        const CurrentBlockedEntry otherEntry{otherCell,otherRow,2};
        if (!blocked.canInsert(otherEntry)) return false;
        blocked.insertPrepared(otherEntry);
        return blocked.top().cell==otherCell&&blocked.find(reusedCell)&&
               blocked.find(otherCell)&&blocked.size()==2;
    } catch (...) {
        return false;
    }
}
} // namespace t3a_detail
#endif

T3AOutcome meshModelRefinedTet4(
    const geo::Model& model,const sizing::MeshSpec& spec,
    const T3ARequest& request,mesh::Mesh& destination,
    const VolumeBoundaryOwnershipProvider& ownershipProvider) noexcept {
    T3AOutcome out;
    try {
        out.report=std::make_unique<T3AReport>();
        out.sourceReceipt=std::make_unique<BoundaryConformingV3Receipt>();
    } catch (const std::bad_alloc&) {
        out.failure={T3AReason::AllocationFailed,T3APhase::ReportEnvelope};
        return out;
    }
    if (!validPolicy(spec,request.refinement)) {
        out.status=T3AStatus::Rejected;
        out.failure={spec.order!=1?T3AReason::UnsupportedElementOrder:
                                  T3AReason::InvalidRefinementPolicy,
                     T3APhase::Preflight};
        out.report->envelopePrepared=true;
        out.report->reportComplete=true;
        out.reportComplete=true;
        return out;
    }

    // Mesh owns standard-library index containers whose default construction
    // can allocate on MSVC.  Construct the private seed inside the structured
    // failure boundary; this route is noexcept and must never let that
    // implementation-dependent allocation reach std::terminate.
    std::unique_ptr<mesh::Mesh> seed;
    try {
        seed=std::make_unique<mesh::Mesh>();
    } catch (const std::bad_alloc&) {
        out.failure={T3AReason::AllocationFailed,T3APhase::SeedTransfer};
        return out;
    }
#if defined(PAACS_TOPO_TESTING)
    auto source=[&]() {
        TestingWallScope timer(t3a_detail::T3AWallPhase::Source);
        return meshModelBoundaryConforming(
            model,spec,request.conformingRequest,*seed,ownershipProvider);
    }();
#else
    auto source=meshModelBoundaryConforming(model,spec,request.conformingRequest,
                                            *seed,ownershipProvider);
#endif
    copyReceipt(*out.sourceReceipt,source);
    if (!source.published() || !source.retainedResult) {
        out.status=sourceStatus(source.status);
        out.failure={sourceReason(source.status),T3APhase::BoundaryConformingV3};
        out.report->envelopePrepared=true;
        out.report->reportComplete=true;
        out.reportComplete=true;
        return out;
    }
#if defined(PAACS_TOPO_TESTING)
    auto& refinementSizingPhase=
        t3a_detail::refinementSizingPhaseForTesting();
    refinementSizingPhase=true;
    struct RefinementSizingPhaseReset {
        bool& value;
        ~RefinementSizingPhaseReset() { value=false; }
    } refinementSizingPhaseReset{refinementSizingPhase};
#endif

    T3APhase activePhase=T3APhase::SeedTransfer;
    try {
#if defined(PAACS_TOPO_TESTING)
        std::optional<TestingWallScope> setupTimer;
        setupTimer.emplace(t3a_detail::T3AWallPhase::Setup);
#endif
        // Keep the private owner local until a successful or intentionally
        // incomplete return.  A TetRefinementTransaction may still be live
        // on any inner failure path; stack unwinding must destroy/abandon the
        // transaction before the workspace owner is released.
        auto workingResult=std::make_unique<T3AResult>();
        auto& result=*workingResult;
        result.authoritative=std::move(source.retainedResult->finalAuthoritative);
        result.constraintState.swap(source.retainedResult->constraintState);
        result.existingGlobalByAuthoritativeNodeSlot=
            std::move(source.retainedResult->globalByAuthoritativeNode);
        result.boundaryAudit=source.retainedResult->boundaryAudit;
        out.sourceReceipt->authoritativeResultConsumed=true;
        if (!out.sourceReceipt->surfaceOutcome||
            !buildImmutableBoundary(*out.sourceReceipt->surfaceOutcome,
                                    result.constraintState,*seed,
                                    result.immutableBoundary))
            throw std::logic_error("v3 result lacks immutable boundary map");
        if (!result.authoritative || !result.authoritative->workspace)
            throw std::logic_error("v3 result lacks authoritative workspace");
        auto& ws=*result.authoritative->workspace;
        if (!bindConstraintStateToAuthoritative(
                result.constraintState,result.immutableBoundary,
                result.existingGlobalByAuthoritativeNodeSlot,ws)) {
            out.status=T3AStatus::Failed;
            out.failure={T3AReason::InvalidConformingResult,
                         T3APhase::SeedTransfer};
            out.retainedResult.reset(); return out;
        }

        activePhase=T3APhase::ReportEnvelope;
        std::uint64_t stops=0,attemptCapacity=0,witnessCapacity=0,
                      topologyCapacity=0;
        const std::uint64_t attemptLimit=(std::max)(
            static_cast<std::uint64_t>(request.refinement.maxCandidateAttemptsPerViolation),
            static_cast<std::uint64_t>(request.refinement.maxSliverCandidateAttempts));
        const std::uint64_t featureCount=
            static_cast<std::uint64_t>(result.constraintState.edgeBindings.size())+
            static_cast<std::uint64_t>(result.constraintState.bindings.size());
        if (!checkedAdd(request.refinement.maxInsertions,1,stops) ||
            !checkedMul(stops,attemptLimit,attemptCapacity) ||
            !checkedMul(attemptCapacity,featureCount,witnessCapacity) ||
            !checkedAdd(request.refinement.maxWorkUnits,
                        request.refinement.maxInsertions,topologyCapacity) ||
            attemptCapacity>request.refinement.maxReportRows ||
            witnessCapacity>request.refinement.maxReportWitnesses ||
            topologyCapacity>request.refinement.maxReportTopologyHandles ||
            !fitsSize(attemptCapacity) || !fitsSize(witnessCapacity) ||
            !fitsSize(request.refinement.maxOutputTets) ||
            !fitsSize(request.refinement.maxInsertions) ||
            request.refinement.maxWorkUnits>
                static_cast<std::uint64_t>(
                    ((std::numeric_limits<std::size_t>::max)()-1)/2) ||
            !fitsSize(topologyCapacity)) {
            out.status=T3AStatus::Rejected;
            out.failure={T3AReason::InvalidRefinementPolicy,T3APhase::ReportEnvelope};
            out.retainedResult.reset();
            return out;
        }
        out.report->attemptStorage.resize(static_cast<std::size_t>(attemptCapacity));
        out.report->protectedWitnessStorage.resize(
            static_cast<std::size_t>(witnessCapacity));
        out.report->blockedStorage.resize(
            static_cast<std::size_t>(attemptCapacity));
        out.report->unattemptedStorage.resize(
            static_cast<std::size_t>(request.refinement.maxOutputTets));
        out.report->committedPrefixStorage.resize(
            static_cast<std::size_t>(request.refinement.maxInsertions));
        out.report->topologyHandleStorage.resize(
            static_cast<std::size_t>(topologyCapacity));
        out.report->envelopePrepared=true;
        out.report->reportComplete=true;
        out.reportComplete=true;

        activePhase=T3APhase::ContextBuild;
        auto& metrics=result.metrics;
        metrics.authoritativeAccountedBytes=workspaceAccountedBytes(ws);
        InsertedMap inserted;
        CoordinateMap coordinates;
        for (const auto node:ws.nodes().activeBySequence())
            coordinates.emplace(node,ws.nodes().get(node).coord);
        CellMap cells;
        EdgeStars edgeStars;
        FaceOwners faceOwners;
        NodeStars nodeStars;
        mesh::Classification region{};
        bool haveRegion=false;
        for (const auto cell:ws.cells().activeBySequence()) {
            const auto data=ws.cells().get(cell);
            if (data.type!=topo::ElementType::Tet4)
                throw std::logic_error("non-Tet4 authoritative cell");
            if (!haveRegion) { region=data.cls; haveRegion=true; }
            else if (!sameClass(region,data.cls)) {
                out.status=T3AStatus::Unsupported;
                out.failure={T3AReason::MultipleRegionsOrInterfaces,
                             T3APhase::InitialAudit};
                out.retainedResult.reset(); return out;
            }
            const std::array<topo::NodeHandle,4> nodes{
                data.nodes[0],data.nodes[1],data.nodes[2],data.nodes[3]};
            cells.emplace(cell,nodes);
            addCellToStars(edgeStars,cell,nodes);
            addCellToFaces(faceOwners,cell,nodes);
            addCellToNodeStars(nodeStars,cell,nodes);
        }
        if (!haveRegion || region.dim!=3 || region.tag<0)
            throw std::logic_error("invalid authoritative region");
        auto baseLocationBuild=bnd::buildTetRecoveryContext(ws,region);
        if (!baseLocationBuild.built()) {
            out.status=T3AStatus::Failed;
            out.failure={baseLocationBuild.reason==
                    bnd::TetRecoveryContextReason::AllocationFailed
                    ? T3AReason::AllocationFailed:T3AReason::InvalidConformingResult,
                T3APhase::ContextBuild};
            return out;
        }
        auto baseLocation=std::move(baseLocationBuild.context);

        ProtectedFeatureIndex protectedIndex;
        const auto protectedBuild=buildProtectedIndex(
            result.constraintState,ws,protectedIndex,metrics);
        if (protectedBuild!=ProtectedIndexBuild::Built) {
            out.status=protectedBuild==ProtectedIndexBuild::Range
                ? T3AStatus::Unsupported:T3AStatus::Failed;
            out.failure={protectedBuild==ProtectedIndexBuild::Range
                    ? T3AReason::ProtectedIndexRangeUnrepresentable
                    : T3AReason::InvalidProtectedFeature,
                T3APhase::ContextBuild};
            out.retainedResult.reset(); return out;
        }
        const std::uint64_t initialNodes=coordinates.size();
        if (request.refinement.maxOutputNodes<
                initialNodes+request.refinement.maxInsertions ||
            request.refinement.maxOutputTets<cells.size()) {
            out.status=T3AStatus::Rejected;
            out.failure={T3AReason::InvalidRefinementPolicy,T3APhase::Preflight};
            out.retainedResult.reset(); return out;
        }

        activePhase=T3APhase::QueueBuild;
        SizingLedger ledger;
        ledger.prepare(static_cast<std::size_t>(
            request.refinement.maxWorkUnits));
        ViolationHeap queue;
        std::uint64_t queueSlotCapacity=0;
        if (!checkedAdd(
                static_cast<std::uint64_t>(ws.cells().stats().materializedSlots),
                request.refinement.maxOutputTets,queueSlotCapacity)||
            !fitsSize(queueSlotCapacity)) {
            out.status=T3AStatus::Rejected;
            out.failure={T3AReason::InvalidRefinementPolicy,
                         T3APhase::QueueBuild};
            out.retainedResult.reset(); return out;
        }
        queue.reserve(static_cast<std::size_t>(request.refinement.maxOutputTets),
                      static_cast<std::size_t>(queueSlotCapacity));
        CurrentBlockedHeap blocked;
        blocked.reserve(static_cast<std::size_t>(attemptCapacity),
                        static_cast<std::size_t>(queueSlotCapacity));
        RefinementLocationIndex locationIndex;
        locationIndex.prepare(
            static_cast<std::size_t>(request.refinement.maxOutputTets),
            static_cast<std::size_t>(queueSlotCapacity));
        WorkBudgetScope workBudget(request.refinement.maxWorkUnits);
        const auto enqueue=[&](QueuedViolation value) {
            recordConstructionWork(metrics);
            queue.push(std::move(value));
        };
        const auto completeQueueReport=[&]() noexcept {
            if (!fillUnattempted(*out.report,queue)) {
                out.report->reportComplete=false;
                out.reportComplete=false;
            }
        };
        bool incomplete=false;
        T3AReason incompleteReason=T3AReason::None;
        T3APhase incompletePhase=T3APhase::None;
        std::uint64_t generation=0;
        activePhase=T3APhase::SizingEvaluation;
        try {
            for (const auto& [cell,nodes]:cells) {
                ++metrics.initialCellsScanned;
                std::array<Vec3,4> p{};
                for (std::size_t i=0;i<4;++i) p[i]=coordinates.at(nodes[i]);
                if (predicates::orient3d(p[0],p[1],p[2],p[3])!=
                        predicates::Sign::Positive)
                    throw std::logic_error("nonpositive initial Tet");
                const double initialVolume=std::abs(
                    mesh::tetVolume(p[0],p[1],p[2],p[3]));
                if (!std::isfinite(initialVolume))
                    throw std::logic_error("nonfinite initial Tet volume");
                if (initialVolume<request.refinement.minimumCreatedTetVolume)
                    ++result.quality.inheritedBelowCreatedFloor;
                auto violation=evaluateViolation(
                    nodes,p,region.tag,
                    result.existingGlobalByAuthoritativeNodeSlot,
                    inserted,spec,request.refinement,ledger,metrics,
                    generation++);
                if (violation) enqueue({*violation,cell,nodes});
            }
        } catch (const WorkBudgetStop&) {
            incomplete=true;
            incompleteReason=T3AReason::WorkBudgetReached;
        }
        metrics.queuePushes=queue.size();
        metrics.queueHighWater=queue.size();
        activePhase=T3APhase::QueueBuild;

        CoordinateIndex coordinateIndex;
        for (const auto& [node,coord]:coordinates)
            if (!coordinateIndex.emplace(coordinateBits(coord),node).second)
                throw std::logic_error("duplicate authoritative coordinate");
        ProtectedQuery protectedHit;
        protectedHit.witnesses.reserve(protectedIndex.rows.size());
        std::vector<std::size_t> protectedStack;
        protectedStack.reserve(protectedIndex.nodes.size());
        std::uint64_t protectedAccounted=0;
        addAccountedBytes(protectedAccounted,
            vectorAccountedBytes(protectedIndex.rows));
        addAccountedBytes(protectedAccounted,
            vectorAccountedBytes(protectedIndex.nodes));
        addAccountedBytes(protectedAccounted,
            vectorAccountedBytes(protectedHit.witnesses));
        addAccountedBytes(protectedAccounted,
            vectorAccountedBytes(protectedStack));
        const std::uint64_t reportAccounted=reportAccountedBytes(*out.report);

        // Accounted transient construction storage.  The convention covers
        // every allocated ordered-map node, every star vector capacity, both
        // queue buffers, the fixed sizing ledger, the dynamic location tree
        // and the complete flat report envelope.  It is a conservative
        // implementation-accounted bound, not a sampled process-RSS claim.
        std::uint64_t currentContextAccounted=completeContextAccountedBytes(
            coordinates,inserted,coordinateIndex,cells,edgeStars,faceOwners,
            nodeStars,queue,blocked,ledger,locationIndex);
        std::uint64_t initialTransient=0;
        addAccountedBytes(initialTransient,protectedAccounted);
        addAccountedBytes(initialTransient,reportAccounted);
        addAccountedBytes(initialTransient,currentContextAccounted);
        metrics.transientAccountedHighWaterBytes=initialTransient;

#if defined(PAACS_TOPO_TESTING)
        setupTimer.reset();
#endif
        while (!incomplete&&!queue.empty()) {
            if (metrics.workUnitsConsumed>=request.refinement.maxWorkUnits) {
                incomplete=true;
                incompleteReason=T3AReason::WorkBudgetReached;
                break;
            }
            recordConstructionWork(metrics);
            const QueuedViolation selectedRoot=queue.top();
            const auto liveSelected=cells.find(selectedRoot.cell);
            if (liveSelected==cells.end()||
                liveSelected->second!=selectedRoot.nodes) {
                (void)queue.pop();
                ++metrics.queuePops;
                ++metrics.staleQueueRows;
                continue;
            }
            if (metrics.committedInsertions>=request.refinement.maxInsertions) {
                incomplete=true; incompleteReason=T3AReason::InsertionBudgetReached;
                break;
            }
            if (metrics.workUnitsConsumed>=request.refinement.maxWorkUnits) {
                incomplete=true; incompleteReason=T3AReason::WorkBudgetReached;
                break;
            }
            if (out.report->attemptsUsed>=out.report->attemptStorage.size()) {
                incomplete=true; incompleteReason=T3AReason::ReportBudgetReached;
                break;
            }
            auto& attempt=out.report->attemptStorage[
                static_cast<std::size_t>(out.report->attemptsUsed++)];
            QueuedViolation selected=queue.pop();
            ++metrics.queuePops;
            T3AAttemptRow* activeAttempt=&attempt;
            try {
            if (selected.row.key.kind==TetViolationKind::Sliver)
                ++result.quality.sliversAttempted;
            attempt.attemptOrdinal=metrics.candidateAttempts++;
            attempt.violation=selected.row.key;
            attempt.violationGeneration=selected.row.generation;
            attempt.workBefore=metrics.workUnitsConsumed;
            activePhase=T3APhase::CandidateGeneration;
            recordConstructionWork(metrics);
            std::array<Vec3,4> selectedPoints{};
            for (std::size_t i=0;i<4;++i)
                selectedPoints[i]=coordinates.at(selected.nodes[i]);
            const auto descriptors=candidateDescriptors(
                selected.row,selected.nodes,
                result.existingGlobalByAuthoritativeNodeSlot,inserted);
            const std::uint32_t requestedCap=selected.row.key.kind==
                    TetViolationKind::Sliver
                ? request.refinement.maxSliverCandidateAttempts
                : request.refinement.maxCandidateAttemptsPerViolation;
            const std::uint8_t descriptorCap=static_cast<std::uint8_t>(
                (std::min<std::uint32_t>)(requestedCap,5));
            if (selected.nextDescriptor>=descriptorCap) {
                attempt.disposition=T3AAttemptDisposition::TerminalFailure;
                attempt.reason=T3AReason::NoAdmissibleCandidate;
                incomplete=true; incompleteReason=T3AReason::NoAdmissibleCandidate;
                out.report->hasTerminalViolation=true;
                out.report->terminalViolation=selected.row.key;
                break;
            }
            const auto& descriptor=descriptors[selected.nextDescriptor];
            attempt.candidateKind=descriptor.kind;
            attempt.hasFace=descriptor.hasFace;
            attempt.face=descriptor.face;
            std::optional<Vec3> candidate;
            std::optional<T3ASizingSample> proposalSizing;
            try {
                candidate=descriptor.kind==T3ACandidateKind::TetCircumcenter
                    ? tetCircumcenter(selectedPoints)
                    : faceOffcenter(descriptor,selectedPoints,spec,ledger,
                                    metrics,&proposalSizing);
            } catch (const SizingFailure& failure) {
                attempt.disposition=T3AAttemptDisposition::TerminalFailure;
                attempt.reason=failure.reason;
                out.status=T3AStatus::Failed;
                out.failure={failure.reason,T3APhase::SizingEvaluation,
                             true,selected.row.key};
                completeQueueReport();
                out.retainedResult.reset(); return out;
            } catch (const std::bad_alloc&) {
                attempt.disposition=T3AAttemptDisposition::TerminalFailure;
                attempt.reason=T3AReason::AllocationFailed;
                out.status=T3AStatus::Failed;
                out.failure={T3AReason::AllocationFailed,
                             T3APhase::SizingEvaluation,true,selected.row.key};
                completeQueueReport();
                out.retainedResult.reset(); return out;
            }
            if (!candidate) {
                attempt.disposition=T3AAttemptDisposition::FormulaUnavailable;
                attempt.reason=T3AReason::NoAdmissibleCandidate;
                ++selected.nextDescriptor;
                if (selected.nextDescriptor<descriptorCap) {
                    enqueue(selected);
                    continue;
                }
                incomplete=true; incompleteReason=T3AReason::NoAdmissibleCandidate;
                out.report->hasTerminalViolation=true;
                out.report->terminalViolation=selected.row.key;
                break;
            }
            attempt.proposedCoordinate=*candidate;
            if (coordinateIndex.count(coordinateBits(*candidate))) {
                attempt.disposition=T3AAttemptDisposition::DuplicateProposal;
                attempt.reason=T3AReason::Nonprogress;
                ++selected.nextDescriptor;
                if (selected.nextDescriptor<descriptorCap) {
                    enqueue(selected);
                    continue;
                }
                incomplete=true; incompleteReason=T3AReason::Nonprogress;
                out.report->hasTerminalViolation=true;
                out.report->terminalViolation=selected.row.key;
                break;
            }
            activePhase=T3APhase::Encroachment;
            try {
#if defined(PAACS_TOPO_TESTING)
                TestingWallScope timer(
                    t3a_detail::T3AWallPhase::ProtectedQuery);
#endif
                protectedEncroachment(*candidate,protectedIndex,ws,metrics,
                                      protectedHit,protectedStack);
            } catch (const std::bad_alloc&) {
                attempt.disposition=T3AAttemptDisposition::TerminalFailure;
                attempt.reason=T3AReason::AllocationFailed;
                out.status=T3AStatus::Failed;
                out.failure={T3AReason::AllocationFailed,T3APhase::Encroachment,
                             true,selected.row.key};
                completeQueueReport();
                out.retainedResult.reset(); return out;
            }
            if (protectedHit.invalid) {
                attempt.disposition=T3AAttemptDisposition::TerminalFailure;
                attempt.reason=T3AReason::InvalidProtectedFeature;
                out.status=T3AStatus::Failed;
                out.failure={T3AReason::InvalidProtectedFeature,T3APhase::Encroachment,
                             true,selected.row.key};
                completeQueueReport();
                out.retainedResult.reset(); return out;
            }
            if (!protectedHit.witnesses.empty()) {
                const auto witnessCount=static_cast<std::uint64_t>(
                    protectedHit.witnesses.size());
                std::uint64_t witnessEnd=0;
#if defined(PAACS_TOPO_TESTING)
                recordBlockerStateDigest(
                    t3a_detail::T3ABlockerDigestCheckpoint::
                        BlockTransitionBefore,
                    ws,queue,blocked,*out.report,metrics);
#endif
#if defined(PAACS_TOPO_TESTING)
                bool forcedReportBudget=false;
                auto& reportBudgetCountdown=
                    t3a_detail::blockerReportBudgetCountdownForTesting();
                if (reportBudgetCountdown>0&&--reportBudgetCountdown==0)
                    forcedReportBudget=true;
#endif
                if (
#if defined(PAACS_TOPO_TESTING)
                    forcedReportBudget||
#endif
                    !checkedAdd(out.report->protectedWitnessesUsed,
                                witnessCount,witnessEnd)||
                    witnessEnd>out.report->protectedWitnessStorage.size()||
                    out.report->blockedUsed>=
                        out.report->blockedStorage.size()) {
                    attempt.disposition=T3AAttemptDisposition::TerminalFailure;
                    attempt.reason=T3AReason::ReportBudgetReached;
                    attempt.workAfter=metrics.workUnitsConsumed;
                    incomplete=true; incompleteReason=T3AReason::ReportBudgetReached;
                    incompletePhase=T3APhase::Encroachment;
                    out.report->hasTerminalViolation=true;
                    out.report->terminalViolation=selected.row.key;
#if defined(PAACS_TOPO_TESTING)
                    recordBlockerStateDigest(
                        t3a_detail::T3ABlockerDigestCheckpoint::
                            BlockTransitionTerminalAfter,
                        ws,queue,blocked,*out.report,metrics);
#endif
                    break;
                }
                const auto blockedRowIndex=out.report->blockedUsed;
                const CurrentBlockedEntry blockedEntry{
                    selected.cell,selected.row,blockedRowIndex};
                if (!blocked.canInsert(blockedEntry)) {
                    attempt.disposition=T3AAttemptDisposition::TerminalFailure;
                    attempt.reason=T3AReason::InternalError;
                    attempt.workAfter=metrics.workUnitsConsumed;
                    out.status=T3AStatus::Failed;
                    out.failure={T3AReason::InternalError,
                                 T3APhase::CoordinatorPreparation,true,
                                 selected.row.key};
                    completeQueueReport();
#if defined(PAACS_TOPO_TESTING)
                    recordBlockerStateDigest(
                        t3a_detail::T3ABlockerDigestCheckpoint::
                            BlockTransitionTerminalAfter,
                        ws,queue,blocked,*out.report,metrics);
#endif
                    out.retainedResult.reset();
                    return out;
                }
                const auto blockDepth=static_cast<std::uint64_t>(
                    std::bit_width(blocked.size()));
                std::uint64_t blockWork=0;
                if (!checkedAdd(3+blockDepth,witnessCount,blockWork))
                    throw std::length_error("blocker work overflow");
#if defined(PAACS_TOPO_TESTING)
                auto& workBudgetCountdown=
                    t3a_detail::blockerWorkBudgetCountdownForTesting();
                if (workBudgetCountdown>0&&--workBudgetCountdown==0)
                    throw WorkBudgetStop{};
#endif
                recordConstructionWork(metrics,blockWork);
#if defined(PAACS_TOPO_TESTING)
                auto& countdown=
                    t3a_detail::reportFaultCountdownForTesting();
                for (std::uint64_t i=0;i<witnessCount;++i)
                        if (countdown>0&&--countdown==0)
                            throw std::bad_alloc{};
#endif
                const auto witnessOffset=out.report->protectedWitnessesUsed;
                for (std::size_t i=0;i<protectedHit.witnesses.size();++i)
                    out.report->protectedWitnessStorage[
                        static_cast<std::size_t>(witnessOffset)+i]=
                        protectedHit.witnesses[i];
                auto& blockedRow=out.report->blockedStorage[
                    static_cast<std::size_t>(blockedRowIndex)];
                blockedRow.violation=selected.row.key;
                blockedRow.violationGeneration=selected.row.generation;
                blockedRow.attemptOrdinal=attempt.attemptOrdinal;
                blockedRow.witnessOffset=witnessOffset;
                blockedRow.witnessCount=witnessCount;
                blockedRow.invalidatedByInsertionOrdinal=
                    kCurrentBlockedInsertion;
                attempt.witnessOffset=witnessOffset;
                attempt.witnessCount=witnessCount;
                attempt.disposition=T3AAttemptDisposition::Encroached;
                attempt.reason=T3AReason::BoundaryRefinementRequired;
                attempt.workAfter=metrics.workUnitsConsumed;
                blocked.insertPrepared(blockedEntry);
                out.report->protectedWitnessesUsed=witnessEnd;
                ++out.report->blockedUsed;
                ++metrics.blockedRowsRecorded;
                ++metrics.blockedRowsCurrent;
                metrics.blockedRowsCurrentHighWater=(std::max)(
                    metrics.blockedRowsCurrentHighWater,
                    metrics.blockedRowsCurrent);
                continue;
            }

            activePhase=T3APhase::CavityDiscovery;
#if defined(PAACS_TOPO_TESTING)
            const auto located=[&]() {
                TestingWallScope timer(t3a_detail::T3AWallPhase::Location);
                return locateExact(
                    *candidate,baseLocation->index,locationIndex,ws,
                    coordinates,metrics);
            }();
#else
            const auto located=locateExact(
                *candidate,baseLocation->index,locationIndex,ws,
                coordinates,metrics);
#endif
            if (!located) {
                attempt.disposition=T3AAttemptDisposition::OutsideDomain;
                attempt.reason=T3AReason::NoAdmissibleCandidate;
                ++selected.nextDescriptor;
                if (selected.nextDescriptor<descriptorCap) {
                    enqueue(selected);
                    continue;
                }
                incomplete=true; incompleteReason=T3AReason::NoAdmissibleCandidate;
                break;
            }
#if defined(PAACS_TOPO_TESTING)
            const auto discovered=[&]() {
                TestingWallScope timer(t3a_detail::T3AWallPhase::Cavity);
                return discoverConflictCavity(
                    located->seed,*candidate,cells,faceOwners,
                    coordinates,metrics);
            }();
#else
            const auto discovered=discoverConflictCavity(
                located->seed,*candidate,cells,faceOwners,coordinates,metrics);
#endif
            if (!discovered) {
                attempt.disposition=T3AAttemptDisposition::InvalidCavity;
                attempt.reason=T3AReason::InvalidCavity;
                ++selected.nextDescriptor;
                if (selected.nextDescriptor<descriptorCap) {
                    enqueue(selected);
                    continue;
                }
                incomplete=true; incompleteReason=T3AReason::InvalidCavity;
                out.report->hasTerminalViolation=true;
                out.report->terminalViolation=selected.row.key;
                break;
            }
            const auto& removed=*discovered;
            auto boundary=cavityBoundary(removed,cells);
            if (!boundary||!pointStrictlyInsideBoundary(
                    *boundary,*candidate,ws,metrics)) {
                attempt.disposition=T3AAttemptDisposition::InvalidCavity;
                attempt.reason=T3AReason::InvalidCavity;
                ++selected.nextDescriptor;
                if (selected.nextDescriptor<descriptorCap) {
                    enqueue(selected);
                    continue;
                }
                incomplete=true; incompleteReason=T3AReason::InvalidCavity;
                out.report->hasTerminalViolation=true;
                out.report->terminalViolation=selected.row.key;
                break;
            }
            metrics.cavityFaceVisits+=boundary->size();
            metrics.cavityHighWater=(std::max<std::uint64_t>)(
                metrics.cavityHighWater,removed.size());
            if (!selectedSizeEdgeAbsent(
                    selected.row,*boundary,
                    result.existingGlobalByAuthoritativeNodeSlot,inserted)) {
                attempt.disposition=T3AAttemptDisposition::NoProgress;
                attempt.reason=T3AReason::Nonprogress;
                ++selected.nextDescriptor;
                if (selected.nextDescriptor<descriptorCap) {
                    enqueue(selected);
                    continue;
                }
                incomplete=true; incompleteReason=T3AReason::Nonprogress;
                break;
            }
            activePhase=T3APhase::CandidateValidation;
#if defined(PAACS_TOPO_TESTING)
            const bool candidateValid=[&]() {
                TestingWallScope timer(
                    t3a_detail::T3AWallPhase::CandidateAudit);
                return validateRefinementCandidate(
                    removed,*boundary,*candidate,cells,edgeStars,nodeStars,
                    coordinates,baseLocation->index,locationIndex,ws,region,
                    metrics);
            }();
#else
            const bool candidateValid=validateRefinementCandidate(
                removed,*boundary,*candidate,cells,edgeStars,nodeStars,
                coordinates,baseLocation->index,locationIndex,ws,region,
                metrics);
#endif
            if (!candidateValid) {
                attempt.disposition=T3AAttemptDisposition::InvalidCandidate;
                attempt.reason=T3AReason::InvalidCandidate;
                out.status=T3AStatus::Failed;
                out.failure={T3AReason::InvalidCandidate,
                             T3APhase::CandidateValidation,true,
                             selected.row.key};
                completeQueueReport();
                return out;
            }

            topo::TetRefinementCavitySpec cavity;
            cavity.stamp=topo::WorkspaceStamp::of(ws);
            cavity.stagedNode={0,*candidate,region};
            cavity.removed=removed;
            for (const auto& face:*boundary) {
                topo::TetRefinementCellSpec cell;
                for (std::size_t i=0;i<3;++i)
                    cell.vertices[i]=topo::TetRefinementNodeRef::existingNode(
                        face.directed[i]);
                cell.vertices[3]=topo::TetRefinementNodeRef::stagedNode();
                cell.classification=region;
                std::array<Vec3,4> p{
                    coordinates.at(face.directed[0]),
                    coordinates.at(face.directed[1]),
                    coordinates.at(face.directed[2]),*candidate};
                const double volume=std::abs(mesh::tetVolume(p[0],p[1],p[2],p[3]));
                if (!std::isfinite(volume) ||
                    volume<request.refinement.minimumCreatedTetVolume) {
                    attempt.disposition=T3AAttemptDisposition::CreatedTetVolumeBlocked;
                    attempt.reason=T3AReason::CreatedTetVolumeBlocked;
                    ++selected.nextDescriptor;
                    if (selected.nextDescriptor<descriptorCap)
                        enqueue(selected);
                    else {
                        incomplete=true;
                        incompleteReason=T3AReason::CreatedTetVolumeBlocked;
                    }
                    break;
                }
                cavity.added.push_back(cell);
            }
            if (cavity.added.size()!=boundary->size()) {
                if (incomplete) break;
                continue;
            }

#if defined(PAACS_TOPO_TESTING)
            bool hasPendingBlockedInvalidation=false;
            for (const auto removedCell:removed)
                hasPendingBlockedInvalidation=
                    hasPendingBlockedInvalidation||
                    blocked.find(removedCell)!=nullptr;
            if (hasPendingBlockedInvalidation)
                recordBlockerStateDigest(
                    t3a_detail::T3ABlockerDigestCheckpoint::
                        InvalidationBeforeTransaction,
                    ws,queue,blocked,*out.report,metrics);
#endif
            auto tx=ws.beginTetRefinementTransaction();
            bool topologyPrepared=false;
            try {
                activePhase=T3APhase::TopologyPreparation;
#if defined(PAACS_TOPO_TESTING)
                std::optional<TestingWallScope> topologyTimer;
                topologyTimer.emplace(
                    t3a_detail::T3AWallPhase::TopologyPreparation);
#endif
                (void)tx.stage(cavity);
                const auto& prediction=tx.prepare();
#if defined(PAACS_TOPO_TESTING)
                topologyTimer.reset();
#endif
                topologyPrepared=true;
                ++metrics.preparedTransactions;
                const std::uint64_t topologyCount=1+prediction.cells.size()+
                    prediction.edges.size()+prediction.faces.size();
                recordConstructionWork(metrics,topologyCount);
                if (out.report->topologyHandlesUsed+topologyCount>
                    out.report->topologyHandleStorage.size()) {
                    tx.abandon();
                    attempt.disposition=T3AAttemptDisposition::TerminalFailure;
                    attempt.reason=T3AReason::ReportBudgetReached;
                    incomplete=true; incompleteReason=T3AReason::ReportBudgetReached;
                    break;
                }
                if (coordinates.size()+1>request.refinement.maxOutputNodes ||
                    cells.size()-removed.size()+boundary->size()>
                        request.refinement.maxOutputTets) {
                    tx.abandon();
                    attempt.disposition=T3AAttemptDisposition::TerminalFailure;
                    attempt.reason=T3AReason::OutputBudgetReached;
                    incomplete=true; incompleteReason=T3AReason::OutputBudgetReached;
                    break;
                }

                const PendingInsertionIdentity pendingIdentity{
                    prediction.stagedNode,metrics.committedInsertions};
                double oldWorstQuality=(std::numeric_limits<double>::infinity)();
                double oldMinimumVolume=(std::numeric_limits<double>::infinity)();
                if (selected.row.key.kind==TetViolationKind::Sliver) {
                    for (const auto removedCell:removed) {
                        const auto& conn=cells.at(removedCell);
                        std::array<Vec3,4> p{};
                        for (std::size_t i=0;i<4;++i)
                            p[i]=coordinates.at(conn[i]);
                        oldWorstQuality=(std::min)(oldWorstQuality,
                            mesh::tetQuality(p[0],p[1],p[2],p[3]));
                        oldMinimumVolume=(std::min)(oldMinimumVolume,
                            std::abs(mesh::tetVolume(p[0],p[1],p[2],p[3])));
                    }
                }
                CreatedCellRows createdRows;
                createdRows.reserve(boundary->size());
                for (const auto& transition:prediction.cells)
                    if (transition.transition==topo::TetRefinementTransition::Created) {
                        std::array<topo::NodeHandle,4> conn{};
                        for (std::size_t i=0;i<4;++i)
                            conn[i]=transition.key[i].kind==
                                    topo::TetRefinementNodeRefKind::Staged
                                ? prediction.stagedNode
                                : transition.key[i].existing;
                        createdRows.push_back({transition.predicted,conn});
                }
                bool candidateRejected=false;
                std::vector<QueuedViolation> createdViolations;
                createdViolations.reserve(createdRows.size());
                double newWorstQuality=(std::numeric_limits<double>::infinity)();
                double newMinimumVolume=(std::numeric_limits<double>::infinity)();
                for (const auto& [newCell,conn]:createdRows) {
                    std::array<Vec3,4> p{};
                    for (std::size_t i=0;i<4;++i)
                        p[i]=conn[i]==prediction.stagedNode
                            ? *candidate:coordinates.at(conn[i]);
                    if (selected.row.key.kind==TetViolationKind::Sliver) {
                        newWorstQuality=(std::min)(newWorstQuality,
                            mesh::tetQuality(p[0],p[1],p[2],p[3]));
                        newMinimumVolume=(std::min)(newMinimumVolume,
                            std::abs(mesh::tetVolume(p[0],p[1],p[2],p[3])));
                    }
                    auto next=evaluateViolation(
                        conn,p,region.tag,
                        result.existingGlobalByAuthoritativeNodeSlot,
                        inserted,spec,request.refinement,ledger,metrics,
                        generation++,&pendingIdentity);
                    if (next) {
                        if (selected.row.key.kind==TetViolationKind::RadiusEdge &&
                            !queueLess(selected.row,*next)) {
                            tx.abandon();
                            attempt.disposition=T3AAttemptDisposition::NoProgress;
                            attempt.reason=T3AReason::Nonprogress;
                            ++selected.nextDescriptor;
                            if (selected.nextDescriptor<descriptorCap)
                                enqueue(selected);
                            else {
                                incomplete=true;
                                incompleteReason=T3AReason::Nonprogress;
                            }
                            candidateRejected=true;
                            break;
                        }
                        if (selected.row.key.kind==TetViolationKind::Sliver&&
                            queueLess(*next,selected.row)) {
                            tx.abandon();
                                attempt.disposition=
                                    T3AAttemptDisposition::SliverImprovementFailed;
                                attempt.reason=T3AReason::Nonprogress;
                                ++selected.nextDescriptor;
                                if (selected.nextDescriptor<descriptorCap)
                                    enqueue(selected);
                                else {
                                    incomplete=true;
                                    incompleteReason=T3AReason::SliverPolicyExhausted;
                                }
                                candidateRejected=true;
                                break;
                        }
                        createdViolations.push_back({*next,newCell,conn});
                    }
                }
                if (!candidateRejected&&
                    selected.row.key.kind==TetViolationKind::Sliver&&
                    (!(newWorstQuality>=oldWorstQuality+
                           request.refinement.minimumSliverImprovement)||
                     !(newMinimumVolume>=oldMinimumVolume))) {
                    tx.abandon();
                    attempt.disposition=
                        T3AAttemptDisposition::SliverImprovementFailed;
                    attempt.reason=T3AReason::Nonprogress;
                    ++selected.nextDescriptor;
                    if (selected.nextDescriptor<descriptorCap)
                        enqueue(selected);
                    else {
                        incomplete=true;
                        incompleteReason=T3AReason::SliverPolicyExhausted;
                    }
                    candidateRejected=true;
                }
                if (candidateRejected) {
                    if (incomplete) break;
                    continue;
                }
                T3AReason sizingReason=T3AReason::None;
                if (!verifyPreparationSizing(
                        spec,selected.row,proposalSizing,metrics,sizingReason)) {
                    tx.abandon();
                    attempt.disposition=T3AAttemptDisposition::TerminalFailure;
                    attempt.reason=sizingReason;
                    out.status=T3AStatus::Failed;
                    out.failure={sizingReason,T3APhase::SizingEvaluation,
                                 true,selected.row.key};
                    completeQueueReport();
                    out.retainedResult.reset(); return out;
                }

                // Prepare only affected coordinator rows.  The authoritative
                // resolver maps, heap and spatial index remain untouched.
                activePhase=T3APhase::CoordinatorPreparation;
#if defined(PAACS_TOPO_TESTING)
                std::optional<TestingWallScope> coordinatorTimer;
                coordinatorTimer.emplace(
                    t3a_detail::T3AWallPhase::CoordinatorPreparation);
#endif
                auto prepared=std::make_unique<
                    PreparedTetRefinementCoordinatorDelta>();
                prepareCoordinatorDelta(
                    *prepared,prediction.stagedNode,*candidate,
                    metrics.committedInsertions,removed,createdRows,
                    createdViolations,coordinates,cells,edgeStars,faceOwners,
                    nodeStars,queue,blocked,*out.report,locationIndex,metrics);

                prepared->topologyRows.reserve(
                    static_cast<std::size_t>(topologyCount));
                T3ATopologyHandleRow nodeTopology;
                nodeTopology.kind=T3ATopologyKind::Node;
                nodeTopology.transition=T3ATopologyTransition::Created;
                nodeTopology.predictedNode=prediction.stagedNode;
                prepared->topologyRows.push_back(nodeTopology);
                for (const auto& row:prediction.cells) {
                    T3ATopologyHandleRow dest;
                    dest.kind=T3ATopologyKind::Cell;
                    dest.transition=publicTransition(row.transition);
                    dest.predictedCell=row.predicted;
                    prepared->topologyRows.push_back(dest);
                }
                for (const auto& row:prediction.edges) {
                    T3ATopologyHandleRow dest;
                    dest.kind=T3ATopologyKind::Edge;
                    dest.transition=publicTransition(row.transition);
                    dest.predictedEdge=row.predicted;
                    prepared->topologyRows.push_back(dest);
                }
                for (const auto& row:prediction.faces) {
                    T3ATopologyHandleRow dest;
                    dest.kind=T3ATopologyKind::Face;
                    dest.transition=publicTransition(row.transition);
                    dest.predictedFace=row.predicted;
                    prepared->topologyRows.push_back(dest);
                }
                if (prepared->topologyRows.size()!=topologyCount)
                    throw std::logic_error("prepared topology report mismatch");
#if defined(PAACS_TOPO_TESTING)
                coordinatorTimer.reset();
#endif

                const std::uint64_t nextContextAccounted=
                    applyAccountedContextDelta(
                        currentContextAccounted,*prepared);
                std::uint64_t preparedAccounted=
                    preparedCoordinatorAccountedBytes(*prepared);
                addAccountedBytes(preparedAccounted,
                    vectorAccountedBytes(createdRows));
                addAccountedBytes(preparedAccounted,
                    vectorAccountedBytes(createdViolations));
                addAccountedBytes(preparedAccounted,
                    vectorAccountedBytes(removed));
                addAccountedBytes(preparedAccounted,
                    vectorAccountedBytes(*boundary));
                addAccountedBytes(preparedAccounted,
                    vectorAccountedBytes(cavity.removed));
                addAccountedBytes(preparedAccounted,
                    vectorAccountedBytes(cavity.added));
                addAccountedBytes(preparedAccounted,
                    vectorAccountedBytes(prediction.cells));
                addAccountedBytes(preparedAccounted,
                    vectorAccountedBytes(prediction.edges));
                addAccountedBytes(preparedAccounted,
                    vectorAccountedBytes(prediction.faces));
                std::uint64_t insertionTransient=0;
                addAccountedBytes(insertionTransient,protectedAccounted);
                addAccountedBytes(insertionTransient,reportAccounted);
                addAccountedBytes(insertionTransient,currentContextAccounted);
                addAccountedBytes(insertionTransient,preparedAccounted);
                metrics.transientAccountedHighWaterBytes=(std::max)(
                    metrics.transientAccountedHighWaterBytes,
                    insertionTransient);

                const std::size_t topologyOffset=static_cast<std::size_t>(
                    out.report->topologyHandlesUsed);

                activePhase=T3APhase::TopologyCommit;
                const auto commit=tx.commit();
#if defined(PAACS_TOPO_TESTING)
                Phase3TailGuard phase3TailGuard;
#endif
                activePhase=T3APhase::CoordinatorCommit;
                applyCoordinatorDelta(
                    *prepared,coordinates,inserted,coordinateIndex,cells,
                    edgeStars,faceOwners,nodeStars,queue,blocked,*out.report,
                    locationIndex,metrics);
                currentContextAccounted=nextContextAccounted;
                for (std::size_t i=0;i<prepared->topologyRows.size();++i)
                    out.report->topologyHandleStorage[topologyOffset+i]=
                        prepared->topologyRows[i];
                std::size_t actualCursor=topologyOffset;
                out.report->topologyHandleStorage[actualCursor++].actualNode=
                    commit.actualNode;
                for (const auto& row:commit.cells)
                    out.report->topologyHandleStorage[actualCursor++].actualCell=
                        row.actual;
                for (const auto& row:commit.edges)
                    out.report->topologyHandleStorage[actualCursor++].actualEdge=
                        row.actual;
                for (const auto& row:commit.faces)
                    out.report->topologyHandleStorage[actualCursor++].actualFace=
                        row.actual;
                out.report->topologyHandlesUsed+=topologyCount;
                result.constraintState.workspaceStamp=commit.stamp;
                attempt.disposition=T3AAttemptDisposition::Committed;
                attempt.reason=T3AReason::None;
                attempt.workAfter=metrics.workUnitsConsumed;
                const std::size_t commitIndex=static_cast<std::size_t>(
                    out.report->committedPrefixUsed++);
                auto& reportCommit=out.report->committedPrefixStorage[commitIndex];
                reportCommit.insertionOrdinal=metrics.committedInsertions;
                reportCommit.violation=selected.row.key;
                reportCommit.predictedNode=commit.predictedNode;
                reportCommit.actualNode=commit.actualNode;
                reportCommit.topologyOffset=topologyOffset;
                reportCommit.topologyCount=topologyCount;
                reportCommit.postCommitRevision=commit.stamp.revision;
                ++metrics.committedInsertions;
                metrics.authoritativeAccountedBytes=workspaceAccountedBytes(ws);
                metrics.contextDeltaUpdates+=1;
                metrics.locallyReevaluatedCells+=createdRows.size();
                metrics.queuePushes+=createdViolations.size();
                metrics.queueHighWater=(std::max<std::uint64_t>)(
                    metrics.queueHighWater,queue.size());
#if defined(PAACS_TOPO_TESTING)
                phase3TailGuard.stop();
#endif
                activePhase=T3APhase::LocalAudit;
                bool localAudit=false;
                try {
#if defined(PAACS_TOPO_TESTING)
                    TestingWallScope timer(
                        t3a_detail::T3AWallPhase::LocalAudit);
#endif
                    localAudit=auditCommittedDelta(commit,ws)&&
                        cells.count(commit.cells.back().actual)>0;
                } catch (const std::bad_alloc&) {
                    out.status=T3AStatus::Failed;
                    out.failure={T3AReason::AllocationFailed,T3APhase::LocalAudit,
                                 true,selected.row.key};
                    completeQueueReport();
                    out.retainedResult.reset(); return out;
                }
#if defined(PAACS_TOPO_TESTING)
                if (t3a_detail::injectPostCommitAuditFailureForTesting())
                    localAudit=false;
#endif
                if (!localAudit) {
                    out.status=T3AStatus::Failed;
                    out.failure={T3AReason::InternalError,T3APhase::LocalAudit,
                                 true,selected.row.key};
                    completeQueueReport();
                    out.retainedResult.reset(); return out;
                }
            } catch (const WorkBudgetStop&) {
                tx.abandon();
                throw;
            } catch (const std::bad_alloc&) {
                tx.abandon();
#if defined(PAACS_TOPO_TESTING)
                if (hasPendingBlockedInvalidation)
                    recordBlockerStateDigest(
                        t3a_detail::T3ABlockerDigestCheckpoint::
                            InvalidationFailureAfterAbandon,
                        ws,queue,blocked,*out.report,metrics);
#endif
                attempt.disposition=T3AAttemptDisposition::TerminalFailure;
                attempt.reason=T3AReason::AllocationFailed;
                out.status=T3AStatus::Failed;
                out.failure={T3AReason::AllocationFailed,
                             topologyPrepared?T3APhase::CoordinatorPreparation:
                                               T3APhase::TopologyPreparation,
                             true,selected.row.key};
                completeQueueReport();
                out.retainedResult.reset(); return out;
            } catch (...) {
                if (topologyPrepared) {
                    tx.abandon();
                    attempt.disposition=T3AAttemptDisposition::TerminalFailure;
                    attempt.reason=T3AReason::CoordinatorPreparationFailed;
                    out.status=T3AStatus::Failed;
                    out.failure={T3AReason::CoordinatorPreparationFailed,
                                 T3APhase::CoordinatorPreparation,true,
                                 selected.row.key};
                    completeQueueReport();
                    out.retainedResult.reset(); return out;
                }
                attempt.disposition=T3AAttemptDisposition::InvalidCandidate;
                attempt.reason=T3AReason::InvalidCandidate;
                incomplete=true; incompleteReason=T3AReason::InvalidCandidate;
                break;
            }
            } catch (const std::bad_alloc&) {
                if (activeAttempt) {
                    activeAttempt->disposition=
                        T3AAttemptDisposition::TerminalFailure;
                    activeAttempt->reason=T3AReason::AllocationFailed;
                    activeAttempt->workAfter=metrics.workUnitsConsumed;
                }
#if defined(PAACS_TOPO_TESTING)
                if (activePhase==T3APhase::Encroachment)
                    recordBlockerStateDigest(
                        t3a_detail::T3ABlockerDigestCheckpoint::
                            BlockTransitionTerminalAfter,
                        ws,queue,blocked,*out.report,metrics);
#endif
                out.status=T3AStatus::Failed;
                out.failure={T3AReason::AllocationFailed,activePhase,
                             true,selected.row.key};
                completeQueueReport();
                return out;
            } catch (const WorkBudgetStop&) {
                if (activeAttempt) {
                    activeAttempt->disposition=
                        T3AAttemptDisposition::TerminalFailure;
                    activeAttempt->reason=T3AReason::WorkBudgetReached;
                    activeAttempt->workAfter=metrics.workUnitsConsumed;
                }
                out.report->hasTerminalViolation=true;
                out.report->terminalViolation=selected.row.key;
                incomplete=true;
                incompleteReason=T3AReason::WorkBudgetReached;
                incompletePhase=activePhase;
#if defined(PAACS_TOPO_TESTING)
                if (activePhase==T3APhase::Encroachment)
                    recordBlockerStateDigest(
                        t3a_detail::T3ABlockerDigestCheckpoint::
                            BlockTransitionTerminalAfter,
                        ws,queue,blocked,*out.report,metrics);
#endif
                break;
            }
        }

        if (!incomplete&&queue.empty()&&!blocked.empty()) {
            incomplete=true;
            incompleteReason=T3AReason::BoundaryRefinementRequired;
            out.report->hasTerminalViolation=true;
            out.report->terminalViolation=blocked.top().row.key;
        }

        workBudget.stopCounting();
        activePhase=T3APhase::FinalQualityAudit;
#if defined(PAACS_TOPO_TESTING)
        TestingWallScope finalTimer(t3a_detail::T3AWallPhase::FinalAudit);
#endif
        T3AReason finalSizing=T3AReason::None;
        if (!verifySizingLedger(spec,ledger,metrics,finalSizing)) {
            out.status=T3AStatus::Failed;
            out.failure={finalSizing,T3APhase::FinalQualityAudit};
            completeQueueReport();
            out.retainedResult.reset(); return out;
        }
        if (!topo::workspaceTopologyReady(
                ws.nodes(),ws.cells(),ws.edges(),ws.faces())) {
            out.status=T3AStatus::Failed;
            out.failure={T3AReason::QualityAuditFailed,
                         T3APhase::FinalQualityAudit};
            completeQueueReport();
            out.retainedResult.reset(); return out;
        }
        if (!finalProtectedAudit(result.constraintState,ws,
                result.existingGlobalByAuthoritativeNodeSlot)) {
            out.status=T3AStatus::Failed;
            out.failure={T3AReason::ProtectedIdentityChanged,T3APhase::LocalAudit};
            completeQueueReport();
            out.retainedResult.reset(); return out;
        }
        if (!fillUnattempted(*out.report,queue)) {
            out.status=T3AStatus::Failed;
            out.failure={T3AReason::InternalError,T3APhase::FinalQualityAudit};
            out.retainedResult.reset(); return out;
        }
#if defined(PAACS_TOPO_TESTING)
        recordBlockerStateDigest(
            t3a_detail::T3ABlockerDigestCheckpoint::Final,
            ws,queue,blocked,*out.report,metrics);
#endif

        activePhase=T3APhase::MeshPreparation;
        mesh::Mesh publication;
        if (!buildPublicationMesh(
                *seed,ws,result.existingGlobalByAuthoritativeNodeSlot,inserted,
                publication,result.publishedGlobalByAuthoritativeNodeSlot,
                result.publishedGlobalByAuthoritativeCellSlot,result.quality)) {
            out.status=T3AStatus::Failed;
            out.failure={T3AReason::QualityAuditFailed,T3APhase::FinalQualityAudit};
            out.retainedResult.reset(); return out;
        }
        activePhase=T3APhase::FinalQualityAudit;
        if (!auditBoundaryBijection(
                publication,result.immutableBoundary,ws,
                result.publishedGlobalByAuthoritativeNodeSlot)) {
            out.status=T3AStatus::Failed;
            out.failure={T3AReason::BoundaryBijectionFailed,
                         T3APhase::FinalQualityAudit};
            out.retainedResult.reset(); return out;
        }
        if (!completeQualityReport(publication,spec,request.refinement,ledger,
                                   metrics,result.quality)) {
            out.status=T3AStatus::Failed;
            out.failure={T3AReason::QualityAuditFailed,
                         T3APhase::FinalQualityAudit};
            out.retainedResult.reset(); return out;
        }
        metrics.finalCellsAudited=cells.size();
        metrics.authoritativeAccountedBytes=workspaceAccountedBytes(ws);
        metrics.publicationAccountedBytes=flatMeshAccountedBytes(publication);
        metrics.qualityComplete=true;
        result.quality.sliverPolicy=request.refinement.sliverPolicy;
        result.quality.sliverGoalComplete=
            request.refinement.sliverPolicy==TetSliverPolicy::ReportOnly ||
            result.quality.sliversRemaining==0;

        if (incomplete) {
            out.status=T3AStatus::Incomplete;
            out.failure={incompleteReason,
                         incompletePhase!=T3APhase::None
                             ? incompletePhase
                             : incompleteReason==
                                   T3AReason::BoundaryRefinementRequired
                                 ? T3APhase::Encroachment
                                 : T3APhase::CandidateValidation};
            if (out.report->hasTerminalViolation) {
                out.failure.hasViolation=true;
                out.failure.violation=out.report->terminalViolation;
            }
            if (incompleteReason==T3AReason::BoundaryRefinementRequired&&
                !blocked.empty()) {
                const auto rowIndex=blocked.top().blockedRowIndex;
                if (rowIndex>=out.report->blockedUsed) std::abort();
                const auto& row=out.report->blockedStorage[
                    static_cast<std::size_t>(rowIndex)];
                if (row.witnessCount==0||
                    row.witnessOffset>=out.report->protectedWitnessesUsed)
                    std::abort();
                out.failure.hasProtectedWitness=true;
                out.failure.protectedWitness=
                    out.report->protectedWitnessStorage[
                        static_cast<std::size_t>(row.witnessOffset)];
            }
            result.publicationReady=false;
            out.retainedResult=std::move(workingResult);
            return out;
        }

        result.publicationReady=true;
        metrics.runComplete=true;
        activePhase=T3APhase::Publication;
        boundary_conforming_detail::publicationTailActive()=true;
        destination.swap(publication);
        boundary_conforming_detail::publicationTailActive()=false;
        out.status=T3AStatus::Published;
        out.destinationUnchanged=false;
        out.retainedResult=std::move(workingResult);
        return out;
    } catch (const SizingFailure& failure) {
        out.status=T3AStatus::Failed;
        out.failure={failure.reason,
                     activePhase==T3APhase::FinalQualityAudit
                         ? T3APhase::FinalQualityAudit
                         : T3APhase::SizingEvaluation};
        out.retainedResult.reset();
        return out;
    } catch (const std::bad_alloc&) {
        out.status=T3AStatus::Failed;
        out.failure={T3AReason::AllocationFailed,activePhase};
        out.retainedResult.reset();
        return out;
    } catch (...) {
        out.status=T3AStatus::Failed;
        out.failure={T3AReason::InternalError,T3APhase::InitialAudit};
        out.retainedResult.reset();
        return out;
    }
}

namespace {

bool validT4Policy(const T4AQualityPolicy& policy) noexcept {
    if (!std::isfinite(policy.minimumRadiusRatio) ||
        policy.minimumRadiusRatio<0.0 || policy.minimumRadiusRatio>1.0 ||
        policy.maxResultNodes==0 || policy.maxResultTets==0 ||
        policy.maxReportRows==0 || policy.maxWorstTailRows==0)
        return false;
    if (policy.mode==T4AQualityMode::ReportOnly)
        return policy.minimumImprovement==0.0 && policy.maxOperations==0 &&
               policy.maxSelectedGenerations==0 &&
               policy.maxCandidateAttemptsPerCell==0 &&
               policy.maxCreatedTetsPerOperation==0 &&
               policy.maxWorkUnits==0 && policy.maxNewNodes==0 &&
               policy.maxDependencyRowsPerAttempt==0 &&
               policy.maxProtectedWitnessesPerAttempt==0 &&
               policy.maxTopologyRowsPerOperation==0;
    return policy.minimumRadiusRatio>0.0 &&
           std::isfinite(policy.minimumImprovement) &&
           policy.minimumImprovement>0.0 &&
           policy.minimumImprovement<=1.0 && policy.maxOperations>0 &&
           policy.maxSelectedGenerations>0 &&
           policy.maxCandidateAttemptsPerCell==1 &&
           policy.maxCreatedTetsPerOperation>0 &&
           policy.maxWorkUnits>0 && policy.maxNewNodes>0 &&
           policy.maxDependencyRowsPerAttempt>0 &&
           policy.maxProtectedWitnessesPerAttempt>0 &&
           policy.maxTopologyRowsPerOperation>0;
}

T4AStatus t4SourceStatus(T3AStatus value) noexcept {
    switch (value) {
        case T3AStatus::Published: return T4AStatus::Published;
        case T3AStatus::Incomplete: return T4AStatus::Incomplete;
        case T3AStatus::Rejected: return T4AStatus::Rejected;
        case T3AStatus::Unsupported: return T4AStatus::Unsupported;
        case T3AStatus::Failed: return T4AStatus::Failed;
    }
    return T4AStatus::Failed;
}

std::uint64_t digestIndexVector(std::span<const Index> values) noexcept {
    std::uint64_t digest=0xcbf29ce484222325ull;
    digestWord(digest,static_cast<std::uint64_t>(values.size()));
    for (const auto value:values)
        digestWord(digest,static_cast<std::uint64_t>(
            static_cast<std::int64_t>(value)));
    return digest;
}

bool buildT4PublicationBase(const mesh::Mesh& source,mesh::Mesh& base) {
    std::size_t cells=0,connectivity=0;
    for (Index i=0;i<source.cellCount();++i) {
        const auto cell=source.cell(i);
        if (cell.type==mesh::CellType::Tet4) continue;
        ++cells;
        connectivity+=static_cast<std::size_t>(mesh::nodesPerCell(cell.type));
    }
    base.reserveAppend(static_cast<std::size_t>(source.nodeCount()),
                       cells,connectivity);
    for (Index i=0;i<source.nodeCount();++i)
        base.addNode(source.node(i),source.nodeClass(i));
    for (Index i=0;i<source.cellCount();++i) {
        const auto cell=source.cell(i);
        if (cell.type==mesh::CellType::Tet4) continue;
        const int count=mesh::nodesPerCell(cell.type);
        std::vector<Index> nodes(static_cast<std::size_t>(count));
        const Index* input=source.cellNodes(i);
        for (int j=0;j<count;++j)
            nodes[static_cast<std::size_t>(j)]=input[j];
        base.addCell(cell.type,nodes,cell.onEntity);
    }
    return true;
}

bool reconstructT4VertexIdentity(
    const T3AResult& source,const T3AReport& report,
    const topo::TopologyWorkspace& workspace,
    std::map<topo::NodeHandle,RefinementVertexKey>& result,
    std::uint64_t& committedCount) {
    committedCount=source.metrics.committedInsertions;
    if (report.committedPrefixUsed!=committedCount ||
        report.committedPrefixUsed>report.committedPrefixStorage.size())
        return false;
    const auto active=workspace.nodes().activeBySequence();
    for (const auto node:active) {
        if (node.slot<0) return false;
        const auto slot=static_cast<std::size_t>(node.slot);
        if (slot<source.existingGlobalByAuthoritativeNodeSlot.size() &&
            source.existingGlobalByAuthoritativeNodeSlot[slot]!=kInvalidIndex) {
            const auto global=source.existingGlobalByAuthoritativeNodeSlot[slot];
            if (global<0 ||
                !result.emplace(node,RefinementVertexKey{
                    RefinementVertexKind::ExistingGlobal,
                    static_cast<std::uint64_t>(global)}).second)
                return false;
        }
    }
    for (std::uint64_t i=0;i<committedCount;++i) {
        const auto& row=report.committedPrefixStorage[
            static_cast<std::size_t>(i)];
        if (row.insertionOrdinal!=i || row.actualNode!=row.predictedNode ||
            !workspace.nodes().isCurrent(row.actualNode) ||
            result.count(row.actualNode)!=0 ||
            !result.emplace(row.actualNode,RefinementVertexKey{
                RefinementVertexKind::InteriorInsertion,i}).second)
            return false;
    }
    if (result.size()!=active.size()) return false;
    std::set<RefinementVertexKey> unique;
    for (const auto& [handle,key]:result) {
        if (!workspace.nodes().isCurrent(handle) || !unique.insert(key).second)
            return false;
    }
    for (const auto node:active) {
        const auto slot=static_cast<std::size_t>(node.slot);
        if (slot>=source.publishedGlobalByAuthoritativeNodeSlot.size())
            return false;
        const auto global=source.publishedGlobalByAuthoritativeNodeSlot[slot];
        if (global<0) return false;
    }
    for (const auto cell:workspace.cells().activeBySequence()) {
        if (cell.slot<0 ||
            static_cast<std::size_t>(cell.slot)>=
                source.publishedGlobalByAuthoritativeCellSlot.size() ||
            source.publishedGlobalByAuthoritativeCellSlot[
                static_cast<std::size_t>(cell.slot)]<0)
            return false;
    }
    return true;
}

bool checkedTetQualityCanonical(
    const std::array<Vec3,4>& coordinates,
    const std::array<RefinementVertexKey,4>& keys,double& quality) noexcept {
    for (const auto& coordinate:coordinates)
        if (!finite(coordinate)) return false;
    for (std::size_t i=0;i<4;++i)
        for (std::size_t j=i+1;j<4;++j)
            if (keys[i]==keys[j]) return false;
    if (predicates::orient3d(coordinates[0],coordinates[1],coordinates[2],
                             coordinates[3])!=predicates::Sign::Positive)
        return false;
    std::array<std::pair<RefinementVertexKey,Vec3>,4> canonical{};
    for (std::size_t i=0;i<4;++i) canonical[i]={keys[i],coordinates[i]};
    std::sort(canonical.begin(),canonical.end(),[](const auto& a,const auto& b) {
        return a.first<b.first;
    });
    quality=mesh::tetQuality(canonical[0].second,canonical[1].second,
                             canonical[2].second,canonical[3].second);
    if (!std::isfinite(quality) || quality<0.0 || quality>1.0) return false;
    if (quality==0.0) quality=0.0;
    return true;
}

T4AQualityKey t4QualityKey(const TetViolationRow& row) noexcept {
    T4AQualityKey result;
    result.radiusRatioBits=
        (std::numeric_limits<std::uint64_t>::max)()-row.severityBits;
    result.tet=row.key.tet;
    result.generation=row.generation;
    return result;
}

double t4QualityValue(const TetViolationRow& row) noexcept {
    return std::bit_cast<double>(t4QualityKey(row).radiusRatioBits);
}

std::optional<TetViolationRow> evaluateT4QualityRow(
    const std::array<topo::NodeHandle,4>& nodes,
    const std::array<Vec3,4>& coordinates,
    const std::vector<Index>& existing,const InsertedMap& inserted,
    double threshold,std::uint64_t generation,T4AMetrics& metrics,
    TetRefinementMetrics& coreMetrics,
    const PendingInsertionIdentity* pending=nullptr) {
    std::array<RefinementVertexKey,4> keys{};
    for (std::size_t i=0;i<4;++i)
        keys[i]=vertexKey(nodes[i],existing,inserted,pending);
    double quality=0.0;
    ++metrics.qualityEvaluations;
    ++coreMetrics.violationEvaluations;
    recordConstructionWork(coreMetrics);
    if (!checkedTetQualityCanonical(coordinates,keys,quality))
        throw std::logic_error("invalid T4 Tet quality input");
    if (!(quality<threshold)) return std::nullopt;
    TetViolationRow row;
    row.key=makeTetKey(nodes,existing,inserted,pending);
    row.key.kind=TetViolationKind::Sliver;
    row.normalizedSeverity=quality;
    row.severityBits=(std::numeric_limits<std::uint64_t>::max)()-
                     std::bit_cast<std::uint64_t>(quality);
    row.generation=generation;
    return row;
}

bool prepareT4Report(T4AReport& report,const T4AQualityPolicy& policy,
                     std::uint64_t initialTets) {
    std::uint64_t attempts=0,dependencies=0,witnesses=0,topology=0,
                  unattempted=0,tailEach=0,tail=0,total=0,next=0;
    if (!checkedMul(policy.maxSelectedGenerations,
                    policy.maxCandidateAttemptsPerCell,attempts) ||
        !checkedMul(attempts,policy.maxDependencyRowsPerAttempt,dependencies) ||
        !checkedMul(attempts,policy.maxProtectedWitnessesPerAttempt,witnesses) ||
        !checkedMul(policy.maxOperations,policy.maxTopologyRowsPerOperation,
                    topology) ||
        !checkedMul(policy.maxOperations,policy.maxCreatedTetsPerOperation,next) ||
        (policy.mode==T4AQualityMode::BoundedImprove &&
         !checkedAdd(initialTets,next,unattempted)))
        return false;
    tailEach=(std::min)(policy.maxWorstTailRows,policy.maxResultTets);
    if (!checkedMul(tailEach,2,tail)) return false;
    for (const auto value:{attempts,policy.maxSelectedGenerations,
             policy.maxOperations,dependencies,witnesses,topology,
             unattempted,tail}) {
        if (!checkedAdd(total,value,next)) return false;
        total=next;
    }
    if (total>policy.maxReportRows || !fitsSize(attempts) ||
        !fitsSize(policy.maxSelectedGenerations) ||
        !fitsSize(policy.maxOperations) || !fitsSize(dependencies) ||
        !fitsSize(witnesses) || !fitsSize(topology) ||
        !fitsSize(unattempted) || !fitsSize(tail))
        return false;
    report.attemptStorage.resize(static_cast<std::size_t>(attempts));
    report.blockedStorage.resize(
        static_cast<std::size_t>(policy.maxSelectedGenerations));
    report.commitStorage.resize(static_cast<std::size_t>(policy.maxOperations));
    report.dependencyStorage.resize(static_cast<std::size_t>(dependencies));
    report.witnessStorage.resize(static_cast<std::size_t>(witnesses));
    report.topologyStorage.resize(static_cast<std::size_t>(topology));
    report.unattemptedStorage.resize(static_cast<std::size_t>(unattempted));
    report.tailStorage.resize(static_cast<std::size_t>(tail));
    report.reportComplete=true;
    return true;
}

bool fillT4Summary(
    const CellMap& cells,const CoordinateMap& coordinates,
    const std::vector<Index>& existing,const InsertedMap& inserted,
    double threshold,std::uint64_t tailOffset,std::uint64_t tailCapacity,
    T4AReport& report,T4AQualitySummary& summary,T4AMetrics& metrics,
    bool finalPass,ViolationHeap* queue=nullptr,
    std::uint64_t* generation=nullptr,
    TetRefinementMetrics* coreMetrics=nullptr) {
    std::vector<T4AQualityKey> rows;
    rows.reserve(cells.size());
    summary={};
    summary.minimum=cells.empty()?0.0:
        (std::numeric_limits<double>::infinity)();
    for (const auto& [cell,nodes]:cells) {
        if (coreMetrics) {
            const auto maximum=(std::numeric_limits<std::uint64_t>::max)();
            if (coreMetrics->workUnitsConsumed<maximum)
                ++coreMetrics->workUnitsConsumed;
        }
        std::array<Vec3,4> p{};
        std::array<RefinementVertexKey,4> keys{};
        for (std::size_t i=0;i<4;++i) {
            p[i]=coordinates.at(nodes[i]);
            keys[i]=vertexKey(nodes[i],existing,inserted);
        }
        if (predicates::orient3d(p[0],p[1],p[2],p[3])!=
                predicates::Sign::Positive)
            return false;
        ++summary.exactPositive;
        double quality=0.0;
        ++metrics.qualityEvaluations;
        if (!checkedTetQualityCanonical(p,keys,quality)) return false;
        T4AQualityKey key;
        key.radiusRatioBits=std::bit_cast<std::uint64_t>(quality);
        key.tet=makeTetKey(nodes,existing,inserted).tet;
        rows.push_back(key);
        if (quality<threshold && queue && generation && coreMetrics) {
            TetViolationRow violation;
            violation.key=makeTetKey(nodes,existing,inserted);
            violation.key.kind=TetViolationKind::Sliver;
            violation.normalizedSeverity=quality;
            violation.severityBits=
                (std::numeric_limits<std::uint64_t>::max)()-
                key.radiusRatioBits;
            violation.generation=(*generation)++;
            recordConstructionWork(*coreMetrics);
            ++coreMetrics->violationEvaluations;
            queue->push({violation,cell,nodes});
            ++metrics.queuePushes;
            metrics.queueHighWater=(std::max<std::uint64_t>)(
                metrics.queueHighWater,queue->size());
        }
        summary.minimum=(std::min)(summary.minimum,quality);
        summary.belowThreshold+=quality<threshold?1U:0U;
        const auto bin=(std::min<std::size_t>)(9,
            static_cast<std::size_t>(quality*10.0));
        ++summary.histogram[bin];
    }
    summary.tetCount=rows.size();
    std::sort(rows.begin(),rows.end());
    const auto percentile=[&](double p) {
        if (rows.empty()) return 0.0;
        const auto index=static_cast<std::size_t>(
            std::ceil(p*static_cast<double>(rows.size())))-1;
        return std::bit_cast<double>(rows[(std::min)(index,rows.size()-1)].
                                     radiusRatioBits);
    };
    summary.p1=percentile(0.01);
    summary.p5=percentile(0.05);
    summary.p50=percentile(0.50);
    summary.tailOffset=tailOffset;
    summary.tailCount=(std::min<std::uint64_t>)(tailCapacity,rows.size());
    summary.tailTruncated=summary.tailCount<rows.size();
    if (tailOffset+summary.tailCount>report.tailStorage.size()) return false;
    for (std::uint64_t i=0;i<summary.tailCount;++i)
        report.tailStorage[static_cast<std::size_t>(tailOffset+i)]=
            rows[static_cast<std::size_t>(i)];
    report.tailUsed=(std::max)(report.tailUsed,tailOffset+summary.tailCount);
    if (finalPass) metrics.finalCellsAudited=rows.size();
    else metrics.initialCellsScanned=rows.size();
    return true;
}

void syncT4CoreMetrics(const TetRefinementMetrics& core,T4AMetrics& metrics) {
    metrics.protectedIndexVisits=core.protectedIndexVisits;
    metrics.locationIndexVisits=core.locationIndexVisits;
    metrics.locationTetTests=core.locationTetTests;
    metrics.cavityCellsVisited=core.cavityCellsVisited;
    metrics.cavityFaceVisits=core.cavityFaceVisits;
    metrics.survivorIndexVisits=core.survivorIndexVisits;
    metrics.candidateSurvivorTests=core.candidateSurvivorTests;
    metrics.localRowsPrepared=core.coordinatorRowsPrepared;
    metrics.localRowsApplied=core.coordinatorRowsPrepared;
    metrics.blockedInvalidations=core.blockedRowsInvalidated;
    metrics.preparedTransactions=core.preparedTransactions;
    metrics.perOperationGlobalTetScans=core.perInsertionGlobalTetScans;
    metrics.contextRebuilds=core.contextRebuilds;
    metrics.phase3Allocations=core.phase3Allocations;
    metrics.publicationTailAllocations=core.publicationTailAllocations;
    metrics.workUnitsConsumed=core.workUnitsConsumed;
    metrics.transientAccountedHighWaterBytes=
        core.transientAccountedHighWaterBytes;
}

bool fillT4Unattempted(T4AReport& report,const ViolationHeap& queue,
                       const CellMap& cells) noexcept {
    report.unattemptedUsed=0;
    for (const auto& queued:queue.rows()) {
        const auto found=cells.find(queued.cell);
        if (found==cells.end()||found->second!=queued.nodes) continue;
        if (report.unattemptedUsed>=report.unattemptedStorage.size())
            return false;
        report.unattemptedStorage[
            static_cast<std::size_t>(report.unattemptedUsed++)]=
            t4QualityKey(queued.row);
    }
    std::sort(report.unattemptedStorage.begin(),
              report.unattemptedStorage.begin()+
                  static_cast<std::ptrdiff_t>(report.unattemptedUsed));
    return true;
}

T4ADependencyKey t4CellDependency(const TetViolationRow& selected) noexcept {
    T4ADependencyKey result;
    result.kind=T4ADependencyKind::Cell;
    result.vertices=selected.key.tet.vertices;
    result.vertexCount=4;
    return result;
}

template<std::size_t N>
T4ADependencyKey t4Dependency(
    T4ADependencyKind kind,const std::array<topo::NodeHandle,N>& nodes,
    const std::vector<Index>& existing,const InsertedMap& inserted,
    const PendingInsertionIdentity* pending=nullptr) {
    static_assert(N>=2&&N<=4);
    T4ADependencyKey result;
    result.kind=kind;
    result.vertexCount=static_cast<std::uint8_t>(N);
    for (std::size_t i=0;i<N;++i)
        result.vertices[i]=vertexKey(nodes[i],existing,inserted,pending);
    std::sort(result.vertices.begin(),result.vertices.begin()+N);
    return result;
}

void canonicalizeT4Dependencies(std::vector<T4ADependencyKey>& values) {
    std::sort(values.begin(),values.end());
    values.erase(std::unique(values.begin(),values.end()),values.end());
}

bool auditT4Identity(
    const topo::TopologyWorkspace& workspace,
    const std::map<topo::NodeHandle,RefinementVertexKey>& identity,
    const std::vector<Index>& existing,const InsertedMap& inserted,
    const CellMap& cells) {
    const auto activeNodes=workspace.nodes().activeBySequence();
    const auto activeCells=workspace.cells().activeBySequence();
    if (identity.size()!=activeNodes.size()||cells.size()!=activeCells.size())
        return false;
    std::set<RefinementVertexKey> unique;
    std::vector<std::uint64_t> insertionOrdinals;
    insertionOrdinals.reserve(inserted.size());
    for (const auto node:activeNodes) {
        const auto found=identity.find(node);
        if (found==identity.end()||!unique.insert(found->second).second)
            return false;
        if (!(found->second==vertexKey(node,existing,inserted))) return false;
        if (found->second.kind==RefinementVertexKind::InteriorInsertion)
            insertionOrdinals.push_back(found->second.value);
    }
    std::sort(insertionOrdinals.begin(),insertionOrdinals.end());
    for (std::size_t i=0;i<insertionOrdinals.size();++i)
        if (insertionOrdinals[i]!=i) return false;
    for (const auto cell:activeCells) {
        const auto found=cells.find(cell);
        if (found==cells.end()) return false;
        const auto data=workspace.cells().get(cell);
        std::array<topo::NodeHandle,4> actual{
            data.nodes[0],data.nodes[1],data.nodes[2],data.nodes[3]};
        auto expected=found->second;
        std::sort(actual.begin(),actual.end());
        std::sort(expected.begin(),expected.end());
        if (actual!=expected) return false;
    }
    return true;
}

} // namespace

T4AOutcome meshModelQualityTet4(
    const geo::Model& model,const sizing::MeshSpec& spec,
    const T4ARequest& request,mesh::Mesh& destination,
    const VolumeBoundaryOwnershipProvider& ownershipProvider) noexcept {
#if defined(PAACS_TOPO_TESTING)
    TestingT4AllocationScope allocationScope;
#endif
    T4AOutcome out;
    try {
        out.report=std::make_unique<T4AReport>();
        out.sourceOutcome=std::make_unique<T3AOutcome>();
        out.report->shellPrepared=true;
    } catch (const std::bad_alloc&) {
        out.failure={T4AReason::AllocationFailed,T4APhase::ReportShell};
        return out;
    }
    if (!validT4Policy(request.quality)) {
        out.status=T4AStatus::Rejected;
        out.failure={T4AReason::InvalidRequest,T4APhase::RequestPreflight};
        out.report->reportComplete=true;
        out.reportComplete=true;
        return out;
    }

    std::unique_ptr<mesh::Mesh> sourceMesh;
    try {
        sourceMesh=std::make_unique<mesh::Mesh>();
    } catch (const std::bad_alloc&) {
        out.failure={T4AReason::AllocationFailed,T4APhase::SourceT3A};
        return out;
    }
#if defined(PAACS_TOPO_TESTING)
    allocationScope.suspend();
#endif
    *out.sourceOutcome=meshModelRefinedTet4(
        model,spec,request.source,*sourceMesh,ownershipProvider);
#if defined(PAACS_TOPO_TESTING)
    allocationScope.resume();
#endif
    if (!out.sourceOutcome->published()||
        !out.sourceOutcome->retainedResult) {
        out.status=t4SourceStatus(out.sourceOutcome->status);
        out.failure={T4AReason::SourceT3AStopped,T4APhase::SourceT3A};
        out.report->reportComplete=true;
        out.reportComplete=true;
        return out;
    }

    T4APhase activePhase=T4APhase::ReportEnvelope;
    try {
        auto& source=*out.sourceOutcome;
        auto& sourceResult=*source.retainedResult;
        if (!sourceResult.authoritative||!sourceResult.authoritative->workspace||
            !source.report||!source.sourceReceipt)
            throw std::logic_error("T3 result is incomplete");
        auto& sourceWorkspace=*sourceResult.authoritative->workspace;
        const auto sourceNodes=static_cast<std::uint64_t>(
            sourceWorkspace.nodes().activeBySequence().size());
        const auto sourceTets=static_cast<std::uint64_t>(
            sourceWorkspace.cells().activeBySequence().size());

        if (!prepareT4Report(*out.report,request.quality,sourceTets)) {
            out.status=T4AStatus::Rejected;
            out.failure={T4AReason::CapacityOverflow,T4APhase::ReportEnvelope};
            return out;
        }
        out.report->reportComplete=true;
        out.reportComplete=true;
        activePhase=T4APhase::SourceT3A;

        std::uint64_t requiredNodes=0,createdAllowance=0,requiredTets=0;
        if (!checkedAdd(sourceNodes,request.quality.maxNewNodes,requiredNodes)||
            !checkedMul(request.quality.maxOperations,
                        request.quality.maxCreatedTetsPerOperation,
                        createdAllowance)||
            !checkedAdd(sourceTets,createdAllowance,requiredTets)) {
            out.status=T4AStatus::Rejected;
            out.failure={T4AReason::CapacityOverflow,T4APhase::ReportEnvelope};
            return out;
        }
        if (request.quality.maxResultNodes<requiredNodes||
            request.quality.maxResultTets<requiredTets) {
            out.status=T4AStatus::Incomplete;
            out.failure={T4AReason::ResultBudgetReached,
                         T4APhase::ReportEnvelope};
            return out;
        }

        auto receipt=std::make_unique<T4AT3SourceReceipt>();
        auto workingResult=std::make_unique<T4AResult>();
        auto artifact=std::make_unique<T4AMutableArtifact>();
        auto publicationBase=std::make_unique<mesh::Mesh>();
        buildT4PublicationBase(*sourceMesh,*publicationBase);

        std::map<topo::NodeHandle,RefinementVertexKey> vertexIdentity;
        std::uint64_t t3InsertionCount=0;
        if (!reconstructT4VertexIdentity(
                sourceResult,*source.report,sourceWorkspace,
                vertexIdentity,t3InsertionCount)) {
            out.status=T4AStatus::Failed;
            out.failure={T4AReason::InvalidAuthoritativeSource,
                         T4APhase::SourceT3A};
            return out;
        }

        artifact->sourceGlobalByAuthoritativeNodeSlot=
            sourceResult.publishedGlobalByAuthoritativeNodeSlot;
        artifact->sourceGlobalByAuthoritativeCellSlot=
            sourceResult.publishedGlobalByAuthoritativeCellSlot;
        receipt->status=source.status;
        receipt->failure=source.failure;
        receipt->reportComplete=source.reportComplete;
        receipt->metrics=sourceResult.metrics;
        receipt->quality=sourceResult.quality;
        receipt->boundaryAudit=sourceResult.boundaryAudit;
        receipt->strictHash=mesh::meshHash(*sourceMesh);
        receipt->semanticHash=mesh::semanticHashEx(*sourceMesh);
        receipt->nodeCount=static_cast<std::uint64_t>(sourceMesh->nodeCount());
        receipt->cellCount=static_cast<std::uint64_t>(sourceMesh->cellCount());
        receipt->nodeMappingDigest=digestIndexVector(
            sourceResult.publishedGlobalByAuthoritativeNodeSlot);
        receipt->cellMappingDigest=digestIndexVector(
            sourceResult.publishedGlobalByAuthoritativeCellSlot);
        workingResult->boundaryAudit=sourceResult.boundaryAudit;
        artifact->t3CommittedInsertionCount=t3InsertionCount;

        // Consuming handoff.  Every destination object and fallible copy above
        // is complete; only nothrow ownership/container swaps follow the first
        // mutation of the accepted T3 result.
        receipt->sourceReceipt=std::move(source.sourceReceipt);
        receipt->sourceReport=std::move(source.report);
        artifact->authoritative=std::move(sourceResult.authoritative);
        artifact->constraintState.swap(sourceResult.constraintState);
        artifact->existingGlobalByAuthoritativeNodeSlot.swap(
            sourceResult.existingGlobalByAuthoritativeNodeSlot);
        artifact->immutableBoundary.swap(sourceResult.immutableBoundary);
        artifact->vertexKeyByHandle.swap(vertexIdentity);
        workingResult->artifact=std::move(artifact);
        out.sourceReceipt=std::move(receipt);
        out.sourceOutcome.reset();
        sourceMesh.reset();

        auto& result=*workingResult;
        auto& mutableArtifact=*result.artifact;
        if (!mutableArtifact.authoritative||
            !mutableArtifact.authoritative->workspace)
            throw std::logic_error("T4 handoff lost authoritative workspace");
        auto& ws=*mutableArtifact.authoritative->workspace;

        activePhase=T4APhase::ContextBuild;
        CoordinateMap coordinates;
        InsertedMap inserted;
        CoordinateIndex coordinateIndex;
        for (const auto node:ws.nodes().activeBySequence()) {
            const auto coord=ws.nodes().get(node).coord;
            coordinates.emplace(node,coord);
            if (!coordinateIndex.emplace(coordinateBits(coord),node).second)
                throw std::logic_error("duplicate authoritative coordinate");
            const auto identity=mutableArtifact.vertexKeyByHandle.find(node);
            if (identity==mutableArtifact.vertexKeyByHandle.end())
                throw std::logic_error("missing T4 vertex identity");
            if (identity->second.kind==RefinementVertexKind::InteriorInsertion)
                inserted.emplace(node,identity->second.value);
        }
        CellMap cells;
        EdgeStars edgeStars;
        FaceOwners faceOwners;
        NodeStars nodeStars;
        mesh::Classification region{};
        bool haveRegion=false;
        for (const auto cell:ws.cells().activeBySequence()) {
            const auto data=ws.cells().get(cell);
            if (data.type!=topo::ElementType::Tet4)
                throw std::logic_error("non-Tet4 T4 source");
            if (!haveRegion) { region=data.cls; haveRegion=true; }
            else if (!sameClass(region,data.cls)) {
                out.status=T4AStatus::Unsupported;
                out.failure={T4AReason::InvalidAuthoritativeSource,
                             T4APhase::ContextBuild};
                return out;
            }
            const std::array<topo::NodeHandle,4> nodes{
                data.nodes[0],data.nodes[1],data.nodes[2],data.nodes[3]};
            cells.emplace(cell,nodes);
            addCellToStars(edgeStars,cell,nodes);
            addCellToFaces(faceOwners,cell,nodes);
            addCellToNodeStars(nodeStars,cell,nodes);
        }
        if (!haveRegion||region.dim!=3||region.tag<0)
            throw std::logic_error("invalid T4 region");

        std::uint64_t nodeSlotCapacity=0,cellSlotCapacity=0;
        if (!checkedAdd(ws.nodes().stats().materializedSlots,
                        request.quality.maxNewNodes,nodeSlotCapacity)||
            !checkedAdd(ws.cells().stats().materializedSlots,
                        createdAllowance,cellSlotCapacity)||
            !fitsSize(nodeSlotCapacity)||!fitsSize(cellSlotCapacity))
            throw std::length_error("T4 mapping slot capacity overflow");
        mutableArtifact.publishedGlobalByAuthoritativeNodeSlot.resize(
            static_cast<std::size_t>(nodeSlotCapacity),kInvalidIndex);
        mutableArtifact.publishedGlobalByAuthoritativeCellSlot.resize(
            static_cast<std::size_t>(cellSlotCapacity),kInvalidIndex);

        T4AMetrics& metrics=result.metrics;
        TetRefinementMetrics coreMetrics;
        ViolationHeap queue;
        CurrentBlockedHeap blocked;
        T4BlockedDependencyIndex blockedDependencies;
        T3AReport coordinatorReport;
        RefinementLocationIndex locationIndex;
        SizingLedger emptyLedger;
        std::uint64_t queueSlotCapacity=0;
        if (!checkedAdd(ws.cells().stats().materializedSlots,
                        requiredTets,queueSlotCapacity)||
            !fitsSize(queueSlotCapacity)||!fitsSize(requiredTets))
            throw std::length_error("T4 queue capacity overflow");
        queue.reserve(static_cast<std::size_t>(requiredTets),
                      static_cast<std::size_t>(queueSlotCapacity));
        blocked.reserve(
            static_cast<std::size_t>(request.quality.maxSelectedGenerations),
            static_cast<std::size_t>(queueSlotCapacity));
        blockedDependencies.prepare(out.report->dependencyStorage.size());
        coordinatorReport.blockedStorage.resize(
            static_cast<std::size_t>(request.quality.maxSelectedGenerations));
        locationIndex.prepare(static_cast<std::size_t>(requiredTets),
                              static_cast<std::size_t>(queueSlotCapacity));

        activePhase=T4APhase::QueueBuild;
        std::uint64_t generation=0;
        const auto tailCapacity=(std::min)(
            request.quality.maxWorstTailRows,request.quality.maxResultTets);
        if (!fillT4Summary(
                cells,coordinates,
                mutableArtifact.existingGlobalByAuthoritativeNodeSlot,
                inserted,request.quality.minimumRadiusRatio,0,tailCapacity,
                *out.report,result.before,metrics,false,
                request.quality.mode==T4AQualityMode::BoundedImprove
                    ? &queue:nullptr,
                request.quality.mode==T4AQualityMode::BoundedImprove
                    ? &generation:nullptr,
                request.quality.mode==T4AQualityMode::BoundedImprove
                    ? &coreMetrics:nullptr))
            throw std::logic_error("invalid initial T4 quality audit");
        out.report->before=result.before;

        auto finishAudited=[&](bool incomplete,T4AReason incompleteReason,
                               T4APhase incompletePhase)->T4AOutcome {
            activePhase=T4APhase::FinalAudit;
            if (!topo::workspaceTopologyReady(
                    ws.nodes(),ws.cells(),ws.edges(),ws.faces())||
                !auditT4Identity(
                    ws,mutableArtifact.vertexKeyByHandle,
                    mutableArtifact.existingGlobalByAuthoritativeNodeSlot,
                    inserted,cells)||
                !finalProtectedAudit(
                    mutableArtifact.constraintState,ws,
                    mutableArtifact.sourceGlobalByAuthoritativeNodeSlot)||
                !fillT4Unattempted(*out.report,queue,cells)||
                !fillT4Summary(
                    cells,coordinates,
                    mutableArtifact.existingGlobalByAuthoritativeNodeSlot,
                    inserted,request.quality.minimumRadiusRatio,tailCapacity,
                    tailCapacity,*out.report,result.after,metrics,true,
                    nullptr,nullptr,&coreMetrics)) {
                out.status=T4AStatus::Failed;
                out.failure={T4AReason::QualityAuditFailed,
                             T4APhase::FinalAudit};
                return std::move(out);
            }
#if defined(PAACS_TOPO_TESTING)
            auto& finalLive=t4a_detail::finalLiveTraceForTesting();
            if (finalLive.enabled) {
                finalLive.used=0;
                finalLive.complete=true;
                for (const auto& queued:queue.rows()) {
                    const auto found=cells.find(queued.cell);
                    if (found==cells.end()||found->second!=queued.nodes)
                        continue;
                    if (finalLive.used>=finalLive.rows.size()) {
                        finalLive.complete=false;
                        break;
                    }
                    finalLive.rows[finalLive.used++]=t4QualityKey(queued.row);
                }
                std::sort(finalLive.rows.begin(),
                          finalLive.rows.begin()+finalLive.used);
            }
#endif
            out.report->after=result.after;
            metrics.qualityComplete=result.after.belowThreshold==0;
            syncT4CoreMetrics(coreMetrics,metrics);
            metrics.runComplete=true;
            if (request.quality.mode==T4AQualityMode::ReportOnly)
                metrics.qualityComplete=result.before.belowThreshold==0;

            activePhase=T4APhase::PublicationPreparation;
            result.publicationReady=false;
            out.retainedResult=std::move(workingResult);
            // MSVC's default unordered-map construction inside Mesh may
            // allocate. Transfer the audited owner before constructing any
            // publication object so every preparation allocation retains it.
            mesh::Mesh publication;
            TetRefinementQualityReport publicationQuality;
            if (!buildPublicationMesh(
                    *publicationBase,ws,
                    mutableArtifact.sourceGlobalByAuthoritativeNodeSlot,
                    inserted,publication,
                    mutableArtifact.publishedGlobalByAuthoritativeNodeSlot,
                    mutableArtifact.publishedGlobalByAuthoritativeCellSlot,
                    publicationQuality)) {
                out.status=T4AStatus::Failed;
                out.failure={T4AReason::PublicationPreparationFailed,
                             T4APhase::PublicationPreparation};
                out.retainedResult.reset();
                return std::move(out);
            }
            activePhase=T4APhase::FinalAudit;
            if (!auditBoundaryBijection(
                    publication,mutableArtifact.immutableBoundary,ws,
                    mutableArtifact.publishedGlobalByAuthoritativeNodeSlot)) {
                out.status=T4AStatus::Failed;
                out.failure={T4AReason::BoundaryBijectionFailed,
                             T4APhase::FinalAudit};
                out.retainedResult.reset();
                return std::move(out);
            }
            result.boundaryAudit=out.sourceReceipt->boundaryAudit;
            result.finalStrictHash=mesh::meshHash(publication);
            result.finalSemanticHash=mesh::semanticHashEx(publication);
            if (incomplete &&
                request.quality.mode==T4AQualityMode::BoundedImprove) {
                result.publicationReady=false;
                out.status=T4AStatus::Incomplete;
                out.failure={incompleteReason,incompletePhase};
                if (incompleteReason==T4AReason::ResidualQuality&&
                    !blocked.empty()) {
                    out.failure.hasCell=true;
                    out.failure.cell=t4QualityKey(blocked.top().row).tet;
                }
                return std::move(out);
            }
            result.publicationReady=true;
            activePhase=T4APhase::Publication;
            boundary_conforming_detail::publicationTailActive()=true;
            destination.swap(publication);
            boundary_conforming_detail::publicationTailActive()=false;
            out.status=T4AStatus::Published;
            out.destinationUnchanged=false;
            return std::move(out);
        };

        if (request.quality.mode==T4AQualityMode::ReportOnly)
            return finishAudited(false,T4AReason::None,T4APhase::FinalAudit);

        activePhase=T4APhase::ContextBuild;
        auto baseLocationBuild=bnd::buildTetRecoveryContext(ws,region);
        if (!baseLocationBuild.built()) {
            out.status=T4AStatus::Failed;
            out.failure={baseLocationBuild.reason==
                    bnd::TetRecoveryContextReason::AllocationFailed
                    ? T4AReason::AllocationFailed:
                      T4AReason::InvalidAuthoritativeSource,
                T4APhase::ContextBuild};
            return out;
        }
        auto baseLocation=std::move(baseLocationBuild.context);
        ProtectedFeatureIndex protectedIndex;
        const auto protectedBuild=buildProtectedIndex(
            mutableArtifact.constraintState,ws,protectedIndex,coreMetrics);
        if (protectedBuild!=ProtectedIndexBuild::Built) {
            out.status=protectedBuild==ProtectedIndexBuild::Range
                ? T4AStatus::Unsupported:T4AStatus::Failed;
            out.failure={T4AReason::InvalidAuthoritativeSource,
                         T4APhase::ContextBuild};
            return out;
        }
        ProtectedQuery protectedHit;
        protectedHit.witnesses.reserve(protectedIndex.rows.size());
        std::vector<std::size_t> protectedStack;
        protectedStack.reserve(protectedIndex.nodes.size());
        std::vector<T4ADependencyKey> dependencyScratch;
        dependencyScratch.reserve(
            static_cast<std::size_t>(
                request.quality.maxDependencyRowsPerAttempt)+1);

        std::uint64_t t4CurrentContextAccounted=
            completeContextAccountedBytes(
                coordinates,inserted,coordinateIndex,cells,edgeStars,
                faceOwners,nodeStars,queue,blocked,emptyLedger,locationIndex);
        std::uint64_t t4FixedAccounted=0;
        addAccountedBytes(t4FixedAccounted,workspaceAccountedBytes(ws));
        addAccountedBytes(t4FixedAccounted,
                          flatMeshAccountedBytes(*publicationBase));
        addAccountedBytes(t4FixedAccounted,reportAccountedBytes(*out.report));
        addAccountedBytes(t4FixedAccounted,t4RetainedAccountedBytes(
            mutableArtifact,*out.sourceReceipt));
        addAccountedBytes(t4FixedAccounted,
                          recoveryContextAccountedBytes(*baseLocation));
        addAccountedBytes(t4FixedAccounted,
                          vectorAccountedBytes(protectedIndex.rows));
        addAccountedBytes(t4FixedAccounted,
                          vectorAccountedBytes(protectedIndex.nodes));
        addAccountedBytes(t4FixedAccounted,
                          vectorAccountedBytes(protectedHit.witnesses));
        addAccountedBytes(t4FixedAccounted,
                          vectorAccountedBytes(protectedStack));
        addAccountedBytes(t4FixedAccounted,
                          blockedDependencies.accountedBytes());
        addAccountedBytes(t4FixedAccounted,
                          vectorAccountedBytes(dependencyScratch));
        coreMetrics.transientAccountedHighWaterBytes=t4FixedAccounted;
        addAccountedBytes(coreMetrics.transientAccountedHighWaterBytes,
                          t4CurrentContextAccounted);
        WorkBudgetScope workBudget(request.quality.maxWorkUnits);
        bool incomplete=false;
        T4AReason incompleteReason=T4AReason::None;
        T4APhase incompletePhase=T4APhase::CandidateGeneration;

        const auto setTerminal=[&](const QueuedViolation& selected,
                                   T4AAttemptRow& attempt,T4AReason reason,
                                   T4APhase phase) noexcept {
            attempt.disposition=T4AAttemptDisposition::TerminalFailure;
            attempt.reason=reason;
            attempt.workAfter=coreMetrics.workUnitsConsumed;
            out.report->hasTerminal=true;
            out.report->terminal={t4QualityKey(selected.row),reason,phase,
                                  attempt.attemptOrdinal};
        };

        const auto blockSelected=[&](const QueuedViolation& selected,
                                     T4AAttemptRow& attempt,T4AReason reason,
                                     std::span<const T4ADependencyKey>
                                         dependencies,
                                     std::span<const T3AProtectedFeatureWitness>
                                         witnesses)->bool {
            std::uint64_t dependencyEnd=0,witnessEnd=0;
            if (dependencies.empty()||
                dependencies.size()>
                    request.quality.maxDependencyRowsPerAttempt||
                witnesses.size()>
                    request.quality.maxProtectedWitnessesPerAttempt||
                !checkedAdd(out.report->dependenciesUsed,
                            dependencies.size(),dependencyEnd)||
                !checkedAdd(out.report->witnessesUsed,witnesses.size(),
                            witnessEnd)||
                dependencyEnd>out.report->dependencyStorage.size()||
                witnessEnd>out.report->witnessStorage.size()||
                out.report->blockedUsed>=out.report->blockedStorage.size()||
                coordinatorReport.blockedUsed>=
                    coordinatorReport.blockedStorage.size())
                return false;
            const auto rowIndex=out.report->blockedUsed;
            const CurrentBlockedEntry entry{
                selected.cell,selected.row,rowIndex};
            if (!blocked.canInsert(entry)||
                !blockedDependencies.canInsert(
                    dependencies,selected.cell,rowIndex)) return false;
            const auto dependencyOffset=out.report->dependenciesUsed;
            for (std::size_t i=0;i<dependencies.size();++i)
                out.report->dependencyStorage[
                    static_cast<std::size_t>(dependencyOffset)+i]=
                    dependencies[i];
            const auto witnessOffset=out.report->witnessesUsed;
            for (std::size_t i=0;i<witnesses.size();++i)
                out.report->witnessStorage[
                    static_cast<std::size_t>(witnessOffset)+i]=witnesses[i];
            auto& row=out.report->blockedStorage[
                static_cast<std::size_t>(rowIndex)];
            row.selected=t4QualityKey(selected.row);
            row.reason=reason;
            row.attemptOrdinal=attempt.attemptOrdinal;
            row.dependencyOffset=dependencyOffset;
            row.dependencyCount=dependencies.size();
            row.witnessOffset=witnessOffset;
            row.witnessCount=witnesses.size();
            row.invalidatedByOperation=kCurrentT4Operation;
            auto& coreRow=coordinatorReport.blockedStorage[
                static_cast<std::size_t>(rowIndex)];
            coreRow.violation=selected.row.key;
            coreRow.violationGeneration=selected.row.generation;
            coreRow.attemptOrdinal=attempt.attemptOrdinal;
            coreRow.invalidatedByInsertionOrdinal=kCurrentBlockedInsertion;
            attempt.dependencyOffset=dependencyOffset;
            attempt.dependencyCount=dependencies.size();
            attempt.witnessOffset=witnessOffset;
            attempt.witnessCount=witnesses.size();
            if (!witnesses.empty())
                attempt.disposition=T4AAttemptDisposition::Encroached;
            else if (attempt.disposition==
                     T4AAttemptDisposition::TerminalFailure)
                attempt.disposition=T4AAttemptDisposition::Blocked;
            attempt.reason=reason;
            attempt.workAfter=coreMetrics.workUnitsConsumed;
            blocked.insertPrepared(entry);
            blockedDependencies.insertPrepared(
                dependencies,selected.cell,rowIndex);
            out.report->dependenciesUsed=dependencyEnd;
            out.report->witnessesUsed=witnessEnd;
            ++out.report->blockedUsed;
            ++coordinatorReport.blockedUsed;
            ++metrics.blockedRows;
            ++coreMetrics.blockedRowsRecorded;
            ++coreMetrics.blockedRowsCurrent;
            return true;
        };

        while (!incomplete&&!queue.empty()) {
#if defined(PAACS_TOPO_TESTING)
            activePhase=T4APhase::CandidateGeneration;
            try {
                recordT4StateDigest(
                    t4a_detail::T4ADigestCheckpoint::PreSelection,false,
                    ws,coordinates,inserted,coordinateIndex,cells,edgeStars,
                    faceOwners,nodeStars,queue,blocked,blockedDependencies,
                    locationIndex,mutableArtifact,*out.report,metrics);
            } catch (const std::bad_alloc&) {
                auto& trace=t4a_detail::digestTraceForTesting();
                trace.faultPending=false;
                recordT4StateDigest(
                    t4a_detail::T4ADigestCheckpoint::PreSelection,true,
                    ws,coordinates,inserted,coordinateIndex,cells,edgeStars,
                    faceOwners,nodeStars,queue,blocked,blockedDependencies,
                    locationIndex,mutableArtifact,*out.report,metrics);
                out.status=T4AStatus::Failed;
                out.failure={T4AReason::AllocationFailed,activePhase};
                return out;
            }
#endif
            const auto& root=queue.top();
            const auto live=cells.find(root.cell);
            if (live==cells.end()||live->second!=root.nodes) {
                (void)queue.pop();
                ++metrics.queuePops;
                ++metrics.staleRows;
                continue;
            }
            if (out.report->attemptsUsed>=out.report->attemptStorage.size()) {
                incomplete=true;
                incompleteReason=T4AReason::ReportBudgetReached;
                break;
            }
            if (coreMetrics.workUnitsConsumed>=request.quality.maxWorkUnits) {
                incomplete=true;
                incompleteReason=T4AReason::WorkBudgetReached;
                break;
            }
            if (metrics.selectedGenerations>=
                    request.quality.maxSelectedGenerations) {
                incomplete=true;
                incompleteReason=T4AReason::SelectedGenerationBudgetReached;
                break;
            }
            if (metrics.committedOperations>=request.quality.maxOperations) {
                incomplete=true;
                incompleteReason=T4AReason::OperationBudgetReached;
                break;
            }
            if (metrics.committedOperations>=request.quality.maxNewNodes) {
                incomplete=true;
                incompleteReason=T4AReason::NodeBudgetReached;
                break;
            }
            QueuedViolation selected=queue.pop();
            ++metrics.queuePops;
            ++metrics.selectedGenerations;
            auto& attempt=out.report->attemptStorage[
                static_cast<std::size_t>(out.report->attemptsUsed++)];
            attempt.attemptOrdinal=metrics.candidateAttempts++;
            attempt.operationOrdinal=metrics.committedOperations;
            attempt.selected=t4QualityKey(selected.row);
            attempt.selectedQuality=t4QualityValue(selected.row);
            attempt.workBefore=coreMetrics.workUnitsConsumed;
            activePhase=T4APhase::CandidateGeneration;

            try {
#if defined(PAACS_TOPO_TESTING)
                recordT4StateDigest(
                    t4a_detail::T4ADigestCheckpoint::SelectedRowSecured,false,
                    ws,coordinates,inserted,coordinateIndex,cells,edgeStars,
                    faceOwners,nodeStars,queue,blocked,blockedDependencies,
                    locationIndex,mutableArtifact,*out.report,metrics);
#endif
            dependencyScratch.clear();
            bool dependencyOverflow=false;
            const auto addDependency=[&](T4ADependencyKey value) {
                if (dependencyScratch.size()>=
                        request.quality.maxDependencyRowsPerAttempt) {
                    dependencyOverflow=true;
                    return;
                }
                dependencyScratch.push_back(value);
            };
            addDependency(t4CellDependency(selected.row));
            const auto addCellDependency=[&](topo::CellHandle cell) {
                const auto found=cells.find(cell);
                if (found!=cells.end())
                    addDependency(t4Dependency(
                        T4ADependencyKind::Cell,found->second,
                        mutableArtifact.existingGlobalByAuthoritativeNodeSlot,
                        inserted));
            };
            const auto addBoundaryDependencies=[&](
                    const std::vector<BoundaryFace>& boundary) {
                for (const auto& face:boundary) {
                    addDependency(t4Dependency(
                        T4ADependencyKind::Face,face.directed,
                        mutableArtifact.existingGlobalByAuthoritativeNodeSlot,
                        inserted));
                    for (std::size_t i=0;i<3;++i) {
                        const std::array<topo::NodeHandle,2> edge{
                            face.directed[i],face.directed[(i+1)%3]};
                        addDependency(t4Dependency(
                            T4ADependencyKind::Edge,edge,
                            mutableArtifact.
                                existingGlobalByAuthoritativeNodeSlot,
                            inserted));
                    }
                }
            };
            const auto blockCurrent=[&](
                    T4AReason reason,
                    std::span<const T3AProtectedFeatureWitness> witnesses) {
                canonicalizeT4Dependencies(dependencyScratch);
                return !dependencyOverflow&&blockSelected(
                    selected,attempt,reason,dependencyScratch,witnesses);
            };

                recordConstructionWork(coreMetrics);
                const double required=attempt.selectedQuality+
                    request.quality.minimumImprovement;
                if (!std::isfinite(required)||required>1.0) {
                    ++metrics.formulaRejects;
                    attempt.disposition=
                        T4AAttemptDisposition::NoQualityProgress;
                    if (!blockCurrent(T4AReason::NoImprovingCavity,{})) {
                        setTerminal(selected,attempt,
                                    T4AReason::ReportBudgetReached,activePhase);
                        incomplete=true;
                        incompleteReason=T4AReason::ReportBudgetReached;
                        break;
                    }
                    continue;
                }
                std::array<Vec3,4> selectedPoints{};
                for (std::size_t i=0;i<4;++i)
                    selectedPoints[i]=coordinates.at(selected.nodes[i]);
                auto candidate=tetCircumcenter(selectedPoints);
#if defined(PAACS_TOPO_TESTING)
                if (t4a_detail::forceCandidateForTesting()) {
                    candidate=t4a_detail::forcedCandidateForTesting();
                    t4a_detail::forceCandidateForTesting()=false;
                }
#endif
                if (!candidate) {
                    ++metrics.formulaRejects;
                    attempt.disposition=T4AAttemptDisposition::FormulaUnavailable;
                    if (!blockCurrent(T4AReason::NoImprovingCavity,{})) {
                        setTerminal(selected,attempt,
                                    T4AReason::ReportBudgetReached,activePhase);
                        incomplete=true;
                        incompleteReason=T4AReason::ReportBudgetReached;
                        break;
                    }
                    continue;
                }
                attempt.proposedCoordinate=*candidate;
                if (coordinateIndex.count(coordinateBits(*candidate))) {
                    ++metrics.duplicateRejects;
                    attempt.disposition=T4AAttemptDisposition::DuplicateProposal;
                    if (!blockCurrent(T4AReason::NoImprovingCavity,{})) {
                        setTerminal(selected,attempt,
                                    T4AReason::ReportBudgetReached,activePhase);
                        incomplete=true;
                        incompleteReason=T4AReason::ReportBudgetReached;
                        break;
                    }
                    continue;
                }

                activePhase=T4APhase::Encroachment;
                protectedEncroachment(*candidate,protectedIndex,ws,coreMetrics,
                                      protectedHit,protectedStack);
                if (protectedHit.invalid) {
                    setTerminal(selected,attempt,T4AReason::InvalidCandidate,
                                activePhase);
                    out.status=T4AStatus::Failed;
                    out.failure={T4AReason::InvalidCandidate,activePhase,true,
                                 selected.row.key.tet};
                    return out;
                }
                if (!protectedHit.witnesses.empty()) {
                    for (const auto& witness:protectedHit.witnesses) {
                        if (witness.kind==T3AProtectedFeatureKind::Edge) {
                            const auto binding=mutableArtifact.constraintState.
                                edgeBindings.find(witness.edge);
                            if (binding==mutableArtifact.constraintState.
                                    edgeBindings.end()||
                                !ws.edges().isCurrent(binding->second.edge)) {
                                dependencyOverflow=true;
                                break;
                            }
                            const auto edge=ws.edges().get(binding->second.edge);
                            addDependency(t4Dependency(
                                T4ADependencyKind::Edge,
                                std::array<topo::NodeHandle,2>{edge.a,edge.b},
                                mutableArtifact.
                                    existingGlobalByAuthoritativeNodeSlot,
                                inserted));
                        } else {
                            const auto binding=mutableArtifact.constraintState.
                                bindings.find(witness.facet);
                            if (binding==mutableArtifact.constraintState.
                                    bindings.end()) {
                                dependencyOverflow=true;
                                break;
                            }
                            addDependency(t4Dependency(
                                T4ADependencyKind::Face,
                                binding->second.canonicalVertices,
                                mutableArtifact.
                                    existingGlobalByAuthoritativeNodeSlot,
                                inserted));
                        }
                    }
                    if (!blockCurrent(T4AReason::ProtectedFeatureBlocked,
                            protectedHit.witnesses)) {
                        setTerminal(selected,attempt,
                                    T4AReason::ReportBudgetReached,activePhase);
                        incomplete=true;
                        incompleteReason=T4AReason::ReportBudgetReached;
                        incompletePhase=activePhase;
                        break;
                    }
                    continue;
                }

                activePhase=T4APhase::CavityDiscovery;
                const auto located=locateExact(
                    *candidate,baseLocation->index,locationIndex,ws,
                    coordinates,coreMetrics);
                if (!located) {
                    attempt.disposition=
                        T4AAttemptDisposition::NoQualityProgress;
                    if (!blockCurrent(T4AReason::NoImprovingCavity,{})) {
                        setTerminal(selected,attempt,
                                    T4AReason::ReportBudgetReached,activePhase);
                        incomplete=true;
                        incompleteReason=T4AReason::ReportBudgetReached;
                        break;
                    }
                    continue;
                }
                for (const auto cell:located->star) addCellDependency(cell);
                const auto discovered=discoverConflictCavity(
                    located->seed,*candidate,cells,faceOwners,coordinates,
                    coreMetrics);
                if (!discovered ||
                    std::find(discovered->begin(),discovered->end(),
                              selected.cell)==discovered->end()) {
                    attempt.disposition=
                        T4AAttemptDisposition::NoQualityProgress;
                    if (!blockCurrent(T4AReason::NoImprovingCavity,{})) {
                        setTerminal(selected,attempt,
                                    T4AReason::ReportBudgetReached,activePhase);
                        incomplete=true;
                        incompleteReason=T4AReason::ReportBudgetReached;
                        break;
                    }
                    continue;
                }
                const auto& removed=*discovered;
                for (const auto cell:removed) addCellDependency(cell);
                auto boundary=cavityBoundary(removed,cells);
                if (boundary) addBoundaryDependencies(*boundary);
                if (!boundary||!pointStrictlyInsideBoundary(
                        *boundary,*candidate,ws,coreMetrics)) {
                    attempt.disposition=
                        T4AAttemptDisposition::NoQualityProgress;
                    if (!blockCurrent(T4AReason::NoImprovingCavity,{})) {
                        setTerminal(selected,attempt,
                                    T4AReason::ReportBudgetReached,activePhase);
                        incomplete=true;
                        incompleteReason=T4AReason::ReportBudgetReached;
                        break;
                    }
                    continue;
                }
                attempt.removedCells=removed.size();
                attempt.createdCells=boundary->size();
                coreMetrics.cavityFaceVisits+=boundary->size();

                activePhase=T4APhase::CandidateValidation;
                if (!validateRefinementCandidate(
                        removed,*boundary,*candidate,cells,edgeStars,nodeStars,
                        coordinates,baseLocation->index,locationIndex,ws,
                        region,coreMetrics)) {
                    attempt.disposition=T4AAttemptDisposition::InvalidCandidate;
                    setTerminal(selected,attempt,T4AReason::InvalidCandidate,
                                activePhase);
                    out.status=T4AStatus::Failed;
                    out.failure={T4AReason::InvalidCandidate,activePhase,true,
                                 selected.row.key.tet};
                    return out;
                }

                std::uint64_t persistentOrdinal=0;
                if (!checkedAdd(mutableArtifact.t3CommittedInsertionCount,
                                metrics.committedOperations,persistentOrdinal)) {
                    setTerminal(selected,attempt,T4AReason::InternalError,
                                activePhase);
                    out.status=T4AStatus::Failed;
                    out.failure={T4AReason::InternalError,activePhase,true,
                                 selected.row.key.tet};
                    return out;
                }
                const RefinementVertexKey pendingKey{
                    RefinementVertexKind::InteriorInsertion,persistentOrdinal};
                // The consuming handoff proves T3 insertion identities are
                // exactly [0,t3CommittedInsertionCount).  Every successful
                // T4 commit appends the checked operation ordinal, so this
                // key is new without a global identity-map scan.

                topo::TetRefinementCavitySpec cavity;
                cavity.stamp=topo::WorkspaceStamp::of(ws);
                cavity.stagedNode={0,*candidate,region};
                cavity.removed=removed;
                double createdMinimum=(std::numeric_limits<double>::infinity)();
                bool progress=true;
                for (const auto& face:*boundary) {
                    topo::TetRefinementCellSpec cell;
                    std::array<Vec3,4> p{};
                    std::array<RefinementVertexKey,4> keys{};
                    for (std::size_t i=0;i<3;++i) {
                        cell.vertices[i]=
                            topo::TetRefinementNodeRef::existingNode(
                                face.directed[i]);
                        p[i]=coordinates.at(face.directed[i]);
                        keys[i]=vertexKey(
                            face.directed[i],
                            mutableArtifact.
                                existingGlobalByAuthoritativeNodeSlot,
                            inserted);
                    }
                    cell.vertices[3]=
                        topo::TetRefinementNodeRef::stagedNode();
                    cell.classification=region;
                    p[3]=*candidate;
                    keys[3]=pendingKey;
                    const auto createdSign=predicates::orient3d(
                        p[0],p[1],p[2],p[3]);
                    ++coreMetrics.orient3dTests;
                    recordConstructionWork(coreMetrics);
                    if (createdSign==predicates::Sign::Zero) {
                        progress=false;
                        break;
                    }
                    if (createdSign==predicates::Sign::Negative) {
                        std::swap(cell.vertices[0],cell.vertices[1]);
                        std::swap(p[0],p[1]);
                        std::swap(keys[0],keys[1]);
                    }
                    double quality=0.0;
                    ++metrics.qualityEvaluations;
                    if (!checkedTetQualityCanonical(p,keys,quality)) {
                        progress=false;
                        break;
                    }
                    createdMinimum=(std::min)(createdMinimum,quality);
                    if (!(quality>=required)) progress=false;
                    cavity.added.push_back(cell);
                }
                attempt.createdMinimumQuality=createdMinimum;
                if (!progress||cavity.added.size()!=boundary->size()) {
                    attempt.disposition=
                        T4AAttemptDisposition::NoQualityProgress;
                    if (!blockCurrent(T4AReason::NoImprovingCavity,{})) {
                        setTerminal(selected,attempt,
                                    T4AReason::ReportBudgetReached,activePhase);
                        incomplete=true;
                        incompleteReason=T4AReason::ReportBudgetReached;
                        break;
                    }
                    continue;
                }

                activePhase=T4APhase::TransactionPreparation;
                if (cavity.added.size()>
                        request.quality.maxCreatedTetsPerOperation) {
                    setTerminal(selected,attempt,
                                T4AReason::ResultBudgetReached,activePhase);
                    incomplete=true;
                    incompleteReason=T4AReason::ResultBudgetReached;
                    break;
                }
                const auto topologyCountResult=
                    t4TopologyTransitionCount(cavity,cells);
                std::uint64_t topologyEnd=0;
                if (!topologyCountResult) {
                    setTerminal(selected,attempt,T4AReason::InternalError,
                                activePhase);
                    out.status=T4AStatus::Failed;
                    out.failure={T4AReason::InternalError,activePhase,true,
                                 selected.row.key.tet};
                    return out;
                }
                const std::uint64_t topologyCount=*topologyCountResult;
                if (topologyCount>
                        request.quality.maxTopologyRowsPerOperation||
                    !checkedAdd(out.report->topologyUsed,topologyCount,
                                topologyEnd)||
                    topologyEnd>out.report->topologyStorage.size()) {
                    setTerminal(selected,attempt,
                                T4AReason::ReportBudgetReached,activePhase);
                    incomplete=true;
                    incompleteReason=T4AReason::ReportBudgetReached;
                    break;
                }
#if defined(PAACS_TOPO_TESTING)
                recordT4StateDigest(
                    t4a_detail::T4ADigestCheckpoint::PrePreparation,false,
                    ws,coordinates,inserted,coordinateIndex,cells,edgeStars,
                    faceOwners,nodeStars,queue,blocked,blockedDependencies,
                    locationIndex,mutableArtifact,*out.report,metrics);
#endif
                auto tx=ws.beginTetRefinementTransaction();
                bool topologyPrepared=false;
                try {
                    activePhase=T4APhase::TransactionPreparation;
                    (void)tx.stage(cavity);
                    const auto& prediction=tx.prepare();
                    topologyPrepared=true;
                    ++coreMetrics.preparedTransactions;
                    const std::uint64_t preparedTopologyCount=1+
                        prediction.cells.size()+prediction.edges.size()+
                        prediction.faces.size();
                    if (preparedTopologyCount!=topologyCount)
                        throw std::logic_error(
                            "T4 local topology transition count mismatch");

                    const PendingInsertionIdentity pending{
                        prediction.stagedNode,persistentOrdinal};
                    CreatedCellRows createdRows;
                    createdRows.reserve(boundary->size());
                    for (const auto& transition:prediction.cells)
                        if (transition.transition==
                                topo::TetRefinementTransition::Created) {
                            std::array<topo::NodeHandle,4> conn{};
                            for (std::size_t i=0;i<4;++i)
                                conn[i]=transition.key[i].kind==
                                        topo::TetRefinementNodeRefKind::Staged
                                    ? prediction.stagedNode:
                                      transition.key[i].existing;
                            std::array<Vec3,4> predictedCoordinates{};
                            for (std::size_t i=0;i<4;++i)
                                predictedCoordinates[i]=
                                    conn[i]==prediction.stagedNode
                                        ? *candidate:coordinates.at(conn[i]);
                            const auto predictedSign=predicates::orient3d(
                                predictedCoordinates[0],predictedCoordinates[1],
                                predictedCoordinates[2],predictedCoordinates[3]);
                            ++coreMetrics.orient3dTests;
                            recordConstructionWork(coreMetrics);
                            if (predictedSign==predicates::Sign::Zero)
                                throw std::logic_error(
                                    "degenerate predicted T4 cell");
                            if (predictedSign==predicates::Sign::Negative)
                                std::swap(conn[0],conn[1]);
                            createdRows.push_back({transition.predicted,conn});
                        }
                    if (createdRows.size()!=cavity.added.size())
                        throw std::logic_error(
                            "T4 prepared created-cell count mismatch");
                    std::vector<QueuedViolation> createdViolations;
                    createdViolations.reserve(createdRows.size());
                    double verifiedMinimum=
                        (std::numeric_limits<double>::infinity)();
                    for (const auto& [createdCell,conn]:createdRows) {
                        std::array<Vec3,4> p{};
                        for (std::size_t i=0;i<4;++i)
                            p[i]=conn[i]==prediction.stagedNode
                                ? *candidate:coordinates.at(conn[i]);
                        std::array<RefinementVertexKey,4> keys{};
                        for (std::size_t i=0;i<4;++i)
                            keys[i]=conn[i]==prediction.stagedNode
                                ? pendingKey:vertexKey(
                                    conn[i],mutableArtifact.
                                        existingGlobalByAuthoritativeNodeSlot,
                                    inserted);
                        double quality=0.0;
                        ++metrics.qualityEvaluations;
                        if (!checkedTetQualityCanonical(p,keys,quality)||
                            !(quality>=required))
                            throw std::logic_error(
                                "prepared T4 candidate lost progress");
                        verifiedMinimum=(std::min)(verifiedMinimum,quality);
                        auto violation=evaluateT4QualityRow(
                            conn,p,
                            mutableArtifact.
                                existingGlobalByAuthoritativeNodeSlot,
                            inserted,request.quality.minimumRadiusRatio,
                            generation++,metrics,coreMetrics,&pending);
                        if (violation)
                            createdViolations.push_back(
                                {*violation,createdCell,conn});
                    }
                    attempt.createdMinimumQuality=verifiedMinimum;

                    activePhase=T4APhase::CoordinatorPreparation;
                    auto prepared=std::make_unique<
                        PreparedTetRefinementCoordinatorDelta>();
                    prepareCoordinatorDelta(
                        *prepared,prediction.stagedNode,*candidate,
                        persistentOrdinal,removed,createdRows,
                        createdViolations,coordinates,cells,edgeStars,
                        faceOwners,nodeStars,queue,blocked,coordinatorReport,
                        locationIndex,coreMetrics);
                    std::vector<T4ADependencyKey> affectedDependencies;
                    affectedDependencies.reserve(
                        removed.size()+prepared->affectedEdges.size()+
                        prepared->affectedFaces.size());
                    for (const auto cell:removed) {
                        const auto found=cells.find(cell);
                        if (found==cells.end())
                            throw std::logic_error(
                                "missing removed T4 dependency cell");
                        affectedDependencies.push_back(t4Dependency(
                            T4ADependencyKind::Cell,found->second,
                            mutableArtifact.
                                existingGlobalByAuthoritativeNodeSlot,
                            inserted));
                    }
                    for (const auto& edge:prepared->affectedEdges)
                        affectedDependencies.push_back(t4Dependency(
                            T4ADependencyKind::Edge,edge,
                            mutableArtifact.
                                existingGlobalByAuthoritativeNodeSlot,
                            inserted,&pending));
                    for (const auto& face:prepared->affectedFaces)
                        affectedDependencies.push_back(t4Dependency(
                            T4ADependencyKind::Face,face,
                            mutableArtifact.
                                existingGlobalByAuthoritativeNodeSlot,
                            inserted,&pending));
                    canonicalizeT4Dependencies(affectedDependencies);
                    prepared->blockedInvalidations.reserve(
                        prepared->blockedInvalidations.size()+blocked.size());
                    for (const auto& dependency:affectedDependencies) {
                        recordConstructionWork(coreMetrics);
                        blockedDependencies.forEach(dependency,
                            [&](const T4BlockedDependencyIndex::Entry& entry) {
                                const bool already=std::any_of(
                                    prepared->blockedInvalidations.begin(),
                                    prepared->blockedInvalidations.end(),
                                    [&](const auto& value) {
                                        return value.cell==entry.selected;
                                    });
                                if (already) return;
                                const auto* current=blocked.find(entry.selected);
                                if (!current||current->blockedRowIndex!=
                                        entry.blockedRowIndex)
                                    throw std::logic_error(
                                        "stale T4 dependency blocker");
                                prepared->blockedInvalidations.push_back({
                                    entry.selected,entry.blockedRowIndex,
                                    persistentOrdinal});
                                ++coreMetrics.coordinatorRowsPrepared;
                                recordConstructionWork(coreMetrics);
                            });
                    }
                    std::sort(prepared->blockedInvalidations.begin(),
                              prepared->blockedInvalidations.end(),
                        [](const auto& a,const auto& b) {
                            return a.cell<b.cell;
                        });
                    prepared->topologyRows.reserve(
                        static_cast<std::size_t>(topologyCount));
                    T3ATopologyHandleRow nodeTopology;
                    nodeTopology.kind=T3ATopologyKind::Node;
                    nodeTopology.transition=T3ATopologyTransition::Created;
                    nodeTopology.predictedNode=prediction.stagedNode;
                    prepared->topologyRows.push_back(nodeTopology);
                    for (const auto& row:prediction.cells) {
                        T3ATopologyHandleRow value;
                        value.kind=T3ATopologyKind::Cell;
                        value.transition=publicTransition(row.transition);
                        value.predictedCell=row.predicted;
                        prepared->topologyRows.push_back(value);
                    }
                    for (const auto& row:prediction.edges) {
                        T3ATopologyHandleRow value;
                        value.kind=T3ATopologyKind::Edge;
                        value.transition=publicTransition(row.transition);
                        value.predictedEdge=row.predicted;
                        prepared->topologyRows.push_back(value);
                    }
                    for (const auto& row:prediction.faces) {
                        T3ATopologyHandleRow value;
                        value.kind=T3ATopologyKind::Face;
                        value.transition=publicTransition(row.transition);
                        value.predictedFace=row.predicted;
                        prepared->topologyRows.push_back(value);
                    }
                    if (prepared->topologyRows.size()!=topologyCount)
                        throw std::logic_error("T4 topology report mismatch");
                    std::map<topo::NodeHandle,RefinementVertexKey>
                        preparedVertexIdentity;
                    if (!preparedVertexIdentity.emplace(
                            prediction.stagedNode,pendingKey).second)
                        throw std::logic_error("duplicate prepared T4 key");
                    const auto topologyOffset=out.report->topologyUsed;
                    const auto nextAccounted=applyAccountedContextDelta(
                        t4CurrentContextAccounted,*prepared);
                    std::uint64_t preparedAccounted=
                        preparedCoordinatorAccountedBytes(*prepared);
                    addAccountedBytes(preparedAccounted,
                                      vectorAccountedBytes(createdRows));
                    addAccountedBytes(preparedAccounted,
                                      vectorAccountedBytes(createdViolations));
                    addAccountedBytes(preparedAccounted,
                                      vectorAccountedBytes(removed));
                    addAccountedBytes(preparedAccounted,
                                      vectorAccountedBytes(*boundary));
                    addAccountedBytes(preparedAccounted,
                                      vectorAccountedBytes(cavity.removed));
                    addAccountedBytes(preparedAccounted,
                                      vectorAccountedBytes(cavity.added));
                    addAccountedBytes(preparedAccounted,
                                      vectorAccountedBytes(prediction.cells));
                    addAccountedBytes(preparedAccounted,
                                      vectorAccountedBytes(prediction.edges));
                    addAccountedBytes(preparedAccounted,
                                      vectorAccountedBytes(prediction.faces));
                    addAccountedBytes(preparedAccounted,
                                      vectorAccountedBytes(
                                          affectedDependencies));
                    addAccountedBytes(preparedAccounted,
                                      orderedMapAccountedBytes(
                                          preparedVertexIdentity));
                    std::uint64_t insertionTransient=t4FixedAccounted;
                    addAccountedBytes(insertionTransient,
                                      t4CurrentContextAccounted);
                    addAccountedBytes(insertionTransient,preparedAccounted);
                    coreMetrics.transientAccountedHighWaterBytes=(std::max)(
                        coreMetrics.transientAccountedHighWaterBytes,
                        insertionTransient);

                    std::uint64_t incidentalRetirements=0;
                    for (const auto cell:removed)
                        if (cell!=selected.cell&&queue.find(cell))
                            ++incidentalRetirements;

#if defined(PAACS_TOPO_TESTING)
                    std::uint64_t preparedPrediction=
                        0xcbf29ce484222325ull;
                    digestHandle(preparedPrediction,prediction.stagedNode);
                    for (const auto& row:prediction.cells)
                        digestHandle(preparedPrediction,row.predicted);
                    for (const auto& row:prediction.edges)
                        digestHandle(preparedPrediction,row.predicted);
                    for (const auto& row:prediction.faces)
                        digestHandle(preparedPrediction,row.predicted);
                    recordT4StateDigest(
                        t4a_detail::T4ADigestCheckpoint::Prepared,false,
                        ws,coordinates,inserted,coordinateIndex,cells,
                        edgeStars,faceOwners,nodeStars,queue,blocked,
                        blockedDependencies,locationIndex,mutableArtifact,
                        *out.report,metrics,preparedPrediction);
#endif
                    const auto commit=tx.commit();
#if defined(PAACS_TOPO_TESTING)
                    Phase3TailGuard phase3TailGuard;
#endif
                    applyCoordinatorDelta(
                        *prepared,coordinates,inserted,coordinateIndex,cells,
                        edgeStars,faceOwners,nodeStars,queue,blocked,
                        coordinatorReport,locationIndex,coreMetrics);
                    t4CurrentContextAccounted=nextAccounted;
                    for (const auto& invalidation:
                            prepared->blockedInvalidations) {
                        if (invalidation.blockedRowIndex>=
                            out.report->blockedUsed) std::abort();
                        auto& row=out.report->blockedStorage[
                            static_cast<std::size_t>(
                                invalidation.blockedRowIndex)];
                        if (row.invalidatedByOperation!=kCurrentT4Operation)
                            std::abort();
                        const auto dependencyBegin=
                            out.report->dependencyStorage.data()+
                            static_cast<std::ptrdiff_t>(row.dependencyOffset);
                        blockedDependencies.erasePrepared(
                            std::span<const T4ADependencyKey>(
                                dependencyBegin,
                                static_cast<std::size_t>(row.dependencyCount)),
                            invalidation.cell,
                            invalidation.blockedRowIndex);
                        row.invalidatedByOperation=metrics.committedOperations;
                    }
                    transferPreparedMapNodes(
                        preparedVertexIdentity,
                        mutableArtifact.vertexKeyByHandle);
                    mutableArtifact.constraintState.workspaceStamp=commit.stamp;
                    for (std::size_t i=0;i<prepared->topologyRows.size();++i)
                        out.report->topologyStorage[
                            static_cast<std::size_t>(topologyOffset)+i]=
                            prepared->topologyRows[i];
                    std::size_t actual=static_cast<std::size_t>(topologyOffset);
                    out.report->topologyStorage[actual++].actualNode=
                        commit.actualNode;
                    for (const auto& row:commit.cells)
                        out.report->topologyStorage[actual++].actualCell=row.actual;
                    for (const auto& row:commit.edges)
                        out.report->topologyStorage[actual++].actualEdge=row.actual;
                    for (const auto& row:commit.faces)
                        out.report->topologyStorage[actual++].actualFace=row.actual;
                    out.report->topologyUsed+=topologyCount;
                    auto& committed=out.report->commitStorage[
                        static_cast<std::size_t>(out.report->commitsUsed++)];
                    committed.operationOrdinal=metrics.committedOperations;
                    committed.attemptOrdinal=attempt.attemptOrdinal;
                    committed.selected=attempt.selected;
                    committed.predictedNode=commit.predictedNode;
                    committed.actualNode=commit.actualNode;
                    committed.topologyOffset=topologyOffset;
                    committed.topologyCount=topologyCount;
                    committed.selectedQuality=attempt.selectedQuality;
                    committed.createdMinimumQuality=
                        attempt.createdMinimumQuality;
                    committed.postCommitStamp=commit.stamp;
                    attempt.disposition=T4AAttemptDisposition::Committed;
                    attempt.reason=T4AReason::None;
                    attempt.workAfter=coreMetrics.workUnitsConsumed;
                    ++metrics.committedOperations;
                    metrics.staleRows+=incidentalRetirements;
                    metrics.blockedInvalidations=
                        coreMetrics.blockedRowsInvalidated;
                    metrics.queuePushes+=createdViolations.size();
                    metrics.queueHighWater=(std::max<std::uint64_t>)(
                        metrics.queueHighWater,queue.size());
#if defined(PAACS_TOPO_TESTING)
                    phase3TailGuard.stop();
#endif
                    activePhase=T4APhase::LocalAudit;
                    bool localAudit=auditCommittedDelta(commit,ws)&&
                        mutableArtifact.vertexKeyByHandle.find(
                            commit.actualNode)!=
                            mutableArtifact.vertexKeyByHandle.end();
#if defined(PAACS_TOPO_TESTING)
                    if (t4a_detail::injectPostCommitAuditFailureForTesting())
                        localAudit=false;
#endif
                    if (!localAudit) {
                        out.status=T4AStatus::Failed;
                        out.failure={T4AReason::InternalError,activePhase,true,
                                     selected.row.key.tet};
                        return out;
                    }
                } catch (const WorkBudgetStop&) {
                    tx.abandon();
                    throw;
                } catch (const std::bad_alloc&) {
                    tx.abandon();
                    setTerminal(selected,attempt,T4AReason::AllocationFailed,
                                activePhase);
#if defined(PAACS_TOPO_TESTING)
                    auto& trace=t4a_detail::digestTraceForTesting();
                    if (trace.faultPending) {
                        trace.faultPending=false;
                        recordT4StateDigest(
                            trace.faultCheckpoint,true,ws,coordinates,
                            inserted,coordinateIndex,cells,edgeStars,
                            faceOwners,nodeStars,queue,blocked,
                            blockedDependencies,locationIndex,mutableArtifact,
                            *out.report,metrics);
                    }
#endif
                    out.status=T4AStatus::Failed;
                    out.failure={T4AReason::AllocationFailed,activePhase,true,
                                 selected.row.key.tet};
                    return out;
                } catch (...) {
                    tx.abandon();
                    setTerminal(selected,attempt,
                        topologyPrepared
                            ? T4AReason::CoordinatorPreparationFailed:
                              T4AReason::TopologyPreparationFailed,
                        activePhase);
                    out.status=T4AStatus::Failed;
                    out.failure={attempt.reason,activePhase,true,
                                 selected.row.key.tet};
                    return out;
                }
            } catch (const WorkBudgetStop&) {
                setTerminal(selected,attempt,T4AReason::WorkBudgetReached,
                            activePhase);
                incomplete=true;
                incompleteReason=T4AReason::WorkBudgetReached;
                incompletePhase=activePhase;
                break;
            } catch (const std::bad_alloc&) {
                setTerminal(selected,attempt,T4AReason::AllocationFailed,
                            activePhase);
#if defined(PAACS_TOPO_TESTING)
                auto& trace=t4a_detail::digestTraceForTesting();
                if (trace.faultPending) {
                    trace.faultPending=false;
                    recordT4StateDigest(
                        trace.faultCheckpoint,true,ws,coordinates,inserted,
                        coordinateIndex,cells,edgeStars,faceOwners,nodeStars,
                        queue,blocked,blockedDependencies,locationIndex,
                        mutableArtifact,*out.report,metrics);
                }
#endif
                out.status=T4AStatus::Failed;
                out.failure={T4AReason::AllocationFailed,activePhase,true,
                             selected.row.key.tet};
                return out;
            }
        }

        if (!incomplete&&queue.empty()&&!blocked.empty()) {
            incomplete=true;
            incompleteReason=T4AReason::ResidualQuality;
            incompletePhase=T4APhase::CandidateGeneration;
        }
        workBudget.stopCounting();
        syncT4CoreMetrics(coreMetrics,metrics);
        return finishAudited(
            incomplete,incompleteReason,incompletePhase);
    } catch (const std::bad_alloc&) {
        out.status=T4AStatus::Failed;
        out.failure={T4AReason::AllocationFailed,activePhase};
        if (activePhase!=T4APhase::PublicationPreparation)
            out.retainedResult.reset();
        return out;
    } catch (...) {
        out.status=T4AStatus::Failed;
        out.failure={T4AReason::InternalError,activePhase};
        out.retainedResult.reset();
        return out;
    }
}

#if defined(PAACS_TOPO_TESTING)
namespace t4a_detail {
bool checkedTetQualityForTesting(
    std::array<Vec3,4> coordinates,
    std::array<RefinementVertexKey,4> keys,double& quality) noexcept {
    return checkedTetQualityCanonical(coordinates,keys,quality);
}
} // namespace t4a_detail
#endif

} // namespace paacs::algo
