// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace paacs::algo::detail {

// Private stop signal. It is caught inside the noexcept candidate boundary;
// no user callback, process-global state, or thread-local state is involved.
struct TetConstructionStopped { bool cancelled; };

struct TetConstructionBudget {
    std::uint64_t maxWork;
    std::uint64_t maxCells;
    std::uint64_t maxBytes;
    const std::atomic_bool* cancellation;
    std::uint64_t work = 0;
    std::uint64_t bytes = 0;
    std::uint64_t peakBytes = 0;

    void step(std::uint64_t count = 1) {
        if (cancellation && cancellation->load(std::memory_order_relaxed))
            throw TetConstructionStopped{true};
        if (count > maxWork - work) throw TetConstructionStopped{false};
        work += count;
    }
    static std::uint64_t product(std::uint64_t count, std::uint64_t size) {
        if (size && count > (std::numeric_limits<std::uint64_t>::max)()/size)
            throw TetConstructionStopped{false};
        return count*size;
    }
    void admit(std::uint64_t count) {
        step(0);
        if (count > maxBytes - bytes) throw TetConstructionStopped{false};
        bytes += count;
        peakBytes = (std::max)(peakBytes,bytes);
    }
    void release(std::uint64_t count) noexcept { bytes -= count; }
    void cells(std::uint64_t count) {
        step();
        if (count > maxCells) throw TetConstructionStopped{false};
    }
};

inline void constructionStep(TetConstructionBudget* budget,
                             std::uint64_t count = 1) {
    if (budget) budget->step(count);
}

// The admitted envelope includes twice each requested vector payload, plus
// old storage during reallocation. It models capacity and growth overlap,
// rather than allocator-exact heap use. Unlimited vectors retain the standard
// library's original growth policy and capacities.
template<class T>
class ConstructionVector {
    using Base = std::vector<T>;
    Base values_;
    TetConstructionBudget* budget_ = nullptr;
    std::uint64_t admitted_ = 0;
public:
    using iterator = typename Base::iterator;
    using const_iterator = typename Base::const_iterator;
    auto begin() noexcept { return values_.begin(); }
    auto end() noexcept { return values_.end(); }
    auto begin() const noexcept { return values_.begin(); }
    auto end() const noexcept { return values_.end(); }
    auto size() const noexcept { return values_.size(); }
    auto capacity() const noexcept { return values_.capacity(); }
    auto max_size() const noexcept { return values_.max_size(); }
    bool empty() const noexcept { return values_.empty(); }
    T& operator[](std::size_t i) noexcept { return values_[i]; }
    const T& operator[](std::size_t i) const noexcept { return values_[i]; }
    T& back() noexcept { return values_.back(); }
    const T& back() const noexcept { return values_.back(); }
    void clear() noexcept { values_.clear(); }
    void pop_back() noexcept { values_.pop_back(); }
    void releaseStorage() noexcept {
        Base{}.swap(values_);
        if (budget_) budget_->release(admitted_);
        admitted_=0;
    }
    ConstructionVector() = default;
    ConstructionVector(const ConstructionVector&) = delete;
    ConstructionVector& operator=(const ConstructionVector&) = delete;
    ConstructionVector(ConstructionVector&& other) noexcept
        : values_(std::move(other.values_)), budget_(other.budget_), admitted_(other.admitted_) {
        other.budget_=nullptr; other.admitted_=0;
    }
    ConstructionVector& operator=(ConstructionVector&& other) noexcept {
        if (budget_) budget_->release(admitted_);
        values_=std::move(other.values_);
        budget_=other.budget_; admitted_=other.admitted_;
        other.budget_=nullptr; other.admitted_=0;
        return *this;
    }
    ~ConstructionVector() { if (budget_) budget_->release(admitted_); }
    void setConstructionBudget(TetConstructionBudget* budget) noexcept { budget_=budget; }
    void reserve(std::size_t count) {
        if (!budget_) { values_.reserve(count); return; }
        budget_->step(0);
        if (count <= values_.capacity()) return;
        budget_->step(values_.size());
        const auto next=TetConstructionBudget::product(count,2U*sizeof(T));
        budget_->admit(next);
        try { values_.reserve(count); }
        catch (...) { budget_->release(next); throw; }
        budget_->release(admitted_);
        admitted_=next;
    }
    void ensure(std::size_t count) {
        if (!budget_ || count<=values_.capacity()) return;
        const auto max=values_.max_size();
        const auto doubled=values_.capacity()>max/2 ? max : values_.capacity()*2;
        reserve((std::max)(count,doubled));
    }
    void resize(std::size_t count) {
        if (budget_) budget_->step(count>values_.size() ? count-values_.size() : 0);
        ensure(count); values_.resize(count);
    }
    void resize(std::size_t count,const T& value) {
        if (budget_) budget_->step(count>values_.size() ? count-values_.size() : 0);
        ensure(count); values_.resize(count,value);
    }
    void assign(std::size_t count,const T& value) {
        if (budget_) budget_->step(count);
        ensure(count); values_.assign(count,value);
    }
    void push_back(const T& value) { ensure(values_.size()+1); values_.push_back(value); }
    void push_back(T&& value) { ensure(values_.size()+1); values_.push_back(std::move(value)); }
    template<class... Args> T& emplace_back(Args&&... args) {
        ensure(values_.size()+1); return values_.emplace_back(std::forward<Args>(args)...);
    }
    void swap(ConstructionVector& other) noexcept {
        values_.swap(other.values_); std::swap(budget_,other.budget_); std::swap(admitted_,other.admitted_);
    }
};

} // namespace paacs::algo::detail
