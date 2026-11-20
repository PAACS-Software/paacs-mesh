// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/core/Parallel.hpp"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <mutex>
#include <stdexcept>

namespace paacs::core {

namespace {

std::atomic<int> gThreads{0}; // 0 = auto

// Persistent worker pool. Spawning/joining std::thread teams per call is
// catastrophic for fine-grained phases (measured externally: smoothing went
// 40 ms -> 1.27 s on a many-core machine because each graph color of each
// iteration paid a full team spawn). Workers are created once on first
// parallel use and parked on a condition variable between jobs.
class Pool {
  public:
    static Pool &instance() {
        static Pool p;
        return p;
    }

    // Blocks until every chunk of the job has run. Chunk index determines
    // WHAT is computed, never in which order results are combined, so pool
    // scheduling cannot affect determinism.
    void run(std::size_t nChunks, const std::function<void(std::size_t)> &fn, int maxWorkers) {
        // One external dispatcher at a time (see Parallel.hpp contract).
        // Nested same-thread dispatch never reaches here (re-entrancy guard
        // degrades it to serial), so a second concurrent dispatch can only
        // be an unsupported cross-thread use: fail loudly.
        // EXCEPTION CONTRACT (review round 3, option 1): body exceptions
        // are CAPTURED (first wins), remaining chunks are CANCELLED
        // (fetched but not executed, so completion accounting stays
        // exact), workers drain, and the exception is RETHROWN on the
        // dispatcher after the pool is quiescent. This replaces the
        // previous behavior where a worker exception terminated the
        // process and a caller exception unwound run() while workers
        // still held the caller-owned function. Partial side effects of
        // cancelled/executed chunks are the CALLER'S responsibility:
        // transactional topology code must validate before any parallel
        // apply (directive 2.3).
        if (dispatching_.exchange(true, std::memory_order_acq_rel))
            throw std::logic_error("core::Parallel: concurrent external dispatch is "
                                   "unsupported (one dispatcher at a time; see Parallel.hpp)");
        struct Release {
            std::atomic<bool> &f;
            ~Release() { f.store(false, std::memory_order_release); }
        } release{dispatching_};
        ensureWorkers(maxWorkers);
        cancel_.store(false, std::memory_order_relaxed);
        jobError_ = nullptr;
        {
            std::unique_lock<std::mutex> lk(mx_);
            job_ = &fn;
            jobChunks_ = nChunks;
            jobWorkerLimit_ = std::max(0, maxWorkers - 1);
            next_.store(0, std::memory_order_relaxed);
            pending_.store(nChunks, std::memory_order_relaxed);
            ++jobId_;
        }
        cv_.notify_all();
        work(&fn, nChunks); // the calling thread participates
        // CORRECTNESS-CRITICAL (fixed race): a straggler worker may have
        // captured job_ under the lock but not yet fetched a chunk. If
        // run() returned on pending_==0 alone, the caller's fn would be
        // destroyed while the straggler still holds a pointer to it, and
        // the NEXT job's reset of next_ would hand that straggler fresh
        // chunk indices to execute against the DANGLING previous function
        // (observed as std::bad_function_call; could also silently consume
        // the new job's chunks). activeWorkers_ counts workers between
        // capture and completion; run() waits for both.
        std::unique_lock<std::mutex> lk(mx_);
        doneCv_.wait(lk, [&] { return pending_.load() == 0 && activeWorkers_ == 0; });
        job_ = nullptr;
        if (jobError_) {
            std::exception_ptr e = jobError_;
            jobError_ = nullptr;
            cancel_.store(false, std::memory_order_relaxed);
            lk.unlock();
            std::rethrow_exception(e);
        }
        cancel_.store(false, std::memory_order_relaxed);
    }

    ~Pool() {
        {
            std::unique_lock<std::mutex> lk(mx_);
            stopping_ = true;
        }
        cv_.notify_all();
        for (auto &t : workers_)
            t.join();
    }

  private:
    Pool() = default;

    void ensureWorkers(int wanted) {
        std::unique_lock<std::mutex> lk(mx_);
        while (static_cast<int>(workers_.size()) < wanted - 1) {
            const int workerIndex = static_cast<int>(workers_.size());
            workers_.emplace_back([this, workerIndex] { workerLoop(workerIndex); });
        }
    }

    void work(const std::function<void(std::size_t)> *fn, std::size_t n) {
        for (;;) {
            const std::size_t c = next_.fetch_add(1);
            if (c >= n)
                break;
            if (!cancel_.load(std::memory_order_relaxed)) {
                try {
                    (*fn)(c);
                } catch (...) {
                    std::lock_guard<std::mutex> lk(mx_);
                    if (!jobError_)
                        jobError_ = std::current_exception();
                    cancel_.store(true, std::memory_order_relaxed);
                }
            }
            // Cancelled chunks are consumed WITHOUT execution so the
            // completion accounting (pending_) stays exact.
            if (pending_.fetch_sub(1) == 1) {
                doneCv_.notify_all();
            }
        }
    }

    void workerLoop(int workerIndex) {
        std::uint64_t seen = 0;
        for (;;) {
            const std::function<void(std::size_t)> *fn = nullptr;
            std::size_t n = 0;
            {
                std::unique_lock<std::mutex> lk(mx_);
                cv_.wait(lk, [&] { return stopping_ || (job_ && jobId_ != seen); });
                if (stopping_)
                    return;
                seen = jobId_;
                if (workerIndex < jobWorkerLimit_) {
                    fn = job_;
                    n = jobChunks_;
                    ++activeWorkers_; // same critical section as capture
                }
            }
            if (fn) {
                work(fn, n);
                {
                    std::unique_lock<std::mutex> lk(mx_);
                    --activeWorkers_;
                }
                doneCv_.notify_all();
            }
        }
    }

    std::mutex mx_;
    std::condition_variable cv_, doneCv_;
    std::vector<std::thread> workers_;
    const std::function<void(std::size_t)> *job_ = nullptr;
    std::size_t jobChunks_ = 0;
    int jobWorkerLimit_ = 0;
    std::uint64_t jobId_ = 0;
    std::atomic<std::size_t> next_{0}, pending_{0};
    int activeWorkers_ = 0; // workers holding a captured job pointer
    std::atomic<bool> dispatching_{false};
    std::atomic<bool> cancel_{false};
    std::exception_ptr jobError_;
    bool stopping_ = false;
};

// Re-entrancy guard: a parallel body must never dispatch to the pool again
// (it would deadlock waiting on itself); nested calls degrade to serial.
thread_local bool tlInPool = false;

// Below this many items dispatch overhead outweighs the work; run serially.
// Deliberately conservative - meshing bodies are heavyweight, and small
// batches (e.g. one graph color of a small mesh) must not pay pool latency.
constexpr std::size_t kMinParallelItems = 1024;

} // namespace

// Test/introspection: is THIS thread currently inside a pool dispatch
// (re-entrancy guard set)? Proves the concurrent-dispatch rejection path
// restores the caller's guard (PM-PAR-001). External linkage.
bool inPoolContext() {
    return tlInPool;
}

int threadCount() {
    return gThreads.load(std::memory_order_relaxed);
}
void setThreadCount(int n) {
    gThreads.store(std::max(0, n), std::memory_order_relaxed);
}

int effectiveThreads() {
    const int configured = threadCount();
    if (configured > 0)
        return configured;
    const unsigned hw = std::thread::hardware_concurrency();
    return hw ? static_cast<int>(hw) : 1;
}

void parallelChunksWithThreads(
    std::size_t begin, std::size_t end, std::size_t chunkSize, int requestedThreads,
    const std::function<void(std::size_t, std::size_t, std::size_t)> &body) {
    if (end <= begin)
        return;
    chunkSize = std::max<std::size_t>(1, chunkSize);
    const std::size_t n = end - begin;
    const std::size_t nChunks = n / chunkSize + (n % chunkSize != 0);
    const int resolved = requestedThreads > 0 ? requestedThreads : effectiveThreads();
    const int nt = static_cast<int>(
        std::min<std::size_t>(static_cast<std::size_t>(std::max(1, resolved)), nChunks));

    if (nt <= 1 || n < kMinParallelItems || tlInPool) {
        for (std::size_t c = 0; c < nChunks; ++c) {
            const std::size_t lo = begin + c * chunkSize;
            body(c, lo, lo + std::min(end - lo, chunkSize));
        }
        return;
    }

    struct TlGuard {
        bool previous = tlInPool;
        TlGuard() { tlInPool = true; }
        ~TlGuard() { tlInPool = previous; }
    } guard;
    std::function<void(std::size_t)> chunkFn = [&](std::size_t c) {
        TlGuard workerGuard;
        const std::size_t lo = begin + c * chunkSize;
        body(c, lo, lo + std::min(end - lo, chunkSize));
    };
    Pool::instance().run(nChunks, chunkFn, nt);
}

void parallelChunks(std::size_t begin, std::size_t end, std::size_t chunkSize,
                    const std::function<void(std::size_t, std::size_t, std::size_t)> &body) {
    parallelChunksWithThreads(begin, end, chunkSize, 0, body);
}

void parallelForWithThreads(std::size_t begin, std::size_t end, int requestedThreads,
                            const std::function<void(std::size_t)> &body) {
    if (end <= begin)
        return;
    const std::size_t n = end - begin;
    const int resolved = requestedThreads > 0 ? requestedThreads : effectiveThreads();
    const std::size_t target = static_cast<std::size_t>(std::max(1, resolved)) * 8U;
    const std::size_t chunk = std::max<std::size_t>(1, n / std::max<std::size_t>(1, target));
    parallelChunksWithThreads(begin, end, chunk, requestedThreads,
                              [&](std::size_t, std::size_t lo, std::size_t hi) {
                                  for (std::size_t i = lo; i < hi; ++i)
                                      body(i);
                              });
}

void parallelFor(std::size_t begin, std::size_t end, const std::function<void(std::size_t)> &body) {
    parallelForWithThreads(begin, end, 0, body);
}

} // namespace paacs::core
