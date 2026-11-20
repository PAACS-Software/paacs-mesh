// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// Minimal parallel primitives (std::thread, no dependencies).
//
// Design rule for the whole library: DETERMINISM ACROSS THREAD COUNTS.
// Reductions are chunked with a fixed chunk size and folded in chunk order,
// so floating-point results are bit-identical whether run on 1 or 64 cores.
// Mutating phases either write to disjoint indices (parallelFor) or are
// restructured as parallel-evaluate + ordered-serial-commit.
//
#include <cstddef>
#include <algorithm>
#include <functional>
#include <thread>
#include <vector>

namespace paacs::core {

// 0 = auto (hardware_concurrency). Atomic process-wide compatibility setting.
// Prefer request-local thread counts for independent callers.
int  threadCount();
void setThreadCount(int n);
int  effectiveThreads();   // resolved value (>= 1)
// Introspection: is the calling thread inside a pool dispatch (re-entrancy
// guard set)? Exists mainly for qualification (proving the concurrent-
// dispatch rejection restores the guard). Cheap and side-effect-free;
// exposed rather than test-gated so the parallel contract test links
// without a special build, but it is not intended for production control
// flow.
bool inPoolContext();

// Parallel for over [begin, end). body(i) must write only to index-i state.
// ---------------------------------------------------------------------------
// CONCURRENCY CONTRACT (PM-PAR-001 closure documentation, per directive):
//   * The persistent pool supports ONE external dispatcher at a time.
//     Concurrent independent calls to parallelFor/parallelReduce from
//     unrelated threads are UNSUPPORTED and rejected loudly (logic_error),
//     not silently corrupted, until serialization, per-session pools, or a
//     queued job system is implemented (M14 scope).
//   * NESTED dispatch from inside a worker body is DEFINED behavior: it
//     degrades to serial execution via the re-entrancy guard (a nested
//     parallel wait would deadlock on the pool's own workers).
//   * EXCEPTIONS FROM BODIES ARE PROPAGATED (round-3 contract): the
//     first exception is captured, remaining chunks are cancelled
//     (fetched, not executed), the pool drains, and the exception is
//     rethrown on the DISPATCHING thread. Some indices may therefore
//     never execute; side-effect atomicity is the caller's
//     responsibility (transactional code validates before any parallel
//     apply). The serial and nested-serial paths propagate naturally
//     with identical observable semantics.
// ---------------------------------------------------------------------------
void parallelFor(std::size_t begin, std::size_t end,
                 const std::function<void(std::size_t)>& body);

// Request-local variant. A positive requestedThreads value limits this one
// dispatch without mutating the process-wide compatibility setting.
void parallelForWithThreads(
    std::size_t begin, std::size_t end, int requestedThreads,
    const std::function<void(std::size_t)>& body);

// Parallel over contiguous chunks: body(chunkIndex, lo, hi).
void parallelChunks(std::size_t begin, std::size_t end, std::size_t chunkSize,
                    const std::function<void(std::size_t, std::size_t,
                                             std::size_t)>& body);
void parallelChunksWithThreads(
    std::size_t begin, std::size_t end, std::size_t chunkSize,
    int requestedThreads,
    const std::function<void(std::size_t, std::size_t, std::size_t)>& body);

// Deterministic map-reduce: each fixed-size chunk accumulates sequentially
// into its own Acc via fold(acc, i); chunk results merge in chunk order.
template <class Acc, class Fold, class Merge>
Acc parallelReduce(std::size_t begin, std::size_t end, std::size_t chunkSize,
                   Acc init, Fold&& fold, Merge&& merge) {
    if (end <= begin) return init;
    chunkSize = std::max<std::size_t>(1, chunkSize);
    const std::size_t n = end - begin;
    const std::size_t nChunks = n / chunkSize + (n % chunkSize != 0);
    std::vector<Acc> partial(nChunks, init);
    parallelChunks(begin, end, chunkSize,
                   [&](std::size_t ci, std::size_t lo, std::size_t hi) {
                       Acc a = init;
                       for (std::size_t i = lo; i < hi; ++i) fold(a, i);
                       partial[ci] = a;
                   });
    Acc out = init;
    for (const Acc& a : partial) merge(out, a);   // fixed order => bitwise
    return out;                                   //   reproducible
}

} // namespace paacs::core
