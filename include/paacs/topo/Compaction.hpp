// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// M4-B.7 — Coordinated workspace compaction. Design authority:
// docs/M4B7_DESIGN.md. M4-B.7 is CLOSED/FROZEN in full (Stage H, C1,
// C2) - closure by managerial risk acceptance 2026-07-16 with the
// external CMake gate DEFERRED as a release/CI obligation. This file is
// FROZEN: no further compaction extensions without a reopening review.
//
// Compaction is a coordinated global identity migration: it builds a NEW
// identity domain (fresh pools, new pool IDs) and connects the domains
// through explicit remap tables. It is NON-DESTRUCTIVE: the source is
// const and remains observably unchanged on success and on failure.
// The destination is temporary until returned; ordinary RAII destruction
// provides the exception guarantee — there is deliberately NO
// transaction-style rollback machinery here.
//
// Scope now implemented: MemoryReport (approximate retained memory with
// a partitioned byte classification), compactWorkspace, remaps for ALL
// FOUR handle types (nodes/cells by live-sequence rank; edges/faces
// reconstructed through canonical logical keys), the V2 equivalence
// gate, deterministic strict rebuild, and the CompactionOptions benefit
// policy (a no-benefit outcome is reported, never an error). The
// destination rebuilds complete, metadata-faithful edge/face topology,
// and the TopologyWorkspace constructor validates topology and builds
// adjacency — so compaction never rebuilds adjacency itself.
//
// Allocation-failure qualification: tests/probe_compaction.cpp sweeps
// EVERY allocation this function makes (global operator-new countdown),
// proving the source is observably unchanged after each one.

#include "paacs/topo/Workspace.hpp"

#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace paacs::topo {

// APPROXIMATE RETAINED MEMORY (design amendment 2026-07-14; classes
// corrected in the C1 review). Aggregated from the pools' OWN narrow
// byte fields — never derived by subtraction. The three classes
// PARTITION the total:
//   totalBytes == recordBytes + allocatorBytes + estimatedOtherBytes
// recordBytes and allocatorBytes are capacity arithmetic over known
// element sizes; estimatedOtherBytes covers hash-table bucket/node and
// incidence-container ESTIMATES. Consistent for before/after
// comparison; NOT process RSS.
struct MemoryReport {
    std::size_t nodeBytes = 0;
    std::size_t cellBytes = 0;
    std::size_t edgeBytes = 0;
    std::size_t faceBytes = 0;
    std::size_t recordBytes = 0;          // KNOWN record capacities
    std::size_t allocatorBytes = 0;       // KNOWN allocator capacities
    std::size_t estimatedOtherBytes = 0;  // hash/incidence ESTIMATES
    std::size_t totalBytes = 0;
    // Invariant (C1 review): the classes PARTITION totalBytes.
    bool consistent() const {
        return totalBytes ==
               recordBytes + allocatorBytes + estimatedOtherBytes;
    }
};

inline MemoryReport measureWorkspaceMemory(const TopologyWorkspace& ws) {
    MemoryReport r;
    const auto n = ws.nodes().stats();
    const auto c = ws.cells().stats();
    const auto e = ws.edges().stats();
    const auto f = ws.faces().stats();
    r.nodeBytes = n.totalBytes;
    r.cellBytes = c.totalBytes;
    r.edgeBytes = e.totalBytes;
    r.faceBytes = f.totalBytes;
    // C1 review blocker: the classes are taken from the pools' own
    // narrow fields (PoolStats::recordBytes / estimatedOtherBytes),
    // never derived by subtraction.
    r.recordBytes =
        n.recordBytes + c.recordBytes + e.recordBytes + f.recordBytes;
    r.allocatorBytes = n.allocatorBytes + c.allocatorBytes +
                       e.allocatorBytes + f.allocatorBytes;
    r.estimatedOtherBytes = n.estimatedOtherBytes + c.estimatedOtherBytes +
                            e.estimatedOtherBytes + f.estimatedOtherBytes;
    r.totalBytes = r.nodeBytes + r.cellBytes + r.edgeBytes + r.faceBytes;
    return r;
}

// Old→new handle map for one handle type. Keys are COMPLETE source
// handles (pool, slot, generation): lookups with free-slot handles,
// retired-slot handles, or stale generations of surviving slots find
// nothing (the source pool never issued a mapping for them).
template <class H>
class HandleMap {
public:
    void insert(const H& oldHandle, const H& newHandle) {
        map_.emplace(pack(oldHandle), newHandle);
    }
    std::optional<H> find(const H& oldHandle) const {
        const auto it = map_.find(pack(oldHandle));
        if (it == map_.end()) return std::nullopt;
        return it->second;
    }
    std::size_t size() const { return map_.size(); }

private:
    struct Key {
        std::uint64_t pool;
        std::int32_t slot;
        std::uint32_t gen;
        bool operator==(const Key&) const = default;
    };
    struct KeyHash {
        std::size_t operator()(const Key& k) const {
            std::uint64_t h = 1469598103934665603ULL;
            for (std::uint64_t v :
                 {k.pool, static_cast<std::uint64_t>(
                              static_cast<std::uint32_t>(k.slot)),
                  static_cast<std::uint64_t>(k.gen)}) {
                h ^= v;
                h *= 1099511628211ULL;
            }
            return static_cast<std::size_t>(h);
        }
    };
    static Key pack(const H& h) { return Key{h.pool, h.slot, h.gen}; }
    std::unordered_map<Key, H, KeyHash> map_;
};

// Benefit policy. Default: beneficial iff the rebuild retains
// strictly fewer approximate bytes. Thresholds let callers refuse
// negligible gains. A no-benefit result is NOT an error - the compacted
// workspace is still returned so the caller decides.
struct CompactionOptions {
    std::size_t minSavedBytes = 0;
    // FINITE value in [0,1]. NaN, infinities, negatives and values > 1
    // are rejected (C2 review: NaN silently made beneficial=false,
    // negatives silently relaxed the threshold, > 1 disabled it).
    double minSavedFraction = 0.0;
};

struct CompactionResult {
    std::unique_ptr<TopologyWorkspace> workspace;  // the NEW identity
                                                   // domain (non-movable
                                                   // type: owned by ptr)
    HandleMap<NodeHandle> nodeMap;
    HandleMap<CellHandle> cellMap;
    HandleMap<EdgeHandle> edgeMap;   // canonical-key reconstruction
    HandleMap<FaceHandle> faceMap;   // canonical-key reconstruction
    MemoryReport before, after;
    bool beneficial = false;
    std::size_t savedBytes = 0;      // 0 when the rebuild is not smaller
};

namespace detail {
#if defined(PAACS_TOPO_TESTING)
// Compaction-owned stage fault seam (this file is NOT frozen): counts
// down across the construction stages so the failure proof can inject
// bad_alloc at every stage boundary. Ordinary RAII destruction of the
// temporaries is the exception guarantee.
inline int& compactionFaultCountdown() {
    static int v = 0;
    return v;
}
inline void compactionMaybeFault() {
    if (compactionFaultCountdown() > 0 && --compactionFaultCountdown() == 0)
        throw std::bad_alloc();
}
#endif
inline void compactionStageSeam() {
#if defined(PAACS_TOPO_TESTING)
    compactionMaybeFault();
#endif
}
}  // namespace detail

// Non-destructive coordinated compaction. Deterministic construction
// order (reviewer sequence): validate options → live nodes in source
// live-sequence order → node map → live cells in source live-sequence
// order with mapped connectivity → cell map → per-cell edge incidence →
// per-3D-cell face incidence → edge/face metadata AND edge/face remaps
// via canonical mapped keys → workspace constructed only once incidence
// is complete (the constructor validates topology and builds adjacency)
// → map totals over all four pools → V2 equivalence gate → memory
// report → benefit evaluation.
inline CompactionResult compactWorkspace(
    const TopologyWorkspace& source, const CompactionOptions& options = {}) {
    // Option-domain validation BEFORE any construction work.
    if (!std::isfinite(options.minSavedFraction) ||
        options.minSavedFraction < 0.0 || options.minSavedFraction > 1.0)
        throw std::invalid_argument(
            "compactWorkspace: CompactionOptions::minSavedFraction must "
            "be a finite value in [0,1]");
    CompactionResult result;
    result.before = measureWorkspaceMemory(source);
    const std::uint64_t sourceSemantic = source.semanticTopologyHash();

    detail::compactionStageSeam();  // stage 1: node construction
    NodePool nodes;
    for (const NodeHandle& oldN : source.nodes().activeBySequence()) {
        const auto d = source.nodes().get(oldN);
        const NodeHandle newN = nodes.add(d.coord, d.cls);
        result.nodeMap.insert(oldN, newN);
    }

    detail::compactionStageSeam();  // stage 2: cell construction
    CellPool cells;
    std::vector<std::pair<CellHandle, CellHandle>> cellPairs;  // old,new
    for (const CellHandle& oldC : source.cells().activeBySequence()) {
        const auto sn = source.cells().get(oldC);
        std::vector<NodeHandle> conn;
        conn.reserve(sn.count);
        for (std::uint8_t i = 0; i < sn.count; ++i) {
            const auto mapped = result.nodeMap.find(sn.nodes[i]);
            if (!mapped)
                throw std::logic_error(
                    "compactWorkspace: source cell references an unmapped "
                    "node (source incidence defect)");
            conn.push_back(*mapped);
        }
        const CellHandle newC = cells.add(sn.type, conn, sn.cls, nodes);
        result.cellMap.insert(oldC, newC);
        cellPairs.emplace_back(oldC, newC);
    }

    detail::compactionStageSeam();  // stage 3: edge incidence
    EdgePool edges;
    for (const auto& [oldC, newC] : cellPairs)
        edges.insertCellEdges(nodes, cells, newC);

    detail::compactionStageSeam();  // stage 4: face incidence
    FacePool faces;
    for (const auto& [oldC, newC] : cellPairs) {
        const auto t = checkedTraits(cells.get(newC).type);
        if (t.faceCount > 0) faces.insertCellFaces(nodes, cells, edges, newC);
    }

    detail::compactionStageSeam();  // stage 5: edge/face metadata
    // Edge metadata AND the edge remap via canonical mapped endpoint
    // keys — never iteration order.
    for (const EdgeHandle& oldE : source.edges().activeBySequence()) {
        const auto d = source.edges().get(oldE);
        const auto na = result.nodeMap.find(d.a);
        const auto nb = result.nodeMap.find(d.b);
        if (!na || !nb)
            throw std::logic_error(
                "compactWorkspace: source edge references an unmapped "
                "node");
        const auto newE = edges.find(nodes, *na, *nb);
        if (!newE)
            throw std::logic_error(
                "compactWorkspace: rebuilt topology is missing a source "
                "edge (source incidence defect)");
        result.edgeMap.insert(oldE, *newE);      // canonical-key map
        edges.setClassification(*newE, d.cls);
        edges.setRequiredConstraint(*newE, d.requiredConstraint);
    }
    for (const FaceHandle& oldF : source.faces().activeBySequence()) {
        const auto d = source.faces().get(oldF);
        std::vector<NodeHandle> corners;
        corners.reserve(d.cornerCount);
        for (std::uint8_t i = 0; i < d.cornerCount; ++i) {
            const auto m = result.nodeMap.find(d.corners[i]);
            if (!m)
                throw std::logic_error(
                    "compactWorkspace: source face references an unmapped "
                    "node");
            corners.push_back(*m);
        }
        const auto newF = faces.find(nodes, corners);
        if (!newF)
            throw std::logic_error(
                "compactWorkspace: rebuilt topology is missing a source "
                "face (source incidence defect)");
        result.faceMap.insert(oldF, *newF);      // canonical-key map
        faces.setClassification(*newF, d.cls);
    }

    detail::compactionStageSeam();  // stage 6: workspace construction
    // Constructed ONLY now, with incidence complete: the constructor
    // enforces the composite readiness gate. Direct construction into
    // the unique_ptr — the workspace type is non-movable.
    result.workspace = std::make_unique<TopologyWorkspace>(
        std::move(nodes), std::move(cells), std::move(edges),
        std::move(faces));

    // Internal map-total validation over ALL FOUR pools, BEFORE the
    // equivalence gate, matching the documented validation order: the
    // library asserts totality itself rather than relying on tests.
    if (result.nodeMap.size() != source.nodes().stats().live ||
        result.cellMap.size() != source.cells().stats().live ||
        result.edgeMap.size() != source.edges().stats().live ||
        result.faceMap.size() != source.faces().stats().live)
        throw std::logic_error(
            "compactWorkspace: remap tables are not total over live "
            "source handles (compaction defect)");

    detail::compactionStageSeam();  // stage 7: equivalence gate
    // NOTE (C1 review): the TopologyWorkspace constructor already builds
    // adjacency and sets adjBuiltFrom_ == revision_, so the workspace is
    // Ready on construction; the previous rebuildAdjacency() here was a
    // redundant second full build (temporary allocation + peak memory).
    if (result.workspace->semanticTopologyHash() != sourceSemantic)
        throw std::logic_error(
            "compactWorkspace: V2 semantic equivalence gate failed "
            "(compaction defect)");
    result.after = measureWorkspaceMemory(*result.workspace);
    // Benefit policy: a no-benefit outcome is reported, never an error;
    // the compacted workspace is returned either way.
    if (result.after.totalBytes < result.before.totalBytes) {
        result.savedBytes =
            result.before.totalBytes - result.after.totalBytes;
        const double frac =
            result.before.totalBytes == 0
                ? 0.0
                : static_cast<double>(result.savedBytes) /
                      static_cast<double>(result.before.totalBytes);
        result.beneficial = result.savedBytes >= options.minSavedBytes &&
                            frac >= options.minSavedFraction;
    } else {
        result.savedBytes = 0;
        result.beneficial = false;
    }
    return result;
}

}  // namespace paacs::topo
