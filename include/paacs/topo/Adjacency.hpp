// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// paacs::topo - M4-B.5: cell ownership / adjacency through faces.
//
// A DERIVED, immutable-after-build view over the frozen M4-B.2/.3/.4
// layers - it introduces no new mutable topology state and never modifies
// the pools. Built by a factory from a fully READY state and validated
// afterwards by semantic recomputation against the FacePool.
//
// GATE (reviewer contract): construction requires
//     storageEdgeReady(nodes, cells)
//  && edgeTopologyReady(nodes, cells, edges)
//  && faceTopologyReady(nodes, cells, edges, faces)
// and throws otherwise.
//
// RELATIONS PROVED:
//   * cell -> faces                (copied from FacePool cell incidence)
//   * face -> incident cells      (via the FacePool; adjacency records the
//                                  derived per-class views below)
//   * cell -> neighboring cells   (through exactly-2-incident interior
//                                  faces; one Neighbor per shared face)
//   * boundary face ownership     (face + its single owner)
//   * interior face ownership     (face + both owners)
//   * nonmanifold faces           (>=3 incident cells: EXPLICITLY
//                                  represented, with all incident cells;
//                                  NO pairwise neighbor relations are
//                                  fabricated through them)
//
// INVALIDATION RULE: adjacency is a snapshot. It is valid only while every
// recorded cell and face handle is still current AND the recomputed
// relations from the live FacePool match the recorded ones exactly.
// checkValid()/auditCellAdjacency() prove this semantically; any pool
// mutation after build (cell removal, face insertion/removal, slot reuse)
// invalidates it and is DETECTED, never silently tolerated. Queries are
// generation-validated so a reused cell slot cannot resolve against stale
// entries (no ghost adjacency).
//
// STORAGE HONESTY: entries are slot-indexed with generation validation
// (same documented model as the pools' incidence tables), bound to the
// originating pool ids. 2D cells (no faces) have no entries.
//
// SCOPE (reviewer): the relations above only. DEFERRED: Delaunay
// migration, public Mesh replacement, boundary recovery, memory reports
// (PM-SOLV-001), parallel construction.
//
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "paacs/diag/Diagnostics.hpp"
#include "paacs/topo/FaceTopology.hpp"

namespace paacs::topo {

struct CellAdjacencyCapacitySnapshot {
    std::size_t entryBucketCount = 0;
    std::size_t entryCount = 0;
    std::size_t entryFaceCapacity = 0;
    std::size_t entryNeighborCapacity = 0;
    std::size_t boundaryCapacity = 0;
    std::size_t interiorCapacity = 0;
    std::size_t nonManifoldCapacity = 0;
    std::size_t nonManifoldCellCapacity = 0;
    bool complete = true;
};

class CellAdjacency {
    struct DirectBuildTag {};

public:
    struct Neighbor {
        CellHandle cell;      // the neighboring cell
        FaceHandle via;       // the shared interior face
    };
    struct BoundaryFace {
        FaceHandle face;
        CellHandle owner;     // the single incident cell
    };
    struct InteriorFace {
        FaceHandle face;
        CellHandle a, b;      // both owners (opposite orientation senses,
                              // guaranteed by the frozen face audit)
    };
    struct NonManifoldFace {
        FaceHandle face;
        std::vector<CellHandle> cells;   // all >=3 incident cells
    };

    CellAdjacency() = default;
    CellAdjacency(CellAdjacency&&) noexcept = default;
    CellAdjacency& operator=(CellAdjacency&&) noexcept = default;
    CellAdjacency(const CellAdjacency&) = delete;
    CellAdjacency& operator=(const CellAdjacency&) = delete;

    // Factory: builds the full adjacency snapshot or throws (readiness gate
    // or allocation failure). A throw produces NO object, so a failed build
    // cannot leave partial adjacency anywhere. Deterministic: faces are
    // walked in sequence order, so neighbor lists have reproducible order.
    static CellAdjacency build(const NodePool& nodes, const CellPool& cells,
                               const EdgePool& edges, const FacePool& faces) {
        return CellAdjacency(DirectBuildTag{}, nodes, cells, edges, faces);
    }

    // Rebuild in place with strong exception safety: build a complete
    // replacement first, then publish it through the verified non-allocating
    // member swap. Any throw (readiness gate or allocation failure) leaves
    // *this unchanged.
    void rebuild(const NodePool& nodes, const CellPool& cells,
                 const EdgePool& edges, const FacePool& faces) {
        CellAdjacency scratch(DirectBuildTag{}, nodes, cells, edges, faces);
        swap(scratch);
    }

    // Publication primitive for a completed scratch snapshot. Every member
    // uses the standard container's allocator-preserving swap; the compile-
    // time checks pin the no-throw contract on every supported standard
    // library. No allocation or element move occurs here.
    void swap(CellAdjacency& other) noexcept {
        static_assert(noexcept(entries_.swap(other.entries_)));
        static_assert(noexcept(boundary_.swap(other.boundary_)));
        static_assert(noexcept(interior_.swap(other.interior_)));
        static_assert(noexcept(nonManifold_.swap(other.nonManifold_)));
        entries_.swap(other.entries_);
        boundary_.swap(other.boundary_);
        interior_.swap(other.interior_);
        nonManifold_.swap(other.nonManifold_);
        using std::swap;
        swap(nodePool_, other.nodePool_);
        swap(cellPool_, other.cellPool_);
        swap(edgePool_, other.edgePool_);
        swap(facePool_, other.facePool_);
    }

    std::uint64_t boundNodePool() const { return nodePool_; }
    std::uint64_t boundCellPool() const { return cellPool_; }
    std::uint64_t boundEdgePool() const { return edgePool_; }
    std::uint64_t boundFacePool() const { return facePool_; }

    // Fully validated queries - no ghost adjacency through ANY route:
    //   * foreign CellPool            -> nothing (pool id check)
    //   * stale handle (cell removed) -> nothing (isCurrent check)
    //   * reused slot (new gen)       -> nothing (recorded-gen check)
    // A stale SNAPSHOT queried with a CURRENT cell still answers (the
    // caller may be mid-update); staleness of the snapshot as a whole is
    // what checkValid()/auditCellAdjacency() detect.
    bool hasCell(const CellPool& cells, const CellHandle& cell) const {
        if (cellPool_ != 0 && cells.poolId() != cellPool_) return false;
        if (!cells.isCurrent(cell)) return false;
        const auto it = entries_.find(cell.slot);
        return it != entries_.end() && it->second.gen == cell.gen;
    }
    std::span<const FaceHandle> cellFaces(const CellPool& cells,
                                          const CellHandle& cell) const {
        if (!hasCell(cells, cell)) return {};
        const auto& e = entries_.at(cell.slot);
        return {e.faces.data(), e.faces.size()};
    }
    std::span<const Neighbor> neighbors(const CellPool& cells,
                                        const CellHandle& cell) const {
        if (!hasCell(cells, cell)) return {};
        const auto& e = entries_.at(cell.slot);
        return {e.neighbors.data(), e.neighbors.size()};
    }
    std::span<const BoundaryFace> boundaryFaces() const {
        return {boundary_.data(), boundary_.size()};
    }
    std::span<const InteriorFace> interiorFaces() const {
        return {interior_.data(), interior_.size()};
    }
    std::span<const NonManifoldFace> nonManifoldFaces() const {
        return {nonManifold_.data(), nonManifold_.size()};
    }
    std::size_t cellCount() const { return entries_.size(); }

    CellAdjacencyCapacitySnapshot capacitySnapshot() const noexcept {
        CellAdjacencyCapacitySnapshot out;
        out.entryBucketCount = entries_.bucket_count();
        out.entryCount = entries_.size();
        out.boundaryCapacity = boundary_.capacity();
        out.interiorCapacity = interior_.capacity();
        out.nonManifoldCapacity = nonManifold_.capacity();
        for (const auto& [slot, entry] : entries_) {
            (void)slot;
            if (!capacityCheckedAdd(out.entryFaceCapacity,
                                    entry.faces.capacity()))
                out.complete = false;
            if (!capacityCheckedAdd(out.entryNeighborCapacity,
                                    entry.neighbors.capacity()))
                out.complete = false;
        }
        for (const auto& row : nonManifold_)
            if (!capacityCheckedAdd(out.nonManifoldCellCapacity,
                                    row.cells.capacity()))
                out.complete = false;
        return out;
    }

    // Semantic validity: every recorded cell/face still current, and the
    // relations recomputed from the LIVE FacePool match the recorded ones
    // exactly (order included - both sides are sequence-deterministic).
    // Any pool mutation after build makes this false.
    bool checkValid(const NodePool& nodes, const CellPool& cells,
                    const EdgePool& edges, const FacePool& faces) const {
        if (cellPool_ == 0) return false;                 // never built
        if (nodes.poolId() != nodePool_ || cells.poolId() != cellPool_ ||
            edges.poolId() != edgePool_ || faces.poolId() != facePool_)
            return false;                                  // foreign pools
        if (!faceTopologyReady(nodes, cells, edges, faces)) return false;
        try {
            // Construct the comparison snapshot directly. Do not route
            // through the public move assignment: on MSVC a standard-
            // container move may allocate despite the retained public
            // noexcept move contract.
            CellAdjacency fresh(DirectBuildTag{}, nodes, cells, edges, faces);
            return sameSnapshot(fresh);
        } catch (...) {
            return false;
        }
    }

#if defined(PAACS_TOPO_TESTING)
    // Corruption seams for adversarial audit probes.
    void testSetEntryGen(std::int32_t slot, std::uint32_t gen) {
        const auto it = entries_.find(slot);
        if (it != entries_.end()) it->second.gen = gen;
    }
    void testAddPhantomNeighbor(std::int32_t slot, const Neighbor& n) {
        const auto it = entries_.find(slot);
        if (it != entries_.end()) it->second.neighbors.push_back(n);
    }
    void testDropEntry(std::int32_t slot) { entries_.erase(slot); }
    static bool testCapacityCheckedAdd(std::size_t& total,
                                       std::size_t value) noexcept {
        return capacityCheckedAdd(total, value);
    }
#endif

private:
    static bool capacityCheckedAdd(std::size_t& total,
                                   std::size_t value) noexcept {
        if (value > (std::numeric_limits<std::size_t>::max)() - total) {
            total = (std::numeric_limits<std::size_t>::max)();
            return false;
        }
        total += value;
        return true;
    }

    struct Entry {
        std::uint32_t gen = 0;
        std::vector<FaceHandle> faces;
        std::vector<Neighbor> neighbors;
    };
    std::unordered_map<std::int32_t, Entry> entries_;
    std::vector<BoundaryFace> boundary_;
    std::vector<InteriorFace> interior_;
    std::vector<NonManifoldFace> nonManifold_;
    std::uint64_t nodePool_ = 0, cellPool_ = 0, edgePool_ = 0, facePool_ = 0;

    CellAdjacency(DirectBuildTag, const NodePool& nodes,
                  const CellPool& cells, const EdgePool& edges,
                  const FacePool& faces) {
        populateDirect(nodes, cells, edges, faces);
    }

    // The single construction implementation. It populates only a new,
    // private scratch/direct object; no published adjacency is ever touched
    // until the caller's final noexcept swap.
    void populateDirect(const NodePool& nodes, const CellPool& cells,
                        const EdgePool& edges, const FacePool& faces) {
        if (!storageEdgeReady(nodes, cells))
            throw std::invalid_argument(
                "CellAdjacency::build: storage is not edge-ready");
        if (!edgeTopologyReady(nodes, cells, edges))
            throw std::invalid_argument(
                "CellAdjacency::build: edge topology is not ready");
        if (!faceTopologyReady(nodes, cells, edges, faces))
            throw std::invalid_argument(
                "CellAdjacency::build: face topology is not ready");
        // COMPLETENESS GATE (review blocker, third recurrence of the
        // pattern): readiness above proves the pools are INTERNALLY
        // consistent, but faceTopologyReady is trivially true for an empty
        // or partial FacePool. Adjacency claims cell->faces / cell->
        // neighbors / ownership for CURRENT cells, so every current cell
        // with faces (3D: traits faceCount > 0) MUST have face incidence
        // recorded before adjacency may build. 2D / no-face cells are a
        // DOCUMENTED SKIP: they have no faces, so they legitimately carry
        // no adjacency entries and do not gate construction.
        // Rule (recorded in KNOWN_ISSUES for every derived layer):
        //   readiness = internal consistency
        //             + complete incidence coverage for all relevant
        //               current cells.
        for (const CellHandle& c : cells.activeBySequence()) {
            const auto snap = cells.get(c);
            if (!isValidElementType(snap.type)) continue;   // storage audits
            const ElementTraits& tr = checkedTraits(snap.type);
            if (tr.faceCount > 0 && !faces.hasCellFaces(cells, c))
                throw std::invalid_argument(
                    "CellAdjacency::build: a current 3D cell has no face "
                    "incidence in the FacePool; insert every cell's faces "
                    "before building adjacency");
        }
        nodePool_ = nodes.poolId();
        cellPool_ = cells.poolId();
        edgePool_ = edges.poolId();
        facePool_ = faces.poolId();
#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        // Per-cell face lists (copied from the FacePool's incidence).
        for (const CellHandle& ch : faces.incidentCells()) {
            const auto fs = faces.cellFaces(cells, ch);
            Entry e;
            e.gen = ch.gen;
            e.faces.assign(fs.begin(), fs.end());
            entries_.emplace(ch.slot, std::move(e));
#if defined(PAACS_TOPO_TESTING)
            maybeFault();
#endif
        }
        // Face classification walk, sequence order (deterministic).
        for (const FaceHandle& fh : faces.activeBySequence()) {
            const auto fd = faces.get(fh);
            if (fd.incident.size() == 1) {
                const auto& in = fd.incident[0];
                boundary_.push_back(BoundaryFace{
                    fh, CellHandle{in.cellSlot, in.cellGen, cellPool_}});
            } else if (fd.incident.size() == 2) {
                const CellHandle a{fd.incident[0].cellSlot,
                                   fd.incident[0].cellGen, cellPool_};
                const CellHandle b{fd.incident[1].cellSlot,
                                   fd.incident[1].cellGen, cellPool_};
                interior_.push_back(InteriorFace{fh, a, b});
                entries_[a.slot].neighbors.push_back(Neighbor{b, fh});
                entries_[b.slot].neighbors.push_back(Neighbor{a, fh});
            } else {
                NonManifoldFace nm;
                nm.face = fh;
                nm.cells.reserve(fd.incident.size());
                for (const auto& in : fd.incident)
                    nm.cells.push_back(
                        CellHandle{in.cellSlot, in.cellGen, cellPool_});
                nonManifold_.push_back(std::move(nm));
                // Explicitly represented; NO pairwise neighbor relations.
            }
#if defined(PAACS_TOPO_TESTING)
            maybeFault();
#endif
        }
    }

    bool sameSnapshot(const CellAdjacency& fresh) const {
        // Compare recorded vs recomputed, exactly.
        if (fresh.entries_.size() != entries_.size()) return false;
        for (const auto& [slot, e] : entries_) {
            const auto it = fresh.entries_.find(slot);
            if (it == fresh.entries_.end()) return false;
            const Entry& f = it->second;
            if (f.gen != e.gen) return false;             // reused slot
            if (f.faces != e.faces) return false;
            if (f.neighbors.size() != e.neighbors.size()) return false;
            for (std::size_t i = 0; i < f.neighbors.size(); ++i)
                if (!(f.neighbors[i].cell == e.neighbors[i].cell) ||
                    !(f.neighbors[i].via == e.neighbors[i].via))
                    return false;
        }
        auto sameB = [](const BoundaryFace& x, const BoundaryFace& y) {
            return x.face == y.face && x.owner == y.owner;
        };
        if (fresh.boundary_.size() != boundary_.size()) return false;
        for (std::size_t i = 0; i < boundary_.size(); ++i)
            if (!sameB(fresh.boundary_[i], boundary_[i])) return false;
        if (fresh.interior_.size() != interior_.size()) return false;
        for (std::size_t i = 0; i < interior_.size(); ++i) {
            const auto& x = fresh.interior_[i];
            const auto& y = interior_[i];
            if (!(x.face == y.face) || !(x.a == y.a) || !(x.b == y.b))
                return false;
        }
        if (fresh.nonManifold_.size() != nonManifold_.size()) return false;
        for (std::size_t i = 0; i < nonManifold_.size(); ++i) {
            const auto& x = fresh.nonManifold_[i];
            const auto& y = nonManifold_[i];
            if (!(x.face == y.face) || x.cells.size() != y.cells.size())
                return false;
            for (std::size_t j = 0; j < x.cells.size(); ++j)
                if (!(x.cells[j] == y.cells[j])) return false;
        }
        return true;
    }
};

// Structured adjacency audit. Ordering: pool bindings, then the face layer
// gate (which itself proves structural safety BEFORE any traversal - a
// corrupt face pool is reported here, never dereferenced), then semantic
// recomputation.
inline diag::Report auditCellAdjacency(const NodePool& nodes,
                                       const CellPool& cells,
                                       const EdgePool& edges,
                                       const FacePool& faces,
                                       const CellAdjacency& adj) {
    diag::Report r;
    auto add = [&](diag::Code code, diag::Severity sev, std::string msg) {
        diag::Diagnostic d;
        d.code = code; d.severity = sev; d.stage = "adjacency-audit";
        d.message = std::move(msg);
        r.add(std::move(d));
    };
    if (adj.boundCellPool() != 0 &&
        (adj.boundNodePool() != nodes.poolId() ||
         adj.boundCellPool() != cells.poolId() ||
         adj.boundEdgePool() != edges.poolId() ||
         adj.boundFacePool() != faces.poolId())) {
        add(diag::Code::AdjacencyPoolBindingMismatch,
            diag::Severity::InvalidOutput,
            "adjacency bound to different pools than audited");
        return r;
    }
    if (!faceTopologyReady(nodes, cells, edges, faces)) {
        add(diag::Code::FaceInvariantViolation, diag::Severity::InvalidOutput,
            "face topology is not ready; adjacency requires valid faces");
        return r;
    }
    // Completeness (review blocker): faceTopologyReady is trivially true
    // for an empty/partial FacePool; every current 3D cell must have face
    // incidence, or the adjacency snapshot is over-approving live cells it
    // silently ignored (2D / no-face cells: documented skip).
    for (const CellHandle& c : cells.activeBySequence()) {
        const auto snap = cells.get(c);
        if (!isValidElementType(snap.type)) continue;
        if (checkedTraits(snap.type).faceCount > 0 &&
            !faces.hasCellFaces(cells, c)) {
            add(diag::Code::AdjacencyIncompleteIncidence,
                diag::Severity::InvalidOutput,
                "a current 3D cell has no face incidence in the FacePool; "
                "adjacency cannot be valid over incomplete incidence");
            return r;
        }
    }
    // Ghost check first for a specific diagnostic: any recorded cell no
    // longer current means the pool mutated after build.
    bool ghost = false;
    for (const auto& bf : adj.boundaryFaces())
        if (!cells.isCurrent(bf.owner)) { ghost = true; break; }
    if (!ghost)
        for (const auto& itf : adj.interiorFaces())
            if (!cells.isCurrent(itf.a) || !cells.isCurrent(itf.b)) {
                ghost = true; break;
            }
    if (ghost) {
        add(diag::Code::GhostAdjacency, diag::Severity::InvalidOutput,
            "adjacency references a removed cell (pool mutated after "
            "build; rebuild required)");
        return r;
    }
    if (!adj.checkValid(nodes, cells, edges, faces))
        add(diag::Code::AdjacencyMismatch, diag::Severity::InvalidOutput,
            "recorded adjacency does not match relations recomputed from "
            "the live face topology (rebuild required)");
    return r;
}

inline bool cellAdjacencyReady(const NodePool& nodes, const CellPool& cells,
                               const EdgePool& edges, const FacePool& faces,
                               const CellAdjacency& adj) {
    return auditCellAdjacency(nodes, cells, edges, faces, adj)
        .acceptable(diag::Severity::Warning);
}

} // namespace paacs::topo
