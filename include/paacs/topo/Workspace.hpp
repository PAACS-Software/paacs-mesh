// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// paacs::topo - M4-B.6A stages 1-2: TopologyWorkspace + TopologyTransaction
// state machine over the frozen pools. Cell mutation operations (stage 3+)
// are NOT yet implemented; this layer establishes the visibility contract,
// single-active-transaction rule, revision-based adjacency invalidation,
// and the workspace hash pair - all per docs/M4B6_DESIGN.md.
//
// VISIBILITY CONTRACT (binding): while a TopologyTransaction is active,
// mutable workspace topology is exclusively owned by that transaction.
// Public workspace topology queries, audits, hashes and derived-view
// access REJECT with std::logic_error. Validators/diagnostic callbacks
// (stage 4) will receive preview data only - never pool references.
//
#include <bit>
#include <functional>
#include <memory>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <optional>
#include <stdexcept>

#include "paacs/topo/Adjacency.hpp"
#include "paacs/topo/TransactionAccess.hpp"

namespace paacs::topo {

// M4-B.8 STAGE A (disclosed frozen-file additive edit #1, managerial-
// waiver closure 2026-07-16): the atomic-cavity primitive lives in the
// NEW unfrozen Cavity.hpp; this frozen file takes only the forward
// declaration, the friend, and the move-only factory declaration, so the
// existing single-active-transaction and workspace-state rules remain
// authoritative for cavity transactions too.
class CavityTransaction;
// BND-2A / ASCR (disclosed frozen-file ADDITIVE edit, rev-3 §0.5.3): the
// atomic Steiner-cavity transaction lives in the NEW unfrozen
// SteinerCavity.hpp; this frozen file takes only the forward declaration,
// the friend and the move-only factory declaration, exactly as it did for
// CavityTransaction, so the single-active-transaction and workspace-state
// rules stay authoritative for Steiner cavities too. No existing signature
// or behaviour changes.
class SteinerCavityTransaction;
// T2-A accepted revision-2 contact: the unfrozen TetCavity.hpp owns the
// complete Tet4 net-delta transaction. This frozen file contributes only the
// forward declaration, friend and factory declaration.
class TetCavityTransaction;
// T3-A accepted revision-2 contact: the unfrozen TetRefinement.hpp owns the
// staged-node Tet4 net-delta transaction. This file contributes only the
// forward declaration, friend and factory declaration.
class TetRefinementTransaction;

enum class WorkspaceState {
    Ready,                    // committed, adjacency current: all reads ok
    TransactionStaging,       // Phase 1 active: preview APIs only
    TransactionPrepared,      // Phase 2 done, placeholders live:
                              //   only commit() or abandon() legal
    CommittedAdjacencyStale   // committed; ordinary adjacency REJECTS
};

// NonManifoldPolicy / TopologyPolicy are defined in TransactionAccess.hpp
// (included above) so Phase-1 planning can consume them directly.

class TopologyTransaction;

// ELEVATED BINDING RULE (review, after the FOURTH recurrence of the
// readiness-vs-coverage defect class): every TopologyWorkspace state -
// initially adopted or transactionally committed - must have complete
// edge incidence for every cell with edges AND complete face incidence
// for every cell with faces. This is a constructor invariant, a
// post-commit invariant (enforced when stage 3 commits mutations), and
// the named audit condition below.
enum class IncidenceCoverage { Complete, MissingEdges, MissingFaces };

inline IncidenceCoverage workspaceIncidenceCoverage(
    const NodePool& nodes, const CellPool& cells, const EdgePool& edges,
    const FacePool& faces) {
    (void)nodes;
    for (const CellHandle& c : cells.activeBySequence()) {
        const auto cell = cells.get(c);
        if (!isValidElementType(cell.type)) continue;   // storage audits
        const ElementTraits& tr = checkedTraits(cell.type);
        if (tr.edgeCount > 0 && !edges.hasCellEdges(cells, c))
            return IncidenceCoverage::MissingEdges;
        if (tr.faceCount > 0 && !faces.hasCellFaces(cells, c))
            return IncidenceCoverage::MissingFaces;
    }
    return IncidenceCoverage::Complete;
}
inline bool workspaceIncidenceComplete(
    const NodePool& nodes, const CellPool& cells, const EdgePool& edges,
    const FacePool& faces) {
    return workspaceIncidenceCoverage(nodes, cells, edges, faces) ==
           IncidenceCoverage::Complete;
}

// THE composite readiness gate (reviewer-required single definition so
// the rules cannot drift apart): a coordinated topology state is ready
// only when it is internally valid, semantically consistent, and
// complete for every relevant current entity. Used for workspace
// construction, post-commit assertion, hash eligibility, and adjacency
// rebuilding (and, later, Phase-1 preview validation and finalization).
inline bool workspaceTopologyReady(
    const NodePool& nodes, const CellPool& cells, const EdgePool& edges,
    const FacePool& faces) {
    return storageEdgeReady(nodes, cells) &&
           edgeTopologyReady(nodes, cells, edges) &&
           faceTopologyReady(nodes, cells, edges, faces) &&
           workspaceIncidenceComplete(nodes, cells, edges, faces);
}

class TopologyWorkspace {
public:
    // Adopt fully READY, incidence-COMPLETE topology (the M4-B.5 rule:
    // readiness = internal consistency + complete incidence coverage).
    TopologyWorkspace(NodePool&& n, CellPool&& c, EdgePool&& e,
                      FacePool&& f)
        : nodes_(std::move(n)), cells_(std::move(c)),
          edges_(std::move(e)), faces_(std::move(f)) {
        // Constructor invariant (elevated rule): complete EDGE incidence
        // for every cell with edges, complete FACE incidence for every
        // cell with faces. CellAdjacency::build below rechecks the 3D
        // face side; the explicit gate here also covers 1D/2D cells,
        // which the adjacency-level rule intentionally skips (review:
        // the workspace previously adopted a live Segment2 with no edge
        // incidence - the fourth readiness-vs-coverage recurrence).
        switch (workspaceIncidenceCoverage(nodes_, cells_, edges_, faces_)) {
            case IncidenceCoverage::Complete: break;
            case IncidenceCoverage::MissingEdges:
                throw std::invalid_argument(
                    "TopologyWorkspace: missing cell-edge incidence for a "
                    "current cell with edges");
            case IncidenceCoverage::MissingFaces:
                throw std::invalid_argument(
                    "TopologyWorkspace: missing cell-face incidence for a "
                    "current cell with faces");
        }
        if (!workspaceTopologyReady(nodes_, cells_, edges_, faces_))
            throw std::invalid_argument(
                "TopologyWorkspace: topology is not ready (storage, edge, "
                "or face audit failed)");
        // BND-2A Increment-5 reviewed frozen amendment (2026-07-31):
        // default-emplace the destination, then let CellAdjacency perform a
        // complete scratch build and noexcept swap. Passing a returned
        // adjacency into optional::emplace can invoke the retained public
        // noexcept move, whose standard-container move allocates on MSVC.
        adj_.emplace();
        adj_->rebuild(nodes_, cells_, edges_, faces_);
        adjBuiltFrom_ = revision_;
    }
    TopologyWorkspace(const TopologyWorkspace&) = delete;
    TopologyWorkspace& operator=(const TopologyWorkspace&) = delete;
    TopologyWorkspace(TopologyWorkspace&&) = delete;
    TopologyWorkspace& operator=(TopologyWorkspace&&) = delete;

    WorkspaceState state() const noexcept { return state_; }
    std::uint64_t topologyRevision() const noexcept { return revision_; }

    // Gated const topology access (visibility contract).
    const NodePool& nodes() const { requireReadable("nodes"); return nodes_; }
    const CellPool& cells() const { requireReadable("cells"); return cells_; }
    const EdgePool& edges() const { requireReadable("edges"); return edges_; }
    const FacePool& faces() const { requireReadable("faces"); return faces_; }

    // Revision-based adjacency semantics (design amendment 1):
    //   adjacencyReady = hasValue && builtFromRevision == topologyRevision
    // Ordinary adjacency access REJECTS while stale or during a
    // transaction; debugStaleAdjacency() is the explicit diagnostic
    // accessor (readable states only).
    std::uint64_t adjacencyBuiltFromRevision() const noexcept {
        return adjBuiltFrom_;
    }
    bool adjacencyReady() const noexcept {
        return readable() && adj_.has_value() && adjBuiltFrom_ == revision_;
    }
    const CellAdjacency& adjacency() const {
        requireReadable("adjacency");
        if (!adj_.has_value() || adjBuiltFrom_ != revision_)
            throw std::logic_error(
                "TopologyWorkspace::adjacency: snapshot is stale; call "
                "rebuildAdjacency() (debugStaleAdjacency() may inspect it)");
        return *adj_;
    }
    const CellAdjacency& debugStaleAdjacency() const {
        requireReadable("debugStaleAdjacency");
        if (!adj_.has_value())
            throw std::logic_error("no adjacency snapshot exists");
        return *adj_;
    }
    void rebuildAdjacency() {
        requireReadable("rebuildAdjacency");
        if (!workspaceTopologyReady(nodes_, cells_, edges_, faces_))
            throw std::logic_error(
                "TopologyWorkspace::rebuildAdjacency: workspace topology "
                "is not ready");
        if (!adj_.has_value())
            throw std::logic_error(
                "TopologyWorkspace::rebuildAdjacency: missing retained "
                "adjacency snapshot");
        // Strong guarantee: CellAdjacency::rebuild leaves the old snapshot
        // unchanged until its complete scratch snapshot is ready. Therefore
        // state_, revision_, and adjBuiltFrom_ also remain untouched on throw.
        adj_->rebuild(nodes_, cells_, edges_, faces_);
        adjBuiltFrom_ = revision_;
        if (state_ == WorkspaceState::CommittedAdjacencyStale)
            state_ = WorkspaceState::Ready;
    }

    // Workspace hash pair (PAACS_WORKSPACE_HASH_V1, per design):
    // semantic = logical live content (node seq+payload+cls; cell seq+
    // type+connectivity-as-node-seqs+cls; edge canonical seq key+cls+
    // constraint+useCount; CELL-EDGE INCIDENCE hashed from the cell side
    // in sequence order - per cell, the ordered cellEdges() relation as
    // canonical endpoint node-seq pairs (review: endpoints+useCount alone
    // could not distinguish reassigned incidence, so the hash could not
    // prove no ghost incidence survived a failed transaction); face
    // sorted corner seqs + (cellSeq,sense) incidence + cls). Excludes
    // slots, generations, pool ids, capacities, derived adjacency.
    // strict = semantic content PLUS per-entity (slot,gen) in iteration
    // (sequence) order, the (slot,gen) of each incident EdgeHandle, AND
    // all four pools' totalCreated counters (review: so consuming a
    // production sequence during failed transaction work is DETECTED even
    // if the entity was removed again). Excludes capacities, free-list
    // layout, allocator addresses, derived adjacency.
    // TRANSITION PERIOD (M4-B.7 Stage H): V1 and V2 are both exposed by
    // explicit name; the public alias switches to V2 only after the
    // rank-invariance + sense-normalization proofs, corruption-
    // sensitivity requalification, the full transaction evidence rerun
    // and baseline requalification (design amendment 2026-07-14).
    std::uint64_t semanticTopologyHashV1() const {
        requireReadable("semanticTopologyHashV1");
        requireHashEligible();
        return hashImpl(false);
    }
    // PAACS_WORKSPACE_HASH_V2 (M4-B.7 hash milestone): a PORTABLE
    // semantic mesh hash over deterministic live-order RANKS - invariant
    // to historical sequence gaps while preserving live sequence order
    // (NOT order-independent). Edges/faces enumerate in canonical-key
    // sorted order; face incidence senses are NORMALIZED to the
    // lowest-cell-rank incidence (stored face orientation is
    // identity-domain data that flips when a different cell recreates
    // the face). Excludes pool ids, slots, generations, capacities,
    // historical sequence values/gaps and totalCreated.
    std::uint64_t semanticTopologyHashV2() const {
        requireReadable("semanticTopologyHashV2");
        requireHashEligible();
        return hashImplV2();
    }
    // ALIAS SWITCHED TO V2 (Stage H, after rank-invariance and
    // sense-normalization proofs passed): the public semantic hash is
    // PAACS_WORKSPACE_HASH_V2; V1 remains available above as an
    // explicitly-named legacy diagnostic pending reviewer removal.
    std::uint64_t semanticTopologyHash() const {
        requireReadable("semanticTopologyHash");
        requireHashEligible();
        return hashImplV2();
    }
    std::uint64_t strictWorkspaceHash() const {
        requireReadable("strictWorkspaceHash");
        requireHashEligible();
        return hashImpl(true);
    }

    TopologyTransaction beginTransaction(const TopologyPolicy& p = {});

#if defined(PAACS_TOPO_TESTING)
    // Instrument-test seams (production-absent, seam-checked): mutable
    // edge-pool access plus UNGATED hash access, so hash-instrument
    // corruption testing can perturb COMPLETE incidence and prove the
    // instrument's sensitivity - while the PUBLIC hash accessors remain
    // gated on workspaceTopologyReady and reject the same state.
    EdgePool& testEdgesMutable() { return edges_; }
    std::vector<std::uint64_t> testV2Tokens() const {
        std::vector<std::uint64_t> v;
        emitV2Stream([&](std::uint64_t t) { v.push_back(t); });
        return v;
    }
    std::uint64_t testHashUngated(bool strict) const {
        // Semantic component follows the public alias (V2).
        return strict ? hashImpl(true) : hashImplV2();
    }
#endif

    // M4-B.8 Stage A: cavity transactions share the state machine
    // (defined in Cavity.hpp; declared here so the factory is a member).
    CavityTransaction beginCavityTransaction();

    // BND-2A / ASCR: the atomic Steiner-cavity transaction shares the same
    // state machine and single-active-transaction rule (defined in
    // SteinerCavity.hpp; declared here so the factory is a member).
    SteinerCavityTransaction beginSteinerCavityTransaction();

    // T2-A: defined in TetCavity.hpp; shares the existing state machine.
    TetCavityTransaction beginTetCavityTransaction();

    // T3-A: defined in TetRefinement.hpp; shares the existing state machine.
    TetRefinementTransaction beginTetRefinementTransaction();

private:
    friend class TopologyTransaction;
    friend class CavityTransaction;   // M4-B.8 Stage A (disclosed edit #1)
    friend class SteinerCavityTransaction;   // BND-2A (disclosed additive)
    friend class TetCavityTransaction;       // T2-A accepted additive contact
    friend class TetRefinementTransaction;   // T3-A accepted additive contact

    bool readable() const noexcept {
        return state_ == WorkspaceState::Ready ||
               state_ == WorkspaceState::CommittedAdjacencyStale;
    }
    // Hash eligibility = the composite gate (reviewer): hashing a
    // non-ready workspace would let a corrupt state masquerade as a
    // reproducible topology. The corruption-SENSITIVITY proof therefore
    // runs through the ungated testing seam below (hash-instrument
    // corruption testing), never through this public path.
    void requireHashEligible() const {
        if (!workspaceTopologyReady(nodes_, cells_, edges_, faces_))
            throw std::logic_error(
                "TopologyWorkspace: hash requested on a workspace whose "
                "topology is not ready");
    }

    void requireReadable(const char* what) const {
        if (!readable())
            throw std::logic_error(std::string(
                "TopologyWorkspace: '") + what +
                "' rejected while a transaction is active (visibility "
                "contract: mutable topology is exclusively owned by the "
                "active TopologyTransaction)");
    }

    static void mix(std::uint64_t& h, std::uint64_t v) noexcept {
        for (int b = 0; b < 8; ++b) {
            h ^= (v >> (8 * b)) & 0xFF;
            h *= 1099511628211ULL;
        }
    }
    static std::uint64_t dbits(double d) noexcept {
        return std::bit_cast<std::uint64_t>(d);
    }
    // PAACS_WORKSPACE_HASH_V2 token-stream schema (FIXED and versioned;
    // review hardening 2026-07-14). Any change to this framing requires
    // an intentional hash-version transition.
    //
    //   "PAACS_WORKSPACE_HASH_V2" (chars)
    //   SECTION_NODES  node_count
    //     { NODE x_bits y_bits z_bits class_dim class_tag } *
    //   SECTION_CELLS  cell_count
    //     { CELL element_type connectivity_count node_rank* class_dim
    //       class_tag } *
    //   SECTION_EDGES  edge_count            (canonical (rlo,rhi) order)
    //     { EDGE rank_lo rank_hi class_dim class_tag
    //       required_constraint } *
    //   SECTION_CELL_EDGE_INCIDENCE  cell_count
    //     { CELL_EDGES cell_rank edge_count { rank_lo rank_hi }* } *
    //       (edge_count is 0xFFFFFFFF when the cell has NO incidence
    //        entry - distinct from an empty one)
    //   SECTION_FACES  face_count            (canonical corner-key order)
    //     { FACE corner_count corner_rank* class_dim class_tag
    //       incidence_count { cell_rank normalized_sense }* } *
    static constexpr std::uint64_t kV2SectionNodes = 0x5EC0000000000001ULL;
    static constexpr std::uint64_t kV2SectionCells = 0x5EC0000000000002ULL;
    static constexpr std::uint64_t kV2SectionEdges = 0x5EC0000000000003ULL;
    static constexpr std::uint64_t kV2SectionCellEdgeInc =
        0x5EC0000000000004ULL;
    static constexpr std::uint64_t kV2SectionFaces = 0x5EC0000000000005ULL;
    static constexpr std::uint64_t kV2RecNode = 0x4E0DE00000000001ULL;
    static constexpr std::uint64_t kV2RecCell = 0x4E0DE00000000002ULL;
    static constexpr std::uint64_t kV2RecEdge = 0x4E0DE00000000003ULL;
    static constexpr std::uint64_t kV2RecCellEdges =
        0x4E0DE00000000004ULL;
    static constexpr std::uint64_t kV2RecFace = 0x4E0DE00000000005ULL;
    static constexpr std::uint64_t kV2NoIncidence = 0xFFFFFFFFULL;
    template <class Sink>
    void emitV2Stream(Sink&& sink) const {
        for (const char* p = "PAACS_WORKSPACE_HASH_V2"; *p; ++p)
            sink(static_cast<std::uint64_t>(*p));
        std::unordered_map<std::int32_t, std::uint64_t> nrank, crank;
        {
            std::uint64_t r = 0;
            for (const NodeHandle& nh : nodes_.activeBySequence())
                nrank.emplace(nh.slot, r++);
        }
        {
            std::uint64_t r = 0;
            for (const CellHandle& ch : cells_.activeBySequence())
                crank.emplace(ch.slot, r++);
        }
        sink(kV2SectionNodes);
        sink(nrank.size());
        for (const NodeHandle& nh : nodes_.activeBySequence()) {
            const auto d = nodes_.get(nh);
            sink(kV2RecNode);
            sink(dbits(d.coord.x)); sink(dbits(d.coord.y));
            sink(dbits(d.coord.z));
            sink(static_cast<std::uint64_t>(d.cls.dim));
            sink(static_cast<std::uint64_t>(d.cls.tag));
        }
        sink(kV2SectionCells);
        sink(crank.size());
        for (const CellHandle& ch : cells_.activeBySequence()) {
            const auto sn = cells_.get(ch);
            sink(kV2RecCell);
            sink(static_cast<std::uint64_t>(sn.type));
            sink(static_cast<std::uint64_t>(sn.count));
            for (std::uint8_t i = 0; i < sn.count; ++i)
                sink(nrank.at(sn.nodes[i].slot));
            sink(static_cast<std::uint64_t>(sn.cls.dim));
            sink(static_cast<std::uint64_t>(sn.cls.tag));
        }
        sink(kV2SectionEdges);
        {
            std::vector<std::array<std::uint64_t, 5>> es;
            for (const EdgeHandle& eh : edges_.activeBySequence()) {
                const auto d = edges_.get(eh);
                es.push_back({nrank.at(d.a.slot), nrank.at(d.b.slot),
                              static_cast<std::uint64_t>(d.cls.dim),
                              static_cast<std::uint64_t>(d.cls.tag),
                              d.requiredConstraint ? 1ULL : 0ULL});
            }
            std::sort(es.begin(), es.end());
            sink(es.size());
            for (const auto& e : es) {
                sink(kV2RecEdge);
                for (const auto v : e) sink(v);
            }
        }
        sink(kV2SectionCellEdgeInc);
        sink(crank.size());
        for (const CellHandle& ch : cells_.activeBySequence()) {
            sink(kV2RecCellEdges);
            sink(crank.at(ch.slot));
            if (!edges_.hasCellEdges(cells_, ch)) {
                sink(kV2NoIncidence);
                continue;
            }
            const auto ce = edges_.cellEdges(cells_, ch);
            sink(ce.size());
            for (const EdgeHandle& eh : ce) {
                const auto ed = edges_.get(eh);
                sink(nrank.at(ed.a.slot));
                sink(nrank.at(ed.b.slot));
            }
        }
        sink(kV2SectionFaces);
        {
            struct FEntry {
                std::array<std::uint64_t, 4> key{~0ULL, ~0ULL, ~0ULL,
                                                 ~0ULL};
                std::uint64_t n = 0, dim = 0, tag = 0;
                std::vector<std::pair<std::uint64_t, std::int64_t>> inc;
                bool operator<(const FEntry& o) const {
                    return std::tie(n, key) < std::tie(o.n, o.key);
                }
            };
            std::vector<FEntry> fs;
            for (const FaceHandle& fh : faces_.activeBySequence()) {
                const auto d = faces_.get(fh);
                FEntry e;
                e.n = d.cornerCount;
                for (std::uint8_t i = 0; i < d.cornerCount; ++i)
                    e.key[i] = nrank.at(d.corners[i].slot);
                std::sort(e.key.begin(), e.key.begin() + d.cornerCount);
                e.dim = static_cast<std::uint64_t>(d.cls.dim);
                e.tag = static_cast<std::uint64_t>(d.cls.tag);
                for (const auto& in : d.incident)
                    e.inc.emplace_back(crank.at(in.cellSlot),
                                       static_cast<std::int64_t>(in.sense));
                std::sort(e.inc.begin(), e.inc.end());
                if (!e.inc.empty()) {
                    const std::int64_t norm = e.inc.front().second;
                    if (norm != 0)
                        for (auto& p : e.inc) p.second *= norm;
                }
                fs.push_back(std::move(e));
            }
            std::sort(fs.begin(), fs.end());
            sink(fs.size());
            for (const auto& f : fs) {
                sink(kV2RecFace);
                sink(f.n);
                for (std::uint64_t i = 0; i < f.n; ++i) sink(f.key[i]);
                sink(f.dim); sink(f.tag);
                sink(f.inc.size());
                for (const auto& [cr, sn] : f.inc) {
                    sink(cr);
                    sink(static_cast<std::uint64_t>(sn));
                }
            }
        }
    }
    std::uint64_t hashImplV2() const {
        std::uint64_t h = 1469598103934665603ULL;
        emitV2Stream([&](std::uint64_t v) { mix(h, v); });
        return h;
    }
    std::uint64_t hashImpl(bool strict) const {
        std::uint64_t h = 1469598103934665603ULL;
        for (const char* p = "PAACS_WORKSPACE_HASH_V1"; *p; ++p)
            mix(h, static_cast<std::uint64_t>(*p));
        mix(h, strict ? 2 : 1);
        for (const NodeHandle& nh : nodes_.activeBySequence()) {
            const auto d = nodes_.get(nh);
            mix(h, nodes_.seqOf(nh));
            mix(h, dbits(d.coord.x)); mix(h, dbits(d.coord.y));
            mix(h, dbits(d.coord.z));
            mix(h, static_cast<std::uint64_t>(d.cls.dim));
            mix(h, static_cast<std::uint64_t>(d.cls.tag));
            if (strict) { mix(h, static_cast<std::uint64_t>(nh.slot));
                          mix(h, nh.gen); }
        }
        for (const CellHandle& ch : cells_.activeBySequence()) {
            const auto s = cells_.get(ch);
            mix(h, cells_.seqOf(ch));
            mix(h, static_cast<std::uint64_t>(s.type));
            mix(h, static_cast<std::uint64_t>(s.cls.dim));
            mix(h, static_cast<std::uint64_t>(s.cls.tag));
            for (std::uint8_t i = 0; i < s.count; ++i)
                mix(h, nodes_.seqOf(s.nodes[i]));
            if (strict) { mix(h, static_cast<std::uint64_t>(ch.slot));
                          mix(h, ch.gen); }
        }
        for (const EdgeHandle& eh : edges_.activeBySequence()) {
            const auto d = edges_.get(eh);
            mix(h, nodes_.seqOf(d.a)); mix(h, nodes_.seqOf(d.b));
            mix(h, static_cast<std::uint64_t>(d.cls.dim));
            mix(h, static_cast<std::uint64_t>(d.cls.tag));
            mix(h, d.requiredConstraint ? 1 : 0);
            mix(h, d.useCount);
            if (strict) { mix(h, static_cast<std::uint64_t>(eh.slot));
                          mix(h, eh.gen); }
        }
        // Cell-edge incidence, from the cell side, sequence order.
        for (const CellHandle& ch : cells_.activeBySequence()) {
            mix(h, cells_.seqOf(ch));
            if (!edges_.hasCellEdges(cells_, ch)) {
                mix(h, 0xEDCE0000ULL);       // no-incidence sentinel
                continue;
            }
            const auto ce = edges_.cellEdges(cells_, ch);
            mix(h, 0xEDCE0001ULL);
            mix(h, ce.size());
            for (const EdgeHandle& eh : ce) {
                const auto ed = edges_.get(eh);
                mix(h, nodes_.seqOf(ed.a));  // canonical lo->hi endpoint
                mix(h, nodes_.seqOf(ed.b));  // sequences = logical identity
                if (strict) {
                    mix(h, static_cast<std::uint64_t>(eh.slot));
                    mix(h, eh.gen);
                }
            }
        }
        for (const FaceHandle& fh : faces_.activeBySequence()) {
            const auto d = faces_.get(fh);
            mix(h, d.cornerCount);
            std::array<std::uint64_t, kMaxFaceNodes> cs{};
            for (std::uint8_t i = 0; i < d.cornerCount; ++i)
                cs[i] = nodes_.seqOf(d.corners[i]);
            std::sort(cs.begin(), cs.begin() + d.cornerCount);
            for (std::uint8_t i = 0; i < d.cornerCount; ++i) mix(h, cs[i]);
            mix(h, static_cast<std::uint64_t>(d.cls.dim));
            mix(h, static_cast<std::uint64_t>(d.cls.tag));
            for (const auto& in : d.incident) {
                mix(h, cells_.seqOf(CellHandle{in.cellSlot, in.cellGen,
                                               cells_.poolId()}));
                mix(h, static_cast<std::uint64_t>(
                           static_cast<std::int64_t>(in.sense)));
            }
            if (strict) { mix(h, static_cast<std::uint64_t>(fh.slot));
                          mix(h, fh.gen); }
        }
        if (strict) {
            // Per-pool creation counters: a consumed production sequence
            // is detected even when the entity was removed again.
            mix(h, nodes_.stats().totalCreated);
            mix(h, cells_.stats().totalCreated);
            mix(h, edges_.stats().totalCreated);
            mix(h, faces_.stats().totalCreated);
        }
        return h;
    }

    NodePool nodes_;
    CellPool cells_;
    EdgePool edges_;
    FacePool faces_;
    std::optional<CellAdjacency> adj_;
    std::uint64_t adjBuiltFrom_ = 0;
    std::uint64_t revision_ = 1;
    WorkspaceState state_ = WorkspaceState::Ready;
};

// Preview VALUES handed to caller validators (Phase 1 only): plain data
// derived from the staged plan - never workspace or pool references, so
// a validator cannot observe placeholder or mid-commit state.
struct TransactionPreview {
    enum class Kind { None, Add, Remove, Replace, Metadata, AddNode,
                      RemoveNode };
    Kind kind = Kind::None;
    CellHandle cell{};                   // Remove/Replace target; Add: {}
    ElementType type{};
    std::uint8_t nodeCount = 0;
    std::array<NodeHandle, 10> nodes{};  // Add/Replace connectivity
    std::size_t newEdgeRecords = 0;      // creates
    std::size_t reusedEdgeIncidence = 0; // incidence added to existing
    std::size_t retainedEdges = 0;       // replace: unchanged incidence
    std::size_t removedEdgeIncidence = 0;
    std::size_t retiringEdges = 0;
    std::size_t newFaceRecords = 0;
    std::size_t reusedFaceIncidence = 0;
    std::size_t retainedFaces = 0;
    std::size_t removedFaceIncidence = 0;
    std::size_t retiringFaces = 0;
    // Metadata (6B): exactly ONE named entity is modified - previews
    // state precisely what will change; nothing propagates implicitly.
    enum class MetaTarget { None, Node, Cell, Edge, Face, EdgeConstraint };
    MetaTarget metaTarget = MetaTarget::None;
    NodeHandle nodeTarget{};
    EdgeHandle edgeTarget{};
    FaceHandle faceTarget{};
    mesh::Classification oldClassification{}, newClassification{};
    bool oldConstraint = false, newConstraint = false;
    bool noOp = false;
    Vec3 nodeCoord{};                    // AddNode
};

// Move-only single-active transaction. Stage-3 scope: at most ONE cell
// operation per transaction (addCell OR removeCell; a second staging
// call is a logic_error - multi-operation transactions arrive with the
// cavity work in a later stage). Phases per docs/M4B6_DESIGN.md:
// staging call = Phase 1 (throwing validation + planning, pools
// untouched); prepare() = Phase 2 (throwing capacity + placeholders;
// state -> TransactionPrepared); commit() = Phase 3 (noexcept prepared
// primitives only; revision++, adjacency stale). abandon()/destructor
// erase placeholders and restore visibility.
class TopologyTransaction {
public:
    TopologyTransaction(const TopologyTransaction&) = delete;
    TopologyTransaction& operator=(const TopologyTransaction&) = delete;
    TopologyTransaction(TopologyTransaction&& o) noexcept
        : ws_(o.ws_), policy_(o.policy_), resume_(o.resume_),
          op_(o.op_), add_(std::move(o.add_)), rem_(std::move(o.rem_)),
          rep_(std::move(o.rep_)), meta_(std::move(o.meta_)),
          addNode_(std::move(o.addNode_)), remNode_(std::move(o.remNode_)),
          validatorRejected_(o.validatorRejected_),
          newCell_(o.newCell_) {
        o.ws_ = nullptr;
    }
    TopologyTransaction& operator=(TopologyTransaction&&) = delete;
    ~TopologyTransaction() { if (ws_) abandonImpl(); }

    bool active() const noexcept { return ws_ != nullptr; }
    const TopologyPolicy& policy() const {
        requireActive("policy");
        return policy_;
    }

    // Phase 1: validate + plan (throwing; pools untouched).
    CellHandle addCell(ElementType type, std::span<const NodeHandle> nodes,
                       const mesh::Classification& cls) {
        requireActive("addCell");
        requireUnstaged();
        add_ = std::make_unique<detail::TransactionAccess::AddPlan>(
            detail::TransactionAccess::planAddCell(
                ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_, type,
                nodes, cls, policy_));
        op_ = Op::Add;
        // The handle identity is only predictable, not allocated, until
        // commit; returned by commit() in the CommitResult.
        return CellHandle{};
    }
    void removeCell(const CellHandle& cell) {
        requireActive("removeCell");
        requireUnstaged();
        rem_ = std::make_unique<detail::TransactionAccess::RemovePlan>(
            detail::TransactionAccess::planRemoveCell(
                ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_, cell,
                policy_));
        op_ = Op::Remove;
    }
    // Phase 1: replace ONE current cell's connectivity, PRESERVING the
    // cell handle (same pool/slot/generation/sequence). NET deltas: the
    // plan retains unchanged edge/face incidence rather than removing
    // and recreating it.
    void replaceCellConnectivity(const CellHandle& cell,
                                 std::span<const NodeHandle> nodes) {
        requireActive("replaceCellConnectivity");
        requireUnstaged();
        rep_ = std::make_unique<detail::TransactionAccess::ReplacePlan>(
            detail::TransactionAccess::planReplaceCell(
                ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_, cell,
                nodes, policy_));
        op_ = Op::Replace;
    }

    // M4-B.6C node lifecycle (Phase 1). Creation is rejected under the
    // default permitOrphanNodes == false (a new node is necessarily
    // isolated in a one-operation transaction); removal structurally
    // rejects any node with a surviving reference (cell connectivity,
    // live edge endpoint, live face corner) and never cascades.
    void addNode(const Vec3& coord, const mesh::Classification& cls) {
        requireActive("addNode");
        requireUnstaged();
        addNode_ = std::make_unique<detail::TransactionAccess::AddNodePlan>(
            detail::TransactionAccess::planAddNode(ws_->nodes_, coord, cls,
                                                   policy_));
        op_ = Op::AddNode;
    }
    void removeNode(const NodeHandle& h) {
        requireActive("removeNode");
        requireUnstaged();
        remNode_ =
            std::make_unique<detail::TransactionAccess::RemoveNodePlan>(
                detail::TransactionAccess::planRemoveNode(
                    ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_, h));
        op_ = Op::RemoveNode;
    }

    // M4-B.6B metadata staging (Phase 1; existing compact classification
    // identifiers only; NO implicit propagation).
    void setNodeClassification(const NodeHandle& h,
                               const mesh::Classification& c) {
        requireActive("setNodeClassification");
        requireUnstaged();
        meta_ = std::make_unique<detail::TransactionAccess::MetadataPlan>(
            detail::TransactionAccess::planSetNodeCls(ws_->nodes_, h, c));
        op_ = Op::Metadata;
    }
    void setCellClassification(const CellHandle& h,
                               const mesh::Classification& c) {
        requireActive("setCellClassification");
        requireUnstaged();
        meta_ = std::make_unique<detail::TransactionAccess::MetadataPlan>(
            detail::TransactionAccess::planSetCellCls(ws_->cells_, h, c));
        op_ = Op::Metadata;
    }
    void setEdgeClassification(const EdgeHandle& h,
                               const mesh::Classification& c) {
        requireActive("setEdgeClassification");
        requireUnstaged();
        meta_ = std::make_unique<detail::TransactionAccess::MetadataPlan>(
            detail::TransactionAccess::planSetEdgeCls(ws_->edges_, h, c));
        op_ = Op::Metadata;
    }
    void setFaceClassification(const FaceHandle& h,
                               const mesh::Classification& c) {
        requireActive("setFaceClassification");
        requireUnstaged();
        meta_ = std::make_unique<detail::TransactionAccess::MetadataPlan>(
            detail::TransactionAccess::planSetFaceCls(ws_->faces_, h, c));
        op_ = Op::Metadata;
    }
    void setEdgeRequiredConstraint(const EdgeHandle& h, bool required) {
        requireActive("setEdgeRequiredConstraint");
        requireUnstaged();
        meta_ = std::make_unique<detail::TransactionAccess::MetadataPlan>(
            detail::TransactionAccess::planSetEdgeConstraint(ws_->edges_, h,
                                                             required));
        op_ = Op::Metadata;
    }

    // Caller validator (Phase 1 only, BEFORE prepare): receives preview
    // VALUES derived from the staged plan - never workspace or pool
    // references. Returning false marks the transaction rejected; commit
    // then refuses. Never invoked by prepare() or commit().
    bool validate(
        const std::function<bool(const TransactionPreview&)>& fn) {
        requireActive("validate");
        if (op_ == Op::None)
            throw std::logic_error(
                "TopologyTransaction::validate: nothing staged");
        if (ws_->state_ != WorkspaceState::TransactionStaging)
            throw std::logic_error(
                "TopologyTransaction::validate: validators run in Phase 1 "
                "only (before prepare)");
        const bool ok = fn(preview());
        if (!ok) validatorRejected_ = true;
        return ok;
    }
    TransactionPreview preview() const {
        requireActive("preview");
        TransactionPreview v;
        using RK = detail::TransactionAccess::RelKind;
        if (op_ == Op::Add) {
            v.kind = TransactionPreview::Kind::Add;
            v.type = add_->type;
            v.nodeCount = add_->nodeCount;
            v.nodes = add_->nodes;
            for (const auto& op : add_->edges)
                (op.create ? v.newEdgeRecords : v.reusedEdgeIncidence)++;
            for (const auto& op : add_->faces)
                (op.create ? v.newFaceRecords : v.reusedFaceIncidence)++;
        } else if (op_ == Op::Remove) {
            v.kind = TransactionPreview::Kind::Remove;
            v.cell = rem_->cell;
            v.removedEdgeIncidence = rem_->edges.size();
            for (const auto& e : rem_->edges)
                v.retiringEdges += e.retire ? 1 : 0;
            v.removedFaceIncidence = rem_->faces.size();
            for (const auto& f : rem_->faces)
                v.retiringFaces += f.retire ? 1 : 0;
        } else if (op_ == Op::AddNode) {
            v.kind = TransactionPreview::Kind::AddNode;
            v.nodeCoord = addNode_->coord;
            v.newClassification = addNode_->cls;
        } else if (op_ == Op::RemoveNode) {
            v.kind = TransactionPreview::Kind::RemoveNode;
            v.nodeTarget = remNode_->node;
        } else if (op_ == Op::Metadata) {
            using MK = detail::TransactionAccess::MetaKind;
            v.kind = TransactionPreview::Kind::Metadata;
            v.noOp = meta_->noOp;
            switch (meta_->kind) {
                case MK::NodeCls:
                    v.metaTarget = TransactionPreview::MetaTarget::Node;
                    v.nodeTarget = meta_->node;
                    v.oldClassification = meta_->oldCls;
                    v.newClassification = meta_->newCls;
                    break;
                case MK::CellCls:
                    v.metaTarget = TransactionPreview::MetaTarget::Cell;
                    v.cell = meta_->cell;
                    v.oldClassification = meta_->oldCls;
                    v.newClassification = meta_->newCls;
                    break;
                case MK::EdgeCls:
                    v.metaTarget = TransactionPreview::MetaTarget::Edge;
                    v.edgeTarget = meta_->edge;
                    v.oldClassification = meta_->oldCls;
                    v.newClassification = meta_->newCls;
                    break;
                case MK::FaceCls:
                    v.metaTarget = TransactionPreview::MetaTarget::Face;
                    v.faceTarget = meta_->face;
                    v.oldClassification = meta_->oldCls;
                    v.newClassification = meta_->newCls;
                    break;
                case MK::EdgeConstraint:
                    v.metaTarget =
                        TransactionPreview::MetaTarget::EdgeConstraint;
                    v.edgeTarget = meta_->edge;
                    v.oldConstraint = meta_->oldConstraint;
                    v.newConstraint = meta_->newConstraint;
                    break;
            }
        } else if (op_ == Op::Replace) {
            v.kind = TransactionPreview::Kind::Replace;
            v.cell = rep_->cell;
            v.type = rep_->type;
            v.nodeCount = rep_->nodeCount;
            v.nodes = rep_->newNodes;
            for (const auto& op : rep_->edges) {
                if (op.kind == RK::AddCreate) v.newEdgeRecords++;
                else if (op.kind == RK::AddReuse) v.reusedEdgeIncidence++;
                else v.retainedEdges++;
            }
            v.removedEdgeIncidence = rep_->removedEdges.size();
            for (const auto& e : rep_->removedEdges)
                v.retiringEdges += e.retire ? 1 : 0;
            for (const auto& op : rep_->faces) {
                if (op.kind == RK::AddCreate) v.newFaceRecords++;
                else if (op.kind == RK::AddReuse) v.reusedFaceIncidence++;
                else v.retainedFaces++;
            }
            v.removedFaceIncidence = rep_->removedFaces.size();
            for (const auto& f : rep_->removedFaces)
                v.retiringFaces += f.retire ? 1 : 0;
        }
        return v;
    }

    // Phase 2: allocation preparation (throwing). On failure, placeholders
    // are erased and the workspace is semantically unchanged (capacity-
    // class residue only).
    void prepare() {
        requireActive("prepare");
        if (ws_->state_ == WorkspaceState::TransactionPrepared) return;
        if (op_ == Op::Add) {
            try {
                detail::TransactionAccess::prepareAddCell(
                    ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_,
                    *add_);
            } catch (...) {
                detail::TransactionAccess::abandonAddPrepared(
                    ws_->edges_, ws_->faces_, *add_);
                throw;
            }
        } else if (op_ == Op::Remove) {
            detail::TransactionAccess::prepareRemoveCell(
                ws_->cells_, ws_->edges_, ws_->faces_, *rem_);
        } else if (op_ == Op::AddNode) {
            detail::TransactionAccess::prepareAddNode(ws_->nodes_,
                                                      *addNode_);
        } else if (op_ == Op::RemoveNode) {
            detail::TransactionAccess::prepareRemoveNode(ws_->nodes_,
                                                         *remNode_);
        } else if (op_ == Op::Metadata) {
            detail::TransactionAccess::prepareMetadata(*meta_);
        } else if (op_ == Op::Replace) {
            try {
                detail::TransactionAccess::prepareReplaceCell(
                    ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_,
                    *rep_);
            } catch (...) {
                detail::TransactionAccess::abandonReplacePrepared(
                    ws_->edges_, ws_->faces_, *rep_);
                throw;
            }
        }
        ws_->state_ = WorkspaceState::TransactionPrepared;
    }

    struct CommitResult {
        bool topologyChanged = false;
        std::uint64_t revision = 0;
        CellHandle cell{};                   // add/replace: the cell handle
        NodeHandle node{};                   // addNode: the new node handle
    };
    // Phase 3 (auto-prepares if needed). Empty transactions restore the
    // prior readable state; mutating transactions increment the revision
    // and leave the workspace CommittedAdjacencyStale.
    CommitResult commit() {
        requireActive("commit");
        CommitResult r;
        if (op_ == Op::None) {
            r.topologyChanged = false;
            r.revision = ws_->revision_;
            abandonImpl();
            return r;
        }
        if (validatorRejected_)
            throw std::logic_error(
                "TopologyTransaction::commit: a validator rejected this "
                "transaction");
        if (op_ == Op::Metadata && meta_->noOp) {
            // No-op metadata (new value == current): no revision churn,
            // adjacency stays ready (reviewer rule).
            r.topologyChanged = false;
            r.revision = ws_->revision_;
            abandonImpl();
            return r;
        }
        prepare();                           // Phase 1+2 if not yet done
#if defined(PAACS_TOPO_TESTING)
        detail::phase3Active() = true;
#endif
        if (op_ == Op::Add) {
            r.cell = detail::TransactionAccess::commitAddCell(
                ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_, *add_);
        } else if (op_ == Op::Remove) {
            detail::TransactionAccess::commitRemoveCell(
                ws_->cells_, ws_->edges_, ws_->faces_, *rem_);
        } else if (op_ == Op::AddNode) {
            r.node = detail::TransactionAccess::commitAddNode(ws_->nodes_,
                                                              *addNode_);
        } else if (op_ == Op::RemoveNode) {
            detail::TransactionAccess::commitRemoveNode(ws_->nodes_,
                                                        *remNode_);
        } else if (op_ == Op::Metadata) {
            detail::TransactionAccess::commitMetadata(
                ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_,
                *meta_);
        } else {
            detail::TransactionAccess::commitReplaceCell(
                ws_->nodes_, ws_->cells_, ws_->edges_, ws_->faces_, *rep_);
            r.cell = rep_->cell;             // PRESERVED handle
        }
        ++ws_->revision_;
        ws_->state_ = WorkspaceState::CommittedAdjacencyStale;
#if defined(PAACS_TOPO_TESTING)
        detail::phase3Active() = false;
        // Post-commit invariant (elevated rule + composite gate).
        if (!workspaceTopologyReady(ws_->nodes_, ws_->cells_, ws_->edges_,
                                    ws_->faces_))
            std::abort();
#endif
        r.topologyChanged = true;
        r.revision = ws_->revision_;
        ws_ = nullptr;                       // consumed; no reuse
        return r;
    }

    void abandon() {
        requireActive("abandon");
        abandonImpl();
    }

private:
    friend class TopologyWorkspace;
    enum class Op { None, Add, Remove, Replace, Metadata, AddNode,
                    RemoveNode };

    TopologyTransaction(TopologyWorkspace& ws, const TopologyPolicy& p)
        : ws_(&ws), policy_(p), resume_(ws.state_) {
        ws_->state_ = WorkspaceState::TransactionStaging;
    }
    void requireActive(const char* what) const {
        if (!ws_)
            throw std::logic_error(std::string("TopologyTransaction: '") +
                what + "' on a finished/abandoned transaction (a "
                "transaction cannot be reused after commit/abandon)");
    }
    void requireUnstaged() const {
        if (op_ != Op::None)
            throw std::logic_error(
                "TopologyTransaction: v1 permits exactly one cell "
                "operation per transaction");
    }
    void abandonImpl() noexcept {
        if (op_ == Op::Add && add_)
            detail::TransactionAccess::abandonAddPrepared(
                ws_->edges_, ws_->faces_, *add_);
        if (op_ == Op::Replace && rep_)
            detail::TransactionAccess::abandonReplacePrepared(
                ws_->edges_, ws_->faces_, *rep_);
        // Remove-prepare leaves only free-list capacity: nothing to erase.
        ws_->state_ = resume_;
        ws_ = nullptr;
    }

    TopologyWorkspace* ws_ = nullptr;
    TopologyPolicy policy_{};
    WorkspaceState resume_ = WorkspaceState::Ready;
    Op op_ = Op::None;
    std::unique_ptr<detail::TransactionAccess::AddPlan> add_;
    std::unique_ptr<detail::TransactionAccess::RemovePlan> rem_;
    std::unique_ptr<detail::TransactionAccess::ReplacePlan> rep_;
    std::unique_ptr<detail::TransactionAccess::MetadataPlan> meta_;
    std::unique_ptr<detail::TransactionAccess::AddNodePlan> addNode_;
    std::unique_ptr<detail::TransactionAccess::RemoveNodePlan> remNode_;
    bool validatorRejected_ = false;
    CellHandle newCell_{};
};

inline TopologyTransaction TopologyWorkspace::beginTransaction(
    const TopologyPolicy& p) {
    if (!readable())
        throw std::logic_error(
            "TopologyWorkspace::beginTransaction: a transaction is already "
            "active (exactly one transaction may be active)");
    return TopologyTransaction(*this, p);
}

} // namespace paacs::topo
