// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// M4-B.8 — Atomic Cavity Replacement.
// Design authority: docs/M4B8_DESIGN.md (rev 3.9, accepted baseline).
//
// STAGE A: ACCEPTED/FROZEN 2026-07-17 - the value types, the move-only
// CavityTransaction, beginCavityTransaction() + shared state machine,
// the full Phase-1 combinatorial staging engine (normative rules incl.
// node-link audit, nonempty boundary, orientation-insensitive
// surviving-duplicate keys, composite gate first, unconditional
// incidence guard, closure-only requirements), canonicalization, the
// plain-value preview and the frozen-pattern sticky validator. The
// Stage-A surface is FROZEN; Stages B/C extend this file with NEW
// members under their own reviews.
//
// M4-B.8 ACCEPTED/CLOSED/FROZEN 2026-07-19 (Stages A/B/C). This file
// is FROZEN in full: the Phase-1 combinatorial staging + validator, the
// Phase-2 plan/preparation/TransactionPrepared, the semantic-no-op and
// substantive commits (revision exactly once; CommittedAdjacencyStale;
// commitCavityPrepared touches prepared pool state only). Changes
// require a reopening review.
//
// M4-B.8 owns the COMBINATORIAL cavity contract only (design §1):
// geometric embedding (nondegeneracy, orientation in the parameter
// domain, T-junctions, crossings) belongs to BND-1, reached through the
// bool-returning validator over the plain-value CavityPreview.

#include "paacs/topo/Workspace.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace paacs::topo {

// ----------------------------------------------------------- value types ---
struct TriCavityCell {                    // Tri3-specific by design (§8):
    std::array<NodeHandle, 3> nodes;      // existing, current nodes only
    mesh::Classification cls;
};

struct CavityEdgeRequirement {            // ENSURE-ONLY (§3): "the final
    NodeHandle a, b;                      // edge must be required and
    mesh::Classification cls;             // carry this classification"
};

struct TriCavitySpec {
    std::vector<CellHandle> removed;
    std::vector<TriCavityCell> added;
    std::vector<CavityEdgeRequirement> requirements;
};

// Plain values ONLY (§1): a validator can never observe workspace/pool
// or mid-commit state — the frozen TransactionPreview rule.
struct CavityPreview {
    std::vector<NodeHandle> boundaryCycle;   // canonical (§2.1)
    std::vector<NodeHandle> interiorNodes;   // sorted by live rank
    std::vector<CellHandle> removedCells;    // canonical (live rank)
    std::vector<std::array<NodeHandle, 3>> removedConnectivity;
    std::vector<std::array<NodeHandle, 3>> addedConnectivity; // canonical
    std::vector<mesh::Classification> addedClass;
    struct EdgeRow {                      // one per canonical cavity edge
        NodeHandle a, b;                  // rank(a) < rank(b)
        std::uint32_t initialUse = 0, finalUse = 0;
        std::int8_t insideSenseRemoved = 0;   // boundary edges only
        std::int8_t insideSenseAdded = 0;
        bool wasRequired = false;
        bool willBeRequired = false;
        // Stage-C ADDITIVE field (disclosed): true iff the SPEC requested
        // this edge - CavityCommit::requiredEdges must list exactly the
        // requested set, which willBeRequired (wasRequired || requested)
        // cannot recover.
        bool requested = false;
        mesh::Classification initialCls;
        mesh::Classification finalCls;
    };
    std::vector<EdgeRow> edges;
    // Plan counts (plain numbers; Phase-2 reservation is Stage B).
    std::size_t edgesToCreate = 0, edgesToRetain = 0, edgesToRetire = 0;
    std::size_t cellsRemoved = 0, cellsAdded = 0;
    bool semanticNoop = false;            // §9 of the design
};

// M4-B.8 STAGE B (under review) — the Phase-2 PLAN: exact per-pool
// accounting with in-transaction slot reuse and the predicted
// allocation order, per M4B8_DESIGN.md §6. PLAIN NUMBERS, computed by
// CONST reads only. Structural note (CORRECTED per Stage-B review):
// for a Tri3 disk cavity, |M| == |R| (Euler) forces cell demand ==
// cell retirements, and the invariant interior-edge count forces edge
// creates == edge retires — so per pool, demand and retirements are
// ALWAYS balanced, and phase2Claims == exhaustedRetirements. Those
// claims are satisfied FIRST from the pre-existing free list, so the
// correct invariant is
//     preExistingConsumed + freshSlots == exhaustedRetirements
// and freshSlots == exhaustedRetirements ONLY when no pre-existing
// free slot satisfies the claim.
struct CavityPoolPlan {
    std::size_t demand = 0;              // slots Phase 3 will allocate
    std::size_t plannedRetirements = 0;  // slots Phase 3 will free
    std::size_t exhaustedRetirements = 0;  // gen == max: retire in place
    std::size_t reusableThisCommit = 0;  // retirements minus exhausted
    std::size_t inCommitReuse = 0;       // min(demand, reusableThisCommit)
    std::size_t phase2Claims = 0;        // demand - inCommitReuse
    std::size_t freeListAvailable = 0;   // pre-existing free slots
    std::size_t preExistingConsumed = 0; // min(phase2Claims, available)
    std::size_t freshSlots = 0;          // phase2Claims - preExisting
    std::size_t finalFreeListSize = 0;
    // The committed allocation order predicted from the validated plan
    // AND the captured strict allocator state (design §4 wording):
    // phase2Claims slots first (free-list LIFO then fresh), then
    // inCommitReuse slots as LIFO pops of the just-freed sequence.
    std::vector<std::int32_t> predictedAllocationOrder;
    std::uint64_t sequenceDemand = 0;    // == demand (one seq per alloc)
    std::uint64_t remainingSequences = 0;  // captured at plan time
    std::vector<std::uint32_t> predictedGenerations;  // parallel to order
    bool balanced() const {              // CORRECTED (Stage-B review):
        return demand == plannedRetirements &&
               preExistingConsumed + freshSlots == exhaustedRetirements;
    }
};
struct CavityPlan {
    CavityPoolPlan cells, edges;
    // Per-pool sequence demand and exact remaining space (via the
    // disclosed remainingSequences static); prepare() validates
    // remaining >= demand per pool BEFORE any preparation.
};

// (Stage C will populate this; declared with the accepted field names.)
struct CavityCommit {
    std::vector<CellHandle> resultCells;      // EXISTING handles on a
    std::vector<EdgeHandle> requiredEdges;    // semantic no-op
    std::uint64_t revision = 0;
};

// ---------------------------------------------------------------------------
// BND-2A / ASCR ADDITIVE combinatorial helpers (rev-3 §0.5.4). PLAIN VALUES
// ONLY: they operate on triangles expressed as RANK triples (each rank a
// unique node id - a live sequence for an existing node, or a predicted
// sequence for a staged node), so the SAME disk / node-link / boundary
// invariants the frozen buildPreview enforces can be reused by the
// staged-aware SteinerCavity orchestration WITHOUT touching the frozen
// buildPreview body or the legacy route. No pool access here.
namespace cavaudit {

using Rank = std::uint64_t;
using RankTri = std::array<Rank, 3>;
using RankEdge = std::array<Rank, 2>;          // {lo, hi}, lo < hi

inline RankEdge canonEdge(Rank a, Rank b) {
    return a < b ? RankEdge{a, b} : RankEdge{b, a};
}

// A rank complex is WELL-FORMED iff it is non-empty, every triangle has three
// DISTINCT ranks (no self-edge), and no two triangles are the same up to
// orientation (no identical or reversed duplicate). Every audit below rejects
// a malformed complex BEFORE any edge counting, so degenerate input can never
// be read as a valid disk / link.
inline bool wellFormedComplex(const std::vector<RankTri>& tris) {
    if (tris.empty()) return false;
    std::set<std::array<Rank, 3>> seen;
    for (const auto& t : tris) {
        if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2]) return false;
        std::array<Rank, 3> key{t[0], t[1], t[2]};
        std::sort(key.begin(), key.end());
        if (!seen.insert(key).second) return false;   // orientation-insensitive
    }
    return true;
}

// Per-edge incidence count across the triangle set (canonical rank pairs).
inline std::map<RankEdge, int> edgeCounts(const std::vector<RankTri>& tris) {
    std::map<RankEdge, int> m;
    for (const auto& t : tris)
        for (int i = 0; i < 3; ++i)
            ++m[canonEdge(t[std::size_t(i)], t[std::size_t((i + 1) % 3)])];
    return m;
}

// Boundary edges = edges incident to exactly ONE triangle. (Interior edges
// are incident to two; an edge incident to >2 is non-manifold and makes the
// set NOT a disk - callers detect that via auditCavityDiskCombinatorial.)
inline std::set<RankEdge> boundaryEdgeSet(const std::vector<RankTri>& tris) {
    std::set<RankEdge> b;
    for (const auto& [e, c] : edgeCounts(tris))
        if (c == 1) b.insert(e);
    return b;
}

// The boundary edges must form EXACTLY ONE simple cycle: every boundary rank
// has degree 2 and the walk from the lowest rank visits them all. On success
// fills 'cycle' with the ordered ring (canonical start = min boundary rank).
inline bool boundaryCycle(const std::vector<RankTri>& tris,
                          std::vector<Rank>& cycle) {
    cycle.clear();
    if (!wellFormedComplex(tris)) return false;
    const auto bedges = boundaryEdgeSet(tris);
    if (bedges.empty()) return false;
    std::map<Rank, std::vector<Rank>> adj;
    for (const auto& e : bedges) {
        adj[e[0]].push_back(e[1]);
        adj[e[1]].push_back(e[0]);
    }
    for (const auto& [r, nb] : adj) {
        (void)r;
        if (nb.size() != 2) return false;               // not a simple cycle
    }
    const Rank start = adj.begin()->first;
    Rank cur = start, prev = 0;
    bool first = true;
    std::set<Rank> visited;
    do {
        if (!visited.insert(cur).second) return false;   // premature revisit
        cycle.push_back(cur);
        const auto& nb = adj[cur];
        const Rank nxt = first ? nb[0] : (nb[0] != prev ? nb[0] : nb[1]);
        first = false;
        prev = cur;
        cur = nxt;
    } while (cur != start);
    return cycle.size() == adj.size();                    // one component
}

// Full disk audit over ranks: non-empty; edge-connected through shared
// edges; no edge used more than twice; Euler V - E + F == 1; exactly one
// boundary cycle.
inline bool auditCavityDiskCombinatorial(const std::vector<RankTri>& tris) {
    if (!wellFormedComplex(tris)) return false;
    const auto ec = edgeCounts(tris);
    for (const auto& [e, c] : ec) { (void)e; if (c > 2) return false; }
    // edge-connectivity via shared (interior) edges
    std::map<RankEdge, std::vector<std::size_t>> byEdge;
    for (std::size_t i = 0; i < tris.size(); ++i)
        for (int k = 0; k < 3; ++k)
            byEdge[canonEdge(tris[i][std::size_t(k)],
                             tris[i][std::size_t((k + 1) % 3)])].push_back(i);
    std::set<std::size_t> seen;
    std::vector<std::size_t> stk{0};
    while (!stk.empty()) {
        const std::size_t c = stk.back();
        stk.pop_back();
        if (!seen.insert(c).second) continue;
        for (int k = 0; k < 3; ++k)
            for (std::size_t nb : byEdge[canonEdge(
                     tris[c][std::size_t(k)], tris[c][std::size_t((k + 1) % 3)])])
                if (!seen.count(nb)) stk.push_back(nb);
    }
    if (seen.size() != tris.size()) return false;
    std::set<Rank> V;
    for (const auto& t : tris)
        for (Rank r : t) V.insert(r);
    const long long v = static_cast<long long>(V.size());
    const long long e = static_cast<long long>(ec.size());
    const long long f = static_cast<long long>(tris.size());
    if (v - e + f != 1) return false;
    std::vector<Rank> cyc;
    return boundaryCycle(tris, cyc);
}

// Node-link audit over ranks: for every node the opposite edges of its
// incident triangles form ONE path or cycle (degree <= 2, single component).
inline bool auditCavityNodeLink(const std::vector<RankTri>& tris) {
    if (!wellFormedComplex(tris)) return false;   // empty -> false (not vacuous)
    std::map<Rank, std::vector<RankEdge>> link;
    for (const auto& t : tris)
        for (int i = 0; i < 3; ++i)
            link[t[std::size_t(i)]].push_back(
                RankEdge{t[std::size_t((i + 1) % 3)],
                         t[std::size_t((i + 2) % 3)]});
    for (const auto& [n, segs] : link) {
        (void)n;
        std::map<Rank, std::vector<Rank>> adj;
        std::map<Rank, int> deg;
        for (const auto& s : segs) {
            if (++deg[s[0]] > 2 || ++deg[s[1]] > 2) return false;
            adj[s[0]].push_back(s[1]);
            adj[s[1]].push_back(s[0]);
        }
        if (adj.empty()) return false;
        std::set<Rank> ls;
        std::vector<Rank> st{segs[0][0]};
        while (!st.empty()) {
            const Rank v = st.back();
            st.pop_back();
            if (!ls.insert(v).second) continue;
            for (Rank w : adj[v]) st.push_back(w);
        }
        if (ls.size() != adj.size()) return false;       // pinched
    }
    return true;
}

}  // namespace cavaudit

// ------------------------------------------------------------ transaction ---
class CavityTransaction {
public:
    CavityTransaction(const CavityTransaction&) = delete;
    CavityTransaction& operator=(const CavityTransaction&) = delete;
    CavityTransaction(CavityTransaction&& o) noexcept
        : ws_(o.ws_), resume_(o.resume_), staged_(std::move(o.staged_)),
          plan_(std::move(o.plan_)), prep_(std::move(o.prep_)),
          resultCellsBuf_(std::move(o.resultCellsBuf_)),
          requiredEdgesBuf_(std::move(o.requiredEdgesBuf_)),
          rejected_(o.rejected_) {
        o.ws_ = nullptr;               // source finished; the captured
        o.plan_.reset();               // allocator-state plan and the
        o.prep_.prepared = false;      // prepared placeholders transfer
    }
    CavityTransaction& operator=(CavityTransaction&&) = delete;
    ~CavityTransaction() {
        if (ws_) abandonImpl();
    }

    // Phase 1: full combinatorial validation (design §2.0/§2.2/§3/§9)
    // and canonicalization; const reads only; pools untouched.
    const CavityPreview& stage(const TriCavitySpec& spec) {
        requireActive("stage");
        if (staged_)
            throw std::logic_error(
                "CavityTransaction: exactly one cavity may be staged per "
                "transaction");
        staged_.emplace(buildPreview(spec));
        return *staged_;
    }

    // Frozen validator pattern: bool-returning callable over plain
    // values; returning false marks the transaction REJECTED (ordinary
    // control flow, not an error). Runs in Phase 1 only.
    bool validate(const std::function<bool(const CavityPreview&)>& fn) {
        requireActive("validate");
        if (ws_->state_ != WorkspaceState::TransactionStaging)
            throw std::logic_error(
                "CavityTransaction::validate: validators run in Phase 1 "
                "only (workspace is TransactionPrepared)");
        if (!staged_)
            throw std::logic_error(
                "CavityTransaction::validate: nothing staged");
        const bool ok = fn(*staged_);        // frozen semantics: return
        if (!ok) rejected_ = true;           // THIS call's result;
        return ok;                           // rejection is sticky in
    }                                        // rejected_, not the return

    bool rejected() const noexcept { return rejected_; }
    const CavityPreview& preview() const {
        requireActive("preview");            // Phase-1-only contract
        if (!staged_)
            throw std::logic_error(
                "CavityTransaction::preview: nothing staged");
        return *staged_;
    }

    void abandon() {
        requireActive("abandon");
        abandonImpl();
    }

    // ---- STAGE B (under review): the Phase-2 plan --------------------------
    // CONST reads only; no reservation and no mutation. Reservation of
    // the pool RECORD arrays / edge table / incidence containers is NOT
    // expressible through the granted friendships (see the Stage-B
    // deviation report in docs/M4B8_DESIGN.md §10b) and is deliberately
    // absent pending review.
    const CavityPlan& computePlan() {
        requireActive("computePlan");
        if (!staged_)
            throw std::logic_error(
                "CavityTransaction::computePlan: nothing staged");
        if (rejected_)
            throw std::logic_error(
                "CavityTransaction::computePlan: transaction was "
                "rejected by a validator");
        if (plan_) return *plan_;
        using TA = detail::TransactionAccess;
        CavityPlan p;
        const CavityPreview& v = *staged_;
        // ---- cells ----
        {
            CavityPoolPlan& q = p.cells;
            const auto& a = TA::cellAllocC(ws_->cells_);
            q.demand = v.semanticNoop ? 0 : v.cellsAdded;
            q.plannedRetirements = v.semanticNoop ? 0 : v.cellsRemoved;
            std::vector<std::int32_t> freedOrder;   // canonical removal
            for (const CellHandle& h : v.removedCells) {
                if (v.semanticNoop) break;
                if (TA::slotGen(a, h.slot) ==
                    std::numeric_limits<std::uint32_t>::max())
                    ++q.exhaustedRetirements;
                else
                    freedOrder.push_back(h.slot);
            }
            fillPoolPlan(q, a, freedOrder);
        }
        // ---- edges ----
        {
            CavityPoolPlan& q = p.edges;
            const auto& a = TA::edgeAllocC(ws_->edges_);
            std::vector<std::int32_t> freedOrder;
            if (!v.semanticNoop) {
                q.demand = v.edgesToCreate;
                q.plannedRetirements = v.edgesToRetire;
                for (const auto& row : v.edges) {
                    if (row.initialUse == 0 || row.finalUse > 0) continue;
                    const auto eh =
                        ws_->edges_.find(ws_->nodes_, row.a, row.b);
                    if (!eh)
                        throw std::logic_error(
                            "CavityTransaction::computePlan: a retiring "
                            "edge vanished between stage and plan");
                    if (TA::slotGen(a, eh->slot) ==
                        std::numeric_limits<std::uint32_t>::max())
                        ++q.exhaustedRetirements;
                    else
                        freedOrder.push_back(eh->slot);
                }
            }
            fillPoolPlan(q, a, freedOrder);
        }
        plan_.emplace(std::move(p));
        return *plan_;
    }
    // ---- STAGE B: Phase-2 preparation (amended option 1) -------------------
    // Batched per-pool sequence validation runs BEFORE any preparation;
    // then the disclosed op-shaped prepareCavity places exact
    // placeholders and capacity. No handle is minted, no sequence
    // consumed, no binding/revision/topology touched. Idempotent.
    void prepare() {
        requireActive("prepare");
        if (!staged_)
            throw std::logic_error(
                "CavityTransaction::prepare: nothing staged");
        if (rejected_)
            throw std::logic_error(
                "CavityTransaction::prepare: transaction was rejected by "
                "a validator");
        if (prep_.prepared) return;
        using TA = detail::TransactionAccess;
        const CavityPlan& p = computePlan();
        // Exact remaining-space validation, per pool, BEFORE anything.
        if (p.cells.remainingSequences < p.cells.sequenceDemand)
            throw std::logic_error(
                "CavityTransaction::prepare: insufficient creation-"
                "sequence space (cells) for the batched demand");
        if (p.edges.remainingSequences < p.edges.sequenceDemand)
            throw std::logic_error(
                "CavityTransaction::prepare: insufficient creation-"
                "sequence space (edges) for the batched demand");
        // Phase-3 result buffers are reserved BEFORE any pool
        // placeholder exists (review blocker 2): if these throw, nothing
        // has been placed and the transaction is NOT prepared.
        resultCellsBuf_.reserve(staged_->cellsAdded);
        std::size_t reqUpper = 0;
        for (const auto& row : staged_->edges)
            if (row.willBeRequired) ++reqUpper;
        requiredEdgesBuf_.reserve(reqUpper);
        TA::CavityPrepareInput in;
        in.cellAllocOrder = p.cells.predictedAllocationOrder;
        in.edgeAllocOrder = p.edges.predictedAllocationOrder;
        in.cellFreePushes = p.cells.reusableThisCommit;
        in.edgeFreePushes = p.edges.reusableThisCommit;
        in.cellFresh = p.cells.freshSlots;
        in.edgeFresh = p.edges.freshSlots;
        for (const auto& row : staged_->edges)
            if (row.initialUse == 0 && row.finalUse > 0)
                in.newEdgeEndpoints.emplace_back(row.a, row.b);
        try {
            TA::prepareCavity(ws_->nodes_, ws_->cells_, ws_->edges_, in,
                              prep_);
        } catch (...) {
            TA::abandonCavityPrepared(ws_->cells_, ws_->edges_, prep_);
            throw;                     // capacity residue only (legal)
        }
        // prep_.prepared became true only after EVERY Phase-3 dependency
        // was ready; the shared frozen state machine records it (review
        // blocker 3) as the FINAL, non-throwing step.
        ws_->state_ = WorkspaceState::TransactionPrepared;
    }
    bool prepared() const noexcept { return prep_.prepared; }

    // ---- STAGE C (under review): commit -----------------------------------
    // The SEMANTIC NO-OP path is complete: it returns EXISTING handles,
    // changes no revision, no handle, no metadata, and restores the
    // exact prior workspace state (Ready or CommittedAdjacencyStale).
    // The SUBSTANTIVE path requires the commit-side counterpart of the
    // disclosed prepared-cavity surface (commitCavityPrepared) and is
    // returned for review as section 10c of the design: the frozen
    // Stage-B placeholders are consumable ONLY by a bespoke noexcept
    // static (public insertCellEdges rejects the pre-placed incidence
    // and table entries; public removeCellEdges retires retained edges
    // through zero). No unreviewed frozen edit is made.
    CavityCommit commit() {
        requireActive("commit");
        if (!staged_)
            throw std::logic_error(
                "CavityTransaction::commit: nothing staged");
        if (rejected_)
            throw std::logic_error(
                "CavityTransaction::commit: transaction was rejected by "
                "a validator");
        if (staged_->semanticNoop) {
            CavityCommit out;
            out.resultCells = staged_->removedCells;   // EXISTING handles
            for (const auto& row : staged_->edges)
                if (row.requested) {
                    const auto eh =
                        ws_->edges_.find(ws_->nodes_, row.a, row.b);
                    out.requiredEdges.push_back(*eh);  // already satisfied
                }
            out.revision = ws_->revision_;             // UNCHANGED
            detail::TransactionAccess::abandonCavityPrepared(
                ws_->cells_, ws_->edges_, prep_);      // no-op plan: empty
            ws_->state_ = resume_;                     // EXACT prior state
            ws_ = nullptr;                             // finished
            return out;
        }
        // ---- SUBSTANTIVE PATH (10c approved, amended) ----
        using TA = detail::TransactionAccess;
        if (!prep_.prepared) prepare();          // may throw (auto-prep)
        const CavityPlan& p = *plan_;
        // EVERYTHING potentially throwing happens HERE, before Phase 3:
        // handle resolution, reference building, scratch sizing.
        TA::CavityCommitInput in;
        in.removedCells = staged_->removedCells;
        in.removedExhausted.resize(in.removedCells.size(), 0);
        {
            const auto& a = TA::cellAllocC(ws_->cells_);
            for (std::size_t k = 0; k < in.removedCells.size(); ++k)
                if (TA::slotGen(a, in.removedCells[k].slot) ==
                    std::numeric_limits<std::uint32_t>::max())
                    in.removedExhausted[k] = 1;
        }
        for (std::size_t j = 0; j < staged_->addedConnectivity.size();
             ++j) {
            TA::CavityCommitInput::AddedCell ac;
            for (int i = 0; i < 3; ++i) {
                const NodeHandle& n =
                    staged_->addedConnectivity[j][std::size_t(i)];
                ac.nodeSlot[std::size_t(i)] = n.slot;
                ac.nodeGen[std::size_t(i)] = n.gen;
            }
            ac.cls = staged_->addedClass[j];
            in.added.push_back(ac);
        }
        in.cellSlots = p.cells.predictedAllocationOrder;
        in.cellGens = p.cells.predictedGenerations;
        in.cellPhase2Claims = p.cells.phase2Claims;
        in.edgeSlots = p.edges.predictedAllocationOrder;
        in.edgeGens = p.edges.predictedGenerations;
        in.edgePhase2Claims = p.edges.phase2Claims;
        // Edge rows -> retains / retires / creates, plus a rank-pair map
        // for incidence and requirement resolution.
        std::map<std::pair<std::uint64_t, std::uint64_t>,
                 TA::CavityCommitInput::ERef>
            byKey;
        auto rankOf = [&](const NodeHandle& n) {
            return ws_->nodes_.seqOf(n);
        };
        for (const auto& row : staged_->edges) {
            const auto key = std::make_pair(rankOf(row.a), rankOf(row.b));
            if (row.initialUse > 0 && row.finalUse > 0) {
                TA::CavityCommitInput::Retain rt;
                rt.h = *ws_->edges_.find(ws_->nodes_, row.a, row.b);
                rt.finalUse = row.finalUse;
                rt.makeRequired = row.willBeRequired ? 1 : 0;
                rt.setCls =
                    (row.finalCls.dim != row.initialCls.dim ||
                     row.finalCls.tag != row.initialCls.tag)
                        ? 1
                        : 0;
                rt.cls = row.finalCls;
                byKey[key] = TA::CavityCommitInput::ERef{0, rt.h, 0};
                in.retains.push_back(rt);
            } else if (row.initialUse > 0 && row.finalUse == 0) {
                in.retires.push_back(
                    *ws_->edges_.find(ws_->nodes_, row.a, row.b));
            } else if (row.initialUse == 0 && row.finalUse > 0) {
                TA::CavityCommitInput::Create cs;
                // EdgeData endpoints are ordered lo->hi by seq; rows are
                // already canonical (rank(a) < rank(b)).
                cs.a = row.a;
                cs.b = row.b;
                cs.finalUse = row.finalUse;
                cs.required = row.willBeRequired ? 1 : 0;
                cs.cls = row.finalCls;
                byKey[key] = TA::CavityCommitInput::ERef{
                    1, EdgeHandle{},
                    static_cast<std::uint32_t>(in.creates.size())};
                in.creates.push_back(cs);
            }
        }
        for (const auto& tri : staged_->addedConnectivity) {
            std::array<TA::CavityCommitInput::ERef, 3> refs;
            for (int i = 0; i < 3; ++i) {
                const NodeHandle& x = tri[std::size_t(i)];
                const NodeHandle& y = tri[std::size_t((i + 1) % 3)];
                const auto rx = rankOf(x), ry = rankOf(y);
                refs[std::size_t(i)] = byKey.at(
                    rx < ry ? std::make_pair(rx, ry)
                            : std::make_pair(ry, rx));
            }
            in.addedEdges.push_back(refs);
        }
        for (const auto& row : staged_->edges)
            if (row.requested)
                in.requested.push_back(
                    byKey.at(std::make_pair(rankOf(row.a),
                                            rankOf(row.b))));
        std::vector<EdgeHandle> createdScratch(in.creates.size());
        resultCellsBuf_.clear();                 // capacity retained
        requiredEdgesBuf_.clear();
        TA::CavityCommitBuffers bufs;
        bufs.resultCells = &resultCellsBuf_;
        bufs.requiredEdges = &requiredEdgesBuf_;
        bufs.createdScratch = &createdScratch;
        // ---- Phase 3: bounded, noexcept, allocation-free ----
#if defined(PAACS_TOPO_TESTING)
        detail::phase3Active() = true;   // the frozen probe instrument
#endif
        TA::commitCavityPrepared(ws_->cells_, ws_->edges_, in, prep_,
                                 bufs);
#if defined(PAACS_TOPO_TESTING)
        detail::phase3Active() = false;  // (testing-only, like Workspace)
#endif
        // Only no-throw operations remain (10c contract): revision
        // increment, stale-state transition, buffer moves, finalization.
        ws_->revision_ += 1;                     // exactly ONCE
        ws_->state_ = WorkspaceState::CommittedAdjacencyStale;
        CavityCommit result;
        result.resultCells = std::move(resultCellsBuf_);
        result.requiredEdges = std::move(requiredEdgesBuf_);
        result.revision = ws_->revision_;
        ws_ = nullptr;                           // finished
        return result;
    }

    // prepare()/commit() are DELIBERATELY ABSENT in Stage A: no
    // reservation and no mutation path exists until Stage B/C review.

private:
    friend class TopologyWorkspace;

    explicit CavityTransaction(TopologyWorkspace& ws)
        : ws_(&ws), resume_(ws.state_) {
        ws_->state_ = WorkspaceState::TransactionStaging;
    }
    void requireActive(const char* what) const {
        if (!ws_)
            throw std::logic_error(std::string("CavityTransaction: '") +
                what + "' on a finished/abandoned transaction");
    }
    void abandonImpl() noexcept {
        // Tracking-based and unconditional: cleans exactly what was
        // placed, whether or not preparation completed.
        detail::TransactionAccess::abandonCavityPrepared(
            ws_->cells_, ws_->edges_, prep_);
        ws_->state_ = resume_;
        ws_ = nullptr;
    }

    // ---- Phase-1 combinatorial engine (const reads only) ----
    struct Key {                          // canonical edge key (§2.0.1)
        std::uint64_t ra, rb;             // ranks, ra < rb
        bool operator<(const Key& o) const {
            return std::tie(ra, rb) < std::tie(o.ra, o.rb);
        }
        bool operator==(const Key& o) const {
            return ra == o.ra && rb == o.rb;
        }
    };
    struct EdgeInfo {
        NodeHandle a, b;                  // handles, rank order
        std::uint32_t initialUse = 0;
        std::uint32_t localRemovedUse = 0, localAddedUse = 0;
        int removedSenseSum = 0, addedSenseSum = 0;   // for opposite check
        std::int8_t removedSense = 0, addedSense = 0; // last seen (boundary)
        bool live = false, wasRequired = false;
        mesh::Classification initialCls;
        std::optional<EdgeHandle> liveHandle;
    };

    static void rej(const std::string& why) {
        throw std::invalid_argument("CavityTransaction::stage: " + why);
    }

    CavityPreview buildPreview(const TriCavitySpec& spec) const {
        const NodePool& np = ws_->nodes_;
        const CellPool& cp = ws_->cells_;
        const EdgePool& ep = ws_->edges_;
        // Composite readiness gate FIRST (Stage-A review blocker 3):
        // staging derives from the pools, so it must refuse to read an
        // inconsistent workspace rather than derive garbage from it.
        if (!workspaceTopologyReady(np, cp, ep, ws_->faces_))
            throw std::logic_error(
                "CavityTransaction::stage: workspace is not "
                "topology-ready (composite gate failed)");

        // --- removed: non-empty, current, distinct, Tri3; canonical ---
        if (spec.removed.empty()) rej("removed set is empty");
        std::vector<CellHandle> removed = spec.removed;
        for (const CellHandle& h : removed) {
            if (!cp.isCurrent(h))
                rej("removed cell is not current (stale or foreign)");
            if (cp.get(h).type != ElementType::Tri3)
                rej("M4-B.8 qualifies Tri3 cavities only");
        }
        std::sort(removed.begin(), removed.end(),
                  [&](const CellHandle& x, const CellHandle& y) {
                      return cp.seqOf(x) < cp.seqOf(y);
                  });
        for (std::size_t i = 1; i < removed.size(); ++i)
            if (removed[i - 1] == removed[i])
                rej("removed set contains duplicates");

        auto rankOf = [&](const NodeHandle& n) { return np.seqOf(n); };

        // --- gather removed connectivity + node sets ---
        std::vector<std::array<NodeHandle, 3>> rConn;
        std::set<std::uint64_t> rNodeRanks, mNodeRanks;
        std::map<std::uint64_t, NodeHandle> rankToNode;
        for (const CellHandle& h : removed) {
            const auto sn = cp.get(h);
            std::array<NodeHandle, 3> tri{sn.nodes[0], sn.nodes[1],
                                          sn.nodes[2]};
            for (const NodeHandle& n : tri) {
                rNodeRanks.insert(rankOf(n));
                rankToNode.emplace(rankOf(n), n);
            }
            rConn.push_back(tri);
        }

        // --- added: current, distinct nodes; canonical rotation ---
        std::vector<std::array<NodeHandle, 3>> aConn;
        std::vector<mesh::Classification> aCls;
        for (const TriCavityCell& c : spec.added) {
            for (const NodeHandle& n : c.nodes)
                if (!np.isCurrent(n))
                    rej("added triangle references a non-current node");
            if (c.nodes[0] == c.nodes[1] || c.nodes[1] == c.nodes[2] ||
                c.nodes[0] == c.nodes[2])
                rej("added triangle repeats a node");
            for (const NodeHandle& n : c.nodes) {
                if (!rNodeRanks.count(rankOf(n)))
                    rej("added triangle pulls in a node from outside the "
                        "cavity (nodes(M) must be a subset of nodes(R))");
                mNodeRanks.insert(rankOf(n));
            }
            // canonical rotation: lowest rank first, cyclic order kept
            std::array<NodeHandle, 3> t = c.nodes;
            int lo = 0;
            for (int i = 1; i < 3; ++i)
                if (rankOf(t[i]) < rankOf(t[lo])) lo = i;
            std::array<NodeHandle, 3> rot{t[lo], t[(lo + 1) % 3],
                                          t[(lo + 2) % 3]};
            aConn.push_back(rot);
            aCls.push_back(c.cls);
        }
        // canonical triangle order + duplicate detection
        {
            std::vector<std::size_t> idx(aConn.size());
            for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = i;
            auto keyOf = [&](const std::array<NodeHandle, 3>& t) {
                return std::array<std::uint64_t, 3>{rankOf(t[0]),
                                                    rankOf(t[1]),
                                                    rankOf(t[2])};
            };
            std::sort(idx.begin(), idx.end(),
                      [&](std::size_t x, std::size_t y) {
                          return keyOf(aConn[x]) < keyOf(aConn[y]);
                      });
            std::vector<std::array<NodeHandle, 3>> a2;
            std::vector<mesh::Classification> c2;
            for (std::size_t i : idx) {
                a2.push_back(aConn[i]);
                c2.push_back(aCls[i]);
            }
            for (std::size_t i = 1; i < a2.size(); ++i)
                if (keyOf(a2[i - 1]) == keyOf(a2[i]))
                    rej("added set contains duplicate triangles");
            aConn.swap(a2);
            aCls.swap(c2);
        }

        // --- edge closure: keys, uses, senses (§2.0) ---
        std::map<Key, EdgeInfo> edges;
        auto touch = [&](const NodeHandle& x, const NodeHandle& y)
            -> EdgeInfo& {
            const std::uint64_t rx = rankOf(x), ry = rankOf(y);
            const Key k = rx < ry ? Key{rx, ry} : Key{ry, rx};
            auto it = edges.find(k);
            if (it == edges.end()) {
                EdgeInfo e;
                if (rx < ry) { e.a = x; e.b = y; }
                else         { e.a = y; e.b = x; }
                const auto eh = ep.find(np, x, y);
                if (eh) {
                    e.live = true;
                    e.liveHandle = *eh;
                    const auto d = ep.get(*eh);
                    e.initialUse = d.useCount;
                    e.wasRequired = d.requiredConstraint;
                    e.initialCls = d.cls;
                }
                it = edges.emplace(k, e).first;
            }
            return it->second;
        };
        auto senseOf = [&](const NodeHandle& from, const NodeHandle& to) {
            return rankOf(from) < rankOf(to) ? std::int8_t(+1)
                                             : std::int8_t(-1);
        };
        for (const auto& t : rConn)
            for (int i = 0; i < 3; ++i) {
                EdgeInfo& e = touch(t[i], t[(i + 1) % 3]);
                e.localRemovedUse += 1;
                e.removedSense = senseOf(t[i], t[(i + 1) % 3]);
                e.removedSenseSum += e.removedSense;
            }
        for (const auto& t : aConn)
            for (int i = 0; i < 3; ++i) {
                EdgeInfo& e = touch(t[i], t[(i + 1) % 3]);
                e.localAddedUse += 1;
                e.addedSense = senseOf(t[i], t[(i + 1) % 3]);
                e.addedSenseSum += e.addedSense;
            }

        // --- global-manifold precondition + per-edge rules ---
        std::set<Key> rBoundary, rInterior, mBoundary, mInterior;
        for (auto& [k, e] : edges) {
            if (e.live && e.initialUse > 2)
                rej("source edge in the cavity closure is already "
                    "nonmanifold (initialUse > 2)");
            if (e.localRemovedUse > 2 || e.localAddedUse > 2)
                rej("a cavity edge is used more than twice within one "
                    "complex (nonmanifold candidate)");
            // UNCONDITIONAL (review blocker 3): also fires when the
            // edge lookup FAILED (live=false, initialUse=0), so the
            // unsigned outsideUse subtraction below can never hide
            // missing live incidence.
            if (e.localRemovedUse > e.initialUse)
                rej("removed cells use an edge more often than its live "
                    "use count (missing or corrupted incidence)");
            const std::uint32_t outsideUse =
                e.initialUse - e.localRemovedUse;
            if (e.localRemovedUse == 1) rBoundary.insert(k);
            if (e.localRemovedUse == 2) {
                if (outsideUse != 0)
                    rej("interior(R) edge has outside use");   // unreachable
                                                               // under <=2,
                                                               // kept exact
                rInterior.insert(k);
                if (e.removedSenseSum != 0)
                    rej("interior edge of R has non-opposite senses "
                        "(inconsistent orientation in the removed set)");
            }
            if (e.localAddedUse == 1) mBoundary.insert(k);
            if (e.localAddedUse == 2) {
                if (outsideUse != 0)
                    rej("candidate interior edge is already used by a "
                        "surviving outside triangle (outsideUse != 0)");
                mInterior.insert(k);
                if (e.addedSenseSum != 0)
                    rej("interior edge of M has non-opposite senses "
                        "(inconsistent candidate orientation)");
            }
        }

        // --- boundary sets and inside senses match (§2.0.4/5) ---
        if (rBoundary != mBoundary)
            rej("boundary(M) != boundary(R) (canonical edge sets differ)");
        for (const Key& k : rBoundary) {
            const EdgeInfo& e = edges.at(k);
            if (e.removedSense != e.addedSense)
                rej("a boundary edge's inside sense differs between R "
                    "and M (outside orientation relation would break)");
        }

        // --- disk invariant for R and for M (§2.2) ---
        checkDisk(rConn, rInterior, rBoundary, rankOf, "removed set R");
        checkDisk(aConn, mInterior, mBoundary, rankOf, "candidate set M");

        // --- boundary cycle (canonical §2.1) + orphan rule ---
        std::vector<NodeHandle> cycle =
            canonicalBoundaryCycle(edges, rBoundary, rankToNode);
        std::set<std::uint64_t> boundaryRanks;
        for (const NodeHandle& n : cycle) boundaryRanks.insert(rankOf(n));
        std::vector<NodeHandle> interior;
        for (const std::uint64_t r : rNodeRanks)
            if (!boundaryRanks.count(r)) {
                if (!mNodeRanks.count(r))
                    rej("cavity interior node would be orphaned (every "
                        "interior node of R must appear in M)");
                interior.push_back(rankToNode.at(r));
            }

        // --- |M| == |R| (Euler consequence, asserted §2.2) ---
        if (aConn.size() != rConn.size())
            rej("|M| != |R| - impossible for a valid disk replacement "
                "(Euler); the candidate cannot satisfy the invariant");

        // --- no candidate duplicates a surviving cell ---
        {
            // Review blocker 2: duplicate detection must be ORIENTATION-
            // INSENSITIVE (sorted rank triple), or a surviving (a,c,b)
            // silently escapes a candidate (a,b,c). Oriented canonical
            // keys remain in use for senses and the semantic no-op.
            std::set<std::array<std::uint64_t, 3>> mset;
            for (const auto& t : aConn) {
                std::array<std::uint64_t, 3> k{rankOf(t[0]), rankOf(t[1]),
                                               rankOf(t[2])};
                std::sort(k.begin(), k.end());
                mset.insert(k);
            }
            std::set<std::int32_t> removedSlots;
            for (const CellHandle& h : removed) removedSlots.insert(h.slot);
            for (const CellHandle& h : cp.activeBySequence()) {
                if (removedSlots.count(h.slot)) continue;
                const auto sn = cp.get(h);
                if (sn.type != ElementType::Tri3) continue;
                std::array<std::uint64_t, 3> k{rankOf(sn.nodes[0]),
                                               rankOf(sn.nodes[1]),
                                               rankOf(sn.nodes[2])};
                std::sort(k.begin(), k.end());
                if (mset.count(k))
                    rej("an added triangle duplicates a surviving cell "
                        "(orientation-insensitive corner set)");
            }
        }

        // --- requirements: canonical, ensure-only (§3) ---
        std::map<Key, mesh::Classification> reqs;
        for (const CavityEdgeRequirement& r : spec.requirements) {
            if (!np.isCurrent(r.a) || !np.isCurrent(r.b) || r.a == r.b)
                rej("requirement endpoints must be distinct current nodes");
            const std::uint64_t ra = rankOf(r.a), rb = rankOf(r.b);
            const Key k = ra < rb ? Key{ra, rb} : Key{rb, ra};
            const auto it = reqs.find(k);
            if (it != reqs.end()) {
                if (it->second.dim != r.cls.dim ||
                    it->second.tag != r.cls.tag)
                    rej("duplicate contradictory requirements for one edge");
                continue;                       // idempotent duplicate
            }
            reqs.emplace(k, r.cls);
        }
        const mesh::Classification defCls{};    // {-1,-1}
        auto isDefault = [&](const mesh::Classification& c) {
            return c.dim == defCls.dim && c.tag == defCls.tag;
        };
        for (const auto& [k, cls] : reqs) {
            const auto it = edges.find(k);
            // CLOSURE-ONLY SCOPE (Stage-A review): requirements must
            // name edges of the CAVITY CLOSURE. An untouched outside
            // edge is out of scope for this operation by contract, not
            // merely "absent" - callers mark outside edges through the
            // ordinary frozen metadata transaction instead.
            if (it == edges.end())
                rej("requirement names an edge outside the cavity "
                    "closure (out of scope for a cavity operation)");
            const bool inFinal =
                (it->second.initialUse - it->second.localRemovedUse +
                 it->second.localAddedUse) > 0;
            if (!inFinal)
                rej("requirement names an edge absent from the final "
                    "topology (it retires in this cavity)");
            const EdgeInfo& e = it->second;
            if (e.live && e.wasRequired) {
                if (!(e.initialCls.dim == cls.dim &&
                      e.initialCls.tag == cls.tag))
                    rej("conflicting classification on an already "
                        "constrained edge");
            } else if (e.live && !isDefault(e.initialCls) &&
                       !(e.initialCls.dim == cls.dim &&
                         e.initialCls.tag == cls.tag)) {
                rej("conflicting non-default classification on an "
                    "existing edge (a geometric association may not be "
                    "silently overwritten)");
            }
        }

        // --- required-edge preservation (§2.0.8) ---
        for (const auto& [k, e] : edges) {
            const std::uint32_t finalUse =
                e.initialUse - e.localRemovedUse + e.localAddedUse;
            if (e.live && e.wasRequired && finalUse == 0)
                rej("an initially required edge would be retired by this "
                    "cavity (finalUse == 0)");
        }

        // --- assemble the plain-value preview ---
        CavityPreview p;
        p.boundaryCycle = std::move(cycle);
        std::sort(interior.begin(), interior.end(),
                  [&](const NodeHandle& x, const NodeHandle& y) {
                      return rankOf(x) < rankOf(y);
                  });
        p.interiorNodes = std::move(interior);
        p.removedCells = std::move(removed);
        p.removedConnectivity = std::move(rConn);
        p.addedConnectivity = aConn;
        p.addedClass = aCls;
        p.cellsRemoved = p.removedConnectivity.size();
        p.cellsAdded = aConn.size();
        for (const auto& [k, e] : edges) {
            CavityPreview::EdgeRow row;
            row.a = e.a; row.b = e.b;
            row.initialUse = e.initialUse;
            row.finalUse =
                e.initialUse - e.localRemovedUse + e.localAddedUse;
            if (e.localRemovedUse == 1) row.insideSenseRemoved = e.removedSense;
            if (e.localAddedUse == 1) row.insideSenseAdded = e.addedSense;
            row.wasRequired = e.wasRequired;
            row.initialCls = e.initialCls;
            const auto rit = reqs.find(k);
            row.willBeRequired =
                e.wasRequired || rit != reqs.end();
            row.requested = rit != reqs.end();   // Stage-C additive
            row.finalCls = rit != reqs.end() ? rit->second : e.initialCls;
            if (!e.live && row.finalUse > 0) ++p.edgesToCreate;
            if (e.live && row.finalUse > 0) ++p.edgesToRetain;
            if (e.live && row.finalUse == 0) ++p.edgesToRetire;
            p.edges.push_back(row);
        }
        // --- semantic no-op (§9): oriented connectivity + cls + reqs ---
        {
            auto canonR = rConnCanonical(p.removedConnectivity, cp,
                                         p.removedCells, rankOf);
            std::vector<std::array<std::uint64_t, 3>> canonM;
            std::vector<mesh::Classification> clsM = p.addedClass;
            for (const auto& t : p.addedConnectivity)
                canonM.push_back(
                    {rankOf(t[0]), rankOf(t[1]), rankOf(t[2])});
            bool same = canonR.first == canonM;
            if (same)
                for (std::size_t i = 0; same && i < clsM.size(); ++i)
                    same = canonR.second[i].dim == clsM[i].dim &&
                           canonR.second[i].tag == clsM[i].tag;
            if (same)
                for (const auto& [k, cls] : reqs) {
                    const EdgeInfo& e = edges.at(k);
                    if (!e.wasRequired ||
                        !(e.initialCls.dim == cls.dim &&
                          e.initialCls.tag == cls.tag)) {
                        same = false;
                        break;
                    }
                }
            p.semanticNoop = same;
        }
        return p;
    }

    template <class RankOf>
    static std::pair<std::vector<std::array<std::uint64_t, 3>>,
                     std::vector<mesh::Classification>>
    rConnCanonical(const std::vector<std::array<NodeHandle, 3>>& conn,
                   const CellPool& cp,
                   const std::vector<CellHandle>& cells, RankOf rankOf) {
        std::vector<std::array<std::uint64_t, 3>> keys;
        std::vector<mesh::Classification> cls;
        std::vector<std::size_t> idx(conn.size());
        for (std::size_t i = 0; i < conn.size(); ++i) idx[i] = i;
        std::vector<std::array<std::uint64_t, 3>> rot(conn.size());
        for (std::size_t i = 0; i < conn.size(); ++i) {
            const auto& t = conn[i];
            int lo = 0;
            for (int j = 1; j < 3; ++j)
                if (rankOf(t[j]) < rankOf(t[lo])) lo = j;
            rot[i] = {rankOf(t[lo]), rankOf(t[(lo + 1) % 3]),
                      rankOf(t[(lo + 2) % 3])};
        }
        std::sort(idx.begin(), idx.end(), [&](std::size_t x, std::size_t y) {
            return rot[x] < rot[y];
        });
        for (std::size_t i : idx) {
            keys.push_back(rot[i]);
            cls.push_back(cp.get(cells[i]).cls);
        }
        return {keys, cls};
    }

    template <class RankOf>
    static void checkDisk(
        const std::vector<std::array<NodeHandle, 3>>& tris,
        const std::set<Key>& interiorKeys, const std::set<Key>& boundaryKeys,
        RankOf rankOf, const char* which) {
        // edge-connectivity through interior edges
        if (!tris.empty()) {
            std::vector<int> comp(tris.size(), -1);
            std::map<Key, std::vector<std::size_t>> byEdge;
            for (std::size_t i = 0; i < tris.size(); ++i)
                for (int j = 0; j < 3; ++j) {
                    const std::uint64_t ra = rankOf(tris[i][j]);
                    const std::uint64_t rb = rankOf(tris[i][(j + 1) % 3]);
                    byEdge[ra < rb ? Key{ra, rb} : Key{rb, ra}].push_back(i);
                }
            std::vector<std::size_t> stack{0};
            comp[0] = 0;
            while (!stack.empty()) {
                const std::size_t i = stack.back();
                stack.pop_back();
                for (int j = 0; j < 3; ++j) {
                    const std::uint64_t ra = rankOf(tris[i][j]);
                    const std::uint64_t rb = rankOf(tris[i][(j + 1) % 3]);
                    const Key k = ra < rb ? Key{ra, rb} : Key{rb, ra};
                    if (!interiorKeys.count(k)) continue;
                    for (const std::size_t o : byEdge[k])
                        if (comp[o] == -1) {
                            comp[o] = 0;
                            stack.push_back(o);
                        }
                }
            }
            for (const int c : comp)
                if (c == -1)
                    rej(std::string(which) +
                        " is not edge-connected (not a single complex)");
        }
        // NODE-LINK MANIFOLD AUDIT (review blocker 1): chi == 1 does NOT
        // exclude a pinched vertex - identifying two nonadjacent
        // vertices of a sphere yields V-1, so chi drops from 2 to 1 with
        // every edge still used twice. Each node's incident triangles,
        // connected through the edges AT THAT NODE, must form ONE
        // component; interior links are single cycles
        // (#edges@v == #tris@v), boundary links single paths
        // (#edges@v == #tris@v + 1, exactly two boundary edges).
        {
            std::map<std::uint64_t, std::vector<std::size_t>> trisAt;
            std::map<std::uint64_t, std::set<Key>> edgesAt;
            for (std::size_t i = 0; i < tris.size(); ++i)
                for (int j = 0; j < 3; ++j) {
                    const std::uint64_t rv = rankOf(tris[i][j]);
                    trisAt[rv].push_back(i);
                    const std::uint64_t ra = rv;
                    const std::uint64_t rb = rankOf(tris[i][(j + 1) % 3]);
                    edgesAt[ra].insert(ra < rb ? Key{ra, rb} : Key{rb, ra});
                    edgesAt[rb].insert(ra < rb ? Key{ra, rb} : Key{rb, ra});
                }
            for (const auto& [rv, tv] : trisAt) {
                const auto& ev = edgesAt.at(rv);
                std::size_t bAtV = 0;
                for (const Key& k : ev)
                    if (boundaryKeys.count(k)) ++bAtV;
                if (bAtV == 0) {
                    if (ev.size() != tv.size())
                        rej(std::string(which) +
                            ": an interior node's link is not a single "
                            "cycle (node link not manifold)");
                } else {
                    if (bAtV != 2 || ev.size() != tv.size() + 1)
                        rej(std::string(which) +
                            ": a boundary node's link is not a single "
                            "path (node link not manifold)");
                }
                // single component through the edges at this node
                std::vector<int> comp(tv.size(), -1);
                std::vector<std::size_t> st{0};
                comp[0] = 0;
                while (!st.empty()) {
                    const std::size_t i = st.back();
                    st.pop_back();
                    for (std::size_t o = 0; o < tv.size(); ++o) {
                        if (comp[o] != -1) continue;
                        // share an edge at v?
                        bool share = false;
                        for (int j = 0; j < 3 && !share; ++j) {
                            const std::uint64_t ra =
                                rankOf(tris[tv[i]][j]);
                            const std::uint64_t rb =
                                rankOf(tris[tv[i]][(j + 1) % 3]);
                            if (ra != rv && rb != rv) continue;
                            const Key k =
                                ra < rb ? Key{ra, rb} : Key{rb, ra};
                            for (int m = 0; m < 3 && !share; ++m) {
                                const std::uint64_t rc =
                                    rankOf(tris[tv[o]][m]);
                                const std::uint64_t rd =
                                    rankOf(tris[tv[o]][(m + 1) % 3]);
                                const Key k2 = rc < rd ? Key{rc, rd}
                                                       : Key{rd, rc};
                                if (k == k2) share = true;
                            }
                        }
                        if (share) {
                            comp[o] = 0;
                            st.push_back(o);
                        }
                    }
                }
                for (const int cc : comp)
                    if (cc == -1)
                        rej(std::string(which) +
                            ": a node's link has multiple components "
                            "(pinched vertex; node link not manifold)");
            }
        }
        // NONEMPTY BOUNDARY (review blocker 1): a closed complex is not
        // a disk, whatever its Euler characteristic says.
        if (boundaryKeys.empty())
            rej(std::string(which) +
                ": boundary is empty (closed complex is not a disk)");
        // one simple boundary cycle: every boundary node has exactly 2
        // boundary edges, and traversal closes a single cycle
        std::map<std::uint64_t, std::vector<Key>> atNode;
        for (const Key& k : boundaryKeys) {
            atNode[k.ra].push_back(k);
            atNode[k.rb].push_back(k);
        }
        for (const auto& [r, ks] : atNode)
            if (ks.size() != 2)
                rej(std::string(which) +
                    ": a boundary node does not have exactly two boundary "
                    "edges (pinched or branching boundary)");
        if (!boundaryKeys.empty()) {
            std::set<Key> seen;
            std::uint64_t start = atNode.begin()->first, cur = start;
            std::uint64_t prev = ~0ULL;
            do {
                const auto& ks = atNode.at(cur);
                const Key* next = nullptr;
                for (const Key& k : ks) {
                    const std::uint64_t other = k.ra == cur ? k.rb : k.ra;
                    if (other != prev || seen.count(k)) {
                        if (!seen.count(k)) { next = &k; break; }
                    }
                }
                if (!next) {
                    for (const Key& k : ks)
                        if (!seen.count(k)) { next = &k; break; }
                }
                if (!next) break;
                seen.insert(*next);
                prev = cur;
                cur = next->ra == cur ? next->rb : next->ra;
            } while (cur != start);
            if (seen.size() != boundaryKeys.size())
                rej(std::string(which) +
                    ": boundary is not one simple cycle (annulus or "
                    "multiple loops)");
        }
        // Euler characteristic V - E + F == 1
        std::set<std::uint64_t> verts;
        std::set<Key> allKeys;
        for (const auto& t : tris)
            for (int j = 0; j < 3; ++j) {
                verts.insert(rankOf(t[j]));
                const std::uint64_t ra = rankOf(t[j]);
                const std::uint64_t rb = rankOf(t[(j + 1) % 3]);
                allKeys.insert(ra < rb ? Key{ra, rb} : Key{rb, ra});
            }
        const long chi = long(verts.size()) - long(allKeys.size()) +
                         long(tris.size());
        if (chi != 1)
            rej(std::string(which) +
                ": Euler characteristic V - E + F != 1 (not a "
                "topological disk)");
    }

    static std::vector<NodeHandle> canonicalBoundaryCycle(
        const std::map<Key, EdgeInfo>& edges, const std::set<Key>& bKeys,
        const std::map<std::uint64_t, NodeHandle>& rankToNode) {
        // Start at the lowest-rank boundary node; traverse in the
        // lexicographically smaller of the two directions (§2.1).
        if (bKeys.empty()) return {};
        std::map<std::uint64_t, std::vector<std::uint64_t>> nbr;
        for (const Key& k : bKeys) {
            nbr[k.ra].push_back(k.rb);
            nbr[k.rb].push_back(k.ra);
        }
        const std::uint64_t start = nbr.begin()->first;
        auto walk = [&](std::uint64_t first) {
            std::vector<std::uint64_t> seq{start, first};
            std::uint64_t prev = start, cur = first;
            while (cur != start) {
                const auto& ns = nbr.at(cur);
                const std::uint64_t nxt = ns[0] == prev ? ns[1] : ns[0];
                seq.push_back(nxt);
                prev = cur;
                cur = nxt;
            }
            seq.pop_back();                      // drop repeated start
            return seq;
        };
        const auto& s = nbr.at(start);
        const auto w1 = walk(std::min(s[0], s[1]));
        const auto w2 = walk(std::max(s[0], s[1]));
        const auto& best = w1 <= w2 ? w1 : w2;
        std::vector<NodeHandle> out;
        for (const std::uint64_t r : best) out.push_back(rankToNode.at(r));
        (void)edges;
        return out;
    }

    template <class Alloc>
    void fillPoolPlan(CavityPoolPlan& q, const Alloc& a,
                      const std::vector<std::int32_t>& freedOrder) const {
        using TA = detail::TransactionAccess;
        q.reusableThisCommit = freedOrder.size();
        q.inCommitReuse = std::min(q.demand, q.reusableThisCommit);
        q.phase2Claims = q.demand - q.inCommitReuse;
        // Blocker-3 fix: the free-list size is PUBLIC - no whole-pool
        // prediction probe, no allocation proportional to pool size.
        q.freeListAvailable = a.reusableFreeSlots();
        const std::size_t materialized =
            static_cast<std::size_t>(a.materializedSlots());
        q.sequenceDemand = q.demand;
        q.remainingSequences = TA::remainingSequences(a);
        const auto claims = TA::predictSlots(a, q.phase2Claims);
        q.preExistingConsumed = std::min(q.phase2Claims,
                                         q.freeListAvailable);
        q.freshSlots = q.phase2Claims - q.preExistingConsumed;
        // Predicted committed order: Phase-2 claims first (captured
        // allocator state - LIFO free list then fresh, exactly what
        // predictSlots models), then the in-commit reuses as LIFO pops
        // of the just-freed sequence (canonical removal order pushes;
        // pops reverse it).
        q.predictedAllocationOrder = claims;
        for (std::size_t i = 0; i < q.inCommitReuse; ++i)
            q.predictedAllocationOrder.push_back(
                freedOrder[freedOrder.size() - 1 - i]);
        // Predicted generations, parallel to the order: an in-commit
        // reuse allocates AFTER its free bumps the generation (+1); a
        // pre-existing free slot already carries its bumped generation;
        // a fresh slot starts at 0.
        for (std::size_t i = 0; i < q.predictedAllocationOrder.size(); ++i) {
            const auto s = q.predictedAllocationOrder[i];
            if (i >= claims.size())
                q.predictedGenerations.push_back(TA::slotGen(a, s) + 1);
            else if (static_cast<std::size_t>(s) < materialized)
                q.predictedGenerations.push_back(TA::slotGen(a, s));
            else
                q.predictedGenerations.push_back(0);
        }
        q.finalFreeListSize = (q.freeListAvailable - q.preExistingConsumed) +
                              (q.reusableThisCommit - q.inCommitReuse);
    }

    TopologyWorkspace* ws_ = nullptr;
    WorkspaceState resume_ = WorkspaceState::Ready;
    std::optional<CavityPreview> staged_;
    std::optional<CavityPlan> plan_;
    detail::TransactionAccess::CavityPrepared prep_;
    std::vector<CellHandle> resultCellsBuf_;
    std::vector<EdgeHandle> requiredEdgesBuf_;
    bool rejected_ = false;
};

inline CavityTransaction TopologyWorkspace::beginCavityTransaction() {
    if (!readable())
        throw std::logic_error(
            "TopologyWorkspace::beginCavityTransaction: a transaction is "
            "already active (exactly one transaction may be active)");
    return CavityTransaction(*this);
}

}  // namespace paacs::topo
