// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Shared testing instruments: the TransactionAccess probe and the FULL
// observable-state digest defined by docs/M4B7_DESIGN.md.
//
// C2 review: tests/probe_compaction.cpp previously computed its OWN
// ad-hoc digest that omitted the five pool bindings and the V1 hash,
// which made its "full source digest preserved" claim an overclaim.
// Both the test binary and the standalone probes now share ONE digest
// implementation, so the claim means the same thing everywhere.
#include <bit>
#include <cstdint>
#include <vector>

#include "paacs/topo/Compaction.hpp"
#include "paacs/topo/Workspace.hpp"

// Test-only probe: the ONLY sanctioned caller of the private prepared
// primitives besides TopologyTransaction (friend gated by
// PAACS_TOPO_TESTING; absent from production builds).
namespace paacs::topo::detail {
struct TransactionAccessTestProbe {
    template <class Tag>
    static Handle<Tag> alloc(HandleAllocator<Tag>& a) noexcept {
        static_assert(noexcept(TransactionAccess::allocPrepared(a)),
                      "prepared alloc must be noexcept");
        return TransactionAccess::allocPrepared(a);
    }
    // M4-B.8 Stage-B fixture passthroughs (testing-only seams).
    static void setRemainingEdgeSequences(EdgePool& p, std::uint64_t r) {
        TransactionAccess::testSetRemainingSequences(
            TransactionAccess::edgeAlloc(p), r);
    }
    // BND-2A: node-sequence boundary fixtures for the Steiner batch.
    static void setRemainingNodeSequences(NodePool& p, std::uint64_t r) {
        TransactionAccess::testSetRemainingSequences(
            TransactionAccess::nodeAlloc(p), r);
    }
    static CellHandle exhaustLiveCell(CellPool& p, const CellHandle& h) {
        return TransactionAccess::testExhaustLiveCell(p, h);
    }
    template <class P>
    static std::size_t slotCap(const P& pool) {
        if constexpr (std::is_same_v<P, NodePool>)
            return TransactionAccess::slotCapacityOf(
                TransactionAccess::nodeAllocC(pool));
        else if constexpr (std::is_same_v<P, CellPool>)
            return TransactionAccess::slotCapacityOf(
                TransactionAccess::cellAllocC(pool));
        else if constexpr (std::is_same_v<P, EdgePool>)
            return TransactionAccess::slotCapacityOf(
                TransactionAccess::edgeAllocC(pool));
        else
            return TransactionAccess::slotCapacityOf(
                TransactionAccess::faceAllocC(pool));
    }
    template <class P>
    static std::size_t materialized(const P& pool) {
        if constexpr (std::is_same_v<P, NodePool>)
            return TransactionAccess::nodeAllocC(pool).materializedSlots();
        else if constexpr (std::is_same_v<P, CellPool>)
            return TransactionAccess::cellAllocC(pool).materializedSlots();
        else if constexpr (std::is_same_v<P, EdgePool>)
            return TransactionAccess::edgeAllocC(pool).materializedSlots();
        else
            return TransactionAccess::faceAllocC(pool).materializedSlots();
    }
    static auto planAdd(const NodePool& np, const CellPool& cp,
                        const EdgePool& ep, const FacePool& fp,
                        ElementType t, std::span<const NodeHandle> nn,
                        const mesh::Classification& cls) {
        return TransactionAccess::planAddCell(np, cp, ep, fp, t, nn, cls,
                                              TopologyPolicy{});
    }
    static void prepareAdd(NodePool& np, CellPool& cp, EdgePool& ep,
                           FacePool& fp, TransactionAccess::AddPlan& p) {
        TransactionAccess::prepareAddCell(np, cp, ep, fp, p);
    }
    static CellHandle commitAdd(NodePool& np, CellPool& cp, EdgePool& ep,
                                FacePool& fp,
                                TransactionAccess::AddPlan& p) noexcept {
        return TransactionAccess::commitAddCell(np, cp, ep, fp, p);
    }
    static void abandonAdd(EdgePool& ep, FacePool& fp,
                           TransactionAccess::AddPlan& p) noexcept {
        TransactionAccess::abandonAddPrepared(ep, fp, p);
    }
    static auto planRemoveNode(const NodePool& np, const CellPool& cp,
                               const EdgePool& ep, const FacePool& fp,
                               const NodeHandle& h) {
        return TransactionAccess::planRemoveNode(np, cp, ep, fp, h);
    }
    static void forceEdgeSeqMax(EdgePool& e) {
        TransactionAccess::edgeAlloc(e).testForceSequenceMax();
    }
    static void forceFaceSeqMax(FacePool& f) {
        TransactionAccess::faceAlloc(f).testForceSequenceMax();
    }
    static TransactionAccess::PoolBindings bindings(
        const CellPool& c, const EdgePool& e, const FacePool& f) {
        return TransactionAccess::debugBindings(c, e, f);
    }
    template <class Tag>
    static void free(HandleAllocator<Tag>& a, const Handle<Tag>& h) noexcept {
        static_assert(noexcept(TransactionAccess::freePrepared(a, h)),
                      "prepared free must be noexcept");
        TransactionAccess::freePrepared(a, h);
    }
};
} // namespace paacs::topo::detail

namespace paacs_tests {
namespace topo = paacs::topo;
namespace mesh = paacs::mesh;
// withPoolIdentity=false drops runtime OWNERSHIP identity (pool ids and
// bindings) while keeping semantic content and strict physical layout.
// Two compactions of one source are deliberately DIFFERENT ownership
// domains (fresh pool ids), so cross-run determinism is asserted modulo
// pool identity - the design's three-concept separation in action.
// M4-B.8 Stage B: the named LOGICAL-state digest - hashes, revision,
// bindings, handles, payloads, metadata, incidence, live/free/retired
// counts - EXCLUDING capacity-derived byte fields, which Phase-2
// abandonment may legally have grown (monotonic residue).
inline std::vector<std::uint64_t> observableLogicalDigest(
    const paacs::topo::TopologyWorkspace& ws);

inline std::vector<std::uint64_t> observableDigest(
    const paacs::topo::TopologyWorkspace& ws,
    bool withPoolIdentity = true, bool includeCapacityClasses = true) {
    // The V2 token STREAM leads the digest, per M4B7_DESIGN.md ("V2
    // token stream + both hashes ..."). The stream IS the semantic
    // content; the folds appended below cannot substitute for it. (C2
    // closure review: extracting this function into the shared header
    // accidentally dropped the tokens the standalone probe used to
    // carry - a regression the shared digest is meant to prevent.)
    std::vector<std::uint64_t> d = ws.testV2Tokens();
    auto p = [&](std::uint64_t v) { d.push_back(v); };
    auto pp = [&](std::uint64_t v) { if (withPoolIdentity) d.push_back(v); };
    auto pd = [&](double v) { p(std::bit_cast<std::uint64_t>(v)); };
    auto pcls = [&](const paacs::mesh::Classification& c) {
        p(static_cast<std::uint64_t>(c.dim));
        p(static_cast<std::uint64_t>(c.tag));
    };
    // hashes / revision / adjacency state / workspace state
    p(ws.semanticTopologyHash());
    p(ws.semanticTopologyHashV1());
    p(ws.strictWorkspaceHash());
    p(ws.topologyRevision());
    p(ws.adjacencyBuiltFromRevision());
    p(ws.adjacencyReady() ? 1 : 0);
    p(static_cast<std::uint64_t>(ws.state()));
    // pool bindings
    {
        const auto b = paacs::topo::detail::TransactionAccessTestProbe::bindings(ws.cells(), ws.edges(), ws.faces());
        pp(b.cellNode); pp(b.edgeNode); pp(b.edgeCell);
        pp(b.faceNode); pp(b.faceCell);
    }
    // per-pool counts; capacity/byte classes only when requested (the
    // LOGICAL digest excludes them: Phase-2 residue is capacity-class)
    for (const auto& st : {ws.nodes().stats(), ws.cells().stats(),
                           ws.edges().stats(), ws.faces().stats()}) {
        p(st.live); p(st.reusableFree); p(st.retired); p(st.totalCreated);
        if (includeCapacityClasses) {
            p(st.materializedSlots); p(st.peakLive); p(st.peakMaterialized);
            p(st.payloadBytes); p(st.recordBytes);
            p(st.estimatedOtherBytes); p(st.allocatorBytes);
            p(st.totalBytes);
        }
    }
    if (includeCapacityClasses) {
        const auto m = paacs::topo::measureWorkspaceMemory(ws);
        p(m.nodeBytes); p(m.cellBytes); p(m.edgeBytes); p(m.faceBytes);
        p(m.recordBytes); p(m.allocatorBytes); p(m.estimatedOtherBytes);
        p(m.totalBytes);
    }
    // nodes: iteration order, identity, currentness, coords, class
    for (const paacs::topo::NodeHandle& h : ws.nodes().activeBySequence()) {
        pp(h.pool); p(static_cast<std::uint64_t>(h.slot)); p(h.gen);
        p(ws.nodes().isCurrent(h) ? 1 : 0);
        p(ws.nodes().seqOf(h));
        const auto v = ws.nodes().get(h);
        pd(v.coord.x); pd(v.coord.y); pd(v.coord.z);
        pcls(v.cls);
    }
    // cells: identity, connectivity, class, and INCIDENCE signatures
    for (const paacs::topo::CellHandle& h : ws.cells().activeBySequence()) {
        pp(h.pool); p(static_cast<std::uint64_t>(h.slot)); p(h.gen);
        p(ws.cells().isCurrent(h) ? 1 : 0);
        p(ws.cells().seqOf(h));
        const auto sn = ws.cells().get(h);
        p(static_cast<std::uint64_t>(sn.type));
        p(sn.count);
        for (std::uint8_t i = 0; i < sn.count; ++i) {
            p(static_cast<std::uint64_t>(sn.nodes[i].slot));
            p(sn.nodes[i].gen);
        }
        pcls(sn.cls);
        const auto ce = ws.edges().cellEdges(ws.cells(), h);
        p(ce.size());
        for (const topo::EdgeHandle& eh : ce) {
            p(static_cast<std::uint64_t>(eh.slot)); p(eh.gen);
        }
        if (ws.faces().hasCellFaces(ws.cells(), h)) {
            const auto cf = ws.faces().cellFaces(ws.cells(), h);
            p(cf.size());
            for (const topo::FaceHandle& fh : cf) {
                p(static_cast<std::uint64_t>(fh.slot)); p(fh.gen);
            }
        } else {
            p(0xFFFFFFFFULL);
        }
    }
    // edges: identity, endpoints, class, constraint, use count
    for (const topo::EdgeHandle& h : ws.edges().activeBySequence()) {
        pp(h.pool); p(static_cast<std::uint64_t>(h.slot)); p(h.gen);
        p(ws.edges().isCurrent(h) ? 1 : 0);
        const auto v = ws.edges().get(h);
        p(static_cast<std::uint64_t>(v.a.slot)); p(v.a.gen);
        p(static_cast<std::uint64_t>(v.b.slot)); p(v.b.gen);
        pcls(v.cls);
        p(v.requiredConstraint ? 1 : 0);
        p(v.useCount);
    }
    // faces: identity, corners, class, full incidence (slot,gen,sense)
    for (const topo::FaceHandle& h : ws.faces().activeBySequence()) {
        pp(h.pool); p(static_cast<std::uint64_t>(h.slot)); p(h.gen);
        p(ws.faces().isCurrent(h) ? 1 : 0);
        const auto v = ws.faces().get(h);
        p(static_cast<std::uint64_t>(v.shape));
        p(v.cornerCount);
        for (std::uint8_t i = 0; i < v.cornerCount; ++i) {
            p(static_cast<std::uint64_t>(v.corners[i].slot));
            p(v.corners[i].gen);
        }
        pcls(v.cls);
        p(v.incident.size());
        for (const auto& in : v.incident) {
            p(static_cast<std::uint64_t>(in.cellSlot));
            p(in.cellGen);
            p(static_cast<std::uint64_t>(
                static_cast<std::int64_t>(in.sense)));
        }
    }
    return d;
}
// M4-B.8 Stage B: the named LOGICAL-state digest - hashes, revision,
// bindings, handles, payloads, metadata, incidence, live/free/retired
// counts - EXCLUDING capacity-derived byte fields, which Phase-2
// abandonment may legally have grown (monotonic residue).
inline std::vector<std::uint64_t> observableLogicalDigest(
    const paacs::topo::TopologyWorkspace& ws) {
    return observableDigest(ws, true, /*includeCapacityClasses=*/false);
}
}  // namespace paacs_tests
