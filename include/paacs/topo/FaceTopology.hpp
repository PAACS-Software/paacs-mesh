// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/core/SmallSort.hpp"
//
// paacs::topo - M4-B.4: explicit reusable face topology.
//
// Faces are first-class entities derived from 3D cell connectivity, for
// boundary extraction, interface tracking, and hex-layer logic. Composed
// ON the frozen M4-B.1 HandleAllocator<FaceTag> and the frozen M4-B.2/.3
// NodePool / CellPool / EdgePool. Consumes ONLY storage+edges that pass
// edgeTopologyReady(nodes, cells, edges).
//
// IDENTITY (directive 2.1: topological identity is node identity, never
// geometry): a face's canonical key is the SET of its corner node
// sequences, sorted ascending - independent of which cell contributed it
// or that cell's local winding. Runtime pool ids / slots / generations
// never enter the key, serialization, or hashes.
//
// ORIENTATION: each incident cell contributes the face in its own OUTWARD
// winding (from ElementTraits face tables). The face record stores one
// canonical winding (the first contributor's outward order); each
// incidence records whether that cell's outward winding matches (+1) or is
// reversed (-1) relative to canonical. For a well-formed 2-cell interior
// face the two senses are opposite.
//
// CLASSIFICATION:
//   * boundary : exactly 1 incident cell.
//   * interior : exactly 2 incident cells (opposite senses expected).
//   * nonmanifold : >= 3 incident cells (detected + flagged, not rejected;
//                   the caller decides - some workflows build nonmanifold
//                   intermediate states).
//
// SCOPE (reviewer): canonical face identity; orientation per incident cell;
// FaceHandle ownership; unique reusable face records; cell<->face
// incidence; boundary/interior classification; traits-driven mixed-element
// face extraction; nonmanifold detection; structured rejection of stale
// cells/nodes, missing edge topology, duplicate insertion; three-state
// accounting + bounded churn. DEFERRED: ownership adjacency beyond
// incidence, transaction manager, Delaunay migration, public Mesh
// replacement, parallel face construction.
//
#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "paacs/core/Types.hpp"
#include "paacs/diag/Diagnostics.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/topo/EdgeTopology.hpp"
#include "paacs/topo/ElementTraits.hpp"
#include "paacs/topo/Handles.hpp"
#include "paacs/topo/Storage.hpp"

namespace paacs::topo {

inline constexpr int kMaxFaceNodes = 4;   // Tri or Quad in this rollout

// Canonical face key: up to 4 corner node sequences, sorted ascending,
// with the count. Excludes pool ids/slots/generations.
struct FaceKey {
    std::array<std::uint64_t, kMaxFaceNodes> s{};
    std::uint8_t n = 0;
    bool operator==(const FaceKey& o) const {
        if (n != o.n) return false;
        for (int i = 0; i < n; ++i) if (s[i] != o.s[i]) return false;
        return true;
    }
};

struct FaceKeyHash {
    std::size_t operator()(const FaceKey& k) const noexcept {
        std::uint64_t h = 1469598103934665603ULL;
        auto mix = [&](std::uint64_t v) {
            for (int b = 0; b < 8; ++b) { h ^= (v >> (8*b)) & 0xFF; h *= 1099511628211ULL; }
        };
        mix(k.n);
        for (int i = 0; i < k.n; ++i) mix(k.s[i]);
        return static_cast<std::size_t>(h);
    }
};

enum class FaceClass : std::uint8_t { Boundary, Interior, NonManifold };

struct FacePoolCompleteBuildMetrics {
    std::uint64_t validationCellVisits = 0;
    std::uint64_t buildCellVisits = 0;
    std::uint64_t faceUses = 0;
    std::uint64_t uniqueFaces = 0;
    std::uint64_t suppliedEdgeAudits = 0;
    std::uint64_t finalFaceAudits = 0;
    std::uint64_t perCellGlobalAudits = 0;
    bool runComplete = false;
};

static_assert(std::is_standard_layout_v<FacePoolCompleteBuildMetrics>);
static_assert(std::is_trivially_copyable_v<FacePoolCompleteBuildMetrics>);

class FacePool;
inline diag::Report auditFaceTopology(
    const NodePool&,const CellPool&,const EdgePool&,const FacePool&);

// ---------------------------------------------------------------------------
// FacePool
// ---------------------------------------------------------------------------
class FacePool {
    friend class detail::TransactionAccess;   // M4-B.6 prepared commit
public:
    struct Incidence {
        std::int32_t cellSlot = -1;
        std::uint32_t cellGen = 0;
        std::int8_t sense = 0;            // +1 canonical, -1 reversed
    };
    struct FaceData {
        FaceShape shape{FaceShape::None};
        std::uint8_t cornerCount = 0;
        std::array<NodeHandle, kMaxFaceNodes> corners{};   // canonical winding
        mesh::Classification cls{};
        std::vector<Incidence> incident;  // 1=boundary, 2=interior, >=3 nonman
    };

    FacePool() = default;
    FacePool(const FacePool&) = delete;
    FacePool& operator=(const FacePool&) = delete;
    FacePool(FacePool&& o)
        : FacePool(TokenlessAdoptionTag{}) {
#if defined(PAACS_TOPO_TESTING)
        if (adoptionBeginHook_ != nullptr) adoptionBeginHook_();
#endif
        adoptFrom(o);
    }
    FacePool& operator=(FacePool&&) = delete;

    static FacePool buildComplete(
        const NodePool& nodes,const CellPool& cells,const EdgePool& edges,
        FacePoolCompleteBuildMetrics* metrics = nullptr) {
        return FacePool(CompleteBuildTag{},nodes,cells,edges,metrics);
    }

    std::uint64_t poolId() const { return alloc_.poolId(); }
    std::uint64_t boundNodePool() const { return boundNodePool_; }
    std::uint64_t boundCellPool() const { return boundCellPool_; }
    bool isCurrent(const FaceHandle& h) const { return alloc_.isCurrent(h); }
    std::size_t liveCount() const { return alloc_.liveCount(); }

    // Insert every face of a 3D `cell` (traits-driven). Creates face records
    // on first use, records this cell's incidence + orientation sense on
    // each. STRONGLY EXCEPTION-SAFE (mirrors EdgePool): all validation and
    // container reservation before mutation; commit rolls back on any throw;
    // free-list pre-reserved so retirement cannot allocate.
    std::vector<FaceHandle> insertCellFaces(const NodePool& nodes,
                                            const CellPool& cells,
                                            const EdgePool& edges,
                                            const CellHandle& cell) {
        // ENFORCE the milestone contract: faces are built on trustworthy
        // edges (review blocker). Reject unless the edge topology is ready
        // AND THIS SPECIFIC cell's edge incidence is recorded in the
        // EdgePool - not merely that the edge RECORDS exist (a geometrically
        // identical cell could have created them; that is not this cell's
        // incidence).
        if (!edgeTopologyReady(nodes, cells, edges))
            throw std::invalid_argument(
                "FacePool::insertCellFaces: edge topology is not ready "
                "(edgeTopologyReady must hold before faces are built)");
        if (!cells.isCurrent(cell))
            throw std::invalid_argument(
                "FacePool::insertCellFaces: stale cell handle");
        if (!edges.hasCellEdges(cells, cell))
            throw std::invalid_argument(
                "FacePool::insertCellFaces: this cell has no edge incidence "
                "in the EdgePool; insert THIS cell's edges before its faces");
        bind(nodes, cells);
        if (!cells.isCurrent(cell))
            throw std::invalid_argument(
                "FacePool::insertCellFaces: stale cell handle");
        {
            const auto it = cellFaces_.find(cell.slot);
            if (it != cellFaces_.end()) {
                if (it->second.gen == cell.gen)
                    throw std::invalid_argument(
                        "FacePool::insertCellFaces: cell faces already "
                        "inserted");
                throw std::invalid_argument(
                    "FacePool::insertCellFaces: stale face incidence for a "
                    "reused cell slot; remove the old cell's faces first");
            }
        }
        const auto snap = cells.get(cell);
        const ElementTraits& tr = checkedTraits(snap.type);
        if (tr.faceCount == 0)
            throw std::invalid_argument(
                "FacePool::insertCellFaces: element type has no faces (not "
                "a 3D cell)");
        // Stage: for each face, gather corner handles (outward order),
        // validate currentness, build the canonical key + sense.
        std::vector<Staged> staged;
        staged.reserve(tr.faceCount);
        for (int f = 0; f < tr.faceCount; ++f) {
            const LocalFace& lf = tr.faces[f];
            Staged st;
            st.shape = lf.shape;
            st.n = lf.n;
            for (int i = 0; i < lf.n; ++i) {
                const NodeHandle nh = snap.nodes[lf.v[static_cast<std::size_t>(i)]];
                if (!nodes.isCurrent(nh))
                    throw std::invalid_argument(
                        "FacePool::insertCellFaces: stale face node");
                st.outward[static_cast<std::size_t>(i)] = nh;
            }
            if (!makeKey(nodes, st.outward, lf.n, st.key))
                throw std::invalid_argument(
                    "FacePool::insertCellFaces: degenerate face (repeated "
                    "node)");
            staged.push_back(st);
        }
        // Reserve everything that may allocate.
        std::vector<FaceHandle> out;
        out.reserve(staged.size());
        std::vector<FaceHandle> cellFaceList;
        cellFaceList.reserve(staged.size());
        table_.reserve(table_.size() + staged.size());
        rec_.reserve(rec_.size() + staged.size());
        alloc_.reserveFreeCapacity(staged.size());
        // Commit with rollback.
        std::vector<FaceHandle> created;
        std::vector<FaceHandle> touched;   // existing faces we pushed to
        created.reserve(staged.size());
        touched.reserve(staged.size());
        auto rollback = [&]() noexcept {
            // Undo incidence pushes on touched (existing) faces.
            for (auto it = touched.rbegin(); it != touched.rend(); ++it) {
                Rec& r = rec_[static_cast<std::size_t>(raw::slotUnchecked(*it))];
                if (!r.data.incident.empty()) r.data.incident.pop_back();
            }
            for (auto it = created.rbegin(); it != created.rend(); ++it)
                retireFace(*it);
        };
        try {
#if defined(PAACS_TOPO_TESTING)
            maybeFault();
#endif
            for (const Staged& st : staged) {
                FaceHandle fh = lookup(st.key);
                Incidence inc;
                inc.cellSlot = cell.slot;
                inc.cellGen = cell.gen;
                if (!alloc_.isCurrent(fh)) {
                    fh = createNew(nodes, st, created);   // sense +1 canonical
                    inc.sense = +1;
                    Rec& r =
                        rec_[static_cast<std::size_t>(raw::slotUnchecked(fh))];
                    r.data.incident.reserve(2);           // interior common
                    r.data.incident.push_back(inc);
                } else {
                    Rec& r =
                        rec_[static_cast<std::size_t>(raw::slotUnchecked(fh))];
                    inc.sense = orientationSense(nodes, r.data, st);
                    r.data.incident.reserve(r.data.incident.size() + 1);
                    r.data.incident.push_back(inc);
                    touched.push_back(fh);
                }
                cellFaceList.push_back(fh);
                out.push_back(fh);
#if defined(PAACS_TOPO_TESTING)
                maybeFault();
#endif
            }
            cellFaces_.emplace(cell.slot,
                               CellFaceEntry{cell.gen, std::move(cellFaceList)});
        } catch (...) {
            rollback();
            throw;
        }
        return out;
    }

    // Remove every face incidence contributed by `cell`; faces whose
    // incidence becomes empty are retired. Validates the cell against its
    // CellPool (current + owned). Pre-reserves free-list capacity so the
    // retirement path cannot allocate.
    void removeCellFaces(const CellPool& cells, const CellHandle& cell) {
        if (boundCellPool_ != 0 && cells.poolId() != boundCellPool_)
            throw std::invalid_argument(
                "FacePool::removeCellFaces: foreign CellPool");
        if (!cells.isCurrent(cell))
            throw std::invalid_argument(
                "FacePool::removeCellFaces: stale cell handle (remove faces "
                "before removing the cell)");
        const auto it = cellFaces_.find(cell.slot);
        if (it == cellFaces_.end() || it->second.gen != cell.gen)
            throw std::invalid_argument(
                "FacePool::removeCellFaces: cell faces were not inserted");
        alloc_.reserveFreeCapacity(it->second.faces.size());
        for (const FaceHandle& fh : it->second.faces) {
            if (!alloc_.isCurrent(fh)) continue;
            Rec& r = rec_[static_cast<std::size_t>(raw::slotUnchecked(fh))];
            // Drop this cell's incidence entry.
            auto& inc = r.data.incident;
            for (auto jt = inc.begin(); jt != inc.end(); ++jt)
                if (jt->cellSlot == cell.slot && jt->cellGen == cell.gen) {
                    inc.erase(jt);
                    break;
                }
            if (inc.empty()) retireFace(fh);
        }
        cellFaces_.erase(it);
    }

    std::optional<FaceHandle> find(const NodePool& nodes,
                                   std::span<const NodeHandle> corners) const {
        if (boundNodePool_ != 0 && nodes.poolId() != boundNodePool_)
            return std::nullopt;
        if (corners.size() < 3 || corners.size() > kMaxFaceNodes)
            return std::nullopt;
        std::array<NodeHandle, kMaxFaceNodes> buf{};
        for (std::size_t i = 0; i < corners.size(); ++i) buf[i] = corners[i];
        FaceKey k;
        if (!makeKey(nodes, buf, static_cast<std::uint8_t>(corners.size()), k))
            return std::nullopt;
        const auto it = table_.find(k);
        if (it == table_.end()) return std::nullopt;
        return it->second;
    }

    FaceData get(const FaceHandle& h) const {
        return rec_[static_cast<std::size_t>(alloc_.deref(h))].data;
    }
    std::uint64_t seqOf(const FaceHandle& h) const { return alloc_.seqOf(h); }

    FaceClass classify(const FaceHandle& h) const {
        const auto& inc =
            rec_[static_cast<std::size_t>(alloc_.deref(h))].data.incident;
        if (inc.size() <= 1) return FaceClass::Boundary;
        if (inc.size() == 2) return FaceClass::Interior;
        return FaceClass::NonManifold;
    }
    bool isBoundary(const FaceHandle& h) const {
        return classify(h) == FaceClass::Boundary;
    }

    void setClassification(const FaceHandle& h,
                           const mesh::Classification& c) {
        rec_[static_cast<std::size_t>(alloc_.deref(h))].data.cls = c;
    }

    std::vector<FaceHandle> activeBySequence() const {
        std::vector<FaceHandle> out;
        out.reserve(liveLinks_);
        for (std::int32_t s = head_; s != detail::kNil; s = rec_[s].next)
            out.push_back(FaceHandle{s, rec_[s].gen, alloc_.poolId()});
        return out;
    }

    // Boundary faces only (exactly one incident cell).
    std::vector<FaceHandle> boundaryFaces() const {
        std::vector<FaceHandle> out;
        for (std::int32_t s = head_; s != detail::kNil; s = rec_[s].next)
            if (rec_[s].data.incident.size() == 1)
                out.push_back(FaceHandle{s, rec_[s].gen, alloc_.poolId()});
        return out;
    }

    // Distinct cells that currently have face incidence recorded here.
    std::vector<CellHandle> incidentCells() const {
        std::vector<CellHandle> out;
        out.reserve(cellFaces_.size());
        for (const auto& [slot, entry] : cellFaces_)
            out.push_back(CellHandle{slot, entry.gen, boundCellPool_});
        return out;
    }

    // Cell-face incidence queries (M4-B.5 adjacency requirement): does THIS
    // specific cell have face incidence recorded (slot AND generation)?
    // Pure const reads of existing incidence state - additive to the frozen
    // M4-B.4 file, no behavior change (same category as EdgePool's
    // hasCellEdges/cellEdges, disclosed in KNOWN_ISSUES).
    bool hasCellFaces(const CellPool& cells, const CellHandle& cell) const {
        if (boundCellPool_ != 0 && cells.poolId() != boundCellPool_)
            return false;
        const auto it = cellFaces_.find(cell.slot);
        return it != cellFaces_.end() && it->second.gen == cell.gen;
    }
    std::span<const FaceHandle> cellFaces(const CellPool& cells,
                                          const CellHandle& cell) const {
        if (boundCellPool_ != 0 && cells.poolId() != boundCellPool_)
            return {};
        const auto it = cellFaces_.find(cell.slot);
        if (it == cellFaces_.end() || it->second.gen != cell.gen) return {};
        return std::span<const FaceHandle>(it->second.faces.data(),
                                           it->second.faces.size());
    }

    PoolStats stats() const {
        PoolStats st;
        st.materializedSlots = alloc_.materializedSlots();
        st.live = alloc_.liveCount();
        st.reusableFree = alloc_.reusableFreeSlots();
        st.retired = alloc_.retiredSlots();
        st.totalCreated = alloc_.totalCreated();
        st.peakLive = peakLive_;
        st.peakMaterialized = peakMaterialized_;
        std::size_t incBytes = 0;
        for (const auto& r : rec_)
            incBytes += r.data.incident.capacity() * sizeof(Incidence);
        std::size_t tableBytes =
            table_.bucket_count() * sizeof(void*) +
            table_.size() * (sizeof(FaceKey) + sizeof(FaceHandle) + sizeof(void*));
        std::size_t cfBytes =
            cellFaces_.bucket_count() * sizeof(void*);
        for (const auto& [slot, e] : cellFaces_)
            cfBytes += e.faces.capacity() * sizeof(FaceHandle) + sizeof(void*);
        st.recordBytes = rec_.capacity() * sizeof(Rec);
        st.estimatedOtherBytes = incBytes + tableBytes + cfBytes;
        st.payloadBytes = st.recordBytes + st.estimatedOtherBytes;
        st.allocatorBytes = alloc_.approxBytes();
        st.totalBytes = st.payloadBytes + st.allocatorBytes;
        return st;
    }

    // Bounds-safe record/list/table audit (nodes form) + incidence audit
    // (nodes, cells form): every incident cell current+bound, every face's
    // incidence non-empty, sense values valid, table<->record round-trip.
    bool checkInvariants(const NodePool& nodes) const {
        return checkListAndTable(nodes);
    }
    bool checkInvariants(const NodePool& nodes, const CellPool& cells) const {
        // Structural (bounds-safe) validation MUST precede the incidence
        // traversal: checkIncidence walks head_/next and cellFaces_ edge
        // handles, which is only safe once checkListAndTable has proven the
        // list and table are well-formed (review blocker: the reversed
        // order segfaulted under ASan on a corrupt head_). Same order as
        // the frozen EdgePool.
        return checkListAndTable(nodes) && checkIncidence(nodes, cells);
    }

    // Reason codes for the cell-aware audit, so auditFaceTopology can emit
    // SPECIFIC diagnostics (review: the enum promised orientation/incidence
    // codes the audit didn't actually use).
    enum class IncidenceFault {
        Ok, BindingMismatch, GhostFace, IncidenceMismatch, OrientationMismatch
    };
    IncidenceFault checkIncidenceReason(const NodePool& nodes,
                                        const CellPool& cells) const {
        // Self-protecting: this traverses the intrusive list and follows
        // recorded handles, which is only safe on a structurally valid
        // pool. Gate on the bounds-safe structural check first so a direct
        // public call on corrupt state reports a fault instead of
        // dereferencing it (review blocker class).
        if (!checkListAndTable(nodes)) return IncidenceFault::IncidenceMismatch;
        if ((boundCellPool_ != 0 && cells.poolId() != boundCellPool_) ||
            (boundNodePool_ != 0 && nodes.poolId() != boundNodePool_))
            return IncidenceFault::BindingMismatch;
        std::unordered_map<std::int32_t, std::size_t> counted;
        for (const auto& [slot, entry] : cellFaces_) {
            const CellHandle ch{slot, entry.gen, cells.poolId()};
            if (!cells.isCurrent(ch)) return IncidenceFault::GhostFace;
            const auto snap = cells.get(ch);
            if (!isValidElementType(snap.type))
                return IncidenceFault::IncidenceMismatch;
            const ElementTraits& tr = checkedTraits(snap.type);
            if (static_cast<int>(entry.faces.size()) != tr.faceCount)
                return IncidenceFault::IncidenceMismatch;
            for (int f = 0; f < tr.faceCount; ++f) {
                const LocalFace& lf = tr.faces[f];
                Staged st; st.shape = lf.shape; st.n = lf.n;
                for (int i = 0; i < lf.n; ++i) {
                    const NodeHandle nh = snap.nodes[lf.v[static_cast<std::size_t>(i)]];
                    if (!nodes.isCurrent(nh)) return IncidenceFault::IncidenceMismatch;
                    st.outward[static_cast<std::size_t>(i)] = nh;
                }
                if (!makeKey(nodes, st.outward, lf.n, st.key))
                    return IncidenceFault::IncidenceMismatch;
                const FaceHandle recorded = entry.faces[static_cast<std::size_t>(f)];
                if (!alloc_.isCurrent(recorded)) return IncidenceFault::IncidenceMismatch;
                const auto tit = table_.find(st.key);
                if (tit == table_.end() ||
                    raw::slotUnchecked(tit->second) != raw::slotUnchecked(recorded))
                    return IncidenceFault::IncidenceMismatch;
                const Rec& fr = rec_[static_cast<std::size_t>(raw::slotUnchecked(recorded))];
                const std::int8_t expectSense = orientationSense(nodes, fr.data, st);
                bool found = false;
                for (const Incidence& in : fr.data.incident)
                    if (in.cellSlot == slot && in.cellGen == entry.gen) {
                        if (in.sense != expectSense)
                            return IncidenceFault::OrientationMismatch;
                        found = true; break;
                    }
                if (!found) return IncidenceFault::IncidenceMismatch;
                ++counted[raw::slotUnchecked(recorded)];
            }
        }
        std::size_t liveSeen = 0;
        for (std::int32_t s = head_; s != detail::kNil; s = rec_[s].next) {
            const Rec& r = rec_[static_cast<std::size_t>(s)];
            const auto it = counted.find(s);
            const std::size_t inc = (it == counted.end()) ? 0 : it->second;
            if (r.data.incident.size() != inc) return IncidenceFault::IncidenceMismatch;
            if (inc == 0) return IncidenceFault::GhostFace;
            for (const Incidence& in : r.data.incident)
                if (in.sense != +1 && in.sense != -1)
                    return IncidenceFault::OrientationMismatch;
            if (r.data.incident.size() == 2 &&
                r.data.incident[0].sense + r.data.incident[1].sense != 0)
                return IncidenceFault::OrientationMismatch;
            ++liveSeen;
        }
        if (liveSeen != liveLinks_) return IncidenceFault::IncidenceMismatch;
        return IncidenceFault::Ok;
    }

#if defined(PAACS_TOPO_TESTING)
    using AdoptionBeginHook = void (*)() noexcept;
    static void testSetAdoptionBeginHook(AdoptionBeginHook hook) noexcept {
        adoptionBeginHook_ = hook;
    }
    std::size_t testIncidenceCount(const FaceHandle& h) const {
        return rec_[static_cast<std::size_t>(alloc_.deref(h))]
            .data.incident.size();
    }
    void testSetTail(std::int32_t s) { tail_ = s; }
    void testSetHead(std::int32_t s) { head_ = s; }
    void testSetNext(std::int32_t slot, std::int32_t next) {
        rec_[static_cast<std::size_t>(slot)].next = next;
    }
    void testSetPrev(std::int32_t slot, std::int32_t prev) {
        rec_[static_cast<std::size_t>(slot)].prev = prev;
    }
    std::int32_t testHead() const { return head_; }
    std::int32_t testMaterialized() const {
        return static_cast<std::int32_t>(alloc_.materializedSlots());
    }
    void testCorruptTableEntry(const FaceHandle& live, const FaceHandle& wrong) {
        const Rec& r = rec_[static_cast<std::size_t>(raw::slotUnchecked(live))];
        table_[r.key] = wrong;
    }
    void testDropTableEntry(const FaceHandle& live) {
        const Rec& r = rec_[static_cast<std::size_t>(raw::slotUnchecked(live))];
        table_.erase(r.key);
    }
    void testFlipSense(const FaceHandle& h, std::size_t idx) {
        auto& inc = rec_[static_cast<std::size_t>(raw::slotUnchecked(h))]
                        .data.incident;
        if (idx < inc.size()) inc[idx].sense = static_cast<std::int8_t>(-inc[idx].sense);
    }
#endif

private:
#if defined(PAACS_TOPO_TESTING)
    inline static AdoptionBeginHook adoptionBeginHook_ = nullptr;
#endif
    struct TokenlessAdoptionTag {};
    struct CompleteBuildTag {};

    struct Staged {
        FaceKey key;
        FaceShape shape;
        std::uint8_t n;
        std::array<NodeHandle, kMaxFaceNodes> outward;  // contributing winding
    };
    struct Rec {
        FaceData data{};
        FaceKey key{};
        std::uint64_t seq = 0;
        std::uint32_t gen = 0;
        std::int32_t prev = detail::kNil;
        std::int32_t next = detail::kNil;
    };
    struct CellFaceEntry {
        std::uint32_t gen = 0;
        std::vector<FaceHandle> faces;
    };

    // Stage 1 prepares ordinary empty containers and may propagate their
    // allocation failure before the source is touched. The embedded
    // allocator is constructed tokenless, so no pool identity is consumed.
    explicit FacePool(TokenlessAdoptionTag)
        : alloc_(HandleAllocator<FaceTag>::TokenlessAdoptionTag{}) {}

    FacePool(CompleteBuildTag,const NodePool& nodes,const CellPool& cells,
             const EdgePool& edges,FacePoolCompleteBuildMetrics* metrics) {
        populateComplete(nodes,cells,edges,metrics);
    }

    struct CompleteKeyLess {
        bool operator()(const FaceKey& a,const FaceKey& b) const noexcept {
            if (a.n != b.n) return a.n < b.n;
            for (std::size_t i = 0; i < a.n; ++i) {
                if (a.s[i] != b.s[i]) return a.s[i] < b.s[i];
            }
            return false;
        }
    };

    static void completeMetricAdd(std::uint64_t& value,
                                  std::uint64_t amount) {
        if (amount > (std::numeric_limits<std::uint64_t>::max)() - value)
            throw std::length_error(
                "FacePool::buildComplete: metric count overflow");
        value += amount;
    }

    Staged completeStage(const NodePool& nodes,const CellPool::CellSnapshot& snap,
                         const LocalFace& localFace) const {
        Staged staged;
        staged.shape = localFace.shape;
        staged.n = localFace.n;
        for (int i = 0; i < localFace.n; ++i) {
            const NodeHandle node =
                snap.nodes[localFace.v[static_cast<std::size_t>(i)]];
            if (!nodes.isCurrent(node))
                throw std::invalid_argument(
                    "FacePool::buildComplete: stale face node");
            staged.outward[static_cast<std::size_t>(i)] = node;
        }
        if (!makeKey(nodes,staged.outward,staged.n,staged.key))
            throw std::invalid_argument(
                "FacePool::buildComplete: degenerate face");
        return staged;
    }

    void populateComplete(const NodePool& nodes,const CellPool& cells,
                          const EdgePool& edges,
                          FacePoolCompleteBuildMetrics* publishedMetrics) {
        FacePoolCompleteBuildMetrics metrics;
        const auto active = cells.activeBySequence();
        std::map<FaceKey,std::uint64_t,CompleteKeyLess> multiplicity;

        // Pass 1 validates cells, their complete edge incidence, all face
        // nodes/keys, and exact incidence multiplicities.
        for (const CellHandle cell : active) {
            ++metrics.validationCellVisits;
            if (!cells.isCurrent(cell))
                throw std::invalid_argument(
                    "FacePool::buildComplete: stale cell handle");
            if (!edges.hasCellEdges(cells,cell))
                throw std::invalid_argument(
                    "FacePool::buildComplete: missing cell edge incidence");
            const auto snap = cells.get(cell);
            const ElementTraits& traits = checkedTraits(snap.type);
            for (int face = 0; face < traits.faceCount; ++face) {
                const Staged staged = completeStage(
                    nodes,snap,traits.faces[face]);
                auto [it,inserted] = multiplicity.emplace(staged.key,0);
                (void)inserted;
                if (it->second ==
                    (std::numeric_limits<std::size_t>::max)())
                    throw std::length_error(
                        "FacePool::buildComplete: face incidence overflow");
                ++it->second;
                completeMetricAdd(metrics.faceUses,1);
            }
        }
        metrics.uniqueFaces = multiplicity.size();
        if (multiplicity.size() > static_cast<std::size_t>(
                (std::numeric_limits<std::int32_t>::max)()))
            throw std::length_error(
                "FacePool::buildComplete: face slot limit exceeded");

#if defined(PAACS_TOPO_TESTING)
        maybeFault();
#endif
        alloc_.reserveMaterializedCapacity(multiplicity.size());
        rec_.resize(multiplicity.size());
        table_.reserve(multiplicity.size());
        cellFaces_.reserve(active.size());
        if (!active.empty()) {
            boundNodePool_ = nodes.poolId();
            boundCellPool_ = cells.poolId();
        }

        // Pass 2 alone determines handles and incidence order.
        for (const CellHandle cell : active) {
            ++metrics.buildCellVisits;
            const auto snap = cells.get(cell);
            const ElementTraits& traits = checkedTraits(snap.type);
            CellFaceEntry entry;
            entry.gen = cell.gen;
            entry.faces.resize(static_cast<std::size_t>(traits.faceCount));
            for (int face = 0; face < traits.faceCount; ++face) {
                const Staged staged = completeStage(
                    nodes,snap,traits.faces[face]);
                auto found = table_.find(staged.key);
                FaceHandle handle;
                if (found == table_.end()) {
                    auto [tableIt,inserted] =
                        table_.emplace(staged.key,FaceHandle{});
                    if (!inserted)
                        throw std::logic_error(
                            "FacePool::buildComplete: duplicate table row");
                    handle = alloc_.alloc();
                    const auto slot = static_cast<std::size_t>(
                        raw::slotUnchecked(handle));
                    if (slot >= rec_.size())
                        throw std::logic_error(
                            "FacePool::buildComplete: unprepared slot");
                    Rec& record = rec_[slot];
                    record.data = FaceData{};
                    record.data.shape = staged.shape;
                    record.data.cornerCount = staged.n;
                    for (int i = 0; i < staged.n; ++i)
                        record.data.corners[static_cast<std::size_t>(i)] =
                            staged.outward[static_cast<std::size_t>(i)];
                    const auto count = multiplicity.find(staged.key);
                    if (count == multiplicity.end() ||
                        count->second > record.data.incident.max_size())
                        throw std::length_error(
                            "FacePool::buildComplete: incidence capacity");
                    record.data.incident.reserve(
                        static_cast<std::size_t>(count->second));
                    record.key = staged.key;
                    record.seq = alloc_.seqOf(handle);
                    record.gen = handle.gen;
                    linkTail(static_cast<std::int32_t>(slot));
                    tableIt->second = handle;
                    bumpPeaks();
                } else {
                    handle = found->second;
                }
                Rec& record = rec_[static_cast<std::size_t>(
                    raw::slotUnchecked(handle))];
                const std::int8_t sense = record.data.incident.empty()
                    ? static_cast<std::int8_t>(1)
                    : orientationSense(nodes,record.data,staged);
                if (sense != 1 && sense != -1)
                    throw std::invalid_argument(
                        "FacePool::buildComplete: invalid face orientation");
                record.data.incident.push_back(
                    Incidence{cell.slot,cell.gen,sense});
                entry.faces[static_cast<std::size_t>(face)] = handle;
            }
            const auto [it,inserted] =
                cellFaces_.emplace(cell.slot,std::move(entry));
            (void)it;
            if (!inserted)
                throw std::logic_error(
                    "FacePool::buildComplete: duplicate cell incidence");
        }
        for (const auto& [key,count] : multiplicity) {
            const auto found = table_.find(key);
            if (found == table_.end() ||
                rec_[static_cast<std::size_t>(raw::slotUnchecked(
                    found->second))].data.incident.size() != count)
                throw std::logic_error(
                    "FacePool::buildComplete: multiplicity mismatch");
        }
        // auditFaceTopology performs the one complete supplied-edge audit
        // and the one complete final face audit.
        ++metrics.suppliedEdgeAudits;
        ++metrics.finalFaceAudits;
        if (!auditFaceTopology(nodes,cells,edges,*this)
                 .acceptable(diag::Severity::Warning))
            throw std::invalid_argument(
                "FacePool::buildComplete: final audit failed");
        metrics.runComplete = true;
        if (publishedMetrics != nullptr) *publishedMetrics = metrics;
    }

    // Stage 2 transfers the exact allocator, records, canonical table,
    // per-cell table, face incidence vectors/senses, sequence list, bindings
    // and counters by compile-time-verified noexcept swaps only.
    void adoptFrom(FacePool& o) noexcept {
        alloc_.adoptFrom(o.alloc_);
        static_assert(noexcept(rec_.swap(o.rec_)));
        static_assert(noexcept(table_.swap(o.table_)));
        static_assert(noexcept(cellFaces_.swap(o.cellFaces_)));
        rec_.swap(o.rec_);
        table_.swap(o.table_);
        cellFaces_.swap(o.cellFaces_);
        using std::swap;
        swap(head_, o.head_);
        swap(tail_, o.tail_);
        swap(liveLinks_, o.liveLinks_);
        swap(boundNodePool_, o.boundNodePool_);
        swap(boundCellPool_, o.boundCellPool_);
        swap(peakLive_, o.peakLive_);
        swap(peakMaterialized_, o.peakMaterialized_);
    }

    void bind(const NodePool& nodes, const CellPool& cells) {
        if (boundNodePool_ == 0) {
            boundNodePool_ = nodes.poolId();
            boundCellPool_ = cells.poolId();
        } else if (boundNodePool_ != nodes.poolId() ||
                   boundCellPool_ != cells.poolId()) {
            throw std::invalid_argument(
                "FacePool: bound to different node/cell pools");
        }
    }
    // Canonical key = sorted corner sequences; false on any repeat.
    bool makeKey(const NodePool& nodes,
                 const std::array<NodeHandle, kMaxFaceNodes>& c,
                 std::uint8_t n, FaceKey& out) const {
        // Explicit bounds proof (reviewer-sanctioned post-freeze cleanup):
        // GCC 14 -Wall -Wextra -Wpedantic emitted array-bounds false
        // positives here because it could not prove n <= kMaxFaceNodes.
        // The guard is also correct defensively: a face has 3 or 4 corners.
        if (n < 3 || n > kMaxFaceNodes) return false;
        const std::size_t count = n;
        out.n = n;
        for (std::size_t i = 0; i < count; ++i) out.s[i] = nodes.seqOf(c[i]);
        if (!core::sortSmallPrefix(out.s, count)) return false;
        for (std::size_t i = 1; i < count; ++i)
            if (out.s[i] == out.s[i-1]) return false;   // repeated node
        return true;
    }
    FaceHandle lookup(const FaceKey& key) const {
        const auto it = table_.find(key);
        return (it == table_.end()) ? FaceHandle{} : it->second;
    }
    // Relative orientation of a contributing cell's outward winding vs the
    // face's stored canonical winding. Same cyclic order (allowing rotation)
    // => +1; reversed cyclic order => -1; otherwise 0 (shouldn't happen for
    // a shared face with identical corner set).
    std::int8_t orientationSense(
        const NodePool& nodes, const FaceData& fd, const auto& staged) const {
        const int n = fd.cornerCount;
        if (staged.n != fd.cornerCount) return 0;
        // Sequences of canonical winding and the candidate winding.
        std::array<std::uint64_t, kMaxFaceNodes> canon{}, cand{};
        for (int i = 0; i < n; ++i) {
            canon[static_cast<std::size_t>(i)] = nodes.seqOf(fd.corners[static_cast<std::size_t>(i)]);
            cand[static_cast<std::size_t>(i)] = nodes.seqOf(staged.outward[static_cast<std::size_t>(i)]);
        }
        // Find candidate start matching canon[0].
        int start = -1;
        for (int i = 0; i < n; ++i) if (cand[static_cast<std::size_t>(i)] == canon[0]) { start = i; break; }
        if (start < 0) return 0;
        bool forward = true, backward = true;
        for (int i = 0; i < n; ++i) {
            if (cand[static_cast<std::size_t>((start + i) % n)] != canon[static_cast<std::size_t>(i)]) forward = false;
            if (cand[static_cast<std::size_t>((start - i + n) % n)] != canon[static_cast<std::size_t>(i)]) backward = false;
        }
        if (forward) return +1;
        if (backward) return -1;
        return 0;
    }
    FaceHandle createNew(const NodePool& nodes, const auto& staged,
                         std::vector<FaceHandle>& createdOut) {
        (void)nodes;
        if (alloc_.reusableFreeSlots() == 0 &&
            HandleAllocator<FaceTag>::wouldExhaustSlots(
                alloc_.materializedSlots()))
            throw std::length_error("FacePool: slot space exhausted");
        const std::size_t worst = alloc_.materializedSlots() + 1;
        if (worst > rec_.size()) {
            detail::reserveFor(rec_, worst);
            rec_.resize(worst);
        }
        auto [tabIt, inserted] = table_.emplace(staged.key, FaceHandle{});
        (void)inserted;
        try {
            const FaceHandle h = alloc_.alloc();
            const auto s = static_cast<std::size_t>(raw::slotUnchecked(h));
            if (s >= rec_.size()) rec_.resize(s + 1);
            Rec& r = rec_[s];
            r.data = FaceData{};
            r.data.shape = staged.shape;
            r.data.cornerCount = staged.n;
            for (int i = 0; i < staged.n; ++i)
                r.data.corners[static_cast<std::size_t>(i)] = staged.outward[static_cast<std::size_t>(i)];
            r.key = staged.key;
            r.seq = alloc_.seqOf(h);
            r.gen = h.gen;
            linkTail(static_cast<std::int32_t>(s));
            tabIt->second = h;
            createdOut.push_back(h);
            bumpPeaks();
            return h;
        } catch (...) {
            table_.erase(tabIt);
            throw;
        }
    }
    void retireFace(const FaceHandle& h) noexcept {
        const auto s = static_cast<std::int32_t>(raw::slotUnchecked(h));
        table_.erase(rec_[static_cast<std::size_t>(s)].key);
        unlink(s);
        alloc_.free(h);       // noexcept: free-list pre-reserved by callers
    }
    void linkTail(std::int32_t s) {
        rec_[s].prev = tail_;
        rec_[s].next = detail::kNil;
        if (tail_ != detail::kNil) rec_[tail_].next = s;
        else head_ = s;
        tail_ = s;
        ++liveLinks_;
    }
    void unlink(std::int32_t s) {
        Rec& r = rec_[static_cast<std::size_t>(s)];
        if (r.prev != detail::kNil) rec_[r.prev].next = r.next;
        else head_ = r.next;
        if (r.next != detail::kNil) rec_[r.next].prev = r.prev;
        else tail_ = r.prev;
        r.prev = r.next = detail::kNil;
        --liveLinks_;
    }
    bool checkListAndTable(const NodePool& nodes) const {
        if (alloc_.materializedSlots() !=
            alloc_.liveCount() + alloc_.reusableFreeSlots() +
                alloc_.retiredSlots())
            return false;
        if (liveLinks_ != alloc_.liveCount()) return false;
        if (rec_.size() < alloc_.materializedSlots()) return false;
        const auto validSlot = [&](std::int32_t s) {
            return s >= 0 && static_cast<std::size_t>(s) < rec_.size();
        };
        if (liveLinks_ == 0)
            return head_ == detail::kNil && tail_ == detail::kNil &&
                   table_.empty();
        if (!validSlot(head_) || !validSlot(tail_)) return false;
        if (rec_[head_].prev != detail::kNil) return false;
        if (rec_[tail_].next != detail::kNil) return false;
        std::vector<std::uint8_t> seen(rec_.size(), 0);
        std::int32_t prevSlot = detail::kNil;
        std::uint64_t prevSeq = 0;
        std::size_t walked = 0;
        for (std::int32_t s = head_; s != detail::kNil; ) {
            if (!validSlot(s)) return false;
            if (seen[static_cast<std::size_t>(s)]) return false;
            seen[static_cast<std::size_t>(s)] = 1;
            const Rec& r = rec_[static_cast<std::size_t>(s)];
            if (r.prev != detail::kNil && !validSlot(r.prev)) return false;
            if (r.prev != prevSlot) return false;
            if (r.next != detail::kNil && !validSlot(r.next)) return false;
            const FaceHandle h{s, r.gen, alloc_.poolId()};
            if (!alloc_.isCurrent(h)) return false;
            if (alloc_.seqOf(h) != r.seq) return false;
            if (r.data.incident.empty()) return false;   // exists iff used
            if (r.data.cornerCount < 3 || r.data.cornerCount > kMaxFaceNodes)
                return false;
            for (std::uint8_t i = 0; i < r.data.cornerCount; ++i)
                if (!nodes.isCurrent(r.data.corners[i])) return false;
            // key must match the sorted corner sequences.
            FaceKey k;
            std::array<NodeHandle, kMaxFaceNodes> c{};
            for (std::uint8_t i = 0; i < r.data.cornerCount; ++i) c[i] = r.data.corners[i];
            if (!makeKey(nodes, c, r.data.cornerCount, k)) return false;
            if (!(k == r.key)) return false;
            const auto tit = table_.find(r.key);
            if (tit == table_.end()) return false;
            if (raw::slotUnchecked(tit->second) != s) return false;
            if (walked > 0 && !(r.seq > prevSeq)) return false;
            prevSeq = r.seq;
            prevSlot = s;
            if (++walked > liveLinks_) return false;
            s = r.next;
        }
        if (prevSlot != tail_) return false;
        if (walked != liveLinks_) return false;
        return table_.size() == liveLinks_;
    }
    bool checkIncidence(const NodePool& nodes, const CellPool& cells) const {
        if (boundCellPool_ != 0 && cells.poolId() != boundCellPool_)
            return false;
        if (boundNodePool_ != 0 && nodes.poolId() != boundNodePool_)
            return false;
        // SEMANTIC recomputation (review blockers 2 & 3): for every current
        // incident cell, recompute its traits-extracted faces, build each
        // key, look it up in the table, and verify the recorded FaceHandle
        // AND orientation sense match. Also count incidence per face slot.
        std::unordered_map<std::int32_t, std::size_t> counted;
        for (const auto& [slot, entry] : cellFaces_) {
            const CellHandle ch{slot, entry.gen, cells.poolId()};
            if (!cells.isCurrent(ch)) return false;      // ghost incidence
            const auto snap = cells.get(ch);
            if (!isValidElementType(snap.type)) return false;
            const ElementTraits& tr = checkedTraits(snap.type);
            if (static_cast<int>(entry.faces.size()) != tr.faceCount)
                return false;
            for (int f = 0; f < tr.faceCount; ++f) {
                const LocalFace& lf = tr.faces[f];
                // Rebuild this cell's outward winding + key for face f.
                Staged st;
                st.shape = lf.shape;
                st.n = lf.n;
                for (int i = 0; i < lf.n; ++i) {
                    const NodeHandle nh =
                        snap.nodes[lf.v[static_cast<std::size_t>(i)]];
                    if (!nodes.isCurrent(nh)) return false;
                    st.outward[static_cast<std::size_t>(i)] = nh;
                }
                if (!makeKey(nodes, st.outward, lf.n, st.key)) return false;
                // The recorded face handle for this local face must be the
                // table entry for the recomputed key.
                const FaceHandle recorded =
                    entry.faces[static_cast<std::size_t>(f)];
                if (!alloc_.isCurrent(recorded)) return false;
                const auto tit = table_.find(st.key);
                if (tit == table_.end()) return false;
                if (raw::slotUnchecked(tit->second) !=
                    raw::slotUnchecked(recorded))
                    return false;                        // wrong face handle
                const Rec& fr = rec_[static_cast<std::size_t>(
                    raw::slotUnchecked(recorded))];
                // The recorded incidence entry for THIS cell must exist with
                // the recomputed orientation sense.
                const std::int8_t expectSense =
                    orientationSense(nodes, fr.data, st);
                bool found = false;
                for (const Incidence& in : fr.data.incident)
                    if (in.cellSlot == slot && in.cellGen == entry.gen) {
                        if (in.sense != expectSense) return false;  // sense
                        found = true;
                        break;
                    }
                if (!found) return false;                // missing incidence
                ++counted[raw::slotUnchecked(recorded)];
            }
        }
        std::size_t liveSeen = 0;
        for (std::int32_t s = head_; s != detail::kNil; s = rec_[s].next) {
            const Rec& r = rec_[static_cast<std::size_t>(s)];
            const auto it = counted.find(s);
            const std::size_t inc = (it == counted.end()) ? 0 : it->second;
            if (r.data.incident.size() != inc) return false;
            if (inc == 0) return false;                  // ghost face
            for (const Incidence& in : r.data.incident)
                if (in.sense != +1 && in.sense != -1) return false;
            // A 2-cell INTERIOR face must have OPPOSITE senses (review
            // blocker 2): a same-sense pair is a real orientation defect.
            if (r.data.incident.size() == 2 &&
                r.data.incident[0].sense + r.data.incident[1].sense != 0)
                return false;
            ++liveSeen;
        }
        return liveSeen == liveLinks_;
    }
    void bumpPeaks() {
        peakLive_ = std::max(peakLive_, alloc_.liveCount());
        peakMaterialized_ = std::max(peakMaterialized_,
                                     alloc_.materializedSlots());
    }

    HandleAllocator<FaceTag> alloc_;
    std::vector<Rec> rec_;
    std::unordered_map<FaceKey, FaceHandle, FaceKeyHash> table_;
    std::unordered_map<std::int32_t, CellFaceEntry> cellFaces_;
    std::int32_t head_ = detail::kNil;
    std::int32_t tail_ = detail::kNil;
    std::size_t liveLinks_ = 0;
    std::uint64_t boundNodePool_ = 0;
    std::uint64_t boundCellPool_ = 0;
    std::size_t peakLive_ = 0;
    std::size_t peakMaterialized_ = 0;

    static_assert(std::is_nothrow_swappable_v<decltype(rec_)>);
    static_assert(std::is_nothrow_swappable_v<decltype(table_)>);
    static_assert(std::is_nothrow_swappable_v<decltype(cellFaces_)>);
};

static_assert(!std::is_nothrow_move_constructible_v<FacePool>);
static_assert(std::is_nothrow_move_constructible_v<
              HandleAllocator<FaceTag>>);

// Structured face-topology audit (mirrors auditEdgeTopology). Requires edge
// topology to be ready first (faces build on trustworthy edges).
inline diag::Report auditFaceTopology(const NodePool& nodes,
                                      const CellPool& cells,
                                      const EdgePool& edges,
                                      const FacePool& faces) {
    diag::Report r;
    auto add = [&](diag::Code code, diag::Severity sev, std::string msg) {
        diag::Diagnostic d;
        d.code = code; d.severity = sev; d.stage = "face-audit";
        d.message = std::move(msg);
        r.add(std::move(d));
    };
    if ((faces.boundNodePool() != 0 && faces.boundNodePool() != nodes.poolId()) ||
        (faces.boundCellPool() != 0 && faces.boundCellPool() != cells.poolId())) {
        add(diag::Code::FacePoolBindingMismatch, diag::Severity::InvalidOutput,
            "FacePool bound to different node/cell pools than audited");
        return r;
    }
    if (!edgeTopologyReady(nodes, cells, edges)) {
        add(diag::Code::EdgeInvariantViolation, diag::Severity::InvalidOutput,
            "edge topology is not ready; faces require valid edges first");
        return r;
    }
    if (!faces.checkInvariants(nodes)) {
        add(diag::Code::FaceInvariantViolation, diag::Severity::InvalidOutput,
            "face records / list / table invariants violated");
        return r;
    }
    // Specific cell-aware diagnostics (review: the enum promised distinct
    // orientation/incidence codes; use them).
    using IF = FacePool::IncidenceFault;
    switch (faces.checkIncidenceReason(nodes, cells)) {
        case IF::Ok: break;
        case IF::BindingMismatch:
            add(diag::Code::FacePoolBindingMismatch, diag::Severity::InvalidOutput,
                "face incidence pool binding mismatch");
            return r;
        case IF::GhostFace:
            add(diag::Code::GhostFace, diag::Severity::InvalidOutput,
                "a live face lacks current-cell backing (ghost face)");
            return r;
        case IF::OrientationMismatch:
            add(diag::Code::FaceOrientationMismatch, diag::Severity::InvalidOutput,
                "a 2-cell interior face has non-opposite senses, or a stored "
                "sense disagrees with the recomputed orientation");
            return r;
        case IF::IncidenceMismatch:
            add(diag::Code::FaceIncidenceMismatch, diag::Severity::InvalidOutput,
                "recorded face/handle does not match the cell's recomputed "
                "traits face");
            return r;
    }
    // Every face-incident cell must itself have current edge incidence in
    // the EdgePool (review blocker: otherwise faces can sit on edge records
    // created by a different but geometrically identical cell).
    for (const CellHandle& ch : faces.incidentCells()) {
        if (!edges.hasCellEdges(cells, ch)) {
            add(diag::Code::FaceIncidenceMismatch, diag::Severity::InvalidOutput,
                "a face-incident cell has no edge incidence in the EdgePool");
            return r;
        }
    }
    // Explicit nonmanifold IDENTIFICATION (M4-B.6 stage 4, disclosed
    // additive diagnosis in this frozen file): faces with >= 3 incident
    // cells are representable by design (M4-B.4: flagged, not rejected),
    // and the reserved NonManifoldFace code now names each one at
    // Severity::Info - purely informational, never gating acceptability,
    // so downstream policy can distinguish "representable" from
    // "acceptable under the current operation policy".
    for (const FaceHandle& fh : faces.activeBySequence()) {
        const auto d = faces.get(fh);
        if (d.incident.size() >= 3)
            add(diag::Code::NonManifoldFace, diag::Severity::Info,
                "face has " + std::to_string(d.incident.size()) +
                    " incident cells (nonmanifold; representable - "
                    "acceptability is an operation-policy decision)");
    }
    return r;
}

inline bool faceTopologyReady(const NodePool& nodes, const CellPool& cells,
                              const EdgePool& edges, const FacePool& faces) {
    return auditFaceTopology(nodes, cells, edges, faces)
        .acceptable(diag::Severity::Warning);
}

} // namespace paacs::topo
