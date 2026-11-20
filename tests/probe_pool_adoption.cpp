// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Increment-5 direct pool-move allocation probe.
//
// Exhaustively injects every allocation while EdgePool/FacePool prepare
// their private tokenless destination. A preparation failure must propagate
// before the source changes. A test-only hook marks the first source mutation;
// the global allocation counter must not advance after that boundary.
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <optional>
#include <vector>

#include "paacs/topo/FaceTopology.hpp"

static long gAdoptionFailCountdown = 0;
static long gAllocationCalls = 0;
static long gCallsAtAdoptionBegin = -1;

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t n) {
    ++gAllocationCalls;
    if (gAdoptionFailCountdown > 0 && --gAdoptionFailCountdown == 0)
        throw std::bad_alloc();
    if (void* p = std::malloc(n ? n : 1)) return p;
    throw std::bad_alloc();
}
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace {
using namespace paacs::topo;

int gBad = 0;
#define REQUIRE(c)                                                        \
    do {                                                                  \
        if (!(c)) {                                                       \
            std::printf("FAIL %d: %s\n", __LINE__, #c);                 \
            ++gBad;                                                       \
        }                                                                 \
    } while (0)

void adoptionBegin() noexcept { gCallsAtAdoptionBegin = gAllocationCalls; }

void appendHandle(std::vector<std::uint64_t>& out, const auto& h) {
    out.push_back(static_cast<std::uint64_t>(h.slot));
    out.push_back(static_cast<std::uint64_t>(h.gen));
    out.push_back(h.pool);
}

void appendStats(std::vector<std::uint64_t>& out, const PoolStats& s) {
    out.push_back(s.materializedSlots);
    out.push_back(s.live);
    out.push_back(s.reusableFree);
    out.push_back(s.retired);
    out.push_back(s.totalCreated);
    out.push_back(s.peakLive);
    out.push_back(s.peakMaterialized);
}

std::vector<std::uint64_t> edgeTokens(
    const EdgePool& edges, const CellPool& cells,
    const std::vector<CellHandle>& fixtureCells) {
    std::vector<std::uint64_t> out;
    out.push_back(edges.poolId());
    out.push_back(edges.boundNodePool());
    out.push_back(edges.boundCellPool());
    appendStats(out, edges.stats());
    const auto active = edges.activeBySequence();
    out.push_back(active.size());
    for (const auto& h : active) {
        appendHandle(out, h);
        out.push_back(edges.seqOf(h));
        const auto d = edges.get(h);
        appendHandle(out, d.a);
        appendHandle(out, d.b);
        out.push_back(static_cast<std::uint64_t>(d.cls.dim));
        out.push_back(static_cast<std::uint64_t>(d.cls.tag));
        out.push_back(d.requiredConstraint ? 1U : 0U);
        out.push_back(d.useCount);
    }
    for (const auto& cell : fixtureCells) {
        appendHandle(out, cell);
        out.push_back(edges.hasCellEdges(cells, cell) ? 1U : 0U);
        const auto incident = edges.cellEdges(cells, cell);
        out.push_back(incident.size());
        for (const auto& h : incident) appendHandle(out, h);
    }
    return out;
}

std::vector<std::uint64_t> faceTokens(
    const FacePool& faces, const CellPool& cells,
    const std::vector<CellHandle>& fixtureCells) {
    std::vector<std::uint64_t> out;
    out.push_back(faces.poolId());
    out.push_back(faces.boundNodePool());
    out.push_back(faces.boundCellPool());
    appendStats(out, faces.stats());
    const auto active = faces.activeBySequence();
    out.push_back(active.size());
    for (const auto& h : active) {
        appendHandle(out, h);
        out.push_back(faces.seqOf(h));
        const auto d = faces.get(h);
        out.push_back(static_cast<std::uint64_t>(d.shape));
        out.push_back(d.cornerCount);
        for (std::uint8_t i = 0; i < d.cornerCount; ++i)
            appendHandle(out, d.corners[i]);
        out.push_back(static_cast<std::uint64_t>(d.cls.dim));
        out.push_back(static_cast<std::uint64_t>(d.cls.tag));
        out.push_back(d.incident.size());
        for (const auto& incidence : d.incident) {
            out.push_back(static_cast<std::uint64_t>(incidence.cellSlot));
            out.push_back(incidence.cellGen);
            out.push_back(static_cast<std::uint64_t>(incidence.sense));
        }
    }
    for (const auto& cell : fixtureCells) {
        appendHandle(out, cell);
        out.push_back(faces.hasCellFaces(cells, cell) ? 1U : 0U);
        const auto incident = faces.cellFaces(cells, cell);
        out.push_back(incident.size());
        for (const auto& h : incident) appendHandle(out, h);
    }
    return out;
}

void makeCells(NodePool& nodes, std::vector<NodeHandle>& handles,
               CellPool& cells, std::vector<CellHandle>& made,
               int count) {
    for (int i = 0; i < 6; ++i)
        handles.push_back(nodes.add(
            {double(i), double(i % 3), double(i % 2)}, {3, 1}));
    made.push_back(cells.add(
        ElementType::Tet4,
        {handles[0], handles[1], handles[2], handles[3]}, {3, 1}, nodes));
    if (count >= 2)
        made.push_back(cells.add(
            ElementType::Tet4,
            {handles[0], handles[2], handles[1], handles[4]}, {3, 2}, nodes));
    if (count >= 3)
        made.push_back(cells.add(
            ElementType::Tet4,
            {handles[0], handles[1], handles[2], handles[5]}, {3, 3}, nodes));
}

long edgeMoveSweep(bool populated) {
    long faults = 0;
    for (long k = 1;; ++k) {
        NodePool nodes;
        std::vector<NodeHandle> handles;
        CellPool cells;
        std::vector<CellHandle> made;
        if (populated) makeCells(nodes, handles, cells, made, 2);
        EdgePool source;
        if (populated) {
            for (const auto& cell : made)
                source.insertCellEdges(nodes, cells, cell);
            const auto required = source.find(nodes, handles[0], handles[1]);
            REQUIRE(required.has_value());
            source.setRequiredConstraint(*required, true);
            source.setClassification(*required, {1, 71});
        }
        const auto before = edgeTokens(source, cells, made);
        std::optional<EdgePool> target;
        gCallsAtAdoptionBegin = -1;
        EdgePool::testSetAdoptionBeginHook(adoptionBegin);
        gAdoptionFailCountdown = k;
        bool badAlloc = false;
        try {
            target.emplace(std::move(source));
        } catch (const std::bad_alloc&) {
            badAlloc = true;
        } catch (...) {
            REQUIRE(false);
        }
        gAdoptionFailCountdown = 0;
        EdgePool::testSetAdoptionBeginHook(nullptr);
        if (!badAlloc) {
            REQUIRE(target.has_value());
            REQUIRE(gCallsAtAdoptionBegin >= 0);
            REQUIRE(gAllocationCalls == gCallsAtAdoptionBegin);
            REQUIRE(edgeTokens(*target, cells, made) == before);
            REQUIRE(target->checkInvariants(nodes, cells));
            REQUIRE(source.poolId() == 0 && source.liveCount() == 0);
            REQUIRE(source.stats().materializedSlots == 0);
            break;
        }
        ++faults;
        REQUIRE(gCallsAtAdoptionBegin < 0);
        REQUIRE(edgeTokens(source, cells, made) == before);
        REQUIRE(source.checkInvariants(nodes, cells));
        EdgePool retry(std::move(source));
        REQUIRE(edgeTokens(retry, cells, made) == before);
        REQUIRE(retry.checkInvariants(nodes, cells));
    }
    std::printf("pool-adoption edge %s: %ld preparation faults / 0 escapes\n",
                populated ? "populated" : "empty", faults);
    return faults;
}

long faceMoveSweep(bool populated) {
    long faults = 0;
    for (long k = 1;; ++k) {
        NodePool nodes;
        std::vector<NodeHandle> handles;
        CellPool cells;
        std::vector<CellHandle> made;
        if (populated) makeCells(nodes, handles, cells, made, 3);
        EdgePool edges;
        if (populated)
            for (const auto& cell : made)
                edges.insertCellEdges(nodes, cells, cell);
        FacePool source;
        if (populated)
            for (const auto& cell : made)
                source.insertCellFaces(nodes, cells, edges, cell);
        const auto before = faceTokens(source, cells, made);
        std::optional<FacePool> target;
        gCallsAtAdoptionBegin = -1;
        FacePool::testSetAdoptionBeginHook(adoptionBegin);
        gAdoptionFailCountdown = k;
        bool badAlloc = false;
        try {
            target.emplace(std::move(source));
        } catch (const std::bad_alloc&) {
            badAlloc = true;
        } catch (...) {
            REQUIRE(false);
        }
        gAdoptionFailCountdown = 0;
        FacePool::testSetAdoptionBeginHook(nullptr);
        if (!badAlloc) {
            REQUIRE(target.has_value());
            REQUIRE(gCallsAtAdoptionBegin >= 0);
            REQUIRE(gAllocationCalls == gCallsAtAdoptionBegin);
            REQUIRE(faceTokens(*target, cells, made) == before);
            REQUIRE(target->checkInvariants(nodes, cells));
            REQUIRE(source.poolId() == 0 && source.liveCount() == 0);
            REQUIRE(source.stats().materializedSlots == 0);
            break;
        }
        ++faults;
        REQUIRE(gCallsAtAdoptionBegin < 0);
        REQUIRE(faceTokens(source, cells, made) == before);
        REQUIRE(source.checkInvariants(nodes, cells));
        FacePool retry(std::move(source));
        REQUIRE(faceTokens(retry, cells, made) == before);
        REQUIRE(retry.checkInvariants(nodes, cells));
    }
    std::printf("pool-adoption face %s: %ld preparation faults / 0 escapes\n",
                populated ? "populated" : "empty", faults);
    return faults;
}

void allocatorMovesAllocateNothing() {
    static_assert(std::is_nothrow_move_constructible_v<
                  HandleAllocator<EdgeTag>>);
    static_assert(std::is_nothrow_move_constructible_v<
                  HandleAllocator<FaceTag>>);
    HandleAllocator<EdgeTag> edgeSource;
    const auto edge = edgeSource.alloc();
    const auto beforeEdgeCalls = gAllocationCalls;
    HandleAllocator<EdgeTag> edgeTarget(std::move(edgeSource));
    REQUIRE(gAllocationCalls == beforeEdgeCalls);
    REQUIRE(edgeTarget.isCurrent(edge));
    REQUIRE(edgeSource.poolId() == 0);

    HandleAllocator<FaceTag> faceSource;
    const auto face = faceSource.alloc();
    const auto beforeFaceCalls = gAllocationCalls;
    HandleAllocator<FaceTag> faceTarget(std::move(faceSource));
    REQUIRE(gAllocationCalls == beforeFaceCalls);
    REQUIRE(faceTarget.isCurrent(face));
    REQUIRE(faceSource.poolId() == 0);
}

}  // namespace

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    allocatorMovesAllocateNothing();
    const long edgeEmpty = edgeMoveSweep(false);
    const long edgePopulated = edgeMoveSweep(true);
    const long faceEmpty = faceMoveSweep(false);
    const long facePopulated = faceMoveSweep(true);
    (void)edgeEmpty;
    (void)edgePopulated;
    (void)faceEmpty;
    (void)facePopulated;
#if defined(_MSC_VER)
    REQUIRE(edgeEmpty > 0 && edgePopulated > 0);
    REQUIRE(faceEmpty > 0 && facePopulated > 0);
#endif
    if (gBad) {
        std::printf("POOL ADOPTION PROBE FAILED (%d)\n", gBad);
        return 1;
    }
    std::printf("POOL ADOPTION PROBE PASSED\n");
    return 0;
}
