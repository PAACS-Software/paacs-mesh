// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Standalone M4-B.5 adjacency allocation-failure probe.
//
// Exhaustively injects every allocation in CellAdjacency::rebuild() and
// TopologyWorkspace::rebuildAdjacency(). Each failed direct scratch build
// must propagate std::bad_alloc, preserve the exact observable adjacency
// token stream, retain valid/readable state, and accept a disarmed retry.
// Run without sanitizers: this probe replaces global operator new.
#include <cstdint>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <new>
#include <vector>

#include "paacs/topo/Workspace.hpp"

static long gAdjacencyFailCountdown = 0;
static long gCurrentAdjacencyInjection = 0;
static const char* gCurrentAdjacencyPhase = "startup";
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t n) {
    if (gAdjacencyFailCountdown > 0 &&
        --gAdjacencyFailCountdown == 0)
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

static int gBad = 0;
#define REQUIRE(c)                                                        \
    do {                                                                  \
        if (!(c)) {                                                       \
            std::printf("FAIL %d: %s\n", __LINE__, #c);                   \
            ++gBad;                                                       \
        }                                                                 \
    } while (0)

void appendHandle(std::vector<std::uint64_t>& out, const auto& handle) {
    out.push_back(static_cast<std::uint64_t>(handle.slot));
    out.push_back(static_cast<std::uint64_t>(handle.gen));
    out.push_back(handle.pool);
}

std::vector<std::uint64_t> adjacencyTokens(const CellAdjacency& adj,
                                           const CellPool& cells) {
    std::vector<std::uint64_t> out;
    out.push_back(adj.boundNodePool());
    out.push_back(adj.boundCellPool());
    out.push_back(adj.boundEdgePool());
    out.push_back(adj.boundFacePool());
    out.push_back(static_cast<std::uint64_t>(adj.cellCount()));
    for (const auto& cell : cells.activeBySequence()) {
        appendHandle(out, cell);
        out.push_back(adj.hasCell(cells, cell) ? 1U : 0U);
        const auto faces = adj.cellFaces(cells, cell);
        out.push_back(static_cast<std::uint64_t>(faces.size()));
        for (const auto& face : faces) appendHandle(out, face);
        const auto neighbors = adj.neighbors(cells, cell);
        out.push_back(static_cast<std::uint64_t>(neighbors.size()));
        for (const auto& neighbor : neighbors) {
            appendHandle(out, neighbor.cell);
            appendHandle(out, neighbor.via);
        }
    }
    out.push_back(static_cast<std::uint64_t>(adj.boundaryFaces().size()));
    for (const auto& boundary : adj.boundaryFaces()) {
        appendHandle(out, boundary.face);
        appendHandle(out, boundary.owner);
    }
    out.push_back(static_cast<std::uint64_t>(adj.interiorFaces().size()));
    for (const auto& interior : adj.interiorFaces()) {
        appendHandle(out, interior.face);
        appendHandle(out, interior.a);
        appendHandle(out, interior.b);
    }
    out.push_back(
        static_cast<std::uint64_t>(adj.nonManifoldFaces().size()));
    for (const auto& nonManifold : adj.nonManifoldFaces()) {
        appendHandle(out, nonManifold.face);
        out.push_back(
            static_cast<std::uint64_t>(nonManifold.cells.size()));
        for (const auto& cell : nonManifold.cells) appendHandle(out, cell);
    }
    return out;
}

void makeTwoTetPools(NodePool& nodes, std::vector<NodeHandle>& handles,
                     CellPool& cells, EdgePool& edges, FacePool& faces,
                     CellHandle& first, CellHandle& second) {
    for (int i = 0; i < 5; ++i)
        handles.push_back(
            nodes.add({double(i), double(i % 2), double(i % 3)}, {3, 1}));
    first = cells.add(
        ElementType::Tet4,
        {handles[0], handles[1], handles[2], handles[3]}, {3, 1}, nodes);
    second = cells.add(
        ElementType::Tet4,
        {handles[0], handles[2], handles[1], handles[4]}, {3, 2}, nodes);
    edges.insertCellEdges(nodes, cells, first);
    edges.insertCellEdges(nodes, cells, second);
    faces.insertCellFaces(nodes, cells, edges, first);
    faces.insertCellFaces(nodes, cells, edges, second);
}

long directRebuildSweep() {
    long faults = 0;
    for (long k = 1;; ++k) {
        gCurrentAdjacencyPhase = "populated direct rebuild";
        gCurrentAdjacencyInjection = k;
        NodePool nodes;
        std::vector<NodeHandle> handles;
        CellPool cells;
        EdgePool edges;
        FacePool faces;
        CellHandle first, second;
        makeTwoTetPools(nodes, handles, cells, edges, faces, first, second);
        auto adjacency =
            CellAdjacency::build(nodes, cells, edges, faces);
        const auto before = adjacencyTokens(adjacency, cells);

        bool badAlloc = false;
        bool escapedOther = false;
        gAdjacencyFailCountdown = k;
        try {
            adjacency.rebuild(nodes, cells, edges, faces);
        } catch (const std::bad_alloc&) {
            badAlloc = true;
        } catch (...) {
            escapedOther = true;
        }
        gAdjacencyFailCountdown = 0;
        REQUIRE(!escapedOther);
        if (!badAlloc) {
            REQUIRE(adjacencyTokens(adjacency, cells) == before);
            REQUIRE(adjacency.checkValid(nodes, cells, edges, faces));
            break;
        }
        ++faults;
        REQUIRE(adjacencyTokens(adjacency, cells) == before);
        REQUIRE(adjacency.checkValid(nodes, cells, edges, faces));
        REQUIRE(adjacency.neighbors(cells, first).size() == 1);
        REQUIRE(adjacency.neighbors(cells, second).size() == 1);

        adjacency.rebuild(nodes, cells, edges, faces);
        REQUIRE(adjacencyTokens(adjacency, cells) == before);
        REQUIRE(adjacency.checkValid(nodes, cells, edges, faces));
    }
    std::printf("adjacency direct rebuild: %ld faults / 0 escapes\n",
                faults);
    return faults;
}

long emptyDirectRebuildSweep() {
    long faults = 0;
    for (long k = 1;; ++k) {
        gCurrentAdjacencyPhase = "empty direct rebuild";
        gCurrentAdjacencyInjection = k;
        NodePool nodes;
        std::vector<NodeHandle> handles;
        for (int i = 0; i < 3; ++i)
            handles.push_back(nodes.add(
                {double(i), double(i % 2), 0.0}, {2, 1}));
        CellPool cells;
        const auto triangle = cells.add(
            ElementType::Tri3,
            {handles[0], handles[1], handles[2]}, {2, 1}, nodes);
        EdgePool edges;
        edges.insertCellEdges(nodes, cells, triangle);
        FacePool faces;
        auto adjacency =
            CellAdjacency::build(nodes, cells, edges, faces);
        const auto before = adjacencyTokens(adjacency, cells);
        REQUIRE(adjacency.cellCount() == 0);
        REQUIRE(adjacency.checkValid(nodes, cells, edges, faces));

        bool badAlloc = false;
        bool escapedOther = false;
        gAdjacencyFailCountdown = k;
        try {
            adjacency.rebuild(nodes, cells, edges, faces);
        } catch (const std::bad_alloc&) {
            badAlloc = true;
        } catch (...) {
            escapedOther = true;
        }
        gAdjacencyFailCountdown = 0;
        REQUIRE(!escapedOther);
        if (!badAlloc) {
            REQUIRE(adjacencyTokens(adjacency, cells) == before);
            REQUIRE(adjacency.checkValid(nodes, cells, edges, faces));
            break;
        }
        ++faults;
        REQUIRE(adjacencyTokens(adjacency, cells) == before);
        REQUIRE(adjacency.checkValid(nodes, cells, edges, faces));

        adjacency.rebuild(nodes, cells, edges, faces);
        REQUIRE(adjacencyTokens(adjacency, cells) == before);
        REQUIRE(adjacency.checkValid(nodes, cells, edges, faces));
    }
    std::printf("adjacency empty direct rebuild: %ld faults / 0 escapes\n",
                faults);
    return faults;
}

long staleWorkspaceSweep() {
    long faults = 0;
    for (long k = 1;; ++k) {
        gCurrentAdjacencyPhase = "stale workspace rebuild";
        gCurrentAdjacencyInjection = k;
        NodePool nodes;
        std::vector<NodeHandle> handles;
        for (int i = 0; i < 5; ++i)
            handles.push_back(nodes.add(
                {double(i), double(i % 2), double(i % 3)}, {3, 1}));
        CellPool cells;
        const auto first = cells.add(
            ElementType::Tet4,
            {handles[0], handles[1], handles[2], handles[3]},
            {3, 1}, nodes);
        EdgePool edges;
        edges.insertCellEdges(nodes, cells, first);
        FacePool faces;
        faces.insertCellFaces(nodes, cells, edges, first);
        TopologyWorkspace workspace(std::move(nodes), std::move(cells),
                                    std::move(edges), std::move(faces));

        const auto oldAttribution =
            workspace.adjacencyBuiltFromRevision();
        auto transaction = workspace.beginTransaction();
        const std::array<NodeHandle, 4> addedNodes{
            handles[0], handles[2], handles[1], handles[4]};
        (void)transaction.addCell(
            ElementType::Tet4, addedNodes, {3, 2});
        const auto committed = transaction.commit();
        REQUIRE(committed.topologyChanged);
        REQUIRE(workspace.state() ==
                WorkspaceState::CommittedAdjacencyStale);
        const auto revision = workspace.topologyRevision();
        const auto semantic = workspace.semanticTopologyHash();
        const auto strict = workspace.strictWorkspaceHash();
        const auto before = adjacencyTokens(
            workspace.debugStaleAdjacency(), workspace.cells());

        bool badAlloc = false;
        bool escapedOther = false;
        gAdjacencyFailCountdown = k;
        try {
            workspace.rebuildAdjacency();
        } catch (const std::bad_alloc&) {
            badAlloc = true;
        } catch (...) {
            escapedOther = true;
        }
        gAdjacencyFailCountdown = 0;
        REQUIRE(!escapedOther);
        if (!badAlloc) {
            REQUIRE(workspace.state() == WorkspaceState::Ready);
            REQUIRE(workspace.adjacencyBuiltFromRevision() == revision);
            REQUIRE(workspace.adjacencyReady());
            break;
        }
        ++faults;
        REQUIRE(workspace.state() ==
                WorkspaceState::CommittedAdjacencyStale);
        REQUIRE(workspace.topologyRevision() == revision);
        REQUIRE(workspace.adjacencyBuiltFromRevision() == oldAttribution);
        REQUIRE(adjacencyTokens(workspace.debugStaleAdjacency(),
                                workspace.cells()) == before);
        REQUIRE(workspace.semanticTopologyHash() == semantic);
        REQUIRE(workspace.strictWorkspaceHash() == strict);

        workspace.rebuildAdjacency();
        REQUIRE(workspace.state() == WorkspaceState::Ready);
        REQUIRE(workspace.topologyRevision() == revision);
        REQUIRE(workspace.adjacencyBuiltFromRevision() == revision);
        REQUIRE(workspace.adjacencyReady());
        REQUIRE(workspace.adjacency().checkValid(
            workspace.nodes(), workspace.cells(), workspace.edges(),
            workspace.faces()));
    }
    std::printf("adjacency stale workspace rebuild: %ld faults / 0 escapes\n",
                faults);
    return faults;
}

}  // namespace

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::set_terminate([] {
        std::printf(
            "TERMINATE in adjacency %s injection %ld "
            "(active exception %s)\n",
            gCurrentAdjacencyPhase, gCurrentAdjacencyInjection,
            std::current_exception() ? "yes" : "no");
        std::_Exit(97);
    });
    const long directFaults = directRebuildSweep();
    const long emptyFaults = emptyDirectRebuildSweep();
    const long staleFaults = staleWorkspaceSweep();
    REQUIRE(directFaults >= 4);
    REQUIRE(emptyFaults >= 1);
    REQUIRE(staleFaults >= 4);
    if (gBad) {
        std::printf("ADJACENCY PROBE FAILED (%d)\n", gBad);
        return 1;
    }
    std::printf("ADJACENCY PROBE PASSED\n");
    return 0;
}
