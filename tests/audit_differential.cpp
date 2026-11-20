// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "audit_reference.hpp"
#include "geometry_checks.hpp"
#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/algo/PreparedSweepMesher.hpp"
#include "paacs/topo/ElementTraits.hpp"
#include <algorithm>
#include <bit>
#include <iostream>
#include <numeric>
#include <random>
#include <tuple>

using namespace geometry_test;
namespace {
int checks = 0;
auto diagnostics(const mesh::LinearVolumeAudit &a) {
    return std::tuple(a.complete, a.valid, a.budgetExceeded, a.volumeCells, a.boundaryCells,
                      a.faceUses, a.exteriorFaces, a.interiorFaces, a.nonmanifoldFaces,
                      a.inconsistentOrientation, a.missingBoundary, a.unexpectedBoundary,
                      a.duplicateBoundary, a.invalidConnectivity, a.nonfiniteNodes,
                      a.unsupportedCells);
}
std::vector<std::uint64_t> snapshot(const mesh::Mesh &m) {
    std::vector<std::uint64_t> out;
    for (Index i = 0; i < m.nodeCount(); ++i) {
        const auto p = m.node(i);
        out.insert(out.end(), {std::bit_cast<std::uint64_t>(p.x), std::bit_cast<std::uint64_t>(p.y),
                               std::bit_cast<std::uint64_t>(p.z),
                               static_cast<std::uint64_t>(m.nodeClass(i).dim),
                               static_cast<std::uint64_t>(m.nodeClass(i).tag)});
    }
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto &c = m.cell(i);
        out.insert(out.end(),
                   {static_cast<std::uint64_t>(c.type), static_cast<std::uint64_t>(c.first),
                    static_cast<std::uint64_t>(c.onEntity.dim),
                    static_cast<std::uint64_t>(c.onEntity.tag)});
        for (int k = 0; k < mesh::nodesPerCell(c.type); ++k)
            out.push_back(static_cast<std::uint64_t>(m.cellNodes(i)[k]));
    }
    return out;
}
void compare(const mesh::Mesh &m) {
    const auto before = snapshot(m);
    const auto fast = mesh::auditLinearVolumeBoundary(m);
    const auto reference = mesh::referenceLinearVolumeBoundary(m);
    require(diagnostics(fast) == diagnostics(reference), "audit/reference diagnostic mismatch");
    require(snapshot(m) == before, "audit mutated input");
    ++checks;
}
void volumeWithBoundary(mesh::Mesh &m, mesh::CellType type, const std::vector<Index> &ids) {
    (void)m.addCell(type, ids, {3, 1});
    topo::ElementType e{};
    require(topo::elementTypeOf(type, e), "element type");
    const auto &t = topo::traits(e);
    for (unsigned j = 0; j < t.faceCount; ++j) {
        const auto &f = t.faces[j];
        std::vector<Index> face;
        for (unsigned k = 0; k < f.n; ++k)
            face.push_back(ids[f.v[k]]);
        (void)m.addCell(f.n == 3 ? mesh::CellType::Tri3 : mesh::CellType::Quad4, face, {2, 1});
    }
}
mesh::Mesh mixed() {
    mesh::Mesh cap, m;
    for (Vec3 p : {Vec3{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}})
        (void)cap.addNode(p, {2, 1});
    (void)cap.addCell(mesh::CellType::Quad4, {0, 1, 4, 3}, {2, 1});
    (void)cap.addCell(mesh::CellType::Tri3, {1, 2, 5}, {2, 1});
    (void)cap.addCell(mesh::CellType::Tri3, {1, 5, 4}, {2, 1});
    algo::PreparedMixedSweepRequest request;
    request.layers = 3;
    request.executionThreads = 1;
    require(algo::generatePreparedMixedSweep(cap, request, m).published(), "mixed sweep");
    const auto base = m.nodeCount();
    for (Vec3 p : {Vec3{4, 0, 0},
                   {5, 0, 0},
                   {5, 1, 0},
                   {4, 1, 0},
                   {4.5, .5, 1},
                   {7, 0, 0},
                   {8, 0, 0},
                   {7, 1, 0},
                   {7, 0, 1}})
        (void)m.addNode(p, {3, 1});
    volumeWithBoundary(m, mesh::CellType::Pyramid5, {base, base + 1, base + 2, base + 3, base + 4});
    volumeWithBoundary(m, mesh::CellType::Tet4, {base + 5, base + 6, base + 7, base + 8});
    return m;
}
struct Row {
    mesh::CellType type;
    std::vector<Index> ids;
};
void randomized(const mesh::Mesh &base) {
    std::mt19937_64 random(0x504141435333ULL);
    for (int trial = 0; trial < 2000; ++trial) {
        std::vector<Index> permutation(static_cast<std::size_t>(base.nodeCount()));
        std::iota(permutation.begin(), permutation.end(), 0);
        std::shuffle(permutation.begin(), permutation.end(), random);
        std::vector<Row> rows;
        for (Index i = 0; i < base.cellCount(); ++i) {
            Row r{base.cell(i).type, {}};
            for (int k = 0; k < mesh::nodesPerCell(r.type); ++k)
                r.ids.push_back(permutation[static_cast<std::size_t>(base.cellNodes(i)[k])]);
            rows.push_back(std::move(r));
        }
        const auto at = static_cast<std::size_t>(random() % rows.size());
        switch (trial % 10) {
        case 0:
            break;
        case 1:
            rows.erase(rows.begin() + static_cast<std::ptrdiff_t>(at));
            break;
        case 2:
            rows.push_back(rows[at]);
            break;
        case 3:
            std::swap(rows[at].ids[0], rows[at].ids[1]);
            break;
        case 4:
            rows[at].ids[0] = rows[at].ids[1];
            break;
        case 5:
            for (auto &r : rows)
                if (r.type == mesh::CellType::Tri3 || r.type == mesh::CellType::Quad4) {
                    std::rotate(r.ids.begin(), r.ids.begin() + 1, r.ids.end());
                    if (random() & 1)
                        std::reverse(r.ids.begin(), r.ids.end());
                }
            break;
        case 6:
            rows[at].ids[0] = static_cast<Index>(random() % permutation.size());
            break;
        case 7:
            rows.push_back(rows[at]);
            rows.push_back(rows[at]);
            break;
        case 8:
            rows.push_back({mesh::CellType::Tri6, {0, 1, 2, 3, 4, 5}});
            break;
        default:
            break;
        }
        std::shuffle(rows.begin(), rows.end(), random);
        mesh::Mesh m;
        // Coordinates do not affect this incidence test, except finite checks.
        for (Index i = 0; i < base.nodeCount() + trial % 7; ++i)
            (void)m.addNode({double(i), 0, 0});
        for (const auto &r : rows)
            (void)m.addCell(r.type, r.ids);
        if (trial % 10 == 9)
            m.setNode(0, {NAN, 0, 0});
        compare(m);
        if (trial % 10 == 0 || trial % 10 == 5)
            require(mesh::auditLinearVolumeBoundary(m).valid, "permutation changed validity");
    }
}
void targeted(const mesh::Mesh &base) {
    compare(base);
    require(mesh::auditLinearVolumeBoundary(base).valid, "mixed fixture validity");
    auto exact = mesh::auditLinearVolumeBoundary(base);
    mesh::LinearVolumeAuditLimits limits;
    limits.maximumTemporaryBytes = exact.temporaryBytes;
    require(mesh::auditLinearVolumeBoundary(base, limits).valid, "exact temporary budget");
    --limits.maximumTemporaryBytes;
    auto result = mesh::auditLinearVolumeBoundary(base, limits);
    require(result.budgetExceeded && !result.complete && !result.valid, "temporary budget refusal");
    limits = {};
    limits.maximumFaceUses = exact.faceUses;
    require(mesh::auditLinearVolumeBoundary(base, limits).valid, "exact face budget");
    --limits.maximumFaceUses;
    require(mesh::auditLinearVolumeBoundary(base, limits).budgetExceeded, "face budget refusal");
    std::atomic<bool> stop{true};
    limits = {};
    limits.cancellation = &stop;
    result = mesh::auditLinearVolumeBoundary(base, limits);
    require(result.cancelled && !result.complete && !result.valid, "pre-cancelled audit");
    checks += 5;
    for (int variant = 0; variant < 5; ++variant) {
        auto m = base;
        if (variant == 0)
            const_cast<mesh::Cell &>(m.cell(0)).first = -1;
        if (variant == 1)
            const_cast<mesh::Cell &>(m.cell(0)).first = m.connectivityCount();
        if (variant == 2)
            const_cast<Index *>(m.cellNodes(0))[0] = -1;
        if (variant == 3)
            const_cast<Index *>(m.cellNodes(0))[0] = m.nodeCount();
        if (variant == 4)
            const_cast<mesh::Cell &>(m.cell(0)).type = static_cast<mesh::CellType>(255);
        // Archive hashing assumes structurally valid storage, so don't call it
        // on deliberately corrupted rows; the audit must preflight these.
        require(diagnostics(mesh::auditLinearVolumeBoundary(m)) ==
                    diagnostics(mesh::referenceLinearVolumeBoundary(m)),
                "invalid metadata mismatch");
        require(!mesh::auditLinearVolumeBoundary(m).valid, "invalid metadata accepted");
        ++checks;
    }
    mesh::Mesh fan;
    (void)fan.addNode({0, 0, 0});
    for (int i = 0; i < 2000; ++i) {
        const auto start = fan.nodeCount();
        (void)fan.addNode({1, double(i), 0});
        (void)fan.addNode({1, double(i) + .5, 0});
        (void)fan.addNode({1, double(i), 1});
        volumeWithBoundary(fan, mesh::CellType::Tet4, {0, start, start + 1, start + 2});
    }
    compare(fan);
    require(mesh::auditLinearVolumeBoundary(fan).maximumBucketSize == 12000,
            "high-valence bucket coverage");
    mesh::Mesh empty;
    compare(empty);
    (void)empty.addNode({0, 0, 0});
    compare(empty);
    algo::MappedHexBlockRequest request;
    request.geometry.cellsU = request.geometry.cellsV = request.geometry.cellsW = 1;
    mesh::Mesh hex;
    require(algo::generateMappedHexBlock(request, hex).published(), "hex fixture");
    for (Index i = 0; i < hex.cellCount(); ++i)
        if (hex.cell(i).type == mesh::CellType::Quad4) {
            auto *ids = const_cast<Index *>(hex.cellNodes(i));
            std::swap(ids[1], ids[2]);
            break;
        }
    const auto crossed = mesh::auditLinearVolumeBoundary(hex);
    require(crossed.missingBoundary == 1 && crossed.unexpectedBoundary == 1,
            "crossed quad must not match a cyclic quad");
    compare(hex);
}
} // namespace
int main() {
    try {
        const auto m = mixed();
        targeted(m);
        randomized(m);
        std::cout << checks << " exact audit differential and resource checks passed\n";
    } catch (const std::exception &e) {
        std::cerr << "check " << checks << ": " << e.what() << '\n';
        return 1;
    }
}
