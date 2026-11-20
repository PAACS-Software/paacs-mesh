// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/core/Parallel.hpp"
#include "paacs/core/SmallSort.hpp"
#include "paacs/core/PointGrid.hpp"
#include "paacs/sdk/CAPI.h"
#include "paacs/sdk/MeshArchive.hpp"
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <thread>

using namespace paacs;
namespace {
void require(bool ok) {
    if (!ok)
        throw std::runtime_error("hardening assertion failed");
}
template <class F> void invalid(F f) {
    bool caught = false;
    try {
        f();
    } catch (const std::invalid_argument &) {
        caught = true;
    }
    require(caught);
}
void put(std::vector<std::uint8_t> &b, std::size_t at, std::uint64_t v, int width = 8) {
    for (int i = 0; i < width; ++i)
        b[at + static_cast<std::size_t>(i)] = static_cast<std::uint8_t>(v >> (8 * i));
}
mesh::Mesh seed() {
    mesh::Mesh m;
    (void)m.addNode({17, 19, 23}, {0, 7});
    return m;
}
void unchanged(const mesh::Mesh &m) {
    require(m.nodeCount() == 1 && m.cellCount() == 0 && m.node(0).x == 17 &&
            m.nodeClass(0).tag == 7);
}
void archives() {
    mesh::Mesh m;
    for (int i = 0; i < 10; ++i)
        (void)m.addNode({double(i), -double(i), 0}, {0, i});
    for (int i = 0; i < 10; ++i) {
        const auto type = static_cast<mesh::CellType>(i);
        std::array<Index, 10> conn{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        (void)m.addCell(type, std::span<const Index>(conn.data(), mesh::nodesPerCell(type)),
                        {3, i});
    }
    std::vector<std::uint8_t> bytes;
    require(sdk::serializeMeshArchive(m, bytes).completed);
    mesh::Mesh out;
    require(sdk::deserializeMeshArchive(bytes, out).completed);
    std::vector<std::uint8_t> roundtrip;
    require(sdk::serializeMeshArchive(out, roundtrip).completed && roundtrip == bytes);
    // Independent wire-format fixture: empty mesh is exactly the 40-byte header.
    mesh::Mesh empty;
    std::vector<std::uint8_t> golden(40, 0), encoded;
    const std::array<std::uint8_t, 8> magic{'P', 'A', 'A', 'C', 'S', 'M', '8', 0};
    std::copy(magic.begin(), magic.end(), golden.begin());
    golden[8] = 1;
    require(sdk::serializeMeshArchive(empty, encoded).completed && encoded == golden);
    for (std::size_t n = 0; n < bytes.size(); ++n) {
        auto held = seed();
        require(!sdk::deserializeMeshArchive(std::span<const std::uint8_t>(bytes.data(), n), held)
                     .completed);
        unchanged(held);
    }
    const auto reject = [&](std::vector<std::uint8_t> bad) {
        auto held = seed();
        require(!sdk::deserializeMeshArchive(bad, held).completed);
        unchanged(held);
    };
    auto bad = bytes;
    put(bad, 16, UINT64_MAX);
    reject(bad);
    bad = bytes;
    put(bad, 16, 1ULL << 40);
    reject(bad);
    bad = bytes;
    bad[12] = 1;
    reject(bad);
    bad = bytes;
    bad.push_back(0);
    reject(bad);
    bad = bytes;
    put(bad, 40, std::bit_cast<std::uint64_t>(std::numeric_limits<double>::quiet_NaN()));
    reject(bad);
    constexpr std::size_t firstCell = 40 + 10 * 32;
    bad = bytes;
    bad[firstCell] = 255;
    reject(bad);
    bad = bytes;
    bad[firstCell + 1] = 1;
    reject(bad);
    bad = bytes;
    put(bad, firstCell + 12, 10);
    reject(bad);
    auto held = seed();
    sdk::MeshArchiveLimits limits;
    limits.maxNodes = 9;
    require(
        sdk::deserializeMeshArchive(std::span<const std::uint8_t>(bytes), held, limits).reason ==
        sdk::MeshArchiveReason::ResourceLimit);
    unchanged(held);
    // Deterministic mutation corpus. Valid mutations must retain their exact
    // canonical byte representation; rejected mutations must be transactional.
    std::mt19937 rng(190929);
    for (int i = 0; i < 2000; ++i) {
        bad = bytes;
        bad[rng() % bad.size()] ^= static_cast<std::uint8_t>(1U << (rng() % 8));
        held = seed();
        if (sdk::deserializeMeshArchive(bad, held).completed) {
            require(sdk::serializeMeshArchive(held, roundtrip).completed && roundtrip == bad);
        } else
            unchanged(held);
    }
}
void smallSort() {
    for (std::size_t count = 0; count <= 4; ++count) {
        std::array<int,4> permutation{0,1,2,3};
        do {
            auto sorted = permutation;
            require(core::sortSmallPrefix(sorted, count));
            require(std::is_sorted(sorted.begin(), sorted.begin()+count));
            for (std::size_t k=count;k<4;++k) require(sorted[k]==permutation[k]);
            int sum=0,original=0;
            for(std::size_t k=0;k<count;++k){sum+=sorted[k];original+=permutation[k];}
            require(sum==original);
        } while(std::next_permutation(permutation.begin(),permutation.end()));
    }
    std::array<int,4> values{3,1,2,0};const auto original=values;
    require(!core::sortSmallPrefix(values,5) && values==original);
}
void meshAppend() {
    mesh::Mesh m;
    for (int i = 0; i < 4; ++i)
        (void)m.addNode({double(i), 0, 0}, {0, i});
    invalid([&] { (void)m.addCell(static_cast<mesh::CellType>(255), {}); });
    require(m.cellCount() == 0 && m.connectivityCount() == 0);
    (void)m.addCell(mesh::CellType::Tet4, {0, 1, 2, 3}, {3, 1});
    (void)m.addCell(mesh::CellType::Tet4, std::span<const Index>(m.cellNodes(0), 4), {3, 1});
    require(m.cellCount() == 2 && m.cellNodes(1)[3] == 3);
    invalid([&] { m.setDedupTolerance(0); });
    invalid([&] { m.setDedupTolerance(std::numeric_limits<double>::quiet_NaN()); });
    invalid([&] { (void)m.addNodeUnique({std::numeric_limits<double>::infinity(), 0, 0}); });
}
void parallel() {
    const auto max = std::numeric_limits<std::size_t>::max();
    int calls = 0;
    core::parallelChunksWithThreads(9, max, max, 4, [&](auto c, auto lo, auto hi) {
        require(c == 0 && lo == 9 && hi == max);
        ++calls;
    });
    require(calls == 1);
    calls = 0;
    core::parallelChunksWithThreads(max - 3, max, 2, 1, [&](auto, auto lo, auto hi) {
        require(hi > lo && hi - lo <= 2);
        ++calls;
    });
    require(calls == 2);
    core::setThreadCount(1);
    require(core::parallelReduce<int>(
                0, 10, 0, 0, [](int &a, std::size_t i) { a += static_cast<int>(i); },
                [](int &a, int b) { a += b; }) == 45);
    require(core::parallelReduce<int>(
                max - 3, max, max, 0, [](int &a, auto) { ++a; }, [](int &a, int b) { a += b; }) ==
            3);
    std::atomic<int> count{0};
    core::parallelForWithThreads(0, 4096, 4, [&](auto) {
        require(core::inPoolContext());
        core::parallelForWithThreads(0, 2, 4, [&](auto) { ++count; });
    });
    require(count == 8192 && !core::inPoolContext());
    bool caught = false;
    try {
        core::parallelForWithThreads(0, 4096, 4, [](auto) { throw std::runtime_error("body"); });
    } catch (const std::runtime_error &) {
        caught = true;
    }
    require(caught && !core::inPoolContext());
    count = 0;
    core::parallelForWithThreads(0, 4096, 4, [&](auto) { ++count; });
    require(count == 4096);
    std::jthread a([] {
        for (int i = 0; i < 10000; ++i)
            core::setThreadCount(i % 4);
    });
    std::jthread b([] {
        for (int i = 0; i < 10000; ++i)
            require(core::effectiveThreads() >= 1);
    });
}
void grid() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    invalid([&] { core::PointGrid g(nan); });
    invalid([] { core::PointGrid g(0); });
    core::PointGrid g(1e-300);
    g.insert({1e200, 0, 0});
    require(g.size() == 1 && g.anyWithin({1e200, 0, 0}, 1e-200));
    require(!g.anyWithin({-1e200, 0, 0}, 1e199));
    require(g.anyWithin({0, 0, 0}, 1e300));
    require(!g.anyWithin({1e200, 0, 0}, 0));
    invalid([&] { g.insert({nan, 0, 0}); });
    invalid([&] { (void)g.anyWithin({0, 0, 0}, -1); });
    core::PointGrid sub(std::numeric_limits<double>::denorm_min());
    sub.insert({0, 0, 0});
    require(sub.anyWithin({0, 0, 0}, std::numeric_limits<double>::denorm_min()));
    core::PointGrid normal(0.1);
    std::vector<Vec3> points;
    std::mt19937 rng(19);
    for (int i = 0; i < 100; ++i) {
        Vec3 p{double(rng() % 100) / 100, double(rng() % 100) / 100, double(rng() % 100) / 100};
        points.push_back(p);
        normal.insert(p);
    }
    for (int i = 0; i < 300; ++i) {
        Vec3 p{double(rng() % 100) / 100, double(rng() % 100) / 100, double(rng() % 100) / 100};
        const double r = double(rng() % 20) / 100;
        bool oracle = false;
        for (const auto &q : points)
            oracle |= (q - p).norm2() < r * r;
        require(normal.anyWithin(p, r) == oracle);
    }
}
void capi() {
    auto *s = paacs_v8_session_create();
    require(s != nullptr);
    paacs_v8_hex_request hex{};
    hex.axis_u.x = 1;
    hex.axis_v.y = 1;
    hex.axis_w.z = 1;
    hex.cells_u = hex.cells_v = hex.cells_w = 1;
    hex.volume_tag = 1;
    hex.boundary_tag_base = 101;
    hex.threads = 1;
    require(paacs_v8_generate_mapped_hex(s, &hex) == PAACS_V8_OK);
    const auto hash = paacs_v8_session_mesh_hash(s);
    hex.boundary_tag_base = INT32_MAX;
    require(paacs_v8_generate_mapped_hex(s, &hex) == PAACS_V8_INVALID_ARGUMENT);
    paacs_v010_box_request box{};
    box.upper = {1, 1, std::numeric_limits<double>::infinity()};
    box.options = {0.5, PAACS_V010_TET4_CONFORMING, 0, 1000, 1000, 1};
    require(paacs_v010_generate_box_tet(s, &box) == PAACS_V8_INVALID_ARGUMENT);
    paacs_v010_cylinder_request cylinder{};
    cylinder.options = box.options;
    cylinder.radius = 1;
    cylinder.height = 1;
    cylinder.base_center.x = std::numeric_limits<double>::quiet_NaN();
    require(paacs_v010_generate_cylinder_tet(s, &cylinder) == PAACS_V8_INVALID_ARGUMENT);
    require(paacs_v8_session_mesh_hash(s) == hash);
    paacs_v8_session_destroy(s);
}
} // namespace
int main() {
    try {
        smallSort();
        archives();
        meshAppend();
        parallel();
        grid();
        capi();
        std::cout << "hardening: PASS (archive truncations + 2000 mutations, mesh, parallel, "
                     "spatial hash, C ABI)\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
