// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "audit_reference.hpp"
#include "geometry_checks.hpp"
#include "paacs/algo/ConvexToTetMesher.hpp"
#include "paacs/algo/LayeredPatchMesher.hpp"
#include "paacs/algo/MappedHexMesher.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
using namespace paacs;
using geometry_test::digest;
using geometry_test::require;
namespace {
int passed = 0, refused = 0;
mesh::Mesh mixed(double scale) {
    mesh::Mesh cap;
    for (Vec3 p :
         std::array<Vec3, 6>{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0, 0}, {2, 1, 0}}})
        cap.addNode(p * scale, {2, 1});
    cap.addCell(mesh::CellType::Quad4, {0, 1, 2, 3}, {2, 1});
    cap.addCell(mesh::CellType::Tri3, {1, 4, 5}, {2, 1});
    cap.addCell(mesh::CellType::Tri3, {1, 5, 2}, {2, 1});
    algo::LayeredPatchRequest r;
    r.schedule.wallLayers = 3;
    r.schedule.firstLayerThickness = scale * .2;
    r.schedule.growthRatio = 1;
    r.schedule.totalDepth = scale * .6;
    r.schedule.maximumCoreThickness = scale * .2;
    r.threads = 1;
    mesh::Mesh m;
    const auto made = algo::meshLayeredPatch(cap, r, m);
    require(made.published,
            "mixed source reason=" + std::to_string(static_cast<int>(made.reason)) +
                " schedule=" + std::to_string(static_cast<int>(made.schedule.reason)) +
                " sweep=" + std::to_string(static_cast<int>(made.sweep.reason)));
    return m;
}
mesh::Mesh renumber(const mesh::Mesh &m, unsigned seed) {
    std::vector<Index> order(static_cast<std::size_t>(m.nodeCount()));
    std::iota(order.begin(), order.end(), 0);
    std::mt19937 engine(seed);
    std::shuffle(order.begin(), order.end(), engine);
    std::vector<Index> ids(order.size());
    mesh::Mesh out;
    for (Index i = 0; i < m.nodeCount(); ++i) {
        ids[static_cast<std::size_t>(order[static_cast<std::size_t>(i)])] = i;
        out.addNode(m.node(order[static_cast<std::size_t>(i)]),
                    m.nodeClass(order[static_cast<std::size_t>(i)]));
    }
    for (Index i = 0; i < m.cellCount(); ++i) {
        std::vector<Index> nodes;
        for (int j = 0; j < mesh::nodesPerCell(m.cell(i).type); ++j)
            nodes.push_back(ids[static_cast<std::size_t>(m.cellNodes(i)[j])]);
        out.addCell(m.cell(i).type, nodes, m.cell(i).onEntity);
    }
    return out;
}
mesh::Mesh alteredBoundary(const mesh::Mesh &m, int kind) {
    mesh::Mesh out;
    for (Index i = 0; i < m.nodeCount(); ++i)
        out.addNode(m.node(i), m.nodeClass(i));
    bool done = false;
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto c = m.cell(i);
        std::vector<Index> n(m.cellNodes(i), m.cellNodes(i) + mesh::nodesPerCell(c.type));
        if (!done && c.onEntity.dim == 2) {
            done = true;
            if (kind == 0)
                std::reverse(n.begin(), n.end());
            if (kind == 1)
                continue;
            if (kind == 2)
                out.addCell(c.type, n, c.onEntity);
        }
        out.addCell(c.type, n, c.onEntity);
    }
    return out;
}
void check(const mesh::Mesh &source, double expected) {
    algo::ConvexToTetOptions o;
    o.warningQuality = .95;
    mesh::Mesh out;
    const auto r = algo::convertConvexCellsToTet4(source, o, out);
    require(r.published, std::string("conversion: ") + algo::convexToTetReasonName(r.reason));
    require(out.nodeCount() == source.nodeCount(), "conversion added nodes");
    for (Index i = 0; i < source.nodeCount(); ++i)
        require((out.node(i) - source.node(i)).norm2() == 0 &&
                    out.nodeClass(i) == source.nodeClass(i),
                "changed node/tag");
    require(std::abs(geometry_test::positiveVolume(out) - expected) < expected * 1e-9,
            "volume preservation");
    require(mesh::referenceLinearVolumeBoundary(out).valid, "independent topology");
    require(r.lowQualityCells > 0 && r.minimumQuality > 0, "warning metrics");
    auto alias = source;
    require(algo::convertConvexCellsToTet4(alias, o, alias).published &&
                digest(alias) == digest(out),
            "source/destination alias");
    auto repeated = out;
    require(algo::convertConvexCellsToTet4(source, o, repeated).published &&
                digest(repeated) == digest(out),
            "repeat determinism");
    o.maximumOutputBytes = r.outputBytes;
    o.maximumWorkingBytes = r.workingBytes;
    require(algo::convertConvexCellsToTet4(source, o, repeated).published, "exact budgets");
    ++passed;
}
void refusals() {
    const auto m = mixed(1);
    algo::ConvexToTetOptions defaults;
    mesh::Mesh good;
    const auto baseline = algo::convertConvexCellsToTet4(m, defaults, good);
    require(baseline.published, "baseline");
    auto reject = [&](const mesh::Mesh &source, const algo::ConvexToTetOptions &o,
                      algo::ConvexToTetReason expected) {
        mesh::Mesh held;
        held.addNode({9, 8, 7}, {0, 9});
        const auto before = digest(held);
        const auto r = algo::convertConvexCellsToTet4(source, o, held);
        require(!r.published && r.destinationUnchanged && before == digest(held) &&
                    r.reason == expected,
                std::string("refusal #") + std::to_string(refused) +
                    " expected=" + algo::convexToTetReasonName(expected) +
                    " got=" + algo::convexToTetReasonName(r.reason));
        ++refused;
    };
    using R = algo::ConvexToTetReason;
    std::atomic<bool> token{true};
    for (int i = 0; i < 13; ++i) {
        auto o = defaults;
        R reason = R::InvalidRequest;
        switch (i) {
        case 0:
            o.relativePlanarityTolerance = 0;
            break;
        case 1:
            o.relativePlanarityTolerance = NAN;
            break;
        case 2:
            o.warningQuality = 2;
            break;
        case 3:
            o.minimumQuality = -1;
            break;
        case 4:
            o.maximumNodes = 1;
            reason = R::OutputBudget;
            break;
        case 5:
            o.maximumCells = 1;
            reason = R::OutputBudget;
            break;
        case 6:
            o.maximumOutputBytes = baseline.outputBytes - 1;
            reason = R::OutputBudget;
            break;
        case 7:
            o.maximumWorkingBytes = baseline.workingBytes - 1;
            reason = R::WorkingBudget;
            break;
        case 8:
            o.auditLimits.maximumFaceUses = 1;
            reason = R::WorkingBudget;
            break;
        case 9:
            o.cancellation = &token;
            reason = R::Cancelled;
            break;
        case 10:
            o.auditLimits.cancellation = &token;
            reason = R::Cancelled;
            break;
        case 11:
            o.minimumQuality = .999;
            reason = R::QualityRequirement;
            break;
        case 12:
            o.auditLimits.maximumTemporaryBytes = 1;
            reason = R::WorkingBudget;
            break;
        }
        reject(m, o, reason);
    }
    for (int i = 0; i < 3; ++i)
        reject(alteredBoundary(m, i), defaults, R::BoundaryAuditFailed);
    mesh::Mesh empty;
    reject(empty, defaults, R::InvalidSource);
    auto unused = m;
    unused.addNode({4, 4, 4}, {3, 1});
    reject(unused, defaults, R::InvalidSource);
    algo::MappedHexBlockRequest block;
    block.executionThreads = 1;
    mesh::Mesh hex;
    require(algo::generateMappedHexBlock(block, hex).published(), "hex source");
    auto warped = hex;
    warped.setNode(7, warped.node(7) + Vec3{0, 0, .2});
    reject(warped, defaults, R::NonplanarFace);
    auto flat = hex;
    for (Index i = 0; i < flat.nodeCount(); ++i) {
        auto p = flat.node(i);
        p.z = 0;
        flat.setNode(i, p);
    }
    reject(flat, defaults, R::NonConvexCell);
    auto nonfinite = hex;
    nonfinite.setNode(0, {NAN, 0, 0});
    reject(nonfinite, defaults, R::InvalidSource);
}
} // namespace
int main() {
    try {
        for (double scale : {1e-6, 1., 1e6}) {
            auto source = mixed(scale);
            for (unsigned seed = 0; seed < 16; ++seed)
                check(renumber(source, seed), 1.2 * scale * scale * scale);
        }
        for (double scale : {1e-9, 1e-6, 1., 1e6, 1e9}) {
            algo::LayerScheduleRequest schedule;
            schedule.wallLayers = 3;
            schedule.firstLayerThickness = .2 * scale;
            schedule.growthRatio = 1;
            schedule.totalDepth = .6 * scale;
            schedule.maximumCoreThickness = scale;
            const auto r = algo::planLayerSchedule(schedule);
            require(r.planned && r.offsets.size() == 4 && r.offsets.back() == schedule.totalDepth,
                    "decimal wall/depth roundoff");
            schedule.totalDepth *= 1 - 1e-8;
            const auto refused = algo::planLayerSchedule(schedule);
            require(!refused.planned &&
                        refused.reason == algo::LayerScheduleReason::WallLayersExceedDepth,
                    "real wall overrun accepted");
        }
        refusals();
        std::cout << "convex to Tet4: PASS (" << passed << " scaled/permuted meshes, " << refused
                  << " refusal/rollback cases)\n";
    } catch (const std::exception &e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
