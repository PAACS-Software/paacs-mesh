// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Smoothing.hpp"
#include "paacs/core/Parallel.hpp"

#include <atomic>
#include "paacs/mesh/Quality.hpp"

#include <algorithm>
#include <unordered_set>
#include <vector>

namespace paacs::algo {

using mesh::CellType;
using mesh::Mesh;

namespace {

// Worst radius-ratio over a set of tets, with `node` hypothetically at `p`.
// Returns -1 if any tet would invert (signed volume <= 0).
double cavityMinQuality(const Mesh& m, const std::vector<Index>& tets,
                        Index node, const Vec3& p) {
    double worst = 1e300;
    for (Index ci : tets) {
        const Index* v = m.cellNodes(ci);
        Vec3 q[4];
        for (int k = 0; k < 4; ++k)
            q[k] = (v[k] == node) ? p : m.node(v[k]);
        if (mesh::tetVolume(q[0], q[1], q[2], q[3]) <= 0) return -1.0;
        worst = std::min(worst, mesh::tetQuality(q[0], q[1], q[2], q[3]));
    }
    return worst;
}

} // namespace

SmoothResult smoothInterior(Mesh& m, const SmoothOptions& opt) {
    SmoothResult res;

    // --- adjacency: node -> incident Tet4 cells; node -> neighbour nodes ---
    const auto nNodes = static_cast<size_t>(m.nodeCount());
    std::vector<std::vector<Index>> nodeTets(nNodes);
    std::vector<std::unordered_set<Index>> nodeNbrs(nNodes);

    for (Index ci = 0; ci < m.cellCount(); ++ci) {
        if (m.cell(ci).type != CellType::Tet4) continue;
        const Index* v = m.cellNodes(ci);
        for (int a = 0; a < 4; ++a) {
            nodeTets[static_cast<size_t>(v[a])].push_back(ci);
            for (int b = 0; b < 4; ++b)
                if (a != b) nodeNbrs[static_cast<size_t>(v[a])].insert(v[b]);
        }
    }

    // Movable nodes: interior (dim 3) with at least one incident tet.
    std::vector<Index> movable;
    for (Index i = 0; i < m.nodeCount(); ++i)
        if (m.nodeClass(i).dim == 3 && !nodeTets[static_cast<size_t>(i)].empty())
            movable.push_back(i);

    // Greedy graph coloring over the node adjacency: two nodes of the same
    // color are never mesh-neighbors, so their smoothing cavities are
    // disjoint. Within a color, moves are therefore INDEPENDENT: they can be
    // applied concurrently, and the result is identical for any thread count
    // and any processing order. This is the standard parallel-smoothing
    // decomposition (colors processed as sequential barriers).
    std::vector<int> color(static_cast<size_t>(m.nodeCount()), -1);
    int nColors = 0;
    {
        std::vector<char> used;
        for (Index i : movable) {
            used.assign(static_cast<size_t>(nColors) + 1, 0);
            for (Index j : nodeNbrs[static_cast<size_t>(i)])
                if (color[static_cast<size_t>(j)] >= 0)
                    used[static_cast<size_t>(color[static_cast<size_t>(j)])] = 1;
            int c = 0;
            while (used[static_cast<size_t>(c)]) ++c;
            color[static_cast<size_t>(i)] = c;
            nColors = std::max(nColors, c + 1);
        }
    }
    std::vector<std::vector<Index>> byColor(static_cast<size_t>(nColors));
    for (Index i : movable)
        byColor[static_cast<size_t>(color[static_cast<size_t>(i)])].push_back(i);

    {
        const auto before = mesh::analyzeTets(m);
        res.minQualityBefore = before.radiusRatio.n ? before.radiusRatio.min : 0;
        res.meanQualityBefore = before.radiusRatio.mean;
    }

    // Characteristic length for the convergence test: mean incident edge.
    auto localH = [&](Index i) {
        double sum = 0; int n = 0;
        for (Index j : nodeNbrs[static_cast<size_t>(i)]) {
            sum += (m.node(i) - m.node(j)).norm();
            ++n;
        }
        return n ? sum / n : 1.0;
    };

    for (int it = 0; it < opt.maxIterations; ++it) {
        std::atomic<Index> accepted{0}, rejected{0};
        double maxRelMove = 0;

        for (const auto& group : byColor) {
            std::vector<double> moves(group.size(), 0.0);
            core::parallelFor(0, group.size(), [&](size_t gi) {
                const Index i = group[gi];
                const auto& nbrs = nodeNbrs[static_cast<size_t>(i)];
                if (nbrs.empty()) return;

                Vec3 centroid{};
                for (Index j : nbrs) centroid = centroid + m.node(j);
                centroid = centroid / static_cast<double>(nbrs.size());

                const Vec3 cur = m.node(i);
                const Vec3 proposal = cur + (centroid - cur) * opt.relaxation;

                const auto& tets = nodeTets[static_cast<size_t>(i)];
                const double qOld = cavityMinQuality(m, tets, i, cur);
                const double qNew = cavityMinQuality(m, tets, i, proposal);

                // Guard: never invert, never degrade the cavity's worst tet.
                if (qNew < 0 || qNew + 1e-12 < qOld) {
                    rejected.fetch_add(1, std::memory_order_relaxed);
                    return;
                }
                m.setNode(i, proposal);   // same-color cavities are disjoint
                accepted.fetch_add(1, std::memory_order_relaxed);
                moves[gi] = (proposal - cur).norm() / localH(i);
            });
            for (double mv : moves) maxRelMove = std::max(maxRelMove, mv);
        }

        res.acceptedMoves += accepted.load();
        res.rejectedMoves += rejected.load();
        ++res.iterations;
        if (maxRelMove < opt.convergenceTol) break;
    }

    {
        const auto after = mesh::analyzeTets(m);
        res.minQualityAfter = after.radiusRatio.n ? after.radiusRatio.min : 0;
        res.meanQualityAfter = after.radiusRatio.mean;
    }
    return res;
}

} // namespace paacs::algo
