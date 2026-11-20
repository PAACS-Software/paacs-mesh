// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/PlanarConstraintRecovery.hpp"
#include "paacs/predicates/Predicates.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <map>
#include <set>

namespace paacs::algo {
namespace {
using Edge = std::pair<Index, Index>;
Edge key(Index a, Index b) {
    return std::minmax(a, b);
}
struct Stop {
    PlanarRecoveryReason reason;
};
struct Slot {
    std::array<std::size_t, 2> cell{};
    unsigned count = 0;
};
} // namespace
PlanarRecoveryOutcome recoverPlanarConstraints(std::span<const Vec2> points,
                                               std::span<const PlanarSegment> segments,
                                               std::vector<Tri> &destination,
                                               std::uint64_t maximumWork) noexcept {
    PlanarRecoveryOutcome out;
    try {
        const auto tick = [&](std::uint64_t n = 1) {
            if (n > maximumWork - out.work)
                throw Stop{PlanarRecoveryReason::WorkBudget};
            out.work += n;
        };
        const auto valid = [&](Index a) {
            return a >= 0 && static_cast<std::uint64_t>(a) < points.size();
        };
        for (const auto p : points) {
            tick();
            if (!std::isfinite(p.u) || !std::isfinite(p.v))
                throw Stop{PlanarRecoveryReason::InvalidInput};
        }
        const auto orient = [&](Index a, Index b, Index c) {
            tick();
            return predicates::orient2d(points[a], points[b], points[c]);
        };
        const auto crosses = [&](Edge a, Edge b) {
            tick();
            if (a.first == b.first || a.first == b.second || a.second == b.first ||
                a.second == b.second)
                return false;
            const auto a0 = orient(a.first, a.second, b.first),
                       a1 = orient(a.first, a.second, b.second);
            if (a0 == predicates::Sign::Zero || a1 == predicates::Sign::Zero || a0 == a1)
                return false;
            const auto b0 = orient(b.first, b.second, a.first),
                       b1 = orient(b.first, b.second, a.second);
            return b0 != predicates::Sign::Zero && b1 != predicates::Sign::Zero && b0 != b1;
        };
        std::vector<Edge> requested;
        if (segments.size() > maximumWork - out.work)
            throw Stop{PlanarRecoveryReason::WorkBudget};
        requested.reserve(segments.size());
        for (const auto s : segments) {
            tick();
            if (!valid(s.first) || !valid(s.second) || s.first == s.second)
                throw Stop{PlanarRecoveryReason::InvalidInput};
            requested.push_back(key(s.first, s.second));
        }
        std::sort(requested.begin(), requested.end());
        requested.erase(std::unique(requested.begin(), requested.end()), requested.end());
        for (std::size_t i = 0; i < requested.size(); ++i)
            for (std::size_t j = 0; j < i; ++j)
                if (crosses(requested[i], requested[j]))
                    throw Stop{PlanarRecoveryReason::CrossingConstraints};
        if (destination.size() > (maximumWork - out.work) / 4)
            throw Stop{PlanarRecoveryReason::WorkBudget};
        auto triangles = destination;
        std::map<Edge, Slot> edges;
        const auto insert = [&](std::size_t i) {
            const auto &t = triangles[i].v;
            for (unsigned k = 0; k < 3; ++k) {
                tick();
                auto &slot = edges[key(t[k], t[(k + 1) % 3])];
                if (slot.count == 2)
                    throw Stop{PlanarRecoveryReason::InvalidTriangulation};
                slot.cell[slot.count++] = i;
            }
        };
        const auto erase = [&](std::size_t i) {
            const auto &t = triangles[i].v;
            for (unsigned k = 0; k < 3; ++k) {
                tick();
                const auto edge = key(t[k], t[(k + 1) % 3]);
                auto it = edges.find(edge);
                if (it == edges.end())
                    throw Stop{PlanarRecoveryReason::InvalidTriangulation};
                auto &s = it->second;
                unsigned j = 0;
                while (j < s.count && s.cell[j] != i)
                    ++j;
                if (j == s.count)
                    throw Stop{PlanarRecoveryReason::InvalidTriangulation};
                s.cell[j] = s.cell[--s.count];
                if (!s.count)
                    edges.erase(it);
            }
        };
        for (std::size_t i = 0; i < triangles.size(); ++i) {
            auto &t = triangles[i].v;
            for (auto a : t)
                if (!valid(a))
                    throw Stop{PlanarRecoveryReason::InvalidInput};
            const auto sign = orient(t[0], t[1], t[2]);
            if (sign == predicates::Sign::Zero)
                throw Stop{PlanarRecoveryReason::InvalidTriangulation};
            if (sign == predicates::Sign::Negative)
                std::swap(t[1], t[2]);
            insert(i);
        }
        std::set<Edge> fixed;
        for (const auto target : requested) {
            if (edges.contains(target)) {
                fixed.insert(target);
                ++out.recovered;
                continue;
            }
            // A valid constrained edge cannot pass through an unsplit vertex.
            for (Index i = 0; i < static_cast<Index>(points.size()); ++i)
                if (i != target.first && i != target.second) {
                    tick();
                    const auto p = points[i], a = points[target.first], b = points[target.second];
                    if (p.u >= std::min(a.u, b.u) && p.u <= std::max(a.u, b.u) &&
                        p.v >= std::min(a.v, b.v) && p.v <= std::max(a.v, b.v) &&
                        orient(target.first, target.second, i) == predicates::Sign::Zero)
                        throw Stop{PlanarRecoveryReason::UnsplitSegment};
                }
            std::deque<Edge> queue;
            for (const auto &[edge, slot] : edges) {
                tick();
                if (crosses(target, edge)) {
                    if (slot.count != 2)
                        throw Stop{PlanarRecoveryReason::InvalidTriangulation};
                    queue.push_back(edge);
                }
            }
            std::size_t stalled = 0;
            while (!queue.empty() && !edges.contains(target)) {
                tick();
                const auto edge = queue.front();
                queue.pop_front();
                const auto it = edges.find(edge);
                if (it == edges.end())
                    continue;
                if (fixed.contains(edge))
                    throw Stop{PlanarRecoveryReason::CrossingConstraints};
                const auto slot = it->second;
                if (slot.count != 2)
                    throw Stop{PlanarRecoveryReason::InvalidTriangulation};
                const auto opposite = [&](std::size_t cell) {
                    for (auto v : triangles[cell].v)
                        if (v != edge.first && v != edge.second)
                            return v;
                    throw Stop{PlanarRecoveryReason::InvalidTriangulation};
                };
                const auto c = opposite(slot.cell[0]), d = opposite(slot.cell[1]);
                const auto diagonal = key(c, d);
                if (c == d || edges.contains(diagonal) || !crosses(edge, diagonal)) {
                    queue.push_back(edge);
                    if (++stalled > queue.size())
                        throw Stop{PlanarRecoveryReason::NotRecoverable};
                    continue;
                }
                erase(slot.cell[0]);
                erase(slot.cell[1]);
                auto a = Tri{{c, d, edge.first}}, b = Tri{{d, c, edge.second}};
                if (orient(a.v[0], a.v[1], a.v[2]) == predicates::Sign::Negative)
                    std::swap(a.v[1], a.v[2]);
                if (orient(b.v[0], b.v[1], b.v[2]) == predicates::Sign::Negative)
                    std::swap(b.v[1], b.v[2]);
                triangles[slot.cell[0]] = a;
                triangles[slot.cell[1]] = b;
                insert(slot.cell[0]);
                insert(slot.cell[1]);
                ++out.flips;
                stalled = 0;
                if (crosses(target, diagonal))
                    queue.push_back(diagonal);
            }
            if (!edges.contains(target))
                throw Stop{PlanarRecoveryReason::NotRecoverable};
            fixed.insert(target);
            ++out.recovered;
        }
        destination.swap(triangles);
        out.completed = true;
        out.reason = PlanarRecoveryReason::None;
    } catch (const Stop &s) {
        out.reason = s.reason;
    } catch (const std::bad_alloc &) {
        out.reason = PlanarRecoveryReason::AllocationFailed;
    } catch (...) {
        out.reason = PlanarRecoveryReason::InternalError;
    }
    return out;
}
} // namespace paacs::algo
