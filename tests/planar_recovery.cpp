// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "geometry_checks.hpp"
#include "paacs/algo/Delaunay2DScalable.hpp"
#include "paacs/algo/PlanarConstraintRecovery.hpp"
#include <iostream>
#include <random>
#include <set>
using namespace geometry_test;
namespace {
auto edges(const std::vector<algo::Tri> &triangles) {
    std::set<std::pair<Index, Index>> result;
    for (const auto &t : triangles)
        for (int j = 0; j < 3; ++j)
            result.insert(std::minmax(t.v[j], t.v[(j + 1) % 3]));
    return result;
}
double area(const std::vector<Vec2> &points, const std::vector<algo::Tri> &triangles) {
    long double total = 0;
    for (const auto &t : triangles) {
        const auto a = points[t.v[0]], b = points[t.v[1]], c = points[t.v[2]];
        const double v = ((b.u - a.u) * (c.v - a.v) - (b.v - a.v) * (c.u - a.u)) / 2;
        require(v > 0, "flipped triangle");
        total += v;
    }
    return static_cast<double>(total);
}
bool same(const std::vector<algo::Tri> &a, const std::vector<algo::Tri> &b) {
    if (a.size() != b.size())
        return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i].v != b[i].v)
            return false;
    return true;
}
} // namespace
int main() {
    try {
        std::mt19937_64 random(0x5041414353ULL);
        std::uint64_t flips = 0;
        for (int run = 0; run < 40; ++run) {
            std::vector<Vec2> p{{0, 0}, {1, 0}, {1, 1}, {0, 1}};
            for (int i = 0; i < 36; ++i)
                p.push_back({std::generate_canonical<double, 53>(random),
                             std::generate_canonical<double, 53>(random)});
            auto t = algo::delaunay2DScalable(p);
            const auto original = t;
            const auto before = area(p, t);
            std::vector<algo::PlanarSegment> segments;
            for (Index i = 1; i < static_cast<Index>(p.size()); ++i)
                segments.push_back({0, i});
            const auto r = algo::recoverPlanarConstraints(p, segments, t);
            require(r.completed, "random radial constraint recovery " + std::to_string(run) +
                                     " reason=" + std::to_string(static_cast<int>(r.reason)));
            const auto e = edges(t);
            for (const auto s : segments)
                require(e.contains(std::minmax(s.first, s.second)), "missing constraint");
            require(t.size() == original.size() && std::abs(area(p, t) - before) < 1e-12,
                    "changed area or cell count");
            auto repeated = original;
            require(algo::recoverPlanarConstraints(p, segments, repeated).completed &&
                        same(t, repeated),
                    "recovery nondeterminism");
            repeated = original;
            require(algo::recoverPlanarConstraints(p, segments, repeated, 1).reason ==
                            algo::PlanarRecoveryReason::WorkBudget &&
                        same(repeated, original),
                    "work-bound rollback");
            flips += r.flips;
        }
        std::vector<Vec2> p{{0, 0}, {1, 0}, {1, 1}, {0, 1}};
        auto t = algo::delaunay2DScalable(p);
        const auto before = t;
        std::vector<algo::PlanarSegment> crossed{{0, 2}, {1, 3}};
        require(algo::recoverPlanarConstraints(p, crossed, t).reason ==
                        algo::PlanarRecoveryReason::CrossingConstraints &&
                    same(t, before),
                "crossing refusal");
        std::vector<algo::PlanarSegment> bad{{0, 99}};
        require(algo::recoverPlanarConstraints(p, bad, t).reason ==
                        algo::PlanarRecoveryReason::InvalidInput &&
                    same(t, before),
                "invalid endpoint");
        p.push_back({.5, .5});
        t = algo::delaunay2DScalable(p);
        const auto mid = t;
        crossed.resize(1);
        require(algo::recoverPlanarConstraints(p, crossed, t).reason ==
                        algo::PlanarRecoveryReason::UnsplitSegment &&
                    same(t, mid),
                "unsplit collinear vertex");
        require(flips > 0, "recovery corpus exercised no flips");
        std::cout << "40 seeded radial PSLGs and 3 refusal cases passed; " << flips << " flips\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
