// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/geometry/PolyhedralModel.hpp"
#include <cmath>
#include <string>

namespace geometry_test {
using namespace paacs;
struct Case {
    std::string name;
    geo::PolyhedralInput input;
    double volume = 0, targetSize = .3;
    algo::GeometryStrategy strategy = algo::GeometryStrategy::Auto;
    algo::GeometryStrategy expected = algo::GeometryStrategy::LinearPrismSweep;
};
inline std::vector<Case> cases() {
    using S = algo::GeometryStrategy;
    std::vector<Case> out;
    const auto add = [&](std::string name, std::vector<std::vector<Vec2>> loops, double area,
                         Vec3 sweep = Vec3{0, 0, 1}, Vec3 origin = Vec3{}, Vec3 u = Vec3{1, 0, 0},
                         Vec3 v = Vec3{0, 1, 0}, double h = .3, S expected = S::LinearPrismSweep) {
        out.push_back({std::move(name), geo::extrudedPolygon(loops, sweep, origin, u, v),
                       area * std::abs(u.cross(v).dot(sweep)), h, S::Auto, expected});
    };
    const std::vector<std::vector<Vec2>> square{{{0, 0}, {1, 0}, {1, 1}, {0, 1}}};
    add("block_unit", square, 1, {0, 0, 1}, {}, {1, 0, 0}, {0, 1, 0}, .2, S::MappedHex8);
    add("block_rectangular", square, 1, {0, 0, 3}, {}, {2, 0, 0}, {0, 1, 0}, .25, S::MappedHex8);
    add("block_rotated", square, 1, {.4, 0, .6928203230275509}, {2, -1, 4},
        {.8660254037844386, 0, -.5}, {0, 1, 0}, .2, S::MappedHex8);
    add("block_skew", square, 1, {.3, -.1, 1.4}, {}, {1, .15, .2}, {.2, 1, .05}, .2, S::MappedHex8);
    add("block_thin", square, 1, {0, 0, .002}, {}, {2, 0, 0}, {0, 1, 0}, .2, S::MappedHex8);
    add("block_slender", square, 1, {0, 0, 20}, {}, {.1, 0, 0}, {0, .1, 0}, .5, S::MappedHex8);
    add("block_translated", square, 1, {0, 0, 1}, {1e6, -1e6, 1e6}, {1, 0, 0}, {0, 1, 0}, .2,
        S::MappedHex8);
    add("block_small", square, 1, {0, 0, 1e-3}, {}, {1e-3, 0, 0}, {0, 1e-3, 0}, 2e-4,
        S::MappedHex8);
    add("block_large", square, 1, {0, 0, 1e3}, {}, {1e3, 0, 0}, {0, 1e3, 0}, 200, S::MappedHex8);
    add("triangle", {{{0, 0}, {2, 0}, {.3, 1}}}, 1);
    add("triangle_oblique", {{{0, 0}, {2, 0}, {.3, 1}}}, 1, {.4, .2, 2});
    add("triangle_thin", {{{0, 0}, {2, 0}, {.1, .03}}}, .03, {0, 0, 1}, {}, {1, 0, 0}, {0, 1, 0},
        .2);
    for (int n : {5, 6, 9, 16}) {
        std::vector<Vec2> ring;
        const double pi = std::acos(-1.0);
        for (int i = 0; i < n; ++i)
            ring.push_back({std::cos(2 * pi * i / n), std::sin(2 * pi * i / n)});
        add("polygon_" + std::to_string(n), {ring}, .5 * n * std::sin(2 * pi / n));
    }
    const std::vector<std::vector<Vec2>> l{{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}};
    add("concave_l", l, 3);
    add("concave_l_oblique", l, 3, {.3, .2, 1.3});
    add("concave_l_rotated", l, 3, {1, 0, 0}, {3, 2, 1}, {0, 0, 1}, {0, -1, 0});
    add("concave_u", {{{0, 0}, {3, 0}, {3, 3}, {2, 3}, {2, 1}, {1, 1}, {1, 3}, {0, 3}}}, 7);
    add("concave_c", {{{0, 0}, {3, 0}, {3, 1}, {1, 1}, {1, 2}, {3, 2}, {3, 3}, {0, 3}}}, 7);
    add("thin_ligament", {{{0, 0}, {3, 0}, {3, 3}, {2, 3}, {2, .02}, {1, .02}, {1, 3}, {0, 3}}},
        6.02);
    auto holed = std::vector<std::vector<Vec2>>{{{0, 0}, {3, 0}, {3, 3}, {0, 3}},
                                                {{1, 1}, {1, 2}, {2, 2}, {2, 1}}};
    add("hole_square", holed, 8);
    add("hole_oblique", holed, 8, {.2, -.1, 2});
    holed = {{{0, 0}, {4, 0}, {4, 3}, {0, 3}},
             {{.5, .5}, {.5, 1.5}, {1.5, 1.5}, {1.5, .5}},
             {{2.5, 1.5}, {2.5, 2.5}, {3.5, 2.5}, {3.5, 1.5}}};
    add("holes_two", holed, 10);
    add("hole_triangle", {{{0, 0}, {3, 0}, {3, 3}, {0, 3}}, {{1, 1}, {1.5, 2}, {2, 1}}}, 8.5);
    // Winding is normalized by the convenience constructor.
    add("clockwise_profile", {{{0, 0}, {0, 2}, {1, 2}, {1, 1}, {2, 1}, {2, 0}}}, 3);
    geo::PolyhedralInput tet;
    tet.vertices = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    tet.faces = {{{{0, 2, 1}}}, {{{0, 1, 3}}}, {{{1, 2, 3}}}, {{{2, 0, 3}}}};
    out.push_back({"tetrahedron", tet, 1. / 6, .35, S::Auto, S::BoundaryConformingTet4});
    geo::PolyhedralInput pyramid;
    pyramid.vertices = {{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {.8, .9, 2}};
    pyramid.faces = {{{{0, 3, 2, 1}}}, {{{0, 1, 4}}}, {{{1, 2, 4}}}, {{{2, 3, 4}}}, {{{3, 0, 4}}}};
    out.push_back({"pyramid", pyramid, 8. / 3, .5, S::Auto, S::BoundaryConformingTet4});
    auto frustum = geo::extrudedPolygon(square, {0, 0, 1});
    for (std::size_t i = 4; i < 8; ++i) {
        frustum.vertices[i].x = .25 + .5 * frustum.vertices[i].x;
        frustum.vertices[i].y = .25 + .5 * frustum.vertices[i].y;
    }
    out.push_back({"frustum", frustum, 7. / 12, .3, S::Auto, S::BoundaryConformingTet4});
    auto tetBox = out.front();
    tetBox.name = "block_forced_tet";
    tetBox.strategy = S::BoundaryConformingTet4;
    tetBox.expected = S::BoundaryConformingTet4;
    tetBox.targetSize = .3;
    out.push_back(std::move(tetBox));
    auto sweepBox = out.front();
    sweepBox.name = "block_forced_sweep";
    sweepBox.strategy = S::LinearPrismSweep;
    sweepBox.expected = S::LinearPrismSweep;
    out.push_back(std::move(sweepBox));
    // Engineering profiles: open channels, web/flange junctions, slots and
    // multiple internal voids. Each exact polygonal volume is checked later.
    add("beam_i",
        {{{0, 0},
          {3, 0},
          {3, .3},
          {1.7, .3},
          {1.7, 2.7},
          {3, 2.7},
          {3, 3},
          {0, 3},
          {0, 2.7},
          {1.3, 2.7},
          {1.3, .3},
          {0, .3}}},
        2.76, {0, 0, 4});
    add("beam_t", {{{0, 0}, {3, 0}, {3, .4}, {1.7, .4}, {1.7, 3}, {1.3, 3}, {1.3, .4}, {0, .4}}},
        2.24, {.2, .1, 2});
    add("channel_thin", {{{0, 0}, {3, 0}, {3, 2}, {2.9, 2}, {2.9, .1}, {.1, .1}, {.1, 2}, {0, 2}}},
        .68, {0, 0, 3});
    add("stepped_bracket", {{{0, 0}, {3, 0}, {3, 1}, {2, 1}, {2, 2}, {1, 2}, {1, 3}, {0, 3}}}, 6,
        {.25, .1, .6});
    add("slotted_plate",
        {{{0, 0}, {6, 0}, {6, 2}, {0, 2}}, {{.5, .5}, {.5, 1.5}, {5.5, 1.5}, {5.5, .5}}}, 7,
        {0, 0, .05});
    add("duct_thin_wall",
        {{{0, 0}, {3, 0}, {3, 2}, {0, 2}}, {{.1, .1}, {.1, 1.9}, {2.9, 1.9}, {2.9, .1}}}, .96,
        {.3, 0, 5});
    std::vector<std::vector<Vec2>> grid{{{0, 0}, {4, 0}, {4, 4}, {0, 4}}};
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i) {
            const double x = .5 + i, y = .5 + j;
            grid.push_back({{x, y}, {x, y + .5}, {x + .5, y + .5}, {x + .5, y}});
        }
    add("perforated_plate_9", grid, 13.75, {0, 0, .15});
    std::vector<std::vector<Vec2>> tube(2);
    for (int i = 0; i < 12; ++i) {
        const double a = 2 * std::acos(-1.) * i / 12;
        tube[0].push_back({2 * std::cos(a), 2 * std::sin(a)});
        tube[1].push_back({std::cos(-a), std::sin(-a)});
    }
    add("polygonal_tube", tube, 9, {.1, .2, 3});
    add("near_hole_wall", {{{0, 0}, {3, 0}, {3, 3}, {0, 3}}, {{.02, 1}, {.02, 2}, {1, 2}, {1, 1}}},
        8.02);
    add("acute_profile", {{{0, 0}, {3, 0}, {.02, .03}}}, .045, {0, 0, 1}, {}, {1, 0, 0}, {0, 1, 0},
        .2);
    return out;
}
} // namespace geometry_test
