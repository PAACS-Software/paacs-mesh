// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/geometry/PolyhedralModel.hpp"
#include <cmath>
#include <string>
namespace complex_test {
using namespace paacs;
struct Shape {
    std::string name;
    std::vector<std::vector<Vec2>> loops;
    Vec3 sweep{0, 0, 3};
    double target = .3;
};
inline double area(const Shape &shape) {
    long double sum = 0;
    for (std::size_t j = 0; j < shape.loops.size(); ++j) {
        long double a = 0;
        const auto &ring = shape.loops[j];
        for (std::size_t i = 0; i < ring.size(); ++i) {
            const auto p = ring[i], q = ring[(i + 1) % ring.size()];
            a += p.u * q.v - p.v * q.u;
        }
        sum += (j == 0 ? 1 : -1) * std::abs(a) / 2;
    }
    return static_cast<double>(sum);
}
inline std::vector<Vec2> polygon(double x, double y, double r, int n = 12) {
    std::vector<Vec2> ring;
    for (int i = 0; i < n; ++i) {
        const double a = 2 * std::acos(-1.) * i / n;
        ring.push_back({x + r * std::cos(a), y + r * std::sin(a)});
    }
    return ring;
}
inline std::vector<Shape> shapes() {
    std::vector<Shape> v;
    v.push_back({"i_beam",
                 {{{0, 0},
                   {3, 0},
                   {3, .35},
                   {1.7, .35},
                   {1.7, 2.65},
                   {3, 2.65},
                   {3, 3},
                   {0, 3},
                   {0, 2.65},
                   {1.3, 2.65},
                   {1.3, .35},
                   {0, .35}}},
                 {0, 0, 4},
                 .28});
    Shape bores{"multi_bore", {{{0, 0}, {4, 0}, {4, 3}, {0, 3}}}, {.25, -.1, 2.8}, .36};
    for (double y : {.8, 2.2})
        for (double x : {.7, 2., 3.3})
            bores.loops.push_back(polygon(x, y, .3, 8));
    v.push_back(std::move(bores));
    Shape star{"star_duct", {}, {0, 0, 3.2}, .3};
    std::vector<Vec2> ring;
    for (int i = 0; i < 20; ++i) {
        const double a = 2 * std::acos(-1.) * i / 20, r = i % 2 ? 1.25 : 1.8;
        ring.push_back({r * std::cos(a), r * std::sin(a)});
    }
    star.loops = {ring, polygon(0, 0, .65, 12)};
    v.push_back(std::move(star));
    Shape channels{"cooling_channels", {{{0, 0}, {4, 0}, {4, 3}, {0, 3}}}, {0, 0, 2.4}, .32};
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 4; ++i) {
            const double x = .35 + .95 * i, y = .35 + .95 * j;
            channels.loops.push_back({{x, y}, {x + .4, y}, {x + .4, y + .4}, {x, y + .4}});
        }
    v.push_back(std::move(channels));
    v.push_back({"notched_frame",
                 {{{0, 0}, {4, 0}, {4, 3}, {3, 3}, {3, 2.4}, {2.6, 2.4}, {2.6, 3}, {0, 3}},
                  {{.45, .45}, {2, .45}, {2, 2.4}, {.45, 2.4}},
                  {{2.45, .45}, {3.55, .45}, {3.55, 1.7}, {2.45, 1.7}}},
                 {.2, .15, 2.5},
                 .32});
    v.push_back({"thin_bracket",
                 {{{0, 0}, {3, 0}, {3, .12}, {.12, .12}, {.12, 2}, {0, 2}}},
                 {0, 0, 2.2},
                 .22});
    v.push_back({"mapped_block", {{{0, 0}, {2, 0}, {2, 1.5}, {0, 1.5}}}, {.2, .1, 3}, .3});
    return v;
}
inline geo::PolyhedralInput input(const Shape &s, double scale = 1, bool rotate = false) {
    auto p = geo::extrudedPolygon(s.loops, s.sweep);
    for (auto &v : p.vertices) {
        v = v * scale;
        if (rotate)
            v = {2 * scale + .8 * v.x + .6 * v.z, -3 * scale + v.y, scale - .6 * v.x + .8 * v.z};
    }
    return p;
}
} // namespace complex_test
