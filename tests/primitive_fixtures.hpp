// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/algo/PrimitiveMesher.hpp"
#include <string>
#include <vector>
namespace primitive_test {
using namespace paacs;
struct Case {
    std::string name;
    algo::PrimitiveShape shape;
    double size = .3;
};
inline std::vector<Case> cases() {
    using K = algo::PrimitiveKind;
    std::vector<Case> out;
    for (auto kind : {K::Cylinder, K::Cone, K::Frustum, K::Tube, K::Torus, K::Elbow}) {
        algo::PrimitiveShape s;
        s.kind = kind;
        if (kind == K::Tube)
            s.innerRadius = .8;
        out.push_back({algo::primitiveKindName(kind), s});
    }
    auto add = [&](const char *name, algo::PrimitiveShape s, double h = .3) {
        out.push_back({name, s, h});
    };
    auto s = out[0].shape;
    s.height = .08;
    add("cylinder_thin", s, .15);
    s = out[0].shape;
    s.radius = .15;
    s.height = 5;
    add("cylinder_slender", s, .18);
    s = out[1].shape;
    s.height = .15;
    add("cone_shallow", s, .2);
    s = out[1].shape;
    s.radius = .3;
    s.height = 5;
    add("cone_sharp", s, .2);
    s = out[2].shape;
    s.endRadius = 2;
    add("frustum_expanding", s, .3);
    s = out[2].shape;
    s.endRadius = .015;
    add("frustum_small_cap", s, .2);
    s = out[3].shape;
    s.innerRadius = .99;
    add("tube_thin_wall", s, .2);
    s.innerRadius = .9999;
    add("tube_film", s, .3);
    s = out[3].shape;
    s.innerRadius = .15;
    add("tube_thick_wall", s, .2);
    s = out[4].shape;
    s.radius = .2;
    add("torus_slender", s, .25);
    s = out[4].shape;
    s.majorRadius = 1.05;
    add("torus_near_horn", s, .25);
    s = out[4].shape;
    s.innerRadius = .8;
    add("torus_hollow", s, .3);
    s.innerRadius = .999;
    add("torus_thin_wall", s, .3);
    s = out[5].shape;
    s.sweepAngle = 3 * std::numbers::pi / 2;
    add("elbow_270", s, .3);
    s = out[5].shape;
    s.innerRadius = .9;
    add("elbow_hollow", s, .3);
    s.innerRadius = .999;
    add("elbow_thin_wall", s, .3);
    s = out[5].shape;
    s.sweepAngle = .08;
    add("elbow_short", s, .2);
    for (auto k : {0U, 1U, 4U}) {
        s = out[k].shape;
        s.origin = {7, -3, 4};
        s.axis = {1, 2, 3};
        s.xDirection = {3, -1, .2};
        add((out[k].name + "_rotated").c_str(), s, .3);
    }
    for (double scale : {.001, 1000.}) {
        s = out[4].shape;
        s.radius *= scale;
        s.majorRadius *= scale;
        s.origin = {10 * scale, -20 * scale, 5 * scale};
        add(scale < 1 ? "torus_small" : "torus_large", s, .3 * scale);
    }
    s = out[1].shape;
    s.origin = {1e6, -2e6, 3e6};
    add("cone_translated", s, .25);
    return out;
}
} // namespace primitive_test
