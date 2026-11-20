// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/sizing/ControlField.hpp"
#include <stdexcept>
#include <string>

namespace spatial_test {
using namespace paacs;
inline std::vector<std::string> cases() {
    return {"point", "segment", "triangle", "box", "ball", "minimum", "maximum", "smooth"};
}
inline sizing::ControlFieldSpecification specification(const std::string &name, double scale = 1) {
    using namespace sizing;
    ControlFieldSpecification s;
    s.minimumSize = .025;
    s.maximumSize = .4;
    if (name == "point")
        s.nodes = {field::distancePoint({.23, .5, .37}), field::threshold(0, .025, .42, .065, .4)};
    else if (name == "segment")
        s.nodes = {field::distanceSegment({.15, .5, .1}, {.85, .5, .9}),
                   field::threshold(0, .015, .30, .07, .4)};
    else if (name == "triangle")
        s.nodes = {field::distanceTriangle({.15, .5, .15}, {.85, .5, .2}, {.5, .5, .85}),
                   field::threshold(0, 0, .24, .075, .4)};
    else if (name == "box")
        s.nodes = {field::distanceBox({.35, .35, .35}, {.65, .65, .65}),
                   field::threshold(0, 0, .3, .075, .4)};
    else if (name == "ball")
        s.nodes = {field::distanceBall({.5, .5, .5}, .15), field::threshold(0, 0, .25, .075, .4)};
    else if (name == "minimum" || name == "maximum") {
        s.nodes = {field::distanceBall({.35, .5, .4}, .12), field::threshold(0, 0, .25, .075, .4),
                   field::distanceBall({.65, .5, .6}, .12), field::threshold(2, 0, .25, .075, .4)};
        s.nodes.push_back(name == "minimum" ? field::minimum(1, 3) : field::maximum(1, 3));
    } else if (name == "smooth") {
        s.nodes = {field::distancePoint({.5, .5, .5}), field::threshold(0, .12, .5, .075, .4, true),
                   field::scale(1, .8), field::clamp(2, .06, .32)};
    } else if (name == "fine")
        s.nodes = {field::constant(.065)};
    else if (name == "coarse")
        s.nodes = {field::constant(.4)};
    else
        throw std::runtime_error("unknown spatial case");
    for (auto &n : s.nodes) {
        n.a = n.a * scale;
        n.b = n.b * scale;
        n.c = n.c * scale;
        switch (n.operation) {
        case FieldOperation::Constant:
        case FieldOperation::DistanceBall:
            n.parameters[0] *= scale;
            break;
        case FieldOperation::ThresholdLinear:
        case FieldOperation::ThresholdSmooth:
            for (auto &p : n.parameters)
                p *= scale;
            break;
        case FieldOperation::Clamp:
            n.parameters[0] *= scale;
            n.parameters[1] *= scale;
            break;
        default:
            break;
        }
    }
    s.minimumSize *= scale;
    s.maximumSize *= scale;
    s.output = static_cast<FieldId>(s.nodes.size() - 1);
    return s;
}
inline sizing::ControlField make(const std::string &name, double scale = 1) {
    auto r = sizing::compileControlField(specification(name, scale));
    if (!r.compiled)
        throw std::runtime_error("spatial fixture compile");
    return std::move(r.field);
}
} // namespace spatial_test
