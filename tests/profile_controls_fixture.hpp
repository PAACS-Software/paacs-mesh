// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/algo/PrimitiveSizeControls.hpp"
#include <string>
namespace profile_test {
using namespace paacs;
inline algo::PrimitiveShape shape(double scale = 1) {
    algo::PrimitiveShape s;
    s.radius = .6 * scale;
    s.height = 4.8 * scale;
    return s;
}
inline sizing::CappedSizeSpecification specification(const std::string &mode, double scale = 1,
                                                     double local = .12) {
    sizing::CappedSizeSpecification s;
    s.length = 4.8 * scale;
    s.globalSize = .36 * scale;
    if (mode == "start" || mode == "both")
        s.controls.push_back({0, 0, local * scale, 1.25});
    if (mode == "end" || mode == "both")
        s.controls.push_back({s.length, s.length, local * scale, 1.25});
    if (mode == "band")
        s.controls.push_back({2 * scale, 2.8 * scale, local * scale, 1.25});
    if (mode == "overlap")
        s.controls = {{1 * scale, 2.8 * scale, local * scale, 1.25},
                      {2.4 * scale, 3.2 * scale, .18 * scale, 1.15}};
    if (mode == "fine")
        s.controls = {{0, 0, local * scale, 1}};
    return s;
}
inline algo::PrimitiveMeshingOptions
options(algo::PrimitiveElements family, algo::PrimitiveProfileStrategy strategy, double scale = 1) {
    algo::PrimitiveMeshingOptions o;
    o.elements = family;
    o.profileStrategy = strategy;
    o.maximumChordError = .015 * scale;
    o.maximumRelativeVolumeError = .03;
    return o;
}
} // namespace profile_test
