// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/topo/ElementTraits.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <mutex>
#include <set>
#include <utility>

namespace paacs::topo {
namespace detail {
namespace {

Vec3 V(const RefCoord& r) { return {r.x, r.y, r.z}; }

double mapDet(const Vec3 b[3], const Vec3 c[3]) {
    const double detB = b[0].dot(b[1].cross(b[2]));
    if (std::abs(detB) < 1e-12) return 0.0;
    const double detC = c[0].dot(c[1].cross(c[2]));
    return detC / detB;
}

bool preservesEdgeAndFaceSets(const ElementTraits& tr,
                              const std::vector<std::uint8_t>& p) {
    std::set<std::pair<int, int>> edges;
    for (int e = 0; e < tr.edgeCount; ++e) {
        const int a = tr.edges[e].a, b = tr.edges[e].b;
        edges.insert({std::min(a, b), std::max(a, b)});
    }
    for (int e = 0; e < tr.edgeCount; ++e) {
        const int a = p[tr.edges[e].a], b = p[tr.edges[e].b];
        if (!edges.count({std::min(a, b), std::max(a, b)})) return false;
    }
    std::set<std::set<int>> faces;
    for (int f = 0; f < tr.faceCount; ++f) {
        std::set<int> fs;
        for (int k = 0; k < tr.faces[f].n; ++k)
            fs.insert(tr.faces[f].v[size_t(k)]);
        faces.insert(fs);
    }
    for (int f = 0; f < tr.faceCount; ++f) {
        std::set<int> fs;
        for (int k = 0; k < tr.faces[f].n; ++k)
            fs.insert(p[tr.faces[f].v[size_t(k)]]);
        if (!faces.count(fs)) return false;
    }
    return true;
}

bool orientationPreserving(const ElementTraits& tr,
                           const std::vector<std::uint8_t>& p) {
    if (tr.dimension == 1) {
        for (std::size_t i = 0; i < p.size(); ++i)
            if (p[i] != i) return false;
        return true;
    }
    if (tr.dimension == 2) {
        const int n = tr.cornerCount;
        const int s = p[0];
        for (int i = 0; i < n; ++i)
            if (p[size_t(i)] != (s + i) % n) return false;
        return true;
    }
    // 3D: induced linear map on reference corners must have det > 0.
    // Basis: three linearly independent edges leaving corner 0. The
    // pyramid apex and other valence issues are covered because every
    // 3D type here has >= 3 independent corner-0 edges.
    Vec3 b[3], c[3];
    int found = 0;
    for (int e = 0; e < tr.edgeCount && found < 3; ++e) {
        const int ea = tr.edges[e].a, eb = tr.edges[e].b;
        if (ea != 0 && eb != 0) continue;
        const int other = ea == 0 ? eb : ea;
        b[found] = V(tr.refCoords[other]) - V(tr.refCoords[0]);
        c[found] = V(tr.refCoords[p[size_t(other)]]) -
                   V(tr.refCoords[p[0]]);
        ++found;
    }
    if (found < 3) return false;
    return mapDet(b, c) > 1e-9;
}

} // namespace

std::vector<std::vector<std::uint8_t>> generateRotations(
    const ElementTraits& tr) {
    std::vector<std::vector<std::uint8_t>> out;
    std::vector<std::uint8_t> p(tr.cornerCount);
    for (std::size_t i = 0; i < p.size(); ++i) p[i] = std::uint8_t(i);
    do {
        if (preservesEdgeAndFaceSets(tr, p) && orientationPreserving(tr, p))
            out.push_back(p);
    } while (std::next_permutation(p.begin(), p.end()));
    return out;
}

} // namespace detail

const std::vector<std::vector<std::uint8_t>>& rotations(ElementType t) {
    if (!isValidElementType(t))
        throw std::out_of_range("rotations: ElementType out of range");
    static std::vector<std::vector<std::uint8_t>>
        cache[static_cast<std::size_t>(ElementType::Count_)];
    static std::once_flag onceFlags[static_cast<std::size_t>(
        ElementType::Count_)];
    const auto i = static_cast<std::size_t>(t);
    std::call_once(onceFlags[i], [&] {
        cache[i] = detail::generateRotations(traits(t));
    });
    return cache[i];
}

} // namespace paacs::topo
