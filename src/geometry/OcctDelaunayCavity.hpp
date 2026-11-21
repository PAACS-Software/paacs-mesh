// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "../bnd/TetBistellarRecovery.hpp"
#include "paacs/mesh/Quality.hpp"
#include <algorithm>
#include <cmath>
#include <optional>

namespace paacs::geo::detail {

// Start from the oversized tetrahedron whose circumcenter was selected,
// or a tetrahedron already known to contain an inward support point.
// The connected in-sphere cavity must contain the point. Extend it across
// nonvisible shell faces only within the same region and away from protected
// faces. An optional dimensionless shape floor also expands across visible
// faces whose cones would be nearly flat. The caller certifies and commits
// the proposed replacement.
inline std::optional<bnd::detail::BistellarComplex::Proposal> buildOcctDelaunayCavity(
    const bnd::detail::BistellarComplex& live, std::size_t start, Vec3 point,
    std::size_t maximumCells, std::size_t& containing, bool growInSphere = true,
    double minimumConeShape = 0.) {
    using Complex = bnd::detail::BistellarComplex;
    using predicates::Sign;
    using bnd::tet_facet_detail::outwardFace;
    if (!live.walkConstraints || !maximumCells || start>=live.alive.size() ||
        !live.alive[start] || live.vertices.size()>=bnd::kInvalidTetFacetVertex ||
        !std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(point.z)) return {};

    Complex::Proposal patch;
    auto& cavity = patch.removed;
    cavity.push_back(start);
    const auto region = live.cells[start].classification;
    const auto position = [&](bnd::TetFacetVertexId id) -> const Vec3& {
        return live.vertices[id].position;
    };
    const auto side = [&](const Complex::Face& face) {
        live.charge();
        return predicates::orient3dFiltered(position(face[0]),position(face[1]),
                                            position(face[2]),point);
    };
    const auto contains = [&](std::size_t cell) {
        return std::ranges::find(cavity,cell)!=cavity.end();
    };
    const auto canCross = [&](std::size_t next, const Complex::Face& face) {
        return next!=Complex::noCell && live.cells[next].classification==region &&
               !live.protectedFaces.contains(Complex::key(face));
    };

    auto owner = Complex::noCell;
    for (std::size_t i=0; i<cavity.size(); ++i) {
        const auto cell = cavity[i];
        bool pointInCell = owner==Complex::noCell;
        for (unsigned f=0; f<4; ++f) {
            live.charge();
            const auto face = outwardFace(live.cells[cell],f);
            if (pointInCell && side(face)==Sign::Positive) pointInCell = false;
            const auto next = live.neighbours[cell][f];
            if (!growInSphere || contains(next) || !canCross(next,face)) continue;
            const auto& v = live.cells[next].vertices;
            if (predicates::insphereFiltered(position(v[0]),position(v[1]),position(v[2]),
                                              position(v[3]),point)!=Sign::Positive) continue;
            if (cavity.size()>=maximumCells) return {};
            cavity.push_back(next);
        }
        if (pointInCell) owner = cell;
    }
    if (owner==Complex::noCell) return {};

    // Every addition removes an obstructing shell face. Faces that become
    // exposed belong to the new cells and are visited later in this loop.
    for (std::size_t i=0; i<cavity.size(); ++i) {
        const auto cell = cavity[i];
        for (unsigned f=0; f<4; ++f) {
            const auto next = live.neighbours[cell][f];
            const auto face = outwardFace(live.cells[cell],f);
            if (contains(next)) continue;
            if (side(face)==Sign::Negative) {
                if (!(minimumConeShape>0.)) continue;
                live.charge(4);
                if (mesh::tetScaledJacobian(point,position(face[0]),position(face[1]),
                                            position(face[2]))>=minimumConeShape) continue;
            }
            if (cavity.size()>=maximumCells || !canCross(next,face)) return {};
            cavity.push_back(next);
        }
    }

    const auto staged = static_cast<bnd::TetFacetVertexId>(live.vertices.size());
    for (const auto cell:cavity) for (unsigned f=0; f<4; ++f) {
        live.charge();
        const auto face = outwardFace(live.cells[cell],f);
        if (contains(live.neighbours[cell][f])) {
            if (live.protectedFaces.contains(Complex::key(face))) return {};
            continue;
        }
        patch.added.push_back({{staged,face[0],face[1],face[2]},0,region});
    }
    if (patch.added.empty()) return {};
    containing = owner;
    return patch;
}

} // namespace paacs::geo::detail
