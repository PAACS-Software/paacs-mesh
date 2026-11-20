// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "../bnd/TetBistellarRecovery.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <span>
#include <vector>

namespace paacs::geo::detail {

// Find a point in the kernel of an already embedded, closed interior star.
// Keeping its outward faces and all positive orientations preserves the
// occupied region. The caller applies the move after all spatial-index use.
template<class Admit, class Charge>
std::optional<Vec3> relocateOcctLiveCavityApex(const bnd::detail::BistellarComplex& live,
    bnd::TetFacetVertexId vertex, std::span<const std::size_t> incident,
    std::size_t target, Admit&& admit, Charge&& charge) {
    using Complex = bnd::detail::BistellarComplex;
    using predicates::Sign;
    if (vertex >= live.vertices.size() || target >= live.cells.size() || !live.alive[target] ||
        incident.size() < 4 || incident.size() > 16 ||
        std::ranges::find(incident, target) == incident.end()) return {};

    struct Plane { Vec3 origin, inward; };
    std::vector<Plane> planes;
    planes.reserve(incident.size());
    const auto origin = live.vertices[vertex].position;
    const auto point = [&](bnd::TetFacetVertexId id) { return live.vertices[id].position; };
    const auto finite = [](Vec3 p) {
        return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
    };
    double scale = 0.;
    for (const auto id : incident) {
        charge(1);
        if (id >= live.cells.size() || !live.alive[id]) return {};
        const auto& cell = live.cells[id];
        const auto apex = std::ranges::find(cell.vertices, vertex);
        if (apex == cell.vertices.end() || cell.classification != live.cells[target].classification)
            return {};
        const auto opposite = static_cast<unsigned>(apex-cell.vertices.begin());
        for (unsigned side = 0; side < 4; ++side) {
            if (side == opposite) continue;
            charge(1);
            const auto key = Complex::key(bnd::tet_facet_detail::outwardFace(cell, side));
            const auto owners = live.faces.find(key);
            if (live.protectedFaces.contains(key) || owners == live.faces.end() ||
                owners->second.size() != 2 ||
                !std::ranges::all_of(owners->second, [&](std::size_t owner) {
                    return std::ranges::find(incident, owner) != incident.end();
                })) return {};
        }
        const auto face = bnd::tet_facet_detail::outwardFace(cell, opposite);
        const auto a = point(face[0]), b = point(face[1]), c = point(face[2]);
        const auto inward = (c-a).cross(b-a).normalized();
        if (!finite(inward) || inward.norm2() == 0.) return {};
        planes.push_back({a, inward});
        scale = std::max({scale, (a-origin).norm(), (b-origin).norm(), (c-origin).norm()});
    }
    if (!std::isfinite(scale) || !(scale > 0.)) return {};

    // Normals move off nearby faces, and their cross products allow motion
    // along narrow portions of the kernel. Rank cheaply before CAD queries.
    std::vector<Vec3> directions;
    directions.reserve(26 + planes.size()*planes.size());
    for (int x : {-1, 0, 1}) for (int y : {-1, 0, 1}) for (int z : {-1, 0, 1})
        if (x || y || z) directions.push_back(Vec3(x,y,z).normalized());
    for (const auto& plane : planes) directions.push_back(plane.inward);
    for (std::size_t i = 0; i < planes.size(); ++i)
        for (std::size_t j = i+1; j < planes.size(); ++j) {
            const auto tangent = planes[i].inward.cross(planes[j].inward);
            const auto length = tangent.norm();
            if (length > 1e-8) {
                directions.push_back(tangent/length);
                directions.push_back(tangent*(-1./length));
            }
        }
    const auto margin = [&](Vec3 candidate) {
        charge(1);
        double distance = std::numeric_limits<double>::infinity();
        for (const auto& plane : planes)
            distance = std::min(distance, plane.inward.dot(candidate-plane.origin));
        return distance;
    };
    struct Candidate { Vec3 point; double margin; };
    std::vector<Candidate> candidates;
    Vec3 current = origin;
    double currentMargin = margin(current);
    double radius = scale*.25;
    for (unsigned level = 0; level < 16; ++level, radius *= .5) {
        for (unsigned step = 0; step < 8; ++step) {
            Vec3 best = current;
            double bestMargin = currentMargin;
            for (const auto direction : directions) {
                const auto candidate = current + direction*radius;
                const auto score = margin(candidate);
                if (score > bestMargin && finite(candidate)) {
                    best = candidate;
                    bestMargin = score;
                }
            }
            if (!(bestMargin > currentMargin)) break;
            current = best;
            currentMargin = bestMargin;
            if (currentMargin > 0.) candidates.push_back({current, currentMargin});
        }
    }
    std::ranges::stable_sort(candidates, std::greater<double>{}, &Candidate::margin);

    const auto centroid = [&](std::size_t id, Vec3 apex) {
        Vec3 result;
        for (const auto v : live.cells[id].vertices)
            result = result + (v == vertex ? apex : point(v))*.25;
        return result;
    };
    unsigned tried = 0;
    for (const auto& candidate : candidates) {
        if (tried++ == 16) break;
        const auto apex = candidate.point;
        const auto positive = [&](std::size_t id) {
            charge(1);
            std::array<Vec3, 4> p;
            const auto& ids = live.cells[id].vertices;
            for (unsigned k = 0; k < 4; ++k) p[k] = ids[k] == vertex ? apex : point(ids[k]);
            return predicates::orient3d(p[0], p[1], p[2], p[3]) == Sign::Positive;
        };
        if (!std::ranges::all_of(incident, positive) || !admit(centroid(target, apex)) || !admit(apex))
            continue;
        if (std::ranges::all_of(incident, [&](std::size_t id) {
            return id == target || admit(centroid(id, apex));
        })) return apex;
    }
    return {};
}

} // namespace paacs::geo::detail
