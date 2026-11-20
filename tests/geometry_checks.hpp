// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/topo/ElementTraits.hpp"
#include <cmath>
#include <set>
#include <stdexcept>

namespace geometry_test {
using namespace paacs;
inline void require(bool value, const std::string &message) {
    if (!value)
        throw std::runtime_error(message);
}
inline std::uint64_t digest(const mesh::Mesh &m) {
    std::vector<std::uint8_t> data;
    require(sdk::serializeMeshArchive(m, data).completed, "archive failed");
    std::uint64_t h = 1469598103934665603ULL;
    for (auto b : data) {
        h ^= b;
        h *= 1099511628211ULL;
    }
    return h;
}
// Independent divergence-theorem oracle. Uses oriented face fans about each
// cell's first corner, rather than the generators' tetrahedral decompositions.
inline double positiveVolume(const mesh::Mesh &m) {
    long double total = 0;
    for (Index i = 0; i < m.cellCount(); ++i) {
        topo::ElementType type{};
        require(topo::elementTypeOf(m.cell(i).type, type), "unsupported cell");
        const auto &traits = topo::traits(type);
        if (traits.dimension != 3)
            continue;
        const auto *c = m.cellNodes(i);
        const auto o = m.node(c[0]);
        long double volume = 0;
        for (unsigned j = 0; j < traits.faceCount; ++j) {
            const auto &f = traits.faces[j];
            for (unsigned k = 1; k + 1 < f.n; ++k) {
                const auto a = m.node(c[f.v[0]]) - o, b = m.node(c[f.v[k]]) - o,
                           d = m.node(c[f.v[k + 1]]) - o;
                volume += static_cast<long double>(a.dot(b.cross(d))) / 6;
            }
        }
        require(volume > 0 && std::isfinite(volume), "non-positive oriented cell volume");
        total += volume;
    }
    return static_cast<double>(total);
}
inline void checkClassification(const geo::Model &model, const mesh::Mesh &m, double tolerance) {
    std::set<int> covered;
    for (Index i = 0; i < m.nodeCount(); ++i) {
        const auto p = m.node(i);
        const auto c = m.nodeClass(i);
        if (c.dim == 0) {
            const auto *v = model.vertex(c.tag);
            require(v, "missing vertex class");
            require((v->position() - p).norm() <= tolerance, "vertex classification distance");
        } else if (c.dim == 1) {
            const auto *e = model.edge(c.tag);
            require(e, "missing edge class");
            require((e->closestPoint(p) - p).norm() <= tolerance,
                    "edge classification distance: tag=" + std::to_string(c.tag) +
                        " node=" + std::to_string(i) +
                        " distance=" + std::to_string((e->closestPoint(p) - p).norm()) +
                        " p=" + std::to_string(p.x) + "," + std::to_string(p.y) + "," +
                        std::to_string(p.z));
        } else if (c.dim == 2) {
            const auto *f = model.face(c.tag);
            require(f, "missing face class");
            require((f->pointAt(f->paramAt(p)) - p).norm() <= tolerance,
                    "face classification distance");
        } else
            require(c.dim == 3 && model.solid(c.tag), "missing solid class");
    }
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto &c = m.cell(i);
        if (c.type == mesh::CellType::Tri3 || c.type == mesh::CellType::Quad4) {
            const auto *f = model.face(c.onEntity.tag);
            require(c.onEntity.dim == 2 && f, "missing boundary owner");
            covered.insert(c.onEntity.tag);
            Vec3 center{};
            const auto n = mesh::nodesPerCell(c.type);
            for (int k = 0; k < n; ++k) {
                const auto p = m.node(m.cellNodes(i)[k]);
                center = center + p * (1. / n);
                require((f->pointAt(f->paramAt(p)) - p).norm() <= tolerance, "wrong boundary face");
            }
            require(f->containsParam(f->paramAt(center)), "boundary outside trimmed face");
        }
    }
    require(covered.size() == model.faces().size(), "CAD face lost from mesh boundary");
}
} // namespace geometry_test
