// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Frozen rc.2 global-sort audit, retained only as a differential oracle.
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/topo/ElementTraits.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace paacs::mesh {
namespace {
struct FaceUse {
    std::array<Index, 4> key{-1, -1, -1, -1};
    std::uint8_t n = 0;
    bool boundary = false;
    std::int8_t orientation = 0;
};
bool sameKey(const FaceUse &a, const FaceUse &b) {
    return a.n == b.n && a.key == b.key;
}
FaceUse makeFace(std::array<Index, 4> nodes, std::uint8_t n, bool boundary) {
    std::uint8_t first = 0;
    for (std::uint8_t i = 1; i < n; ++i)
        if (nodes[i] < nodes[first])
            first = i;
    const auto direction =
        static_cast<std::int8_t>(nodes[(first + 1) % n] < nodes[(first + n - 1) % n] ? 1 : -1);
    // Preserve the cyclic adjacency of quad corners: equal vertex sets alone
    // do not identify the same quadrilateral face.
    std::array<Index, 4> canonical{-1, -1, -1, -1};
    for (std::uint8_t i = 0; i < n; ++i)
        canonical[i] = nodes[(first + n + direction * i) % n];
    return {canonical, n, boundary, direction};
}
const topo::ElementTraits &cellTraits(CellType t) {
    topo::ElementType e{};
    (void)topo::elementTypeOf(t, e);
    return topo::traits(e);
}
bool supported(CellType t) {
    return t == CellType::Tet4 || t == CellType::Hex8 || t == CellType::Prism6 ||
           t == CellType::Pyramid5 || t == CellType::Tri3 || t == CellType::Quad4 ||
           t == CellType::Line2;
}
} // namespace
LinearVolumeAudit referenceLinearVolumeBoundary(const Mesh &value,
                                            const LinearVolumeAuditLimits &limits) noexcept {
    LinearVolumeAudit out;
    try {
        for (Index i = 0; i < value.nodeCount(); ++i) {
            const auto p = value.node(i);
            if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
                ++out.nonfiniteNodes;
        }
        std::uint64_t total = 0;
        for (Index i = 0; i < value.cellCount(); ++i) {
            const auto &c = value.cell(i);
            if (!supported(c.type)) {
                ++out.unsupportedCells;
                continue;
            }
            if (c.type == CellType::Line2)
                continue;
            const auto &traits = cellTraits(c.type);
            const bool surface = traits.dimension == 2;
            const auto add = surface ? 1U : traits.faceCount;
            if (total > limits.maximumFaceUses || add > limits.maximumFaceUses - total) {
                out.budgetExceeded = true;
                return out;
            }
            total += add;
            if (surface)
                ++out.boundaryCells;
            else
                ++out.volumeCells;
            const auto n = static_cast<Index>(traits.nodeCount);
            if (c.first < 0 || c.first > value.connectivityCount() ||
                n > value.connectivityCount() - c.first) {
                ++out.invalidConnectivity;
                continue;
            }
            const auto *conn = value.cellNodes(i);
            bool invalid = false;
            for (Index j = 0; j < n; ++j) {
                if (conn[j] < 0 || conn[j] >= value.nodeCount())
                    invalid = true;
                for (Index k = 0; k < j; ++k)
                    if (conn[k] == conn[j])
                        invalid = true;
            }
            if (invalid)
                ++out.invalidConnectivity;
        }
        if (out.invalidConnectivity || out.unsupportedCells || out.nonfiniteNodes) {
            out.complete = true;
            return out;
        }
        if (total > std::numeric_limits<std::size_t>::max() / sizeof(FaceUse)) {
            out.budgetExceeded = true;
            return out;
        }
        std::vector<FaceUse> faces;
        faces.reserve(static_cast<std::size_t>(total));
        out.temporaryBytes = total * sizeof(FaceUse); // benchmark telemetry only
        for (Index i = 0; i < value.cellCount(); ++i) {
            const auto &c = value.cell(i);
            if (c.type == CellType::Line2)
                continue;
            const auto &traits = cellTraits(c.type);
            const auto *conn = value.cellNodes(i);
            if (traits.dimension == 2) {
                std::array<Index, 4> ids{-1, -1, -1, -1};
                for (std::uint8_t k = 0; k < traits.cornerCount; ++k)
                    ids[k] = conn[k];
                faces.push_back(makeFace(ids, traits.cornerCount, true));
            } else
                for (std::uint8_t j = 0; j < traits.faceCount; ++j) {
                    const auto &f = traits.faces[j];
                    std::array<Index, 4> ids{-1, -1, -1, -1};
                    for (std::uint8_t k = 0; k < f.n; ++k)
                        ids[k] = conn[f.v[k]];
                    faces.push_back(makeFace(ids, f.n, false));
                }
        }
        std::sort(faces.begin(), faces.end(), [](const auto &a, const auto &b) {
            return a.n != b.n ? a.n < b.n : a.key < b.key;
        });
        for (std::size_t i = 0; i < faces.size();) {
            std::size_t end = i + 1;
            while (end < faces.size() && sameKey(faces[i], faces[end]))
                ++end;
            std::uint64_t volumes = 0, boundaries = 0;
            int orientation = 0, boundaryOrientation = 0;
            for (std::size_t j = i; j < end; ++j)
                if (faces[j].boundary) {
                    ++boundaries;
                    boundaryOrientation += faces[j].orientation;
                } else {
                    ++volumes;
                    orientation += faces[j].orientation;
                }
            if (volumes > 2)
                ++out.nonmanifoldFaces;
            if (boundaries > 1)
                ++out.duplicateBoundary;
            if (volumes == 1) {
                ++out.exteriorFaces;
                if (boundaries == 0)
                    ++out.missingBoundary;
                // Surface mesh adapters may orient a face to its geometric
                // chart, independently of the owning solid's outward sense.
                // Boundary bijection is required; boundary normal sense is
                // checked by geometry-specific ownership providers.
                (void)boundaryOrientation;
            } else if (volumes == 2) {
                ++out.interiorFaces;
                if (orientation != 0)
                    ++out.inconsistentOrientation;
                if (boundaries)
                    ++out.unexpectedBoundary;
            } else if (boundaries)
                ++out.unexpectedBoundary;
            i = end;
        }
        out.faceUses = total;
        out.complete = true;
        out.valid = out.volumeCells > 0 && !out.nonmanifoldFaces && !out.inconsistentOrientation &&
                    !out.missingBoundary && !out.unexpectedBoundary && !out.duplicateBoundary;
    } catch (...) {
        out.complete = false;
        out.valid = false;
    }
    return out;
}
} // namespace paacs::mesh
