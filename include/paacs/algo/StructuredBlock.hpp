// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/core/Types.hpp"
#include "paacs/mesh/Mesh.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace paacs::algo {

// Canonical face order shared by the mapped Hex8 and structured-polyhedral
// generators.  It is identical to topo::ElementTraits Hex8 local-face order:
// bottom/top in W, then min-V, max-U, max-V and min-U.
enum class StructuredBlockVertex : std::uint8_t {
    UMinVMinWMin = 0,
    UMaxVMinWMin,
    UMaxVMaxWMin,
    UMinVMaxWMin,
    UMinVMinWMax,
    UMaxVMinWMax,
    UMaxVMaxWMax,
    UMinVMaxWMax,
    Count_
};

// Canonical edge order is identical to topo::ElementTraits Hex8 edges.
enum class StructuredBlockEdge : std::uint8_t {
    VMinWMin = 0,
    UMaxWMin,
    VMaxWMin,
    UMinWMin,
    VMinWMax,
    UMaxWMax,
    VMaxWMax,
    UMinWMax,
    UMinVMin,
    UMaxVMin,
    UMaxVMax,
    UMinVMax,
    Count_
};

enum class StructuredBlockFace : std::uint8_t {
    WMin = 0,
    WMax,
    VMin,
    UMax,
    VMax,
    UMin,
    Count_
};

// An affine block. axisU/axisV/axisW span the complete physical block, not one
// cell. A positive scalar triple product is required so generated Hex8 cells
// follow PAACS's frozen positive-orientation convention.
struct StructuredBlockGeometry {
    Vec3 origin{};
    Vec3 axisU{1.0, 0.0, 0.0};
    Vec3 axisV{0.0, 1.0, 0.0};
    Vec3 axisW{0.0, 0.0, 1.0};
    std::uint32_t cellsU = 1;
    std::uint32_t cellsV = 1;
    std::uint32_t cellsW = 1;
    mesh::Classification volumeClassification{3, 1};
    std::array<mesh::Classification, 8> vertexClassifications{
        mesh::Classification{0, 1}, mesh::Classification{0, 2},
        mesh::Classification{0, 3}, mesh::Classification{0, 4},
        mesh::Classification{0, 5}, mesh::Classification{0, 6},
        mesh::Classification{0, 7}, mesh::Classification{0, 8}};
    std::array<mesh::Classification, 12> edgeClassifications{
        mesh::Classification{1, 1}, mesh::Classification{1, 2},
        mesh::Classification{1, 3}, mesh::Classification{1, 4},
        mesh::Classification{1, 5}, mesh::Classification{1, 6},
        mesh::Classification{1, 7}, mesh::Classification{1, 8},
        mesh::Classification{1, 9}, mesh::Classification{1, 10},
        mesh::Classification{1, 11}, mesh::Classification{1, 12}};
    std::array<mesh::Classification, 6> boundaryClassifications{
        mesh::Classification{2, 1}, mesh::Classification{2, 2},
        mesh::Classification{2, 3}, mesh::Classification{2, 4},
        mesh::Classification{2, 5}, mesh::Classification{2, 6}};

    // Hard, caller-visible budgets. The defaults are nonbinding for ordinary
    // builds but still make all count arithmetic explicit and overflow-safe.
    std::uint64_t maxNodes = (std::numeric_limits<std::uint64_t>::max)();
    std::uint64_t maxCells = (std::numeric_limits<std::uint64_t>::max)();
    std::uint64_t maxFaces = (std::numeric_limits<std::uint64_t>::max)();
};

// Returns the owning topological entity for a structured node. Corners and
// edges never inherit an arbitrary incident face identity; all 0D/1D/2D/3D
// classifications are explicit in StructuredBlockGeometry.
inline mesh::Classification structuredBlockNodeClassification(
        const StructuredBlockGeometry& g,
        std::uint64_t i, std::uint64_t j, std::uint64_t k) noexcept {
    const bool u0 = i == 0;
    const bool u1 = i == g.cellsU;
    const bool v0 = j == 0;
    const bool v1 = j == g.cellsV;
    const bool w0 = k == 0;
    const bool w1 = k == g.cellsW;
    const unsigned boundaryCount = static_cast<unsigned>(u0 || u1) +
        static_cast<unsigned>(v0 || v1) +
        static_cast<unsigned>(w0 || w1);

    if (boundaryCount == 3) {
        std::size_t vertex = 0;
        if (w1) vertex += 4;
        if (v1) vertex += u1 ? 2 : 3;
        else if (u1) vertex += 1;
        return g.vertexClassifications[vertex];
    }
    if (boundaryCount == 2) {
        std::size_t edge = 0;
        if (w0 && v0) edge = 0;
        else if (w0 && u1) edge = 1;
        else if (w0 && v1) edge = 2;
        else if (w0 && u0) edge = 3;
        else if (w1 && v0) edge = 4;
        else if (w1 && u1) edge = 5;
        else if (w1 && v1) edge = 6;
        else if (w1 && u0) edge = 7;
        else if (u0 && v0) edge = 8;
        else if (u1 && v0) edge = 9;
        else if (u1 && v1) edge = 10;
        else edge = 11;
        return g.edgeClassifications[edge];
    }
    if (boundaryCount == 1) {
        std::size_t face = 0;
        if (w0) face = static_cast<std::size_t>(StructuredBlockFace::WMin);
        else if (w1)
            face = static_cast<std::size_t>(StructuredBlockFace::WMax);
        else if (v0)
            face = static_cast<std::size_t>(StructuredBlockFace::VMin);
        else if (u1)
            face = static_cast<std::size_t>(StructuredBlockFace::UMax);
        else if (v1)
            face = static_cast<std::size_t>(StructuredBlockFace::VMax);
        else
            face = static_cast<std::size_t>(StructuredBlockFace::UMin);
        return g.boundaryClassifications[face];
    }
    return g.volumeClassification;
}

}  // namespace paacs::algo
