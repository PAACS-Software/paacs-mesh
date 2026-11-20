// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/mesh/Mesh.hpp"

#include <array>
#include <cstddef>

namespace paacs::mesh::detail {

// Internal exact-capacity builder for already-validated deterministic engines.
// It removes the per-cell arity/range scan from hot structured generation,
// while preserving vector ownership and ordinary public Mesh output.
class CertifiedMeshAppender {
public:
    CertifiedMeshAppender(Mesh& mesh, std::size_t nodes,
                          std::size_t cells,
                          std::size_t connectivity)
        : mesh_(mesh) {
        mesh_.reserveAppend(nodes, cells, connectivity);
    }

    Index addNode(const Vec3& point, Classification classification) {
        return mesh_.addNode(point, classification);
    }

    template <std::size_t N>
    Index addCell(CellType type, const std::array<Index, N>& connectivity,
                  Classification classification) {
        return mesh_.appendCertifiedCell(type, connectivity, classification);
    }

private:
    Mesh& mesh_;
};

}  // namespace paacs::mesh::detail
