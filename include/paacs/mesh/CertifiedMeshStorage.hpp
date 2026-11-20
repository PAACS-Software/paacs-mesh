// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/mesh/Mesh.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

namespace paacs::mesh::detail {

// Exact-size, fixed-range construction storage for deterministic engines.
// All vectors are resized once before any worker starts. Distinct workers may
// write disjoint node/cell/connectivity ordinals without locks. Publication is
// a single noexcept ownership exchange into the ordinary public Mesh format.
// Each storage object is single-use: call publishTo exactly once.
class CertifiedMeshStorage {
public:
    CertifiedMeshStorage(std::size_t nodes, std::size_t cells,
                         std::size_t connectivity)
        : nodes_(nodes), nodeClass_(nodes), cells_(cells),
          connectivity_(connectivity) {}

    void setNode(std::size_t ordinal, const Vec3& point,
                 Classification classification) noexcept {
        nodes_[ordinal] = point;
        nodeClass_[ordinal] = classification;
    }

    template <std::size_t N>
    void setCell(std::size_t cellOrdinal, std::size_t firstConnectivity,
                 CellType type, const std::array<Index, N>& connectivity,
                 Classification classification) noexcept {
        cells_[cellOrdinal] = {
            type, static_cast<Index>(firstConnectivity), classification};
        std::copy(connectivity.begin(), connectivity.end(),
                  connectivity_.begin() +
                      static_cast<std::ptrdiff_t>(firstConnectivity));
    }

    std::size_t nodeCount() const noexcept { return nodes_.size(); }
    // Read only after the node-writing phase has joined.
    const Vec3 &node(std::size_t ordinal) const noexcept { return nodes_[ordinal]; }
    std::size_t cellCount() const noexcept { return cells_.size(); }
    std::size_t connectivityCount() const noexcept {
        return connectivity_.size();
    }

    void publishTo(Mesh& destination) noexcept {
        replacement_.nodes_.swap(nodes_);
        replacement_.nodeClass_.swap(nodeClass_);
        replacement_.cells_.swap(cells_);
        replacement_.connectivity_.swap(connectivity_);
        destination.swap(replacement_);
    }

private:
    // An empty Mesh may allocate for its deduplication table on some standard
    // libraries. Construct it during fallible preparation, before publication.
    Mesh replacement_;
    std::vector<Vec3> nodes_;
    std::vector<Classification> nodeClass_;
    std::vector<Cell> cells_;
    std::vector<Index> connectivity_;
};

}  // namespace paacs::mesh::detail
