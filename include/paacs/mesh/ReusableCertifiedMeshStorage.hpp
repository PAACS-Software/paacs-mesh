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

// Reusable fixed-range construction storage for deterministic generators.
// It owns the same four flat vectors as Mesh and swaps them directly at
// publication. The old destination buffers become the next private candidate,
// providing a true allocation-free double buffer after one priming exchange.
class ReusableCertifiedMeshStorage {
public:
    bool prepare(std::size_t nodes, std::size_t cells,
                 std::size_t connectivity) {
        const bool reused = nodes_.size() == nodes &&
            nodeClass_.size() == nodes && cells_.size() == cells &&
            connectivity_.size() == connectivity;
        if (!reused) {
            // Construct exact-size replacements so the cold path retains the
            // optimized vector(count) behavior of V8.
            std::vector<Vec3> newNodes(nodes);
            std::vector<Classification> newNodeClass(nodes);
            std::vector<Cell> newCells(cells);
            std::vector<Index> newConnectivity(connectivity);
            nodes_.swap(newNodes);
            nodeClass_.swap(newNodeClass);
            cells_.swap(newCells);
            connectivity_.swap(newConnectivity);
        }
        return reused;
    }

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
    std::size_t cellCount() const noexcept { return cells_.size(); }
    std::size_t connectivityCount() const noexcept {
        return connectivity_.size();
    }

    std::uint64_t logicalBytes() const noexcept {
        return static_cast<std::uint64_t>(nodes_.size()) * sizeof(Vec3) +
            static_cast<std::uint64_t>(nodeClass_.size()) *
                sizeof(Classification) +
            static_cast<std::uint64_t>(cells_.size()) * sizeof(Cell) +
            static_cast<std::uint64_t>(connectivity_.size()) * sizeof(Index);
    }

    void release() noexcept {
        std::vector<Vec3>().swap(nodes_);
        std::vector<Classification>().swap(nodeClass_);
        std::vector<Cell>().swap(cells_);
        std::vector<Index>().swap(connectivity_);
    }

    // One ownership exchange. Destination's former vectors become the next
    // private workspace buffers; no candidate data is copied.
    void publishTo(Mesh& destination) noexcept {
        destination.nodes_.swap(nodes_);
        destination.nodeClass_.swap(nodeClass_);
        destination.cells_.swap(cells_);
        destination.connectivity_.swap(connectivity_);
        destination.dedupGrid_.clear();
        destination.dedupTol_ = -1.0;
        destination.dedupDirty_ = false;
    }

    // Takes ownership of a completed/retired ordinary Mesh buffer set. Used
    // after full audit so the workspace still retains the old destination
    // buffers for the next request.
    void acquireFrom(Mesh& source) noexcept {
        source.nodes_.swap(nodes_);
        source.nodeClass_.swap(nodeClass_);
        source.cells_.swap(cells_);
        source.connectivity_.swap(connectivity_);
        source.dedupGrid_.clear();
        source.dedupTol_ = -1.0;
        source.dedupDirty_ = false;
    }

private:
    std::vector<Vec3> nodes_;
    std::vector<Classification> nodeClass_;
    std::vector<Cell> cells_;
    std::vector<Index> connectivity_;
};

}  // namespace paacs::mesh::detail
