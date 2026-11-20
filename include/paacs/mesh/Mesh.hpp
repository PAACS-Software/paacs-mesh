// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// Cache-friendly array-of-structs mesh container. Every node and element is
// classified onto the geometry entity that owns it (dim, tag) so downstream
// consumers (solvers, UI, boundary-condition assignment) can query by CAD
// entity. Connectivity is flat with per-type offsets.
//
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <initializer_list>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>
#include "paacs/core/Types.hpp"

namespace paacs::mesh {

namespace detail {
class CertifiedMeshAppender;
class CertifiedMeshStorage;
class ReusableCertifiedMeshStorage;
}

enum class CellType : std::uint8_t {
    Line2 = 0, Line3, Tri3, Tri6, Quad4, Tet4, Tet10, Prism6, Hex8, Pyramid5
};

inline int nodesPerCell(CellType t) {
    switch (t) {
        case CellType::Line2:  return 2;  case CellType::Line3: return 3;
        case CellType::Tri3:   return 3;  case CellType::Tri6:  return 6;
        case CellType::Quad4:  return 4;  case CellType::Tet4:  return 4;
        case CellType::Tet10:  return 10; case CellType::Prism6: return 6;
        case CellType::Hex8:   return 8; case CellType::Pyramid5: return 5;
    }
    return 0;
}

struct Classification {
    int dim = -1;
    int tag = -1;

    friend bool operator==(const Classification&,
                           const Classification&) = default;
};

struct Cell {
    CellType type;
    Index    first;   // offset into Mesh::connectivity
    Classification onEntity;
};

struct MeshCapacitySnapshot {
    std::size_t nodeCoordinateCapacity = 0;
    std::size_t nodeClassificationCapacity = 0;
    std::size_t cellCapacity = 0;
    std::size_t connectivityCapacity = 0;
    std::size_t dedupBucketCount = 0;
    std::size_t dedupEntryCount = 0;
    std::size_t dedupIndexCapacity = 0;
    bool complete = true;
};

// Shared only so the contact's otherwise-unreachable overflow branch can be
// pinned directly. It performs one saturating checked addition and never
// allocates or mutates mesh state.
inline bool meshCapacityCheckedAdd(std::size_t& total,
                                   std::size_t value) noexcept {
    if (value > (std::numeric_limits<std::size_t>::max)() - total) {
        total = (std::numeric_limits<std::size_t>::max)();
        return false;
    }
    total += value;
    return true;
}

class Mesh {
public:
    // ---- nodes ----
    Index addNode(const Vec3& p, Classification c = {}) {
        if (nodes_.size() >= static_cast<std::size_t>((std::numeric_limits<Index>::max)()))
            throw std::length_error("Mesh::addNode: index overflow");
        nodes_.push_back(p);
        try {
            nodeClass_.push_back(c);
        } catch (...) {
            nodes_.pop_back();
            throw;
        }
        return static_cast<Index>(nodes_.size()) - 1;
    }
    // Pre-latch the dedup tolerance (normally derived from model scale by
    // the pipeline: ~1e-9 x bbox diagonal). Must be called before the first
    // addNodeUnique; changing it afterwards throws.
    void setDedupTolerance(double tol) {
        if (!(tol > 0.0) || !std::isfinite(tol))
            throw std::invalid_argument("Mesh::setDedupTolerance: expected finite positive tolerance");
        // Exactly immutable once latched: an epsilon comparison here would
        // let tiny models change the spatial-index cell size after nodes
        // were registered (reviewer finding).
        if (dedupTol_ >= 0 && tol != dedupTol_)
            throw std::logic_error(
                "Mesh::setDedupTolerance: tolerance already latched");
        dedupTol_ = tol;
    }
    double dedupTolerance() const { return dedupTol_; }

    // Deduplicating add: nodes within the mesh's dedup tolerance snap to an
    // existing node. Grid-hash lookup, O(1) expected.
    //
    // Lifecycle rules (each was a latent bug before):
    // - One FIXED tolerance per mesh: the first call latches it; later calls
    //   passing a different tolerance throw, because buckets built with one
    //   cell size cannot answer queries at another.
    // - setNode() marks the spatial index dirty; the next addNodeUnique()
    //   rebuilds it from current coordinates, so moved nodes are found at
    //   their new positions, not their insertion positions.
    Index addNodeUnique(const Vec3& p, Classification c = {}) {
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
            throw std::invalid_argument("Mesh::addNodeUnique: non-finite position");
        if (dedupTol_ < 0) dedupTol_ = 1e-9;   // default; see setDedupTolerance
        if (dedupDirty_) rebuildSpatialIndex();

        const double cellSize = dedupTol_ * 4.0;
        const long long ix = cellCoord(p.x, cellSize);
        const long long iy = cellCoord(p.y, cellSize);
        const long long iz = cellCoord(p.z, cellSize);
        for (long long dx = -1; dx <= 1; ++dx)
            for (long long dy = -1; dy <= 1; ++dy)
                for (long long dz = -1; dz <= 1; ++dz) {
                    auto it = dedupGrid_.find(
                        dedupKey(ix + dx, iy + dy, iz + dz));
                    if (it == dedupGrid_.end()) continue;
                    for (Index id : it->second)
                        if ((nodes_[static_cast<size_t>(id)] - p).norm() <=
                            dedupTol_)
                            return id;
                }
        Index id = addNode(p, c);
        try {
            dedupGrid_[dedupKey(ix, iy, iz)].push_back(id);
        } catch (...) {
            nodes_.pop_back();
            nodeClass_.pop_back();
            throw;
        }
        return id;
    }

    // Rebuild the dedup index from current node coordinates. Indexes only
    // nodes previously registered through addNodeUnique, preserving the
    // invariant that plain addNode() nodes never participate in dedup.
    void rebuildSpatialIndex() {
        std::vector<Index> registered;
        for (const auto& [k, ids] : dedupGrid_)
            registered.insert(registered.end(), ids.begin(), ids.end());
        std::sort(registered.begin(), registered.end());
        decltype(dedupGrid_) rebuilt;
        rebuilt.reserve(dedupGrid_.size());
        const double cellSize = dedupTol_ * 4.0;
        for (Index id : registered) {
            const Vec3& q = nodes_[static_cast<size_t>(id)];
            if (!std::isfinite(q.x) || !std::isfinite(q.y) || !std::isfinite(q.z))
                throw std::invalid_argument("Mesh::rebuildSpatialIndex: non-finite position");
            rebuilt[dedupKey(cellCoord(q.x, cellSize),
                                cellCoord(q.y, cellSize),
                                cellCoord(q.z, cellSize))].push_back(id);
        }
        dedupGrid_.swap(rebuilt);
        dedupDirty_ = false;
    }

    Index nodeCount() const { return static_cast<Index>(nodes_.size()); }
    const Vec3& node(Index i) const { return nodes_[static_cast<size_t>(i)]; }
    void setNode(Index i, const Vec3& p) {
        nodes_[static_cast<size_t>(i)] = p;
        dedupDirty_ = true;   // spatial index no longer matches coordinates
    }
    const std::vector<Vec3>& nodes() const { return nodes_; }
    const Classification& nodeClass(Index i) const { return nodeClass_[static_cast<size_t>(i)]; }

    // ---- cells ----
    Index addCell(CellType t, std::span<const Index> conn, Classification c = {}) {
        const int expected = nodesPerCell(t);
        if (expected <= 0 || conn.size() != static_cast<std::size_t>(expected))
            throw std::invalid_argument(
                "Mesh::addCell: connectivity size does not match cell type");
        for (Index v : conn)
            if (v < 0 || v >= nodeCount())
                throw std::invalid_argument(
                    "Mesh::addCell: node index out of range");
        constexpr auto limit = static_cast<std::size_t>((std::numeric_limits<Index>::max)());
        if (cells_.size() >= limit || conn.size() > limit - connectivity_.size())
            throw std::length_error("Mesh::addCell: index overflow");
        const auto oldSize = connectivity_.size();
        Cell cell{t, static_cast<Index>(oldSize), c};
        // A span may refer to this mesh's own connectivity. Copy the bounded
        // input before vector growth can invalidate it.
        std::array<Index, 10> stable{};
        std::copy(conn.begin(), conn.end(), stable.begin());
        connectivity_.insert(connectivity_.end(), stable.begin(), stable.begin() + expected);
        try {
            cells_.push_back(cell);
        } catch (...) {
            connectivity_.resize(oldSize);
            throw;
        }
        return static_cast<Index>(cells_.size()) - 1;
    }
    Index addCell(CellType t, const std::vector<Index>& conn, Classification c = {}) {
        return addCell(t, std::span<const Index>(conn), c);
    }
    Index addCell(CellType t, std::initializer_list<Index> conn, Classification c = {}) {
        return addCell(t, std::span<const Index>(conn.begin(), conn.size()), c);
    }
    Index cellCount() const { return static_cast<Index>(cells_.size()); }
    Index connectivityCount() const {
        return static_cast<Index>(connectivity_.size());
    }
    const Cell& cell(Index i) const { return cells_[static_cast<size_t>(i)]; }
    const Index* cellNodes(Index i) const {
        return connectivity_.data() + cells_[static_cast<size_t>(i)].first;
    }

    Index countCells(CellType t) const {
        Index n = 0;
        for (const auto& c : cells_) if (c.type == t) ++n;
        return n;
    }

    MeshCapacitySnapshot capacitySnapshot() const noexcept {
        MeshCapacitySnapshot out;
        out.nodeCoordinateCapacity = nodes_.capacity();
        out.nodeClassificationCapacity = nodeClass_.capacity();
        out.cellCapacity = cells_.capacity();
        out.connectivityCapacity = connectivity_.capacity();
        out.dedupBucketCount = dedupGrid_.bucket_count();
        out.dedupEntryCount = dedupGrid_.size();
        for (const auto& [key, indices] : dedupGrid_) {
            (void)key;
            if (!meshCapacityCheckedAdd(out.dedupIndexCapacity,
                                        indices.capacity())) {
                out.complete = false;
                break;
            }
        }
        return out;
    }

    // Increment-5 publication preparation. Reserves every flat-mesh buffer
    // before the first logical append; a reserve failure may leave capacity
    // residue but cannot expose a node or cell. The subsequent bounded
    // publisher uses addNode()/addCell() only.
    void reserveAppend(std::size_t nodeCount, std::size_t cellCount,
                       std::size_t connectivityCount) {
        constexpr auto limit = static_cast<std::size_t>((std::numeric_limits<Index>::max)());
        if (nodeCount > limit - nodes_.size() ||
            cellCount > limit - cells_.size() ||
            connectivityCount > limit - connectivity_.size() ||
            nodeCount > nodes_.max_size() - nodes_.size() ||
            nodeCount > nodeClass_.max_size() - nodeClass_.size() ||
            cellCount > cells_.max_size() - cells_.size() ||
            connectivityCount >
                connectivity_.max_size() - connectivity_.size())
            throw std::length_error("Mesh::reserveAppend: size overflow");
        nodes_.reserve(nodes_.size() + nodeCount);
        nodeClass_.reserve(nodeClass_.size() + nodeCount);
        cells_.reserve(cells_.size() + cellCount);
        connectivity_.reserve(connectivity_.size() + connectivityCount);
    }

    // CK-2 whole-pipeline publication. A completed local mesh is made global
    // by one allocation-free ownership exchange. Default allocators and the
    // exact hash/equality members make every container swap non-throwing;
    // keep these assertions beside the operation so a future member change
    // cannot silently weaken atomic publication.
    void swap(Mesh& other) noexcept {
        static_assert(noexcept(nodes_.swap(other.nodes_)));
        static_assert(noexcept(nodeClass_.swap(other.nodeClass_)));
        static_assert(noexcept(cells_.swap(other.cells_)));
        static_assert(noexcept(connectivity_.swap(other.connectivity_)));
        static_assert(noexcept(dedupGrid_.swap(other.dedupGrid_)));
        nodes_.swap(other.nodes_);
        nodeClass_.swap(other.nodeClass_);
        cells_.swap(other.cells_);
        connectivity_.swap(other.connectivity_);
        dedupGrid_.swap(other.dedupGrid_);
        using std::swap;
        swap(dedupTol_, other.dedupTol_);
        swap(dedupDirty_, other.dedupDirty_);
    }

private:
    friend class detail::CertifiedMeshAppender;
    friend class detail::CertifiedMeshStorage;
    friend class detail::ReusableCertifiedMeshStorage;

    template <std::size_t N>
    Index appendCertifiedCell(CellType type,
                              const std::array<Index, N>& conn,
                              Classification classification) {
        Cell cell{type, static_cast<Index>(connectivity_.size()),
                  classification};
        connectivity_.insert(connectivity_.end(), conn.begin(), conn.end());
        cells_.push_back(cell);
        return static_cast<Index>(cells_.size()) - 1;
    }

    std::vector<Vec3>           nodes_;
    std::vector<Classification> nodeClass_;
    std::vector<Cell>           cells_;
    std::vector<Index>          connectivity_;
    // Overflow-safe (PM-NUM-001): tiny tolerances or far-translated models
    // produce cell indices beyond 32-bit (and `long` is 32-bit on Windows).
    static long long cellCoord(double x, double cellSize) {
        const double c = std::floor(x / cellSize);
        constexpr double kLim = 4.0e18;
        return static_cast<long long>(std::clamp(c, -kLim, kLim));
    }
    static std::uint64_t dedupKey(long long ix, long long iy, long long iz) {
        return (static_cast<std::uint64_t>(ix) * 73856093ULL) ^
               (static_cast<std::uint64_t>(iy) * 19349663ULL) ^
               (static_cast<std::uint64_t>(iz) * 83492791ULL);
    }
    std::unordered_map<std::uint64_t, std::vector<Index>> dedupGrid_;
    double dedupTol_ = -1.0;    // latched by first addNodeUnique
    bool   dedupDirty_ = false;
};

} // namespace paacs::mesh
