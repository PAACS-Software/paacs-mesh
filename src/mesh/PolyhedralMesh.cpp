// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/mesh/PolyhedralMesh.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>

namespace paacs::mesh {
namespace {

bool finite(const Vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

bool current(Index value, std::size_t size) noexcept {
    return value >= 0 && static_cast<std::uint64_t>(value) < size;
}

}  // namespace

Index PolyhedralMesh::addNode(const Vec3& p, Classification c) {
    if (nodes_.size() >= static_cast<std::size_t>(
            (std::numeric_limits<Index>::max)()))
        throw std::length_error("PolyhedralMesh::addNode overflow");
    nodes_.push_back(p);
    try {
        nodeClass_.push_back(c);
    } catch (...) {
        nodes_.pop_back();
        throw;
    }
    return static_cast<Index>(nodes_.size() - 1);
}

Index PolyhedralMesh::addFace(std::span<const Index> nodes, Index owner,
                              Index neighbour, Classification c) {
    if (nodes.size() < 3 ||
        nodes.size() > (std::numeric_limits<std::uint32_t>::max)())
        throw std::invalid_argument("PolyhedralMesh::addFace invalid arity");
    if (faces_.size() >= static_cast<std::size_t>(
            (std::numeric_limits<Index>::max)()) ||
        faceNodes_.size() > static_cast<std::size_t>(
            (std::numeric_limits<Index>::max)()) - nodes.size())
        throw std::length_error("PolyhedralMesh::addFace overflow");
    for (const Index node : nodes)
        if (!current(node, nodes_.size()))
            throw std::invalid_argument(
                "PolyhedralMesh::addFace node out of range");
    const Index first = static_cast<Index>(faceNodes_.size());
    faceNodes_.insert(faceNodes_.end(), nodes.begin(), nodes.end());
    try {
        faces_.push_back({first, static_cast<std::uint32_t>(nodes.size()),
                          owner, neighbour, c});
    } catch (...) {
        faceNodes_.resize(static_cast<std::size_t>(first));
        throw;
    }
    return static_cast<Index>(faces_.size() - 1);
}

Index PolyhedralMesh::addCell(std::span<const OrientedPolyFace> faces,
                              Classification c) {
    if (faces.size() < 4 ||
        faces.size() > (std::numeric_limits<std::uint32_t>::max)())
        throw std::invalid_argument("PolyhedralMesh::addCell invalid arity");
    if (cells_.size() >= static_cast<std::size_t>(
            (std::numeric_limits<Index>::max)()) ||
        cellFaces_.size() > static_cast<std::size_t>(
            (std::numeric_limits<Index>::max)()) - faces.size())
        throw std::length_error("PolyhedralMesh::addCell overflow");
    for (const auto& ref : faces)
        if (!current(ref.face, faces_.size()) ||
            (ref.sense != 1 && ref.sense != -1))
            throw std::invalid_argument(
                "PolyhedralMesh::addCell face reference invalid");
    const Index first = static_cast<Index>(cellFaces_.size());
    cellFaces_.insert(cellFaces_.end(), faces.begin(), faces.end());
    try {
        cells_.push_back({first, static_cast<std::uint32_t>(faces.size()), c});
    } catch (...) {
        cellFaces_.resize(static_cast<std::size_t>(first));
        throw;
    }
    return static_cast<Index>(cells_.size() - 1);
}

void PolyhedralMesh::setFaceNeighbour(Index faceIndex, Index neighbour) {
    faces_.at(static_cast<std::size_t>(faceIndex)).neighbour = neighbour;
}

void PolyhedralMesh::setFaceClassification(Index faceIndex,
                                           Classification c) {
    faces_.at(static_cast<std::size_t>(faceIndex)).onEntity = c;
}

void PolyhedralMesh::reserve(std::size_t nodes, std::size_t faces,
                             std::size_t faceNodes, std::size_t cells,
                             std::size_t cellFaces) {
    nodes_.reserve(nodes);
    nodeClass_.reserve(nodes);
    faces_.reserve(faces);
    faceNodes_.reserve(faceNodes);
    cells_.reserve(cells);
    cellFaces_.reserve(cellFaces);
}

void PolyhedralMesh::swap(PolyhedralMesh& other) noexcept {
    static_assert(noexcept(nodes_.swap(other.nodes_)));
    static_assert(noexcept(nodeClass_.swap(other.nodeClass_)));
    static_assert(noexcept(faces_.swap(other.faces_)));
    static_assert(noexcept(faceNodes_.swap(other.faceNodes_)));
    static_assert(noexcept(cells_.swap(other.cells_)));
    static_assert(noexcept(cellFaces_.swap(other.cellFaces_)));
    nodes_.swap(other.nodes_);
    nodeClass_.swap(other.nodeClass_);
    faces_.swap(other.faces_);
    faceNodes_.swap(other.faceNodes_);
    cells_.swap(other.cells_);
    cellFaces_.swap(other.cellFaces_);
}

Index PolyhedralMesh::nodeCount() const noexcept {
    return static_cast<Index>(nodes_.size());
}
Index PolyhedralMesh::faceCount() const noexcept {
    return static_cast<Index>(faces_.size());
}
Index PolyhedralMesh::cellCount() const noexcept {
    return static_cast<Index>(cells_.size());
}
const Vec3& PolyhedralMesh::node(Index i) const {
    return nodes_.at(static_cast<std::size_t>(i));
}
const Classification& PolyhedralMesh::nodeClass(Index i) const {
    return nodeClass_.at(static_cast<std::size_t>(i));
}
const PolyFace& PolyhedralMesh::face(Index i) const {
    return faces_.at(static_cast<std::size_t>(i));
}
const PolyCell& PolyhedralMesh::cell(Index i) const {
    return cells_.at(static_cast<std::size_t>(i));
}
std::span<const Index> PolyhedralMesh::faceNodes(Index i) const {
    const auto& row = face(i);
    return std::span<const Index>(
        faceNodes_.data() + static_cast<std::size_t>(row.firstNode),
        row.nodeCount);
}
std::span<const OrientedPolyFace> PolyhedralMesh::cellFaces(Index i) const {
    const auto& row = cell(i);
    return std::span<const OrientedPolyFace>(
        cellFaces_.data() + static_cast<std::size_t>(row.firstFace),
        row.faceCount);
}

PolyhedralCapacitySnapshot PolyhedralMesh::capacitySnapshot() const noexcept {
    return {nodes_.capacity(), nodeClass_.capacity(), faces_.capacity(),
            faceNodes_.capacity(), cells_.capacity(), cellFaces_.capacity()};
}


std::uint64_t polyhedralStrictHash(const PolyhedralMesh& value) noexcept {
    std::uint64_t hash = UINT64_C(1469598103934665603);
    auto mix = [&](std::uint64_t bits) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            hash ^= (bits >> (byte * 8U)) & UINT64_C(0xff);
            hash *= UINT64_C(1099511628211);
        }
    };
    mix(static_cast<std::uint64_t>(value.nodeCount()));
    mix(static_cast<std::uint64_t>(value.faceCount()));
    mix(static_cast<std::uint64_t>(value.cellCount()));
    for (Index i = 0; i < value.nodeCount(); ++i) {
        const auto& p = value.node(i);
        mix(std::bit_cast<std::uint64_t>(p.x));
        mix(std::bit_cast<std::uint64_t>(p.y));
        mix(std::bit_cast<std::uint64_t>(p.z));
        mix(static_cast<std::uint32_t>(value.nodeClass(i).dim));
        mix(static_cast<std::uint32_t>(value.nodeClass(i).tag));
    }
    for (Index i = 0; i < value.faceCount(); ++i) {
        const auto& face = value.face(i);
        mix(static_cast<std::uint64_t>(face.nodeCount));
        mix(static_cast<std::uint64_t>(face.owner));
        mix(static_cast<std::uint64_t>(face.neighbour));
        mix(static_cast<std::uint32_t>(face.onEntity.dim));
        mix(static_cast<std::uint32_t>(face.onEntity.tag));
        for (const Index node : value.faceNodes(i))
            mix(static_cast<std::uint64_t>(node));
    }
    for (Index i = 0; i < value.cellCount(); ++i) {
        const auto& cell = value.cell(i);
        mix(static_cast<std::uint64_t>(cell.faceCount));
        mix(static_cast<std::uint32_t>(cell.onEntity.dim));
        mix(static_cast<std::uint32_t>(cell.onEntity.tag));
        for (const auto& ref : value.cellFaces(i)) {
            mix(static_cast<std::uint64_t>(ref.face));
            mix(static_cast<std::uint8_t>(ref.sense));
        }
    }
    return hash;
}

PolyhedralAuditReport auditPolyhedralMesh(
        const PolyhedralMesh& value) noexcept {
    PolyhedralAuditReport out;
    out.nodes = static_cast<std::uint64_t>(value.nodeCount());
    out.faces = static_cast<std::uint64_t>(value.faceCount());
    out.cells = static_cast<std::uint64_t>(value.cellCount());
    try {
        for (Index i = 0; i < value.nodeCount(); ++i) {
            if (!finite(value.node(i))) ++out.nonfiniteNodes;
            const auto cls = value.nodeClass(i);
            if (cls.dim < 0 || cls.dim > 3 || cls.tag < 0)
                ++out.invalidClassifications;
        }

        std::vector<std::uint32_t> referenceCount(
            static_cast<std::size_t>(value.faceCount()), 0);
        std::vector<std::uint32_t> ownerCount(
            static_cast<std::size_t>(value.faceCount()), 0);
        std::vector<std::uint32_t> neighbourCount(
            static_cast<std::size_t>(value.faceCount()), 0);
        std::vector<std::uint8_t> faceValid(
            static_cast<std::size_t>(value.faceCount()), 1);

        for (Index fi = 0; fi < value.faceCount(); ++fi) {
            const auto& face = value.face(fi);
            const auto nodes = value.faceNodes(fi);
            bool bad = nodes.size() < 3;
            std::vector<Index> unique(nodes.begin(), nodes.end());
            for (const Index node : nodes)
                if (node < 0 || node >= value.nodeCount()) bad = true;
            std::sort(unique.begin(), unique.end());
            if (std::adjacent_find(unique.begin(), unique.end()) !=
                unique.end()) {
                ++out.duplicateFaceNodes;
                bad = true;
            }
            const bool ownerValid = face.owner >= 0 &&
                face.owner < value.cellCount();
            const bool neighbourValid = face.neighbour == kInvalidIndex ||
                (face.neighbour >= 0 &&
                 face.neighbour < value.cellCount() &&
                 face.neighbour != face.owner);
            if (!ownerValid || !neighbourValid) {
                ++out.invalidOwnerNeighbour;
                bad = true;
            }
            if (face.neighbour == kInvalidIndex) {
                ++out.boundaryFaces;
                if (face.onEntity.dim != 2 || face.onEntity.tag < 0)
                    ++out.invalidClassifications;
            } else {
                ++out.interiorFaces;
                if (face.onEntity.dim != 3 || face.onEntity.tag < 0)
                    ++out.invalidClassifications;
            }
            if (bad) {
                ++out.invalidFaceConnectivity;
                faceValid[static_cast<std::size_t>(fi)] = 0;
            }
        }

        double minVolume = (std::numeric_limits<double>::infinity)();
        double maxVolume = 0.0;
        long double totalVolume = 0.0L;
        for (Index ci = 0; ci < value.cellCount(); ++ci) {
            const auto& cell = value.cell(ci);
            const auto refs = value.cellFaces(ci);
            bool cellBad = refs.size() < 4;
            if (cell.onEntity.dim != 3 || cell.onEntity.tag < 0)
                ++out.invalidClassifications;
            std::vector<Index> uniqueFaces;
            std::vector<Index> uniqueNodes;
            uniqueFaces.reserve(refs.size());
            for (const auto& ref : refs) {
                if (ref.face < 0 || ref.face >= value.faceCount() ||
                    (ref.sense != 1 && ref.sense != -1)) {
                    cellBad = true;
                    continue;
                }
                const auto idx = static_cast<std::size_t>(ref.face);
                ++referenceCount[idx];
                const auto& face = value.face(ref.face);
                if (face.owner == ci && ref.sense == 1)
                    ++ownerCount[idx];
                else if (face.neighbour == ci && ref.sense == -1)
                    ++neighbourCount[idx];
                else
                    cellBad = true;
                uniqueFaces.push_back(ref.face);
                const auto nodes = value.faceNodes(ref.face);
                uniqueNodes.insert(uniqueNodes.end(),
                                   nodes.begin(), nodes.end());
            }
            std::sort(uniqueFaces.begin(), uniqueFaces.end());
            if (std::adjacent_find(uniqueFaces.begin(), uniqueFaces.end()) !=
                uniqueFaces.end()) {
                ++out.duplicateCellFaces;
                cellBad = true;
            }
            std::sort(uniqueNodes.begin(), uniqueNodes.end());
            uniqueNodes.erase(std::unique(uniqueNodes.begin(),
                                          uniqueNodes.end()),
                              uniqueNodes.end());
            if (uniqueNodes.size() < 4) cellBad = true;
            if (cellBad) ++out.invalidCellConnectivity;
            if (uniqueNodes.empty()) continue;

            Vec3 centre{};
            bool finiteCell = true;
            for (const Index node : uniqueNodes) {
                if (node < 0 || node >= value.nodeCount()) {
                    finiteCell = false;
                    break;
                }
                centre = centre + value.node(node);
            }
            if (!finiteCell) continue;
            centre = centre / static_cast<double>(uniqueNodes.size());
            long double volume = 0.0L;
            for (const auto& ref : refs) {
                if (ref.face < 0 || ref.face >= value.faceCount() ||
                    !faceValid[static_cast<std::size_t>(ref.face)])
                    continue;
                const auto nodes = value.faceNodes(ref.face);
                const auto oriented = [&](std::size_t ordinal) -> Index {
                    if (ref.sense == 1) return nodes[ordinal];
                    return ordinal == 0 ? nodes[0]
                                        : nodes[nodes.size() - ordinal];
                };
                const Vec3 a = value.node(oriented(0));
                for (std::size_t n = 1; n + 1 < nodes.size(); ++n) {
                    const Vec3 b = value.node(oriented(n));
                    const Vec3 c = value.node(oriented(n + 1));
                    volume += static_cast<long double>(
                        (a - centre).dot(
                            (b - centre).cross(c - centre))) / 6.0L;
                }
            }
            if (!std::isfinite(static_cast<double>(volume)) ||
                volume <= 0.0L) {
                ++out.nonPositiveCellVolumes;
            } else {
                const double cellVolume = static_cast<double>(volume);
                minVolume = std::min(minVolume, cellVolume);
                maxVolume = std::max(maxVolume, cellVolume);
                totalVolume += volume;
            }
        }

        for (Index fi = 0; fi < value.faceCount(); ++fi) {
            const auto idx = static_cast<std::size_t>(fi);
            const auto& face = value.face(fi);
            const std::uint32_t expected =
                face.neighbour == kInvalidIndex ? 1U : 2U;
            if (referenceCount[idx] != expected || ownerCount[idx] != 1U)
                ++out.ownerReferenceMismatch;
            if ((face.neighbour == kInvalidIndex &&
                 neighbourCount[idx] != 0U) ||
                (face.neighbour != kInvalidIndex &&
                 neighbourCount[idx] != 1U))
                ++out.neighbourReferenceMismatch;
        }

        if (value.cellCount() == 0) minVolume = 0.0;
        out.minimumCellVolume = minVolume;
        out.maximumCellVolume = maxVolume;
        out.totalVolume = static_cast<double>(totalVolume);
        out.valid = value.cellCount() != 0 && out.nonfiniteNodes == 0 &&
            out.invalidFaceConnectivity == 0 &&
            out.duplicateFaceNodes == 0 &&
            out.invalidOwnerNeighbour == 0 &&
            out.invalidCellConnectivity == 0 &&
            out.duplicateCellFaces == 0 &&
            out.ownerReferenceMismatch == 0 &&
            out.neighbourReferenceMismatch == 0 &&
            out.invalidClassifications == 0 &&
            out.nonPositiveCellVolumes == 0;
    } catch (const std::bad_alloc&) {
        out.complete = false;
        out.valid = false;
        out.allocationFailed = true;
    } catch (...) {
        out.complete = false;
        out.valid = false;
    }
    return out;
}

}  // namespace paacs::mesh
