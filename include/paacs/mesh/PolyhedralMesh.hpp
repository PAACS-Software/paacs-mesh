// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/core/Types.hpp"
#include "paacs/mesh/Mesh.hpp"

#include <cstdint>
#include <span>
#include <vector>

namespace paacs::mesh {

struct PolyFace {
    Index firstNode = 0;
    std::uint32_t nodeCount = 0;
    Index owner = kInvalidIndex;
    Index neighbour = kInvalidIndex;
    Classification onEntity{};
};

struct OrientedPolyFace {
    Index face = kInvalidIndex;
    // +1 means the stored face loop is outward for this cell; -1 means the
    // reverse loop is outward. Zero is never valid.
    std::int8_t sense = 0;
};

struct PolyCell {
    Index firstFace = 0;
    std::uint32_t faceCount = 0;
    Classification onEntity{};
};

struct PolyhedralCapacitySnapshot {
    std::size_t nodeCapacity = 0;
    std::size_t nodeClassificationCapacity = 0;
    std::size_t faceCapacity = 0;
    std::size_t faceNodeCapacity = 0;
    std::size_t cellCapacity = 0;
    std::size_t cellFaceCapacity = 0;
};

// Face-based variable-arity volume mesh. It is intentionally separate from
// mesh::Mesh, whose fixed-width connectivity is the correct representation for
// standard FEM cells but cannot honestly encode arbitrary polyhedra.
class PolyhedralMesh {
public:
    Index addNode(const Vec3& p, Classification c = {});
    Index addFace(std::span<const Index> nodes, Index owner,
                  Index neighbour = kInvalidIndex,
                  Classification c = {});
    Index addCell(std::span<const OrientedPolyFace> faces,
                  Classification c = {});
    void setFaceNeighbour(Index face, Index neighbour);
    void setFaceClassification(Index face, Classification c);

    void reserve(std::size_t nodes, std::size_t faces,
                 std::size_t faceNodes, std::size_t cells,
                 std::size_t cellFaces);
    void swap(PolyhedralMesh& other) noexcept;

    Index nodeCount() const noexcept;
    Index faceCount() const noexcept;
    Index cellCount() const noexcept;

    const Vec3& node(Index i) const;
    const Classification& nodeClass(Index i) const;
    const PolyFace& face(Index i) const;
    const PolyCell& cell(Index i) const;
    std::span<const Index> faceNodes(Index i) const;
    std::span<const OrientedPolyFace> cellFaces(Index i) const;

    PolyhedralCapacitySnapshot capacitySnapshot() const noexcept;

private:
    std::vector<Vec3> nodes_;
    std::vector<Classification> nodeClass_;
    std::vector<PolyFace> faces_;
    std::vector<Index> faceNodes_;
    std::vector<PolyCell> cells_;
    std::vector<OrientedPolyFace> cellFaces_;
};

struct PolyhedralAuditReport {
    bool complete = true;
    bool valid = false;
    bool allocationFailed = false;
    std::uint64_t nodes = 0;
    std::uint64_t faces = 0;
    std::uint64_t cells = 0;
    std::uint64_t boundaryFaces = 0;
    std::uint64_t interiorFaces = 0;
    std::uint64_t nonfiniteNodes = 0;
    std::uint64_t invalidFaceConnectivity = 0;
    std::uint64_t duplicateFaceNodes = 0;
    std::uint64_t invalidOwnerNeighbour = 0;
    std::uint64_t invalidCellConnectivity = 0;
    std::uint64_t duplicateCellFaces = 0;
    std::uint64_t ownerReferenceMismatch = 0;
    std::uint64_t neighbourReferenceMismatch = 0;
    std::uint64_t invalidClassifications = 0;
    std::uint64_t nonPositiveCellVolumes = 0;
    double minimumCellVolume = 0.0;
    double maximumCellVolume = 0.0;
    double totalVolume = 0.0;
};

// Audits owner/neighbour consistency, oriented cell-face references, finite
// geometry and positive polyhedral volumes. Faces may contain any number of
// vertices >= 3; cells may contain any number of faces >= 4.
PolyhedralAuditReport auditPolyhedralMesh(
    const PolyhedralMesh& value) noexcept;

// Bitwise deterministic diagnostic hash of the complete face-based storage.
// This is a strict representation hash, not a topology-canonical semantic hash.
std::uint64_t polyhedralStrictHash(const PolyhedralMesh& value) noexcept;

}  // namespace paacs::mesh
