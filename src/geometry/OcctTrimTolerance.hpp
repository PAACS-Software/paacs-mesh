// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private native trimming evidence shared by decomposition adapters.
#include "OcctPartition.hpp"
#include <BRepBuilderAPI_Copy.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS_Shape.hxx>
#include <span>

namespace paacs::geo::detail {
using OcctWitnessShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;
struct OcctTrimToleranceInput {
    const TopoDS_Shape &originalShape;
    const BRepBuilderAPI_Copy &originalCopy;
    const OcctWitnessShapeMap &sourceSolids, &sourceFaces, &sourceEdges, &sourceVertices;
    const OcctWitnessShapeMap &derivedFaces, &derivedEdges, &derivedVertices;
    std::span<const int> selectedSourceBodies, parentBodyByPiece;
    std::span<const mesh::Classification> parentByFace, parentByEdge, parentByVertex;
    // Actual shared native face identities. Single owner: original exterior.
    // Two owners: opposite orientations in pieces of the same original body.
    std::span<const std::vector<std::pair<int, int>>> faceOwners;
    double sourceVolume = 0, derivedVolume = 0;
};
struct OcctTrimToleranceOptions {
    double maximumPhysicalDrift = 0;
    std::uint64_t maximumWork = 2'000'000;
    std::uint64_t maximumLogicalBytes = 128ULL * 1024 * 1024;
    std::uint64_t maximumSnapshotBytes = 64ULL * 1024 * 1024;
    const std::atomic_bool *cancellation = nullptr;
};
// Requires complete copy/split parent mappings, retained original topology,
// unchanged copied surface handles/locations and native-tolerance point/trim
// samples. It never establishes identity by proximity or mutates geometry.
// The volume envelope is sampled trim evidence, not a rigorous continuous
// conservation proof. Full original CAD inclusion/fidelity, exact mesh audits,
// body ownership and final volume limits remain mandatory caller checks.
OcctTrimToleranceReceipt auditOcctTrimTolerance(const OcctTrimToleranceInput &,
                                                const OcctTrimToleranceOptions &) noexcept;
const char *occtTrimToleranceReasonName(OcctTrimToleranceReason) noexcept;
} // namespace paacs::geo::detail
