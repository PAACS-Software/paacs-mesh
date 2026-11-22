// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Internal v4 handoff seam. This deliberately does not live in the installed
// include tree: callers must already own the private, audited surface-artifact
// rows assembled by BoundaryConformingMesher. Public callers cannot supply an
// unstamped endpoint/classification list to the partial-import policy.

#include "paacs/algo/TetCandidate.hpp"
#include "paacs/algo/detail/TetConstructionBudget.hpp"

namespace paacs::algo::tet_candidate_detail {

// Serial constructor whose predicate coordinates use only exact, reversible
// powers of two, with no rounded translation. Original coordinates and Morton
// ordering stay intact. This is a bounded retry of a rejected candidate.
std::vector<Tet> delaunay3DExactInputCoordinates(
    const std::vector<Vec3>&, Delaunay3DStats*);
std::vector<Tet> delaunay3DBounded(
    const std::vector<Vec3>&, Delaunay3DStats*,
    detail::TetConstructionBudget&, bool exactCoordinates);
std::uint64_t measuredDelaunayConstructionWork(const Delaunay3DStats&) noexcept;
void retainDiscardedDelaunayAttempt(Delaunay3DStats& accepted,
                                    const Delaunay3DStats& discarded) noexcept;

TetHandoffOutcome importTetCandidateWithPartialProtectedBoundary(
    const TetCandidate& candidate,
    std::span<const Index> selectedCandidateCells,
    std::span<const Index> existingGlobalByPoint,
    std::span<const mesh::Classification> expectedClassByPoint,
    const mesh::Mesh& globalMesh,
    mesh::Classification cellClassification,
    std::span<const TetImportProtectedEdge> protectedEdges,
    std::span<const TetImportProtectedFace> protectedFaces) noexcept;

}  // namespace paacs::algo::tet_candidate_detail
