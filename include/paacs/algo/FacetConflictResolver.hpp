// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <array>
#include <cstdint>
#include <span>

#include "paacs/algo/FacetBarycentre.hpp"
#include "paacs/bnd/TetFacetSetRecovery.hpp"

namespace paacs::algo {

// Exact retained-artifact resolution for one normalized volume-facet key.
// The resolver returns copied scalar authority plus a read-only span into the
// retained immutable certificate boundary.  The retained SurfaceSetResult
// must therefore outlive the result; no attempt-local handle or flat index is
// exposed and the production overload performs no allocation.
enum class FacetConflictResolveStatus : std::uint8_t {
    Resolved,
    Unsupported,
    Invalid
};

enum class FacetConflictResolveReason : std::uint8_t {
    None,
    MissingRetainedArtifact,
    InvalidFaceIdentity,
    FaceNotFound,
    AmbiguousFace,
    InvalidFaceClassification,
    IncompleteFacePayload,
    InvalidLocalCell,
    InvalidPublishedNode,
    DuplicateCanonicalFacet,
    FacetOrdinalOutOfRange,
    InvalidChartCertificate,
    ScratchCapacityReached,
    WorkBudgetReached,
    WorkAccountingOverflow,
    AllocationFailed,
    UnexpectedFailure
};

// Caller-owned resolver scratch.  Production callers reserve this storage
// before entering the facet-conflict preparation checkpoint; resolution then
// performs no allocation.  Only the first requiredScratchRows entries are
// touched.
struct FacetConflictRankedFacetScratch {
    std::array<Index, 3> canonicalGlobalNodes{};
    std::size_t localCell = 0;
};

struct FacetConflictResolveResult {
    FacetConflictResolveStatus status = FacetConflictResolveStatus::Invalid;
    FacetConflictResolveReason reason =
        FacetConflictResolveReason::MissingRetainedArtifact;
    ArtifactStamp faceStamp{};
    int faceTag = -1;
    bnd::ChartId chart = 0;
    FacetBarycentreCertificateView chartCertificate{};
    std::array<FacetBarycentreVertex, 3> vertices{};
    // Exact resolver-local work: one unit per retained face row, local-cell
    // row and published-node/mapping row inspected. This is a counter, not a
    // byte/capacity gauge.
    std::uint64_t workUnits = 0;

    bool resolved() const noexcept {
        return status == FacetConflictResolveStatus::Resolved;
    }
};

// normalizedLocalCellOrdinal is the rank after sorting every retained Tri3 by
// its canonical published-global-node triple. It is never interpreted as an
// index into SurfaceFaceResult::localCells.
FacetConflictResolveResult resolveFacetConflict(
    const SurfaceSetResult& retained,
    const bnd::VolumeFacetKey& key) noexcept;

FacetConflictResolveResult resolveFacetConflict(
    const SurfaceSetResult& retained,
    const bnd::VolumeFacetKey& key,
    std::span<FacetConflictRankedFacetScratch> scratch,
    std::uint64_t maxWorkUnits) noexcept;

FacetConflictResolveResult resolveFacetConflict(
    const SurfaceSetOutcome& outcome,
    const bnd::VolumeFacetKey& key) noexcept;

FacetConflictResolveResult resolveFacetConflict(
    const SurfaceSetOutcome& outcome,
    const bnd::VolumeFacetKey& key,
    std::span<FacetConflictRankedFacetScratch> scratch,
    std::uint64_t maxWorkUnits) noexcept;

}  // namespace paacs::algo
