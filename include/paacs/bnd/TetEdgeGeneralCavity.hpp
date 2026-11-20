// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// T2-B.1: bounded plain-value existing-node tetrahedralization of one
// validated Tet edge cavity. No workspace or topology mutation is visible
// through this layer.

#include "paacs/bnd/TetFacetGeometry.hpp"

#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace paacs::bnd {

struct TetEdgeGeneralCavityLimits {
    std::uint64_t maxCavityCells = 0;
    std::uint64_t maxCavityVertices = 0;
    std::uint64_t maxEnumeratedTets = 0;
    std::uint64_t maxSearchStates = 0;
    std::uint64_t maxWorkUnits = 0;
    // Opt-in indexed recovery fallback after a certified existing-node
    // refusal. Existing and staged-node solves share maxWorkUnits.
    std::uint64_t maximumSteinerCandidates = 0;

    bool positive() const noexcept {
        return maxCavityCells != 0 && maxCavityVertices != 0 &&
               maxEnumeratedTets != 0 && maxSearchStates != 0 &&
               maxWorkUnits != 0;
    }
};

// Explicitly enabled recovery extension. Work is cumulative across exact
// general-cavity solves in one set-scheduler invocation (including repeated
// assessments); initial discovery and other meshing phases are outside it.
struct TetEdgeGeneralRecoveryPolicy {
    bool enabled = false;
    TetEdgeGeneralCavityLimits limits{64,16,5000,10000,500000};
    std::uint64_t maximumTotalWorkUnits = 5000000;
    // Canonical single/pair proposals from the original unprotected dual
    // one-ring. Zero disables expansion; no proposal adds more than two cells.
    std::uint64_t maximumAdjacentExpansionAttempts = 80;

    bool valid() const noexcept {
        return !enabled || (limits.positive() && maximumTotalWorkUnits != 0);
    }
};

enum class TetEdgeGeneralCavityStatus : std::uint8_t {
    Solved,
    Unavailable,
    BudgetReached,
    InvalidInput,
    InvalidCandidate,
    AllocationFailed
};

struct TetEdgeGeneralCavityMetrics {
    std::uint64_t cavityCells = 0;
    std::uint64_t cavityVertices = 0;
    std::uint64_t combinationCount = 0;
    std::uint64_t enumeratedTets = 0;
    std::uint64_t positiveTets = 0;
    std::uint64_t searchStates = 0;
    std::uint64_t faceLedgerVisits = 0;
    std::uint64_t pairContactTests = 0;
    std::uint64_t validatorVisits = 0;
    std::uint64_t workUnits = 0;
    std::uint64_t stackHighWater = 0;
    bool runComplete = false;
};

struct TetEdgeGeneralCavityInput {
    std::span<const TetFacetVertex> vertices;
    std::span<const TetFacetCell> removedCells;
    TetFacetEdgeKey target{};
    std::span<const TetFacetEdgeKey> protectedEdges;
    std::span<const TetFacetProtectedFace> protectedFaces;
    mesh::Classification region{3,0};
};

struct TetEdgeGeneralCavityOutcome {
    TetEdgeGeneralCavityStatus status =
        TetEdgeGeneralCavityStatus::InvalidInput;
    TetFacetGeometryReason geometryReason =
        TetFacetGeometryReason::InvalidInput;
    std::vector<std::array<TetFacetVertexId,4>> enumeratedTetKeys;
    std::vector<std::array<TetFacetVertexId,4>> orderedTetKeys;
    std::vector<TetFacetCell> cells;
    TetEdgeGeneralCavityMetrics metrics{};

    bool solved() const noexcept {
        return status == TetEdgeGeneralCavityStatus::Solved;
    }
};

TetEdgeGeneralCavityOutcome solveTetEdgeGeneralCavity(
    const TetEdgeGeneralCavityInput& input,
    const TetEdgeGeneralCavityLimits& limits) noexcept;

// Adapted from the retained V0.10 implementation. Coordinates are proposals;
// the exact general-cavity solver certifies the full replacement. The source
// is immutable and the result refers to its vertices plus one staged vertex.
struct TetEdgeOneSteinerOutcome {
    TetEdgeGeneralCavityStatus status = TetEdgeGeneralCavityStatus::InvalidInput;
    TetFacetGeometryReason geometryReason = TetFacetGeometryReason::InvalidInput;
    Vec3 steinerPosition{};
    std::uint64_t candidateAttempts = 0;
    std::uint64_t selectedCandidateOrdinal =
        (std::numeric_limits<std::uint64_t>::max)();
    std::uint64_t cumulativeWorkUnits = 0;
    TetEdgeGeneralCavityOutcome cavity{};
    bool solved() const noexcept { return status == TetEdgeGeneralCavityStatus::Solved; }
};

TetEdgeOneSteinerOutcome solveTetEdgeGeneralCavityOneSteiner(
    const TetEdgeGeneralCavityInput& input,
    const TetEdgeGeneralCavityLimits& perCandidateLimits,
    std::uint64_t maximumCandidateAttempts,
    std::uint64_t maximumCumulativeWorkUnits =
        (std::numeric_limits<std::uint64_t>::max)()) noexcept;

}  // namespace paacs::bnd
