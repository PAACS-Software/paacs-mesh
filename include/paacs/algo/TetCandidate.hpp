// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// CK-3 Increment 1: deterministic, audited Tet candidate and private
// authoritative-topology handoff. The existing delaunay3D predicate and tie
// policy remain unchanged; this layer canonicalizes its finite output,
// validates the complete Tet complex, and imports a selected subset without
// mutating caller-owned mesh state.

#include <array>
#include <atomic>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <type_traits>
#include <vector>

#include "paacs/algo/Delaunay3D.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/topo/Handles.hpp"

namespace paacs::topo {
class TopologyWorkspace;
}

namespace paacs::algo {

enum class TetCandidateStatus : std::uint8_t { Built, Rejected, Failed };

enum class TetCandidateReason : std::uint8_t {
    None,
    NonfiniteInput,
    NonVolumetricInput,
    InvalidPointMapping,
    MissingRequiredPoint,
    InvalidCellIndex,
    RepeatedCellVertex,
    NonPositiveOrientation,
    DuplicateCell,
    NonManifoldFace,
    AdjacencyMismatch,
    DisconnectedComplex,
    BoundaryExtractionMismatch,
    BoundaryTopologyInvalid,
    NonCanonicalOrdering,
    InvalidImportSelection,
    InvalidGlobalMapping,
    InvalidClassification,
    WorkspaceImportFailed,
    AllocationFailed,
    ConstructionFailed,
    BudgetExceeded,
    Cancelled
};

struct TetCandidateCell {
    std::array<Index, 4> vertices{};
    std::array<Index, 4> neighbours{
        kInvalidIndex, kInvalidIndex, kInvalidIndex, kInvalidIndex};
    Index sourceOrdinal = kInvalidIndex;
};

struct TetCandidateBoundaryFace {
    std::array<Index, 3> vertices{};  // outward directed cycle
    Index owner = kInvalidIndex;
    std::uint8_t ownerLocalFace = 0;
};

struct TetCandidateMetrics {
    Delaunay3DStats construction{};
    std::uint64_t inputPoints = 0;
    std::uint64_t uniquePoints = 0;
    std::uint64_t candidateCells = 0;
    std::uint64_t interiorFaces = 0;
    std::uint64_t boundaryFaces = 0;
    std::uint64_t validationFaceUses = 0;
    std::uint64_t validationHullTests = 0;
    std::uint64_t candidateAccountedBytes = 0;
    // Conservative accounted construction envelope: final candidate plus
    // Delaunay arena plus raw returned-Tet capacity. The arena is destroyed
    // before final candidate canonicalization, so this sum is deliberately
    // not described as simultaneous live overlap. Excludes caller input.
    std::uint64_t constructionCandidateConservativeBytes = 0;
    bool identityComplete = false;
    bool adjacencyComplete = false;
    bool boundaryComplete = false;
    bool runComplete = false;
    TetCandidateReason initialCoordinateFailure = TetCandidateReason::None;
    bool originalCoordinateRetryUsed = false;
    // Known canonicalization/validation visits reported by a discarded
    // candidate attempt. Its Delaunay work is retained separately above.
    std::uint64_t discardedCandidateWork = 0;
    // Present only for the additive bounded constructor. Work spans every
    // discarded/retry attempt; payload is a conservative modeled peak.
    std::uint64_t boundedWork = 0;
    std::uint64_t boundedPeakTemporaryBytes = 0;
};

static_assert(std::is_standard_layout_v<TetCandidateMetrics>);
static_assert(std::is_trivially_copyable_v<TetCandidateMetrics>);

struct TetCandidate {
    // Canonical lexicographic point order. pointToInput retains the first
    // caller identity for each unique coordinate; inputToPoint maps every
    // caller identity (including duplicates) back to the canonical point.
    std::vector<Vec3> points;
    std::vector<Index> pointToInput;
    std::vector<Index> inputToPoint;
    std::vector<TetCandidateCell> cells;
    std::vector<TetCandidateBoundaryFace> boundaryFaces;
    TetCandidateMetrics metrics{};
};

struct TetCandidateValidation {
    bool valid = false;
    TetCandidateReason reason = TetCandidateReason::ConstructionFailed;
    Index entity = kInvalidIndex;
};

struct TetCandidateBuildOutcome {
    TetCandidateStatus status = TetCandidateStatus::Failed;
    TetCandidateReason reason = TetCandidateReason::ConstructionFailed;
    TetCandidateMetrics metrics{};
    std::unique_ptr<TetCandidate> candidate;

    bool built() const noexcept {
        return status == TetCandidateStatus::Built && candidate != nullptr;
    }
};

static_assert(std::is_nothrow_move_constructible_v<TetCandidateBuildOutcome>);

// Complete candidate validation. Allocation failures are reported, never
// escaped. The candidate is not mutated.
TetCandidateValidation validateTetCandidate(
    const TetCandidate& candidate) noexcept;

// Calls the accepted exact delaunay3D constructor, then canonicalizes and
// validates its complete finite output. No caller-owned state is mutated.
TetCandidateBuildOutcome buildTetCandidate(
    const std::vector<Vec3>& points) noexcept;

struct TetCandidateBuildLimits {
    std::uint64_t maxWork = (std::numeric_limits<std::uint64_t>::max)();
    // Live finite construction cells and final candidate cells. Symbolic hull
    // cells are admitted against temporary payload, not this finite-cell cap.
    std::uint64_t maxCells = (std::numeric_limits<std::uint64_t>::max)();
    // Includes a fixed 4 MiB exact-predicate/workspace envelope, vector
    // growth overlap, raw output, canonicalization and validation payload.
    std::uint64_t maxTemporaryBytes = (std::numeric_limits<std::uint64_t>::max)();
    const std::atomic_bool* cancellation = nullptr;
};

// Serial, online-bounded construction. Stop outcomes never publish a partial
// candidate. The caller's input storage is excluded from the payload model.
TetCandidateBuildOutcome buildTetCandidate(
    const std::vector<Vec3>& points, const TetCandidateBuildLimits& limits) noexcept;

std::uint64_t tetCandidateAccountedBytes(
    const TetCandidate& candidate) noexcept;

struct TetWorkspaceNodeBinding {
    Index candidatePoint = kInvalidIndex;
    Index existingGlobalIndex = kInvalidIndex;
    topo::NodeHandle local{};
};

struct TetWorkspaceCellBinding {
    Index candidateCell = kInvalidIndex;
    topo::CellHandle local{};
};

struct TetHandoffMetrics {
    std::uint64_t selectedCells = 0;
    std::uint64_t importedNodes = 0;
    std::uint64_t importedCells = 0;
    std::uint64_t candidateAccountedBytes = 0;
    std::uint64_t constructionCandidateConservativeAccountedBytes = 0;
    // Modeled StorageStats/vector payload categories, not allocator-exact
    // heap usage. In particular, scratch excludes implementation-specific
    // std::set node overhead; its operation traffic is counted below.
    std::uint64_t authoritativeModeledStorageBytes = 0;
    std::uint64_t importScratchModeledPayloadBytes = 0;
    // Complete filtered connectivity exists only before authoritative import.
    std::uint64_t filterConnectivityStageAccountedBytes = 0;
    // Compact selected-cell and point/global/classification maps coexist with
    // candidate and authoritative workspace during import.
    std::uint64_t callerMappingAccountedBytes = 0;
    std::uint64_t publishedPayloadAccountedBytes = 0;
    // Sum of the disclosed modeled categories at their import overlap. This
    // is not described as a complete process heap high-water measurement.
    std::uint64_t modeledOverlapPayloadBytes = 0;
    std::uint64_t completeTetRepresentationsAtOverlap = 0;
    std::uint64_t completeTetRepresentationsBeforePublication = 0;
    std::uint64_t importEdgeValidationCellVisits = 0;
    std::uint64_t importEdgeBuildCellVisits = 0;
    std::uint64_t importEdgeUses = 0;
    std::uint64_t importUniqueEdges = 0;
    std::uint64_t importFaceValidationCellVisits = 0;
    std::uint64_t importFaceBuildCellVisits = 0;
    std::uint64_t importFaceUses = 0;
    std::uint64_t importUniqueFaces = 0;
    std::uint64_t importEdgeAudits = 0;
    std::uint64_t importFaceAudits = 0;
    std::uint64_t importPerCellGlobalAudits = 0;
    // The import deliberately revalidates the caller-owned candidate. These
    // counters charge its complete semantic row/predicate traversal again;
    // they are not carried from the earlier candidate-construction metrics.
    std::uint64_t importCandidateValidationCalls = 0;
    std::uint64_t importCandidateValidationPointRows = 0;
    std::uint64_t importCandidateValidationCellRows = 0;
    std::uint64_t importCandidateValidationFaceUses = 0;
    std::uint64_t importCandidateValidationFaceRows = 0;
    std::uint64_t importCandidateValidationBoundaryRows = 0;
    std::uint64_t importCandidateValidationDualNeighbourVisits = 0;
    std::uint64_t importCandidateValidationBoundaryEdgeUses = 0;
    std::uint64_t importCandidateValidationBoundaryNeighbourVisits = 0;
    std::uint64_t importCandidateValidationHullTests = 0;
    std::uint64_t importCandidateValidationSlotsInitialized = 0;
    std::uint64_t importSelectionRowsCopied = 0;
    std::uint64_t importSelectionSortComparisons = 0;
    std::uint64_t importSelectionDuplicateComparisons = 0;
    std::uint64_t importSelectionCellVisits = 0;
    std::uint64_t importSelectionVertexVisits = 0;
    std::uint64_t importSelectionPointSlotsInitialized = 0;
    std::uint64_t importGlobalMapRows = 0;
    std::uint64_t importGlobalMapComparisons = 0;
    std::uint64_t importGlobalMapInsertions = 0;
    std::uint64_t importLocalPointSlotsInitialized = 0;
    std::uint64_t importLocalGlobalSlotsInitialized = 0;
    std::uint64_t importLocalCellSlotsInitialized = 0;
    std::uint64_t importPersistentNodeSlotsInitialized = 0;
    // Exact-present v4 import preflight and resolution work. These remain
    // zero for the established ordinary and complete-boundary import paths.
    std::uint64_t partialImportRowsCopied = 0;
    std::uint64_t partialImportGlobalSlotsInitialized = 0;
    std::uint64_t partialImportCandidateMappingVisits = 0;
    std::uint64_t partialImportEdgeValidationRows = 0;
    std::uint64_t partialImportFaceValidationRows = 0;
    std::uint64_t partialImportProtectedEndpointVisits = 0;
    std::uint64_t partialImportAuthorityPointValidations = 0;
    std::uint64_t partialImportDuplicateKeyComparisons = 0;
    std::uint64_t partialImportResolutionNodeVisits = 0;
    std::uint64_t partialImportEdgeResolutionLookups = 0;
    std::uint64_t partialImportFaceResolutionLookups = 0;
    std::uint64_t partialImportEdgeDataReads = 0;
    std::uint64_t partialImportFaceDataCopies = 0;
    std::uint64_t partialImportFaceIncidenceRowsCopied = 0;
    std::uint64_t partialImportSortComparisons = 0;
    bool workspaceValid = false;
    bool candidateReleasedBeforePublication = false;
    bool filterRepresentationReleasedBeforePublication = false;
    bool workAccountingComplete = false;
    bool runComplete = false;
    friend bool operator==(const TetHandoffMetrics&,
                           const TetHandoffMetrics&) = default;
};

static_assert(std::is_standard_layout_v<TetHandoffMetrics>);
static_assert(std::is_trivially_copyable_v<TetHandoffMetrics>);

struct TetAuthoritativeResult {
    TetAuthoritativeResult();
    ~TetAuthoritativeResult();
    TetAuthoritativeResult(const TetAuthoritativeResult&) = delete;
    TetAuthoritativeResult& operator=(const TetAuthoritativeResult&) = delete;
    TetAuthoritativeResult(TetAuthoritativeResult&&) noexcept;
    TetAuthoritativeResult& operator=(TetAuthoritativeResult&&) noexcept;

    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::vector<TetWorkspaceNodeBinding> nodeBindings;
    // Indexed by canonical candidate-cell id; unselected entries are stale
    // default handles. selectedCellBindings gives deterministic live order.
    std::vector<topo::CellHandle> localCellByCandidate;
    std::vector<TetWorkspaceCellBinding> selectedCellBindings;
    std::vector<Index> existingGlobalByLocalNodeSlot;
    TetHandoffMetrics metrics{};
};

struct TetHandoffOutcome {
    TetCandidateStatus status = TetCandidateStatus::Failed;
    TetCandidateReason reason = TetCandidateReason::WorkspaceImportFailed;
    std::unique_ptr<TetAuthoritativeResult> result;

    bool imported() const noexcept {
        return status == TetCandidateStatus::Built && result != nullptr &&
               result->workspace != nullptr;
    }
};

struct TetImportProtectedEdge {
    std::array<Index,2> existingGlobalNodes{kInvalidIndex,kInvalidIndex};
    mesh::Classification classification{};
};

struct TetImportProtectedFace {
    std::array<Index,3> existingGlobalNodes{
        kInvalidIndex,kInvalidIndex,kInvalidIndex};
    mesh::Classification classification{};
};
static_assert(std::is_standard_layout_v<TetImportProtectedEdge>);
static_assert(std::is_trivially_copyable_v<TetImportProtectedEdge>);
static_assert(std::is_standard_layout_v<TetImportProtectedFace>);
static_assert(std::is_trivially_copyable_v<TetImportProtectedFace>);

static_assert(std::is_nothrow_move_constructible_v<TetHandoffOutcome>);

// Import the selected canonical candidate-cell ids into a new private
// authoritative workspace. selectedCandidateCells must be nonempty.
// existingGlobalByPoint and expectedClassByPoint are indexed by
// candidate.points. Every selected point must map injectively to a current
// globalMesh node whose exact coordinate and classification match the supplied
// candidate/expected values. Imported classification is always derived from
// globalMesh. The source candidate and global mesh remain observably unchanged
// on every outcome.
TetHandoffOutcome importTetCandidate(
    const TetCandidate& candidate,
    std::span<const Index> selectedCandidateCells,
    std::span<const Index> existingGlobalByPoint,
    std::span<const mesh::Classification> expectedClassByPoint,
    const mesh::Mesh& globalMesh,
    mesh::Classification cellClassification) noexcept;

// Additive artifact-aware handoff used by boundary-conforming-v3. The
// original import entry above is unchanged. Required edge/face metadata is
// resolved through the already-validated global identity mapping and applied
// while the destination pools are still private; no topology transaction or
// post-construction mutation is involved.
TetHandoffOutcome importTetCandidateWithProtectedBoundary(
    const TetCandidate& candidate,
    std::span<const Index> selectedCandidateCells,
    std::span<const Index> existingGlobalByPoint,
    std::span<const mesh::Classification> expectedClassByPoint,
    const mesh::Mesh& globalMesh,
    mesh::Classification cellClassification,
    std::span<const TetImportProtectedEdge> protectedEdges,
    std::span<const TetImportProtectedFace> protectedFaces) noexcept;

}  // namespace paacs::algo
