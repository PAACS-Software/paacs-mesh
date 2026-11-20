// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// T2-B: bounded existing-node recovery of one protected Tet edge.
// The implementation is added with the accepted target-edge mode on the
// existing TetCavityTransaction; this payload is the closed public taxonomy.

#include "paacs/algo/SurfaceMesherAdapter.hpp"
#include "paacs/bnd/TetFacetGeometry.hpp"
#include "paacs/bnd/TetEdgeGeneralCavity.hpp"
#include "paacs/topo/Handles.hpp"

#include <cstdint>
#include <array>
#include <algorithm>
#include <limits>
#include <span>
#include <type_traits>
#include <vector>

namespace paacs::topo {
struct TetEdgeCavityPreview;
}

namespace paacs::bnd {

struct VolumeArtifactEdgeKey {
    std::array<Index, 2> canonicalGlobalNodes{kInvalidIndex, kInvalidIndex};
    static VolumeArtifactEdgeKey of(Index a, Index b) noexcept {
        if (b < a) std::swap(a, b);
        return {{a, b}};
    }
    friend bool operator==(const VolumeArtifactEdgeKey&,
                           const VolumeArtifactEdgeKey&) = default;
    friend bool operator<(const VolumeArtifactEdgeKey& a,
                          const VolumeArtifactEdgeKey& b) noexcept {
        return a.canonicalGlobalNodes < b.canonicalGlobalNodes;
    }
};

struct TetEdgeRecoveryRequest {
    topo::WorkspaceStamp stamp{};
    algo::ArtifactStamp artifactStamp{};
    VolumeArtifactEdgeKey key{};
    std::array<topo::NodeHandle, 2> vertices{};
    mesh::Classification classification{};
    mesh::Classification region{};
};

enum class TetEdgeRecoveryStatus : std::uint8_t {
    AlreadyPresent,
    Marked,
    Recovered,
    Refused,
    Rejected,
    Failed
};

enum class TetEdgeRecoveryPhase : std::uint8_t {
    None,
    Preflight,
    Routing,
    Discovery,
    Candidate,
    Preparation,
    StateCopy
};

enum class TetEdgeRecoveryReason : std::uint8_t {
    None,
    StaleWorkspaceStamp,
    StaleRecoveryContext,
    InvalidIdentity,
    InvalidClassification,
    NonfiniteCoordinate,
    DegenerateTarget,
    UnsupportedCell,
    NonPositiveTet,
    MixedRegion,
    StateKeyMismatch,
    StaleBinding,
    WrongBindingEdge,
    WrongBindingClassification,
    IdRebinding,
    EdgeAlias,
    ProtectedConstraintConflict,
    BoundaryStar,
    NonManifoldStar,
    DisconnectedStar,
    MixedRegionStar,
    TangentialTouch,
    CoplanarContact,
    IdentityDistinctContact,
    NonCycleLink,
    NonBallClosure,
    InvalidCandidate,
    PreparationConflict,
    AllocationFailed,
    PendingConstraint,
    RecoveryConflictCycle,
    ResourceBudgetReached
};

struct TetEdgeRecoveryOutcome {
    TetEdgeRecoveryStatus status = TetEdgeRecoveryStatus::Failed;
    TetEdgeRecoveryReason reason = TetEdgeRecoveryReason::None;
    TetEdgeRecoveryPhase phase = TetEdgeRecoveryPhase::None;
    topo::EdgeHandle edge{};
    TetFacetGeometryReason geometryReason = TetFacetGeometryReason::None;
    TetFacetGeometryCounters counters{};
    bool capacityResiduePermitted = false;
    bool stagedNodeCreated = false;

    bool succeeded() const noexcept {
        return status == TetEdgeRecoveryStatus::AlreadyPresent ||
               status == TetEdgeRecoveryStatus::Marked ||
               status == TetEdgeRecoveryStatus::Recovered;
    }
};

struct TetEdgeIndexedOutcome {
    TetEdgeRecoveryOutcome recovery{};
    std::vector<topo::CellHandle> retiredCells;
    std::vector<topo::CellHandle> createdCells;
};

// Evidence for the optional plain-value solve; no topology or state is owned
// here. Limits bound this selected cavity's enumeration/search/validation,
// not import, the initial segment discovery, or an entire meshing job.
struct TetEdgeGeneralRecoveryDiagnostics {
    bool attempted = false;
    TetEdgeGeneralCavityStatus status = TetEdgeGeneralCavityStatus::InvalidInput;
    TetFacetGeometryReason geometryReason = TetFacetGeometryReason::InvalidInput;
    TetEdgeGeneralCavityMetrics metrics{};
    std::uint64_t adjacentExpansionAttempts = 0;
    std::uint64_t adjacentCellsAdded = 0;
    // Boundary preparation proposals, including refused replacements.
    std::uint64_t manifoldClosurePasses = 0;
    std::uint64_t manifoldCellsAdded = 0;
    std::uint64_t steinerCandidateAttempts = 0;
    bool stagedCandidateAvailable = false;
};

// Additive opt-in accounting for topology reads and persistent-state copies
// performed inside indexed edge/facet routing.  The accepted v3 callers pass
// null and therefore retain their exact existing result and execution path.
// Counts are completed operations: a FaceData copy contributes one row plus
// the exact number of copied incidence entries.
struct TetRecoveryBindingWorkCounters {
    std::uint64_t edgeBindingAuditCalls = 0;
    std::uint64_t edgeBindingRowsVisited = 0;
    std::uint64_t edgeAliasRowsVisited = 0;
    std::uint64_t facetBindingAuditCalls = 0;
    std::uint64_t facetBindingRowsVisited = 0;
    std::uint64_t facetBindingRequestRowsVisited = 0;
    std::uint64_t facetBindingDuplicateRowsCompared = 0;
    std::uint64_t facetBindingReceiptBuilds = 0;
    std::uint64_t facetBindingReceiptReuses = 0;
    std::uint64_t facetBindingReceiptRowsCopied = 0;
    std::uint64_t facetBindingReceiptLookups = 0;
    std::uint64_t facetBindingReceiptReserveCalls = 0;
    std::uint64_t edgeDataCopyRows = 0;
    std::uint64_t faceDataCopyRows = 0;
    std::uint64_t faceIncidenceRowsCopied = 0;
    std::uint64_t stateCopies = 0;
    std::uint64_t stateCopyRows = 0;
    bool complete = true;

    void reset() noexcept { *this = {}; }
    friend bool operator==(const TetRecoveryBindingWorkCounters&,
                           const TetRecoveryBindingWorkCounters&) = default;
};

#if defined(PAACS_TOPO_TESTING)
// Non-mutating T2-B.1 discovery evidence. This testing-only snapshot records
// the exact cavity selected by the accepted edge route; it does not alter
// routing, candidate construction, preparation or persistent state.
struct TetEdgeBoundaryFaceDiagnostic {
    std::array<std::uint32_t,3> denseNodes{};
    std::array<std::uint32_t,3> directedDenseNodes{};
    std::uint8_t targetEndpointCount = 0;
    std::uint8_t currentOwnerCount = 0;
    bool protectedFace = false;
    bool regionBoundary = false;
};

struct TetEdgeCavityVertexDiagnostic {
    std::uint32_t denseVertex = (std::numeric_limits<std::uint32_t>::max)();
    std::uint64_t sequence = 0;
    Vec3 position{};
};

struct TetEdgeCavityCellDiagnostic {
    topo::CellHandle cell{};
    std::uint64_t sequence = 0;
    std::array<std::uint32_t,4> denseNodes{};
};

struct TetEdgeRingEdgeDiagnostic {
    std::array<std::uint32_t,2> denseNodes{};
    std::size_t multiplicity = 0;
};

struct TetEdgeLinkVertexDiagnostic {
    std::uint32_t denseVertex = (std::numeric_limits<std::uint32_t>::max)();
    std::size_t degree = 0;
    std::size_t component = 0;
};

enum class TetEdgeTransactionDiagnosticCheckpoint : std::uint8_t {
    None,
    InputCaptured,
    StageEntered,
    StageCompleted,
    PrepareEntered,
    PrepareCompleted
};

struct TetEdgeCavityDiagnostics {
    bool enabled = false;
    bool complete = false;
    bool allocationFailed = false;
    VolumeArtifactEdgeKey targetGlobal{};
    std::array<std::uint32_t,2> targetDense{};
    std::vector<topo::CellHandle> initiallySelected;
    std::vector<topo::CellHandle> closureAdded;
    std::vector<topo::CellHandle> selected;
    std::vector<TetEdgeCavityVertexDiagnostic> initialVertices;
    std::vector<TetEdgeCavityCellDiagnostic> initialCells;
    std::array<std::size_t,3> boundaryFacesByTargetEndpoints{};
    std::vector<TetEdgeBoundaryFaceDiagnostic> boundaryFaces;
    std::vector<TetEdgeRingEdgeDiagnostic> ringEdges;
    std::vector<TetEdgeLinkVertexDiagnostic> linkVertices;
    std::vector<std::array<std::uint32_t,2>> protectedContactEdges;
    std::vector<std::array<std::uint32_t,3>> protectedBoundaryFaces;
    std::vector<std::array<std::uint32_t,3>> regionBoundaryFaces;
    std::size_t linkComponents = 0;
    bool linkIsSingleCycle = false;
    bool ballValid = false;
    TetFacetGeometryReason ballReason = TetFacetGeometryReason::InvalidInput;
    std::vector<topo::CellHandle> completionAdded;
    std::vector<topo::CellHandle> completionSelected;
    std::array<std::size_t,3> completionBoundaryFacesByTargetEndpoints{};
    std::size_t completionPasses = 0;
    std::size_t completionWorkUnits = 0;
    std::size_t completionLinkComponents = 0;
    bool completionAvailable = false;
    bool completionBlockedByBoundary = false;
    bool completionBlockedByProtectedFace = false;
    bool completionBlockedByRegion = false;
    bool completionLinkIsSingleCycle = false;
    bool completionBallValid = false;
    TetFacetGeometryReason completionBallReason =
        TetFacetGeometryReason::InvalidInput;

};

// Allocation-free T2-C preparation-conflict evidence. This intentionally
// lives outside the vector-heavy cavity diagnostic so enabling it cannot
// trigger discovery snapshots or any move/allocation boundary.
struct TetEdgeTransactionDiagnostics {
    bool enabled = false;
    VolumeArtifactEdgeKey targetGlobal{};
    topo::WorkspaceStamp transactionRequestStamp{};
    topo::WorkspaceStamp transactionCurrentStamp{};
    mesh::Classification transactionClassification{};
    mesh::Classification transactionRegion{};
    bool transactionRequiredConstraint = false;
    std::size_t transactionRemovedCount = 0;
    std::size_t transactionCandidateCount = 0;
    std::uint64_t transactionRemovedDigest = 0;
    std::uint64_t transactionCandidateDigest = 0;
    TetEdgeTransactionDiagnosticCheckpoint transactionCheckpoint =
        TetEdgeTransactionDiagnosticCheckpoint::None;
};

void resetTetEdgeCavityDiagnosticsForTesting(
    bool enabled, VolumeArtifactEdgeKey target = {}) noexcept;
const TetEdgeCavityDiagnostics& tetEdgeCavityDiagnosticsForTesting() noexcept;
void resetTetEdgeTransactionDiagnosticsForTesting(
    bool enabled,VolumeArtifactEdgeKey target={}) noexcept;
TetEdgeTransactionDiagnostics tetEdgeTransactionDiagnosticsForTesting()
    noexcept;
#endif

static_assert(std::is_standard_layout_v<TetEdgeRecoveryOutcome>);
static_assert(std::is_trivially_copyable_v<TetEdgeRecoveryOutcome>);
static_assert(std::is_nothrow_move_constructible_v<TetEdgeRecoveryOutcome>);
static_assert(std::is_standard_layout_v<VolumeArtifactEdgeKey>);
static_assert(std::is_trivially_copyable_v<VolumeArtifactEdgeKey>);

struct TetFacetSetState;
struct TetRecoveryContext;

#if defined(PAACS_TOPO_TESTING)
// Testing-only exact prepared-transaction guard.  The callback runs after
// stage()+prepare() have completed for the one authoritative transaction and
// before the fallible metadata-state copy that precedes commitEdge().  The
// workspace is deliberately not exposed: it is unreadable while prepared.
// All topology inputs are the transaction-owned plain preview plus the
// stamp/context captured before staging.  Defer abandons that same prepared
// transaction; Allow proceeds through the state copy and commits it.
enum class TetEdgePreparedGuardDecisionForTesting : std::uint8_t {
    Allow,
    Defer
};

using TetEdgePreparedGuardForTesting =
    TetEdgePreparedGuardDecisionForTesting (*)(
        const TetEdgeRecoveryRequest&,topo::WorkspaceStamp,
        topo::EdgeHandle,const topo::TetEdgeCavityPreview&,
        const TetRecoveryContext&,void*) noexcept;
#endif

TetEdgeRecoveryOutcome recoverTetEdgeIndexed(
    const TetEdgeRecoveryRequest& request, TetFacetSetState& state,
    topo::TopologyWorkspace& workspace, TetRecoveryContext& context) noexcept;

// Pure current-stamp assessment used by the deterministic set scheduler.
// It performs routing/discovery/candidate validation but never prepares or
// commits topology and never changes state/context.
TetEdgeRecoveryOutcome assessTetEdgeIndexed(
    const TetEdgeRecoveryRequest& request, const TetFacetSetState& state,
    const topo::TopologyWorkspace& workspace,
    TetRecoveryContext& context,
    TetRecoveryBindingWorkCounters* bindingWork = nullptr) noexcept;

TetEdgeIndexedOutcome assessTetEdgeIndexedDetailed(
    const TetEdgeRecoveryRequest& request, const TetFacetSetState& state,
    const topo::TopologyWorkspace& workspace,
    TetRecoveryContext& context,
    TetRecoveryBindingWorkCounters* bindingWork = nullptr) noexcept;

// V4-only non-allocating discovery witness. The accepted detailed overload
// above retains its existing allocation behavior and result shape.
TetEdgeRecoveryOutcome assessTetEdgeIndexedFrontier(
    const TetEdgeRecoveryRequest& request, const TetFacetSetState& state,
    const topo::TopologyWorkspace& workspace, TetRecoveryContext& context,
    std::span<topo::CellHandle> selectedCells,
    std::size_t& selectedCellsUsed, bool& witnessComplete,
    TetRecoveryBindingWorkCounters* bindingWork = nullptr) noexcept;

TetEdgeIndexedOutcome recoverTetEdgeIndexedDetailed(
    const TetEdgeRecoveryRequest& request, TetFacetSetState& state,
    topo::TopologyWorkspace& workspace, TetRecoveryContext& context,
    TetRecoveryBindingWorkCounters* bindingWork = nullptr) noexcept;

// Opt-in existing-node general-cavity fallback on an exact NonCycleLink stop.
// Every limit must be positive. The original overloads retain their routing.
// Solved candidates use the same atomic TetCavityTransaction and bindings;
// refusal, budget exhaustion and preparation failure publish no candidate.
TetEdgeIndexedOutcome assessTetEdgeIndexedDetailed(
    const TetEdgeRecoveryRequest& request, const TetFacetSetState& state,
    const topo::TopologyWorkspace& workspace, TetRecoveryContext& context,
    const TetEdgeGeneralCavityLimits& generalLimits,
    TetEdgeGeneralRecoveryDiagnostics* diagnostics = nullptr,
    TetRecoveryBindingWorkCounters* bindingWork = nullptr,
    std::uint64_t maximumAdjacentExpansionAttempts = 0) noexcept;

TetEdgeIndexedOutcome recoverTetEdgeIndexedDetailed(
    const TetEdgeRecoveryRequest& request, TetFacetSetState& state,
    topo::TopologyWorkspace& workspace, TetRecoveryContext& context,
    const TetEdgeGeneralCavityLimits& generalLimits,
    TetEdgeGeneralRecoveryDiagnostics* diagnostics = nullptr,
    TetRecoveryBindingWorkCounters* bindingWork = nullptr,
    std::uint64_t maximumAdjacentExpansionAttempts = 0) noexcept;

TetEdgeRecoveryOutcome assessTetEdgeIndexedFrontier(
    const TetEdgeRecoveryRequest& request, const TetFacetSetState& state,
    const topo::TopologyWorkspace& workspace, TetRecoveryContext& context,
    std::span<topo::CellHandle> selectedCells,
    std::size_t& selectedCellsUsed, bool& witnessComplete,
    const TetEdgeGeneralCavityLimits& generalLimits,
    TetEdgeGeneralRecoveryDiagnostics* diagnostics = nullptr,
    TetRecoveryBindingWorkCounters* bindingWork = nullptr,
    std::uint64_t maximumAdjacentExpansionAttempts = 0) noexcept;

#if defined(PAACS_TOPO_TESTING)
// Executes the ordinary indexed recovery, differing only at the exact
// prepared/commit seam documented above.  Null is mechanically identical to
// recoverTetEdgeIndexedDetailed().
TetEdgeIndexedOutcome recoverTetEdgeIndexedPreparedGuardedForTesting(
    const TetEdgeRecoveryRequest& request,TetFacetSetState& state,
    topo::TopologyWorkspace& workspace,TetRecoveryContext& context,
    TetEdgePreparedGuardForTesting guard,void* guardState,
    TetRecoveryBindingWorkCounters* bindingWork=nullptr) noexcept;

TetEdgeIndexedOutcome recoverTetEdgeIndexedPreparedGuardedForTesting(
    const TetEdgeRecoveryRequest& request,TetFacetSetState& state,
    topo::TopologyWorkspace& workspace,TetRecoveryContext& context,
    TetEdgePreparedGuardForTesting guard,void* guardState,
    const TetEdgeGeneralCavityLimits& generalLimits,
    TetEdgeGeneralRecoveryDiagnostics* diagnostics=nullptr,
    TetRecoveryBindingWorkCounters* bindingWork=nullptr,
    std::uint64_t maximumAdjacentExpansionAttempts=0) noexcept;
#endif

}  // namespace paacs::bnd
