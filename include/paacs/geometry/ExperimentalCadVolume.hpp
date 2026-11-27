// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0
#pragma once

#include "paacs/geometry/OcctMesher.hpp"
#include <span>

namespace paacs::geo {

// Experimental, optional-OCCT, single-solid uniform Tet4 construction only.
// Classical selects ordinary initial placement. Supplied replaces that stage;
// an empty span deliberately requests boundary-only construction.
enum class ExperimentalInteriorMode : std::uint8_t { Classical, Supplied };

struct ExperimentalCadBoundary {
    mesh::Mesh boundary;
    int solidTag = -1;
    double targetSize = 0, maximumChordError = 0, maximumRelativeVolumeError = .01;
    // Diagnostic import-instance token, never portable CAD identity authority.
    std::uint64_t modelIdentity = 0;
};

struct ExperimentalCadBoundaryOutcome {
    bool prepared = false, destinationUnchanged = true, rebound = false;
    NativeCadMeshReason reason = NativeCadMeshReason::InternalError;
    const char* diagnostic = "not_started";
    std::uint64_t work = 0, nodes = 0, triangles = 0;
    double seconds = 0, maximumNodeError = 0, maximumFacetError = 0;
    NativeCadSourceShellReport sourceShell;
};

struct ExperimentalCadVolumeReport {
    const char* diagnostic = "not_started";
    ExperimentalInteriorMode mode = ExperimentalInteriorMode::Classical;
    bool boundaryRevalidated = false, boundaryUnchanged = false;
    bool allSuppliedPointsRetained = false;
    bool initialPointsAvailable = false;
    std::uint64_t suppliedPoints = 0, admittedPoints = 0, retainedPoints = 0;
    std::uint64_t sourceWork = 0, placementWork = 0, recoveryWork = 0;
    double sourceSeconds = 0, placementSeconds = 0, recoverySeconds = 0;
    // Measured named spans, not a decomposition of totalSeconds. Placement
    // includes its auxiliary numerical construction; recovery is initial
    // constraint recovery; fill includes material construction/postprocessing.
    // Validation combines front CAD checks and final assembly checks. Final
    // global audit/quality times remain separately in NativeCadMeshOutcome.
    // Unmeasured setup and initial-topology finalization remain in total time.
    double fillSeconds = 0, validationSeconds = 0, totalSeconds = 0;
    // Exact staged initial coordinates immediately before constraint recovery,
    // excluding source-boundary nodes. Bounded by the same request resources.
    std::vector<Vec3> initialPoints;
};

struct ExperimentalCadVolumeOutcome {
    NativeCadMeshOutcome native;
    ExperimentalCadVolumeReport experimental;
};

// A bounded placement checkpoint, not a qualified volume mesh. Its owned
// ordered cloud excludes boundary nodes; a generated empty cloud is distinct
// from an unavailable/refused checkpoint. Placement includes its auxiliary
// tetrahedral construction and topology finalization, but no recovery/fill.
struct ExperimentalCadInitialPointsOutcome {
    bool generated = false;
    NativeCadMeshReason reason = NativeCadMeshReason::InternalError;
    const char* diagnostic = "not_started";
    std::uint64_t work = 0;
    double seconds = 0;
    ExperimentalCadVolumeReport stages;
    NativeCadSourceShellReport sourceShell;
    NativeCadInitialFrontReport initialFront;
    NativeCadInitialTopologyQualityReport initialTopologyQuality;
};

// Preparation publishes only an exact-certified tagged Tri3 shell, never a
// volume mesh. Every failure preserves the caller's previous snapshot.
ExperimentalCadBoundaryOutcome prepareExperimentalCadBoundary(
    const OcctModel&, const NativeCadMeshOptions&, ExperimentalCadBoundary&) noexcept;

// Untrusted archive data must re-enter through this function after import.
// It rechecks original CAD node/face incidence, geometric error and exact shell
// embedding. A saved modelIdentity/node-key ordinal is never sufficient.
ExperimentalCadBoundaryOutcome rebindExperimentalCadBoundary(
    const OcctModel&, const mesh::Mesh& boundary, const NativeCadMeshOptions&,
    ExperimentalCadBoundary&) noexcept;

// Run the ordinary initial placement against the exact revalidated snapshot
// and stop before constraint recovery/material fill. The request uses the same
// frozen uniform sizing/fidelity and explicit resource contract as replay.
// A generated cloud does not certify any volume connectivity or solver quality.
// Failure never publishes a partial cloud or changes the prepared boundary.
ExperimentalCadInitialPointsOutcome generateExperimentalCadInitialPoints(
    const OcctModel&, const ExperimentalCadBoundary&, const NativeCadMeshOptions&) noexcept;

// Uniform sizing/fidelity options must match the prepared snapshot. Supplied
// sites are admitted all-or-nothing against original material, without an
// automatic placement/boundary-only fallback. Every supplied coordinate must
// remain fixed and participate in a published Tet4. Recovery may add points.
// Exact tagged boundary equality and all ordinary native publication gates
// remain mandatory. Failure never modifies destination.
ExperimentalCadVolumeOutcome replayExperimentalCadVolume(
    const OcctModel&, const ExperimentalCadBoundary&, const NativeCadMeshOptions&,
    ExperimentalInteriorMode, std::span<const Vec3> initialPoints,
    mesh::Mesh& destination) noexcept;

} // namespace paacs::geo
