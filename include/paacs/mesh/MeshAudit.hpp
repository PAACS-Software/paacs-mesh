// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// paacs::mesh::MeshAudit - M5's mandatory external truth layer.
//
// Product invariants #2/#3/#10: no silent invalidity, auditable boundary
// conformity, claims require evidence. auditMesh() inspects a finished (or
// intermediate) mesh with no access to the geometry kernel and reports
// every violation as a structured diagnostic; callers gate on
// result.passed(maxAcceptableSeverity).
//
// TOPOLOGY IS CORNER-BASED: Tet10 participates in boundary extraction,
// manifoldness, components and region volumes through its four corner
// nodes (stored first by convention), and Tri6 prescribes its corner
// triangle - a Tet10-only mesh audits identically to its Tet4 shadow
// (review finding: it previously vanished from every topology check).
//
// Checks (M5-A scope):
//   * non-finite node coordinates
//   * inverted Tet4 / invalid Tet10 (sampled Jacobian, via Quality)
//   * FULL boundary comparison: prescribed Tri3 facets vs the tet
//     complex's actual boundary faces - matched, orphaned (prescribed but
//     absent) and unexpected (present but unprescribed) both counted
//   * face manifoldness (no tet face shared by more than two tets)
//   * unused nodes / point accounting
//   * classification sanity (cells and nodes carry plausible dims/tags)
//   * per-region (dim-3 tag) volume, optional comparison to expectations
//   * connected components of the tet complex (face-connected)
//   * reproducibility hash over full mesh content
//
// NOT yet here (M5-B / later): per-pipeline-stage enforcement hooks, edge
// manifoldness, Tri/Line boundary chains, geometry-aware checks.
//
#include <cstdint>
#include <map>

#include "paacs/diag/Diagnostics.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/topo/ElementTraits.hpp"

namespace paacs::mesh {

// Boundary-comparison engagement policy (review finding: comparing only
// "when both populations exist" let a missing entire surface or volume
// population masquerade as a clean boundary audit).
enum class BoundaryComparisonMode {
    Disabled,             // no boundary diagnostics at all
    WhenSurfacePresent,   // default: engaged iff prescribed facets exist;
                          //   prescribed + no volume => ALL facets orphaned
    Required,             // always engaged: volume + no prescribed => ALL
                          //   hull faces unexpected; empty mesh => failure
};

struct AuditOptions {
    BoundaryComparisonMode boundaryMode =
        BoundaryComparisonMode::WhenSurfacePresent;
    // Known-gap allowance for prescribed-boundary orphans (PM-BND-001):
    // orphans within the budget are Warnings citing the open issue;
    // beyond it, RecoverableFailure. Constrained recovery (M10) sets this
    // to zero everywhere.
    std::uint64_t boundaryOrphanBudget = 0;
    // SEPARATE budget for unexpected hull facets (review: missing
    // prescribed facets and unexpected hull facets are different defects
    // and must not share a policy).
    std::uint64_t boundaryUnexpectedBudget = 0;
    // Expected face-connected component count of the tet complex; 0 = no
    // expectation. A single solid silently splitting in two must not pass.
    int expectedComponents = 0;
    bool unusedNodesAreErrors = false;
    // Strict export gates may escalate duplicate prescribed facets from
    // Warning to RecoverableFailure.
    bool duplicateFacetsAreErrors = false;
    // Expected |volume| per dim-3 region tag; empty = no expectation.
    std::map<int, double> expectedRegionVolumes;
    double regionVolumeRelTol = 1e-9;
};

struct AuditResult {
    diag::Report report;

    Index nNodes = 0, nTet4 = 0, nTet10 = 0, nTri3 = 0, nTri6 = 0;
    Index malformedCells = 0;      // invalid connectivity (defense-in-depth)
    // Cells whose TYPE the topology audit does not yet cover (Hex8,
    // Wedge6 volume cells; Quad4 surface facets): counted and diagnosed
    // loudly instead of silently passing (review: a one-Hex8 mesh
    // previously audited "clean" with 0 components and empty volumes).
    // Coverage arrives with shape-aware faces (M4-B.4) + generic audit
    // (M4-B.8).
    Index deferredAuditCells = 0;
    Index boundaryDuplicates = 0;  // prescribed facets sharing a key
    Index nInverted = 0, nNonFinite = 0, nInvalidTet10 = 0;
    Index unusedNodes = 0;
    Index boundaryPrescribed = 0, boundaryMatched = 0;
    Index boundaryOrphans = 0, boundaryUnexpected = 0;
    Index nonManifoldFaces = 0;
    int connectedComponents = 0;
    std::map<int, double> regionVolumes;   // dim-3 tag -> sum |V|
    std::uint64_t meshHash = 0;

    bool passed(diag::Severity maxAcceptable =
                    diag::Severity::Warning) const {
        return report.acceptable(maxAcceptable);
    }
};

// STRICT EXECUTION HASH (directive 5.1, stream V2 - structured): the
// stream begins with a domain+version marker, then explicit node/cell
// counts, section tags and per-cell connectivity lengths, so
// structurally different serializations cannot alias. Exact coordinate
// bit patterns,
// full entity/connectivity/classification ordering. Purpose: same-build
// regression, thread-count determinism, transaction/storage migration
// verification. May legitimately differ across compilers/platforms/FP
// modes.
std::uint64_t meshHash(const Mesh& m);

// CANONICAL SEMANTIC HASH v3 (directive 5.2, second review redesign).
// Physical: origin (bbox center) and scale (bbox diagonal) hashed
// explicitly at relative precision; node coordinates encoded LOCALLY
// relative to the decoded origin/scale on an absolute 2^-40 grid, so
// local detail survives arbitrary translation. Node-relabeling invariant
// via canonical relabeling; AMBIGUOUS relabelings (distinct nodes with
// identical quantized keys) are DETECTED and reported via
// SemanticHash::canonical=false, never silently input-order dependent.
// Cells take the lexicographic minimum over their orientation-preserving
// rotation group (topo::rotations) with high-order nodes re-slotted -
// generic across all ten trait types. HONEST LIMIT: quantization absorbs
// sub-quantum noise only in bucket INTERIORS; a value straddling a
// rounding boundary can flip on an epsilon change, so cross-platform
// semantic equality holds when computations agree well within a quantum
// (the exact-predicate pipeline guarantees this for topology; for
// coordinates it holds with high probability, not certainty). Recorded consideration: a
// tolerance-aware mesh COMPARISON may eventually complement or replace
// single-hash cross-platform equivalence (review note).
struct SemanticHash {
    std::uint64_t value = 0;
    bool canonical = true;   // false: hash is input-order dependent
};
SemanticHash semanticHashEx(const Mesh& m);
// THROWS std::logic_error when the relabeling is ambiguous - a
// non-canonical hash can never be consumed silently as canonical.
std::uint64_t semanticHash(const Mesh& m);

AuditResult auditMesh(const Mesh& m, const AuditOptions& opt = {});

} // namespace paacs::mesh
