// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// paacs::diag - M5 structured diagnostics.
//
// Product invariant #9: "Diagnostics are part of the API. Failure stages,
// locations, affected entities, attempted recoveries, and suggested
// remedies are machine-readable." This module is the vocabulary; MeshAudit
// (mesh/MeshAudit.hpp) is its principal producer, and pipeline stages
// adopt it incrementally (full per-stage enforcement is M5-B).
//
// Severity ladder (explicit, per review):
//   Info               - observation, no action required.
//   Warning            - suspicious but output remains usable.
//   UnsupportedInput   - the INPUT is outside the supported envelope; the
//                        engine declined rather than produced garbage.
//   RecoverableFailure - a stage failed in a way a caller or later
//                        milestone can remedy (e.g. boundary orphans that
//                        constrained recovery will eliminate).
//   InvalidOutput      - the produced data violates validity guarantees
//                        and MUST NOT be consumed by a solver.
//
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include "paacs/core/Types.hpp"
#include "paacs/mesh/Mesh.hpp"

namespace paacs::diag {

enum class Severity : int {
    Info = 0,
    Warning = 1,
    UnsupportedInput = 2,
    RecoverableFailure = 3,
    InvalidOutput = 4,
};

enum class Code : int {
    EmptyMesh,
    NonFiniteCoordinates,
    InvertedElement,
    InvalidTet10,
    NonManifoldFace,
    BoundaryOrphanFacet,
    BoundaryUnexpectedFacet,
    UnusedNodes,
    ClassificationInconsistent,
    RegionVolumeMismatch,
    MultipleComponents,
    ComponentCountMismatch,
    MalformedConnectivity,
    DuplicateBoundaryFacet,
    AuditCoverageDeferred,
    // M4-B.2 storage audit (the mandatory validity gate for M4-B.3).
    StorageInvariantViolation,     // live-list/allocator invariant broken
    DanglingCellConnectivity,      // cell references a non-current node
    CellPoolNodePoolMismatch,      // cell pool bound to a different node pool
    // M4-B.3 edge topology audit.
    EdgeInvariantViolation,        // edge-record/list/table invariant broken
    GhostEdge,                     // live edge with no current-cell backing
    EdgeIncidenceMismatch,         // use count != incidence, or stale cell
    EdgePoolBindingMismatch,       // edges bound to different node/cell pool
    FaceInvariantViolation,        // face-record/list/table invariant broken
    GhostFace,                     // live face with no current-cell backing
    FacePoolBindingMismatch,       // faces bound to different node/cell pool
    FaceOrientationMismatch,       // 2-cell interior face without opposite senses
    FaceIncidenceMismatch,         // recorded face/sense != recomputed from cell
    GhostAdjacency,                // adjacency references a removed cell
    AdjacencyMismatch,             // recorded adjacency != recomputed relations
    AdjacencyPoolBindingMismatch,  // adjacency bound to different pools
    AdjacencyIncompleteIncidence,  // a current 3D cell lacks face incidence
};

inline const char* toString(Severity s) {
    switch (s) {
        case Severity::Info: return "info";
        case Severity::Warning: return "warning";
        case Severity::UnsupportedInput: return "unsupported-input";
        case Severity::RecoverableFailure: return "recoverable-failure";
        case Severity::InvalidOutput: return "INVALID-OUTPUT";
    }
    return "?";
}

inline const char* toString(Code c) {
    switch (c) {
        case Code::EmptyMesh: return "empty-mesh";
        case Code::NonFiniteCoordinates: return "non-finite-coordinates";
        case Code::InvertedElement: return "inverted-element";
        case Code::InvalidTet10: return "invalid-tet10";
        case Code::NonManifoldFace: return "non-manifold-face";
        case Code::BoundaryOrphanFacet: return "boundary-orphan-facet";
        case Code::BoundaryUnexpectedFacet: return "boundary-unexpected-facet";
        case Code::UnusedNodes: return "unused-nodes";
        case Code::ClassificationInconsistent: return "classification-inconsistent";
        case Code::RegionVolumeMismatch: return "region-volume-mismatch";
        case Code::MultipleComponents: return "multiple-components";
        case Code::ComponentCountMismatch: return "component-count-mismatch";
        case Code::MalformedConnectivity: return "malformed-connectivity";
        case Code::DuplicateBoundaryFacet: return "duplicate-boundary-facet";
        case Code::AuditCoverageDeferred: return "audit-coverage-deferred";
        case Code::StorageInvariantViolation:
            return "storage-invariant-violation";
        case Code::DanglingCellConnectivity:
            return "dangling-cell-connectivity";
        case Code::CellPoolNodePoolMismatch:
            return "cell-pool-node-pool-mismatch";
        case Code::EdgeInvariantViolation: return "edge-invariant-violation";
        case Code::GhostEdge: return "ghost-edge";
        case Code::EdgeIncidenceMismatch: return "edge-incidence-mismatch";
        case Code::EdgePoolBindingMismatch:
            return "edge-pool-binding-mismatch";
        case Code::FaceInvariantViolation: return "face-invariant-violation";
        case Code::GhostFace: return "ghost-face";
        case Code::FacePoolBindingMismatch:
            return "face-pool-binding-mismatch";
        case Code::FaceOrientationMismatch:
            return "face-orientation-mismatch";
        case Code::FaceIncidenceMismatch:
            return "face-incidence-mismatch";
        case Code::GhostAdjacency: return "ghost-adjacency";
        case Code::AdjacencyMismatch: return "adjacency-mismatch";
        case Code::AdjacencyPoolBindingMismatch:
            return "adjacency-pool-binding-mismatch";
        case Code::AdjacencyIncompleteIncidence:
            return "adjacency-incomplete-incidence";
    }
    return "?";
}

struct Diagnostic {
    Code code{};
    Severity severity{};
    std::string stage;                    // producing stage, e.g. "audit"
    std::string message;                  // human-readable specifics
    std::string suggestion;               // machine-readable-ish remedy hint
    mesh::Classification entity{};        // attributed geometry entity
    Index cellIndex = kInvalidIndex;      // attributed mesh cell, if any
    Index nodeIndex = kInvalidIndex;      // attributed mesh node, if any
    Vec3 location{0, 0, 0};
    bool hasLocation = false;
};

struct Report {
    std::vector<Diagnostic> items;

    void add(Diagnostic d) { items.push_back(std::move(d)); }

    Severity worst() const {
        Severity w = Severity::Info;
        for (const auto& d : items)
            if (static_cast<int>(d.severity) > static_cast<int>(w))
                w = d.severity;
        return w;
    }
    std::size_t count(Severity s) const {
        std::size_t n = 0;
        for (const auto& d : items) n += d.severity == s ? 1 : 0;
        return n;
    }
    // True when nothing exceeds maxAcceptable.
    bool acceptable(Severity maxAcceptable = Severity::Warning) const {
        return static_cast<int>(worst()) <= static_cast<int>(maxAcceptable);
    }

    std::string toString() const {
        std::ostringstream os;
        if (items.empty()) return "no diagnostics\n";
        for (const auto& d : items) {
            os << '[' << diag::toString(d.severity) << "] "
               << diag::toString(d.code);
            if (!d.stage.empty()) os << " (" << d.stage << ')';
            os << ": " << d.message;
            if (d.entity.dim >= 0)
                os << " [entity dim " << d.entity.dim << " tag "
                   << d.entity.tag << ']';
            if (d.cellIndex != kInvalidIndex) os << " [cell " << d.cellIndex << ']';
            if (d.nodeIndex != kInvalidIndex) os << " [node " << d.nodeIndex << ']';
            if (!d.suggestion.empty()) os << " -> " << d.suggestion;
            os << '\n';
        }
        return os.str();
    }
};

} // namespace paacs::diag
