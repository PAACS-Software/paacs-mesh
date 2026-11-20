// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/algo/Meshers.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/sizing/SizeField.hpp"

#include <cstdint>
#include <span>
#include <vector>

namespace paacs::algo {

struct CanonicalCurveOutcome;

struct CanonicalCurveRecord {
    int edgeTag = -1;
    int startVertexTag = -1;
    int endVertexTag = -1;
    std::vector<Index> naturalNodeOrder;
    std::vector<Index> naturalLineOrder;
};

class CanonicalCurveSet {
public:
    const mesh::Mesh& mesh() const noexcept { return mesh_; }
    std::span<const CanonicalCurveRecord> curves() const noexcept {
        return curves_;
    }
    const CanonicalCurveRecord* curve(int edgeTag) const noexcept;
    std::vector<Index> orientedNodes(int edgeTag, bool reversed) const;
    std::uint64_t deterministicHash() const noexcept { return hash_; }
    bool empty() const noexcept { return curves_.empty(); }

private:
    friend struct CanonicalCurveOutcome;
    friend CanonicalCurveOutcome buildCanonicalCurves(
        const geo::Model&, const sizing::MeshSpec&) noexcept;

    mesh::Mesh mesh_;
    std::vector<CanonicalCurveRecord> curves_;
    std::uint64_t hash_ = 0;
};

enum class CanonicalCurveStatus : std::uint8_t {
    Built,
    Rejected,
    Failed
};

enum class CanonicalCurveReason : std::uint8_t {
    None,
    EmptyModel,
    InvalidEdge,
    MissingLineCells,
    DisconnectedCurve,
    BranchingCurve,
    EndpointMismatch,
    InvalidClassification,
    AllocationFailed,
    InternalError
};

struct CanonicalCurveOutcome {
    CanonicalCurveStatus status = CanonicalCurveStatus::Failed;
    CanonicalCurveReason reason = CanonicalCurveReason::InternalError;
    CanonicalCurveSet result{};
    bool built() const noexcept {
        return status == CanonicalCurveStatus::Built;
    }
};

// Builds each topological edge exactly once using the production edge mesher,
// then records its natural node order.  Faces can consume the same record
// forward or reversed; no coordinate-based resampling is required.
CanonicalCurveOutcome buildCanonicalCurves(
    const geo::Model& model, const sizing::MeshSpec& spec) noexcept;

// Atomically appends a complete canonical curve artifact to destination.
// Node identity inside the artifact is preserved through one deterministic
// artifact-node mapping.  Existing destination content remains unchanged on
// failure.
bool appendCanonicalCurves(const CanonicalCurveSet& curves,
                           mesh::Mesh& destination) noexcept;

}  // namespace paacs::algo
