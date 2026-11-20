// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/mesh/Mesh.hpp"
#include "paacs/sizing/BoundedSizeField.hpp"
#include <atomic>

namespace paacs::algo {
enum class MeshSizeAuditReason : std::uint8_t {
    None,
    InvalidRequest,
    InvalidMesh,
    UnsupportedCell,
    FieldEvaluationFailed,
    BudgetExceeded,
    Cancelled,
    SizeRequirement
};
struct MeshSizeAuditOptions {
    double sizeCap = 1e100;
    // Zero disables each gate. Positive ratios must be >= 1.
    double warningEdgeToTargetRatio = 0, maximumEdgeToTargetRatio = 0;
    std::uint64_t maximumCells = 5'000'000, maximumSizingWork = 100'000'000;
    const std::atomic<bool> *cancellation = nullptr;
};
struct MeshSizeAudit {
    bool complete = false, requirementMet = false;
    MeshSizeAuditReason reason = MeshSizeAuditReason::InvalidRequest;
    std::uint64_t volumeCells = 0, sizingWork = 0, cellsAboveWarning = 0, cellsAboveLimit = 0;
    double minimumTarget = 0, maximumEdgeLength = 0;
    double maximumConservativeRatio = 0, maximumSampledRatio = 0;
    Index worstCell = kInvalidIndex;
};
// Linear Tet4/Hex8/Prism6/Pyramid5 only. Full cell-box bounds provide a
// conservative maximum-edge/target ratio; it can overestimate the true value.
// Vertex/edge-midpoint samples give a separate diagnostic, never the strict
// acceptance certificate. This checks size, not topology, Jacobians or overlap.
// No allocation, mutation, or percentile sorting. Boundary cells are skipped.
MeshSizeAudit auditMeshSize(const mesh::Mesh &, const sizing::BoundedSizeField &,
                            const MeshSizeAuditOptions & = {}) noexcept;
} // namespace paacs::algo
