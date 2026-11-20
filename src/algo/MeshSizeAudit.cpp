// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MeshSizeAudit.hpp"
#include "paacs/topo/ElementTraits.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace paacs::algo {
MeshSizeAudit auditMeshSize(const mesh::Mesh &mesh, const sizing::BoundedSizeField &field,
                            const MeshSizeAuditOptions &o) noexcept {
    MeshSizeAudit out;
    auto fail = [&](MeshSizeAuditReason r) {
        out.reason = r;
        return out;
    };
    const auto ratio = [](double x) { return std::isfinite(x) && (x == 0 || x >= 1); };
    if (!field.valid() || !(o.sizeCap > 0) || !std::isfinite(o.sizeCap) ||
        !ratio(o.warningEdgeToTargetRatio) || !ratio(o.maximumEdgeToTargetRatio))
        return out;
    const auto cost = field.boxQueryWorkBound();
    if (!cost)
        return out;
    if (static_cast<std::uint64_t>(mesh.cellCount()) > o.maximumCells)
        return fail(MeshSizeAuditReason::BudgetExceeded);
    out.minimumTarget = INFINITY;
    try {
        for (Index i = 0; i < mesh.cellCount(); ++i) {
            if (o.cancellation && o.cancellation->load(std::memory_order_relaxed))
                return fail(MeshSizeAuditReason::Cancelled);
            const auto &cell = mesh.cell(i);
            topo::ElementType type;
            if (!topo::elementTypeOf(cell.type, type))
                return fail(MeshSizeAuditReason::UnsupportedCell);
            const auto &tr = topo::traits(type);
            if (cell.first < 0 || cell.first > mesh.connectivityCount() ||
                tr.nodeCount > mesh.connectivityCount() - cell.first)
                return fail(MeshSizeAuditReason::InvalidMesh);
            const auto *ids = mesh.cellNodes(i);
            if (!topo::validConnectivity(type, {ids, tr.nodeCount}, mesh.nodeCount()))
                return fail(MeshSizeAuditReason::InvalidMesh);
            if (tr.dimension != 3)
                continue;
            if (tr.nodeCount != tr.cornerCount)
                return fail(MeshSizeAuditReason::UnsupportedCell);
            const std::uint64_t queries = 1ULL + tr.nodeCount + tr.edgeCount;
            if (cost > UINT64_MAX / queries ||
                cost * queries > o.maximumSizingWork - out.sizingWork)
                return fail(MeshSizeAuditReason::BudgetExceeded);
            out.sizingWork += cost * queries;
            Vec3 lo = mesh.node(ids[0]), hi = lo;
            double sampled = INFINITY, largest = 0;
            for (unsigned j = 0; j < tr.nodeCount; ++j) {
                const Vec3 p = mesh.node(ids[j]);
                if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
                    return fail(MeshSizeAuditReason::InvalidMesh);
                lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
                hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
                const double h = field.eval(p);
                if (!(h > 0) || !std::isfinite(h))
                    return fail(MeshSizeAuditReason::FieldEvaluationFailed);
                sampled = std::min(sampled, h);
            }
            for (unsigned j = 0; j < tr.edgeCount; ++j) {
                const Vec3 a = mesh.node(ids[tr.edges[j].a]), b = mesh.node(ids[tr.edges[j].b]);
                const auto d = b - a;
                largest = std::max(largest, std::hypot(d.x, d.y, d.z));
                const double h = field.eval(a + d * .5);
                if (!(h > 0) || !std::isfinite(h))
                    return fail(MeshSizeAuditReason::FieldEvaluationFailed);
                sampled = std::min(sampled, h);
            }
            const double bound = field.minimumOnBox(lo, hi);
            if (!(bound > 0) || !std::isfinite(bound) || !std::isfinite(largest))
                return fail(MeshSizeAuditReason::FieldEvaluationFailed);
            const double target = std::min(bound, o.sizeCap), r = largest / target;
            if (!std::isfinite(r))
                return fail(MeshSizeAuditReason::FieldEvaluationFailed);
            out.minimumTarget = std::min(out.minimumTarget, target);
            out.maximumEdgeLength = std::max(out.maximumEdgeLength, largest);
            if (r > out.maximumConservativeRatio) {
                out.maximumConservativeRatio = r;
                out.worstCell = i;
            }
            out.maximumSampledRatio =
                std::max(out.maximumSampledRatio, largest / std::min(sampled, o.sizeCap));
            const double slack = 1 + 128 * std::numeric_limits<double>::epsilon();
            out.cellsAboveWarning +=
                o.warningEdgeToTargetRatio > 0 && r > o.warningEdgeToTargetRatio * slack;
            out.cellsAboveLimit +=
                o.maximumEdgeToTargetRatio > 0 && r > o.maximumEdgeToTargetRatio * slack;
            ++out.volumeCells;
        }
    } catch (...) {
        return fail(MeshSizeAuditReason::FieldEvaluationFailed);
    }
    if (!out.volumeCells)
        return fail(MeshSizeAuditReason::InvalidMesh);
    out.complete = true;
    out.requirementMet = out.cellsAboveLimit == 0;
    out.reason =
        out.requirementMet ? MeshSizeAuditReason::None : MeshSizeAuditReason::SizeRequirement;
    return out;
}
} // namespace paacs::algo
