// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/SurfaceSizingPlanner.hpp"

#include "paacs/core/Parallel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <new>

namespace paacs::algo {
namespace {
struct EdgeKey {
    Index a = -1, b = -1;
    friend bool operator<(const EdgeKey &x, const EdgeKey &y) noexcept {
        return x.a < y.a || (x.a == y.a && x.b < y.b);
    }
};
struct EdgeRow {
    EdgeKey key{};
    int faceTag = -1;
};
bool finite(const Vec3 &p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
} // namespace

SurfaceSizingPlanOutcome planSurfaceSizingPass(const mesh::Mesh &surface,
                                               const SurfaceSizingPlanRequest &request) noexcept {
    SurfaceSizingPlanOutcome out;
    if (surface.nodeCount() == 0 || surface.cellCount() == 0) {
        out.reason = SurfaceSizingPlanReason::EmptySurface;
        return out;
    }
    if (!request.field) {
        out.reason = SurfaceSizingPlanReason::MissingField;
        return out;
    }
    if (!std::isfinite(request.maximumEdgeToTargetRatio) ||
        request.maximumEdgeToTargetRatio <= 0.0) {
        out.reason = SurfaceSizingPlanReason::InvalidTargetSize;
        return out;
    }
    try {
        std::vector<EdgeRow> raw;
        raw.reserve(static_cast<std::size_t>(surface.cellCount()) * 4U);
        for (Index ci = 0; ci < surface.cellCount(); ++ci) {
            const auto &cell = surface.cell(ci);
            if (cell.type != mesh::CellType::Tri3 && cell.type != mesh::CellType::Quad4) {
                out.reason = SurfaceSizingPlanReason::UnsupportedCell;
                return out;
            }
            if (cell.onEntity.dim != 2 || cell.onEntity.tag < 0) {
                out.reason = SurfaceSizingPlanReason::InvalidClassification;
                return out;
            }
            const int count = mesh::nodesPerCell(cell.type);
            const Index *conn = surface.cellNodes(ci);
            for (int e = 0; e < count; ++e) {
                Index a = conn[e], b = conn[(e + 1) % count];
                if (a < 0 || b < 0 || a >= surface.nodeCount() || b >= surface.nodeCount() ||
                    a == b || !finite(surface.node(a)) || !finite(surface.node(b))) {
                    out.reason = SurfaceSizingPlanReason::InvalidConnectivity;
                    return out;
                }
                if (a > b)
                    std::swap(a, b);
                raw.push_back({{a, b}, cell.onEntity.tag});
            }
        }
        std::sort(raw.begin(), raw.end(), [](const EdgeRow &a, const EdgeRow &b) {
            if (a.key < b.key)
                return true;
            if (b.key < a.key)
                return false;
            return a.faceTag < b.faceTag;
        });
        std::vector<EdgeRow> edges;
        edges.reserve(raw.size());
        for (const auto &row : raw) {
            if (edges.empty() || edges.back().key.a != row.key.a || edges.back().key.b != row.key.b)
                edges.push_back(row);
            else
                edges.back().faceTag = std::min(edges.back().faceTag, row.faceTag);
        }
        std::vector<SurfaceSizingPlacement> candidate(edges.size());
        std::vector<std::uint8_t> active(edges.size(), 0), projectionFailed(edges.size(), 0),
            invalid(edges.size(), 0);
        core::parallelForWithThreads(0, edges.size(), request.executionThreads, [&](std::size_t i) {
            const auto &edge = edges[i];
            const Vec3 a = surface.node(edge.key.a), b = surface.node(edge.key.b);
            const Vec3 midpoint = (a + b) * 0.5;
            const double target = request.field->eval(midpoint);
            if (!std::isfinite(target) || !(target > 0.0)) {
                invalid[i] = 1;
                return;
            }
            const double length = (b - a).norm();
            const double ratio = length / target;
            if (ratio <= request.maximumEdgeToTargetRatio)
                return;
            Vec3 point = midpoint;
            if (request.project && !request.project(edge.faceTag, midpoint, point)) {
                projectionFailed[i] = 1;
                return;
            }
            if (!finite(point)) {
                projectionFailed[i] = 1;
                return;
            }
            candidate[i] = {edge.faceTag, edge.key.a, edge.key.b, point, length, target, ratio};
            active[i] = 1;
        });
        if (std::find(invalid.begin(), invalid.end(), 1) != invalid.end()) {
            out.reason = SurfaceSizingPlanReason::InvalidTargetSize;
            return out;
        }
        if (std::find(projectionFailed.begin(), projectionFailed.end(), 1) !=
            projectionFailed.end()) {
            out.reason = SurfaceSizingPlanReason::ProjectionFailed;
            return out;
        }
        for (std::size_t i = 0; i < candidate.size(); ++i)
            if (active[i])
                out.placements.push_back(candidate[i]);
        std::sort(out.placements.begin(), out.placements.end(), [](const auto &a, const auto &b) {
            if (a.violationRatio != b.violationRatio)
                return a.violationRatio > b.violationRatio;
            if (a.faceTag != b.faceTag)
                return a.faceTag < b.faceTag;
            if (a.firstNode != b.firstNode)
                return a.firstNode < b.firstNode;
            return a.secondNode < b.secondNode;
        });
        if (out.placements.size() > request.maximumPlacements) {
            out.reason = SurfaceSizingPlanReason::BudgetExceeded;
            out.placements.clear();
            return out;
        }
        out.metrics.sourceEdges = edges.size();
        out.metrics.violatingEdges = out.placements.size();
        out.metrics.projectedPlacements = request.project ? out.placements.size() : 0;
        for (const auto &p : out.placements)
            out.metrics.maximumViolationRatio =
                std::max(out.metrics.maximumViolationRatio, p.violationRatio);
        out.planned = true;
        out.reason = SurfaceSizingPlanReason::None;
        return out;
    } catch (const std::bad_alloc &) {
        out.reason = SurfaceSizingPlanReason::AllocationFailed;
    } catch (...) {
        out.reason = SurfaceSizingPlanReason::InternalError;
    }
    return out;
}

} // namespace paacs::algo
