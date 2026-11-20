// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/SurfaceFeatures.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>

namespace paacs::algo {
namespace {
using Reason = SurfaceFeatureReason;
struct Stop {
    Reason reason;
};
void need(bool ok, Reason reason) {
    if (!ok)
        throw Stop{reason};
}
void poll(const SurfaceFeatureOptions &o) {
    need(!o.cancellation || !o.cancellation->load(std::memory_order_relaxed), Reason::Cancelled);
}
std::uint64_t add(std::uint64_t a, std::uint64_t b) {
    need(b <= UINT64_MAX - a, Reason::WorkingBudget);
    return a + b;
}
std::uint64_t mul(std::uint64_t a, std::uint64_t b) {
    need(!a || b <= UINT64_MAX / a, Reason::WorkingBudget);
    return a * b;
}
bool finite(Vec3 p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
unsigned corners(mesh::CellType type) {
    switch (type) {
    case mesh::CellType::Tri3:
    case mesh::CellType::Tri6:
        return 3;
    case mesh::CellType::Quad4:
        return 4;
    default:
        return 0;
    }
}
struct Face {
    Vec3 normal;
    Index cell;
    int tag;
};
struct Use {
    std::uint64_t key, faceDirection;
};
constexpr std::uint8_t flag(SurfaceFeatureFlag value) {
    return static_cast<std::uint8_t>(value);
}
} // namespace

SurfaceFeatureResult detectSurfaceFeatures(const mesh::Mesh &m,
                                           const SurfaceFeatureOptions &o) noexcept {
    SurfaceFeatureResult r;
    if (!o.enabled) {
        r.completed = true;
        r.reason = Reason::Disabled;
        return r;
    }
    try {
        poll(o);
        need(std::isfinite(o.angleDegrees) && o.angleDegrees >= 0 && o.angleDegrees <= 180 &&
                 std::isfinite(o.angularToleranceDegrees) && o.angularToleranceDegrees >= 0 &&
                 o.angularToleranceDegrees <= 180 && std::isfinite(o.minimumLength) &&
                 o.minimumLength >= 0,
             Reason::InvalidRequest);
        need(m.nodeCount() <= UINT32_MAX, Reason::WorkingBudget);
        for (Index i = 0; i < m.cellCount(); ++i) {
            if (i % 1024 == 0)
                poll(o);
            const auto &cell = m.cell(i);
            need(mesh::nodesPerCell(cell.type) != 0, Reason::InvalidConnectivity);
            if (const auto n = corners(cell.type)) {
                r.metrics.surfaceFaces = add(r.metrics.surfaceFaces, 1);
                r.metrics.highOrderFaces += cell.type == mesh::CellType::Tri6;
                r.metrics.edgeUses = add(r.metrics.edgeUses, n);
            }
        }
        need(r.metrics.surfaceFaces > 0, Reason::EmptySurface);
        r.metrics.workingBytes =
            add(mul(r.metrics.surfaceFaces, sizeof(Face)), mul(r.metrics.edgeUses, sizeof(Use)));
        need(r.metrics.edgeUses <= o.maximumEdgeUses &&
                 r.metrics.workingBytes <= o.maximumWorkingBytes &&
                 r.metrics.workingBytes <= SIZE_MAX,
             Reason::WorkingBudget);
        std::vector<Face> faces(r.metrics.surfaceFaces);
        std::vector<Use> uses(r.metrics.edgeUses);
        std::uint64_t face = 0, at = 0;
        for (Index i = 0; i < m.cellCount(); ++i) {
            if (i % 1024 == 0)
                poll(o);
            const auto &cell = m.cell(i);
            const auto n = corners(cell.type);
            if (!n)
                continue;
            const auto arity = mesh::nodesPerCell(cell.type);
            need(cell.first >= 0 && cell.first <= m.connectivityCount() &&
                     arity <= m.connectivityCount() - cell.first,
                 Reason::InvalidConnectivity);
            const auto *v = m.cellNodes(i);
            for (int j = 0; j < arity; ++j) {
                need(v[j] >= 0 && v[j] < m.nodeCount(), Reason::InvalidConnectivity);
                need(finite(m.node(v[j])), Reason::NumericallyUnresolved);
                for (int k = 0; k < j; ++k)
                    need(v[j] != v[k], Reason::InvalidConnectivity);
            }
            std::array<Vec3, 4> p{};
            double scale = 0;
            for (unsigned j = 1; j < n; ++j) {
                p[j] = m.node(v[j]) - m.node(v[0]);
                need(finite(p[j]), Reason::NumericallyUnresolved);
                scale = std::max({scale, std::abs(p[j].x), std::abs(p[j].y), std::abs(p[j].z)});
            }
            need(scale > 0, Reason::NumericallyUnresolved);
            for (unsigned j = 1; j < n; ++j)
                p[j] = p[j] / scale;
            Vec3 normal{};
            for (unsigned j = 1; j + 1 < n; ++j)
                normal = normal + p[j].cross(p[j + 1]);
            const double length = std::hypot(normal.x, normal.y, normal.z);
            need(length > 128 * std::numeric_limits<double>::epsilon() && std::isfinite(length),
                 Reason::NumericallyUnresolved);
            faces[face] = {normal / length, i, cell.onEntity.dim == 2 ? cell.onEntity.tag : -1};
            for (unsigned j = 0; j < n; ++j) {
                const auto a = v[j], b = v[(j + 1) % n];
                const auto key =
                    (std::uint64_t(std::min(a, b)) << 32) | std::uint64_t(std::max(a, b));
                uses[at++] = {key, 2 * face + (a > b ? 1 : 0)};
            }
            ++face;
        }
        poll(o);
        std::sort(uses.begin(), uses.end(), [](Use a, Use b) {
            return a.key < b.key || (a.key == b.key && a.faceDirection < b.faceDirection);
        });
        poll(o);
        // Two linear scans: count/filter first, then allocate the exact output.
        // Never reserve an unbounded worst-case edge list or a node hash table.
        auto scan = [&](bool emit) {
            std::uint64_t output = 0, nextPoll = 0;
            for (std::size_t begin = 0; begin < uses.size();) {
                if (begin >= nextPoll) {
                    poll(o);
                    nextPoll = begin + 1024;
                }
                auto end = begin + 1;
                while (end < uses.size() && uses[end].key == uses[begin].key) {
                    if (end % 1024 == 0)
                        poll(o);
                    ++end;
                }
                const auto count = end - begin;
                if (!emit)
                    ++r.metrics.uniqueEdges;
                if (count > 2) {
                    if (!emit)
                        ++r.metrics.nonManifoldEdges;
                    begin = end;
                    continue;
                }
                const auto &a = faces[uses[begin].faceDirection / 2];
                SurfaceFeatureEdge edge;
                edge.first = static_cast<Index>(uses[begin].key >> 32);
                edge.second = static_cast<Index>(uses[begin].key & UINT32_MAX);
                edge.firstFace = a.cell;
                if (count == 1) {
                    if (!emit)
                        ++r.metrics.openBoundaryEdges;
                    if (o.includeOpenBoundaries)
                        edge.flags |= flag(SurfaceFeatureFlag::OpenBoundary);
                } else {
                    const auto &b = faces[uses[begin + 1].faceDirection / 2];
                    edge.secondFace = b.cell;
                    if ((uses[begin].faceDirection & 1) == (uses[begin + 1].faceDirection & 1)) {
                        if (!emit)
                            ++r.metrics.inconsistentWindingEdges;
                        begin = end;
                        continue;
                    }
                    const auto cross = a.normal.cross(b.normal);
                    edge.angleDegrees = std::atan2(std::hypot(cross.x, cross.y, cross.z),
                                                   std::clamp(a.normal.dot(b.normal), -1., 1.)) *
                                        180 / std::numbers::pi;
                    if (edge.angleDegrees > o.angleDegrees + o.angularToleranceDegrees) {
                        edge.flags |= flag(SurfaceFeatureFlag::Sharp);
                        if (!emit)
                            ++r.metrics.sharpEdges;
                    }
                    if (a.tag >= 0 && b.tag >= 0 && a.tag != b.tag) {
                        if (!emit)
                            ++r.metrics.tagBoundaryEdges;
                        if (o.includeTagBoundaries)
                            edge.flags |= flag(SurfaceFeatureFlag::TagBoundary);
                    }
                }
                const auto d = m.node(edge.second) - m.node(edge.first);
                edge.length = std::hypot(d.x, d.y, d.z);
                need(edge.length > 0 && std::isfinite(edge.length), Reason::NumericallyUnresolved);
                if (edge.flags) {
                    if (edge.length < o.minimumLength) {
                        if (!emit)
                            ++r.metrics.filteredShortEdges;
                    } else {
                        if (emit)
                            r.edges[output] = edge;
                        ++output;
                    }
                }
                begin = end;
            }
            return output;
        };
        r.metrics.selectedEdges = scan(false);
        need(r.metrics.nonManifoldEdges == 0, Reason::NonManifoldSurface);
        need(r.metrics.inconsistentWindingEdges == 0, Reason::InconsistentWinding);
        r.metrics.outputBytes = mul(r.metrics.selectedEdges, sizeof(SurfaceFeatureEdge));
        need(r.metrics.selectedEdges <= o.maximumEdges &&
                 r.metrics.outputBytes <= o.maximumOutputBytes && r.metrics.outputBytes <= SIZE_MAX,
             Reason::OutputBudget);
        r.edges.resize(r.metrics.selectedEdges);
        need(scan(true) == r.metrics.selectedEdges, Reason::InternalError);
        poll(o);
        r.completed = true;
        r.reason = Reason::None;
    } catch (const Stop &e) {
        r.reason = e.reason;
    } catch (const std::bad_alloc &) {
        r.reason = Reason::AllocationFailed;
    } catch (...) {
        r.reason = Reason::InternalError;
    }
    if (!r.completed)
        std::vector<SurfaceFeatureEdge>().swap(r.edges);
    return r;
}
const char *surfaceFeatureReasonName(SurfaceFeatureReason reason) noexcept {
    switch (reason) {
    case Reason::None:
        return "none";
    case Reason::Disabled:
        return "disabled";
    case Reason::InvalidRequest:
        return "invalid-request";
    case Reason::EmptySurface:
        return "empty-surface";
    case Reason::InvalidConnectivity:
        return "invalid-connectivity";
    case Reason::NumericallyUnresolved:
        return "numerically-unresolved";
    case Reason::NonManifoldSurface:
        return "nonmanifold-surface";
    case Reason::InconsistentWinding:
        return "inconsistent-winding";
    case Reason::WorkingBudget:
        return "working-budget";
    case Reason::OutputBudget:
        return "output-budget";
    case Reason::Cancelled:
        return "cancelled";
    case Reason::AllocationFailed:
        return "allocation-failed";
    case Reason::InternalError:
        return "internal-error";
    }
    return "unknown";
}
} // namespace paacs::algo
