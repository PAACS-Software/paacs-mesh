// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/LayeredPatchMesher.hpp"
#include "paacs/predicates/Predicates.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <new>
#include <numeric>

namespace paacs::algo {
namespace {
bool positive(double x) noexcept {
    return std::isfinite(x) && x > 0;
}
bool finite(Vec3 p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
double component(Vec3 p, unsigned a) noexcept {
    return a == 0 ? p.x : a == 1 ? p.y : p.z;
}
constexpr double roundoff = 128 * std::numeric_limits<double>::epsilon();
bool charge(std::uint64_t &sum, std::uint64_t count, std::uint64_t width,
            std::uint64_t limit) noexcept {
    if (sum > limit || count > (limit - sum) / width)
        return false;
    sum += count * width;
    return true;
}
} // namespace
LayerSchedule planLayerSchedule(const LayerScheduleRequest &r) noexcept {
    LayerSchedule out;
    auto fail = [&](LayerScheduleReason why) {
        out.reason = why;
        return std::move(out);
    };
    if (!r.wallLayers || !positive(r.firstLayerThickness) || !positive(r.growthRatio) ||
        r.growthRatio < 1 || !positive(r.maximumAllowedGrowth) ||
        r.maximumAllowedGrowth < r.growthRatio || !std::isfinite(r.totalDepth) ||
        r.totalDepth < 0 || !positive(r.maximumCoreThickness))
        return fail(LayerScheduleReason::InvalidRequest);
    if (r.wallLayers > r.maximumLayers)
        return fail(LayerScheduleReason::LayerBudget);
    if (static_cast<std::uint64_t>(r.wallLayers) + 1 > r.maximumBytes / sizeof(double))
        return fail(LayerScheduleReason::ByteBudget);
    long double wall = 0, width = r.firstLayerThickness, last = width;
    for (std::uint32_t i = 0; i < r.wallLayers; ++i) {
        wall += width;
        last = width;
        width *= r.growthRatio;
    }
    if (!std::isfinite(static_cast<double>(wall)))
        return fail(LayerScheduleReason::NumericallyUnresolved);
    out.wallLayers = r.wallLayers;
    out.wallThickness = static_cast<double>(wall);
    out.totalDepth = r.totalDepth == 0 ? out.wallThickness : r.totalDepth;
    // Decimal wall widths can sum a few ulps above an exactly matching
    // requested depth (e.g. three 0.2 layers in depth 0.6). Only the final
    // boundary absorbs roundoff; a meaningful overrun still refuses.
    if (out.totalDepth < out.wallThickness && wall - out.totalDepth > roundoff * out.totalDepth)
        return fail(LayerScheduleReason::WallLayersExceedDepth);
    out.wallThickness = std::min(out.wallThickness, out.totalDepth);
    const long double remaining = std::max(0.0L, static_cast<long double>(out.totalDepth) - wall);
    std::uint32_t core = 0;
    long double coreSum = 0,
                firstCore = std::min(last * r.growthRatio,
                                     static_cast<long double>(r.maximumCoreThickness)),
                lastCore = firstCore, scale = 1;
    if (remaining > roundoff * out.totalDepth) {
        width = firstCore;
        while (coreSum < remaining) {
            if (static_cast<std::uint64_t>(core) + r.wallLayers >= r.maximumLayers)
                return fail(LayerScheduleReason::LayerBudget);
            if (static_cast<std::uint64_t>(core) + r.wallLayers + 2 >
                r.maximumBytes / sizeof(double))
                return fail(LayerScheduleReason::ByteBudget);
            coreSum += width;
            lastCore = width;
            ++core;
            width =
                std::min(width * r.growthRatio, static_cast<long double>(r.maximumCoreThickness));
        }
        scale = remaining / coreSum;
        auto fits = [&] {
            const long double first = firstCore * scale;
            return first > 0 &&
                   std::max(first / last, last / first) <=
                       r.maximumAllowedGrowth * (1 + roundoff) &&
                   lastCore * scale <= r.maximumCoreThickness * (1 + roundoff);
        };
        if (!fits() && core > 1) {
            coreSum -= lastCore;
            --core;
            lastCore = firstCore;
            for (std::uint32_t i = 1; i < core; ++i)
                lastCore = std::min(lastCore * r.growthRatio,
                                    static_cast<long double>(r.maximumCoreThickness));
            scale = remaining / coreSum;
        }
        if (!fits())
            return fail(LayerScheduleReason::TransitionDoesNotFit);
    }
    out.coreLayers = core;
    const auto count = static_cast<std::uint64_t>(r.wallLayers) + core + 1;
    if (count > r.maximumBytes / sizeof(double))
        return fail(LayerScheduleReason::ByteBudget);
    try {
        out.offsets.resize(static_cast<std::size_t>(count));
        long double position = 0;
        width = r.firstLayerThickness;
        for (std::uint32_t i = 0; i < r.wallLayers; ++i) {
            position += width;
            out.offsets[i + 1] = static_cast<double>(position);
            width *= r.growthRatio;
        }
        width = firstCore;
        for (std::uint32_t i = 0; i < core; ++i) {
            position += width * scale;
            out.offsets[static_cast<std::size_t>(r.wallLayers) + i + 1] =
                static_cast<double>(position);
            width =
                std::min(width * r.growthRatio, static_cast<long double>(r.maximumCoreThickness));
        }
        out.offsets.back() = out.totalDepth;
        out.minimumThickness = std::numeric_limits<double>::infinity();
        double previous = 0;
        for (std::size_t i = 1; i < out.offsets.size(); ++i) {
            const double h = out.offsets[i] - out.offsets[i - 1];
            if (!positive(h))
                return fail(LayerScheduleReason::NumericallyUnresolved);
            out.minimumThickness = std::min(out.minimumThickness, h);
            out.maximumThickness = std::max(out.maximumThickness, h);
            if (previous > 0)
                out.maximumAdjacentGrowth =
                    std::max({out.maximumAdjacentGrowth, h / previous, previous / h});
            previous = h;
        }
        if (out.maximumAdjacentGrowth > r.maximumAllowedGrowth * (1 + roundoff))
            return fail(LayerScheduleReason::NumericallyUnresolved);
        out.planned = true;
        out.reason = LayerScheduleReason::None;
    } catch (const std::bad_alloc &) {
        out.reason = LayerScheduleReason::AllocationFailed;
    } catch (...) {
        out.reason = LayerScheduleReason::InternalError;
    }
    return out;
}

namespace {
struct Triangle {
    std::array<Index, 3> ids;
    std::array<Vec2, 3> p;
    double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
};
int orient(Vec2 a, Vec2 b, Vec2 c) {
    return static_cast<int>(predicates::orient2d(a, b, c));
}
bool onSegment(Vec2 p, Vec2 a, Vec2 b) noexcept {
    return p.u >= std::min(a.u, b.u) && p.u <= std::max(a.u, b.u) && p.v >= std::min(a.v, b.v) &&
           p.v <= std::max(a.v, b.v);
}
bool conflict(const Triangle &a, const Triangle &b) {
    unsigned shared = 0;
    for (auto x : a.ids)
        for (auto y : b.ids)
            if (x == y)
                ++shared;
    if (shared == 3)
        return true;
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            const unsigned ii = (i + 1) % 3, jj = (j + 1) % 3;
            const int s0 = orient(a.p[i], a.p[ii], b.p[j]), s1 = orient(a.p[i], a.p[ii], b.p[jj]);
            const int t0 = orient(b.p[j], b.p[jj], a.p[i]), t1 = orient(b.p[j], b.p[jj], a.p[ii]);
            if (s0 * s1 < 0 && t0 * t1 < 0)
                return true;
            if (!s0 && onSegment(b.p[j], a.p[i], a.p[ii]) && b.ids[j] != a.ids[i] &&
                b.ids[j] != a.ids[ii])
                return true;
            if (!s1 && onSegment(b.p[jj], a.p[i], a.p[ii]) && b.ids[jj] != a.ids[i] &&
                b.ids[jj] != a.ids[ii])
                return true;
            if (!t0 && onSegment(a.p[i], b.p[j], b.p[jj]) && a.ids[i] != b.ids[j] &&
                a.ids[i] != b.ids[jj])
                return true;
            if (!t1 && onSegment(a.p[ii], b.p[j], b.p[jj]) && a.ids[ii] != b.ids[j] &&
                a.ids[ii] != b.ids[jj])
                return true;
        }
    auto inside = [](Vec2 p, const Triangle &t) {
        const int a0 = orient(t.p[0], t.p[1], p), a1 = orient(t.p[1], t.p[2], p),
                  a2 = orient(t.p[2], t.p[0], p);
        return a0 * a1 > 0 && a0 * a2 > 0;
    };
    for (auto p : a.p)
        if (inside(p, b))
            return true;
    for (auto p : b.p)
        if (inside(p, a))
            return true;
    return false;
}
LayeredPatchReason sourceGeometry(const mesh::Mesh &m, const LayeredPatchRequest &r, Vec3 normal,
                                  LayeredPatchResult &out) {
    const auto origin = m.node(0);
    Vec3 lo = origin, hi = origin;
    for (Index i = 0; i < m.nodeCount(); ++i) {
        const auto p = m.node(i);
        if (!finite(p))
            return LayeredPatchReason::InvalidSource;
        const auto cls = m.nodeClass(i);
        if (cls.dim < 0 || cls.dim > 2 || cls.tag < 0)
            return LayeredPatchReason::InvalidSource;
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
    const double scale = std::max({hi.x - lo.x, hi.y - lo.y, hi.z - lo.z});
    if (!positive(scale))
        return LayeredPatchReason::InvalidSource;
    unsigned axis = 0;
    for (unsigned a = 1; a < 3; ++a)
        if (std::abs(component(normal, a)) > std::abs(component(normal, axis)))
            axis = a;
    const unsigned u = (axis + 1) % 3, v = (axis + 2) % 3;
    for (Index i = 0; i < m.nodeCount(); ++i)
        if (std::abs((m.node(i) - origin).dot(normal)) > 1e-10 * scale)
            return LayeredPatchReason::NonplanarSource;
    std::vector<Triangle> triangles;
    triangles.reserve(static_cast<std::size_t>(m.cellCount()) * 2);
    std::vector<std::uint8_t> used(static_cast<std::size_t>(m.nodeCount()), 0);
    auto project = [&](Index id) {
        const auto p = (m.node(id) - lo) / scale;
        return Vec2{component(p, u), component(p, v)};
    };
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto &cell = m.cell(i);
        if (cell.type != mesh::CellType::Tri3 && cell.type != mesh::CellType::Quad4)
            return LayeredPatchReason::InvalidSource;
        const auto n = mesh::nodesPerCell(cell.type);
        if (cell.first < 0 || cell.first > m.connectivityCount() ||
            n > m.connectivityCount() - cell.first)
            return LayeredPatchReason::InvalidSource;
        const auto *ids = m.cellNodes(i);
        for (int k = 0; k < n; ++k) {
            if (ids[k] < 0 || ids[k] >= m.nodeCount())
                return LayeredPatchReason::InvalidSource;
            used[static_cast<std::size_t>(ids[k])] = 1;
            for (int j = 0; j < k; ++j)
                if (ids[j] == ids[k])
                    return LayeredPatchReason::InvalidSource;
        }
        // A convex, correctly oriented planar source gives a positive full
        // straight-extrusion Jacobian, including bilinear Quad4 columns.
        for (int j = 0; j < n; ++j)
            if (orient(project(ids[j]), project(ids[(j + 1) % n]), project(ids[(j + 2) % n])) *
                    component(normal, axis) <=
                0)
                return LayeredPatchReason::InvalidSource;
        for (int j = 1; j < n - 1; ++j) {
            Triangle t;
            t.ids = {ids[0], ids[j], ids[j + 1]};
            for (unsigned k = 0; k < 3; ++k)
                t.p[k] = project(t.ids[k]);
            t.xmin = std::min({t.p[0].u, t.p[1].u, t.p[2].u});
            t.xmax = std::max({t.p[0].u, t.p[1].u, t.p[2].u});
            t.ymin = std::min({t.p[0].v, t.p[1].v, t.p[2].v});
            t.ymax = std::max({t.p[0].v, t.p[1].v, t.p[2].v});
            triangles.push_back(t);
        }
    }
    if (std::find(used.begin(), used.end(), 0) != used.end())
        return LayeredPatchReason::InvalidSource;
    std::sort(triangles.begin(), triangles.end(), [](const auto &a, const auto &b) {
        return a.xmin != b.xmin ? a.xmin < b.xmin : a.ids < b.ids;
    });
    for (std::size_t i = 0; i < triangles.size(); ++i) {
        if (r.cancellation && r.cancellation->load(std::memory_order_relaxed))
            return LayeredPatchReason::Cancelled;
        for (std::size_t j = i + 1; j < triangles.size() && triangles[j].xmin <= triangles[i].xmax;
             ++j) {
            if (out.sourcePairChecks >= r.maximumSourcePairChecks)
                return LayeredPatchReason::SourceBudget;
            ++out.sourcePairChecks;
            if (triangles[j].ymin > triangles[i].ymax || triangles[i].ymin > triangles[j].ymax)
                continue;
            if (conflict(triangles[i], triangles[j]))
                return LayeredPatchReason::IntersectingSource;
        }
    }
    return LayeredPatchReason::None;
}
} // namespace
LayeredPatchResult meshLayeredPatch(const mesh::Mesh &source, const LayeredPatchRequest &r,
                                    mesh::Mesh &destination) noexcept {
    LayeredPatchResult out;
    auto fail = [&](LayeredPatchReason why) {
        out.reason = why;
        return std::move(out);
    };
    auto cancelled = [&] {
        return r.cancellation && r.cancellation->load(std::memory_order_relaxed);
    };
    if (cancelled())
        return fail(LayeredPatchReason::Cancelled);
    if (source.nodeCount() == 0 || source.cellCount() == 0)
        return fail(LayeredPatchReason::InvalidSource);
    if (!finite(r.direction) ||
        !positive(std::hypot(r.direction.x, r.direction.y, r.direction.z)) || r.threads < 0)
        return fail(LayeredPatchReason::InvalidRequest);
    if (static_cast<std::uint64_t>(source.nodeCount()) > r.maximumSourceNodes ||
        static_cast<std::uint64_t>(source.cellCount()) > r.maximumSourceCells)
        return fail(LayeredPatchReason::SourceBudget);
    try {
        double allowedGrowth = r.schedule.maximumAllowedGrowth;
        std::uint64_t retainedExcess = 0;
        if (r.sizeProfile) {
            auto limits = r.profileLimits;
            limits.maximumCells = std::min<std::uint64_t>(limits.maximumCells, UINT32_MAX - 1);
            limits.maximumBytes = std::min(limits.maximumBytes, r.maximumWorkingBytes);
            if (!limits.cancellation)
                limits.cancellation = r.cancellation;
            auto profile = sizing::planProfileDistribution(*r.sizeProfile, limits);
            out.usedSizeProfile = true;
            out.profileWorkingBytes = profile.workingBytes;
            out.plateauCellsRemoved = profile.plateauCellsRemoved;
            if (!profile.planned) {
                if (profile.reason == sizing::ProfileReason::Cancelled)
                    return fail(LayeredPatchReason::Cancelled);
                if (profile.reason == sizing::ProfileReason::AllocationFailed)
                    return fail(LayeredPatchReason::AllocationFailed);
                if (profile.reason == sizing::ProfileReason::BudgetExceeded)
                    return fail(LayeredPatchReason::WorkingBudget);
                return fail(LayeredPatchReason::ScheduleRejected);
            }
            allowedGrowth = limits.maximumAdjacentGrowth == 0 ? r.sizeProfile->maximumGrowth()
                                                              : limits.maximumAdjacentGrowth;
            out.schedule.planned = true;
            out.schedule.reason = LayerScheduleReason::None;
            out.schedule.totalDepth = r.sizeProfile->length();
            out.schedule.coreLayers = static_cast<std::uint32_t>(profile.positions.size() - 1);
            out.schedule.minimumThickness = profile.minimumWidth;
            out.schedule.maximumThickness = profile.maximumWidth;
            out.schedule.maximumAdjacentGrowth = profile.maximumAdjacentGrowth;
            retainedExcess = profile.workingBytes - profile.positions.size() * sizeof(double);
            out.schedule.offsets = std::move(profile.positions);
        } else
            out.schedule = planLayerSchedule(r.schedule);
        if (!out.schedule.planned)
            return fail(out.schedule.reason == LayerScheduleReason::AllocationFailed
                            ? LayeredPatchReason::AllocationFailed
                            : LayeredPatchReason::ScheduleRejected);
        const auto layers = static_cast<std::uint64_t>(out.schedule.offsets.size() - 1);
        // Conservatively bound both patch validation and the legacy source
        // analyzer's node-based maps on this reference ABI, plus schedules.
        // Array payload estimates exclude allocator bookkeeping and process RSS.
        out.plannedWorkingBytes = retainedExcess;
        if (!charge(out.plannedWorkingBytes, static_cast<std::uint64_t>(source.cellCount()), 2048,
                    r.maximumWorkingBytes) ||
            !charge(out.plannedWorkingBytes, static_cast<std::uint64_t>(source.nodeCount()), 64,
                    r.maximumWorkingBytes) ||
            !charge(out.plannedWorkingBytes, layers + 1, 24, r.maximumWorkingBytes))
            return fail(LayeredPatchReason::WorkingBudget);
        const auto normal = r.direction / std::hypot(r.direction.x, r.direction.y, r.direction.z);
        const auto reason = sourceGeometry(source, r, normal, out);
        if (reason != LayeredPatchReason::None)
            return fail(reason);
        if (cancelled())
            return fail(LayeredPatchReason::Cancelled);
        const auto sourceCells = static_cast<std::uint64_t>(source.cellCount());
        std::uint64_t nodes = 0, volumeCells = 0, maxCells = 0, maxConnectivity = 0;
        const auto max = std::numeric_limits<std::uint64_t>::max();
        if (!charge(nodes, static_cast<std::uint64_t>(source.nodeCount()), layers + 1,
                    r.maximumNodes) ||
            !charge(volumeCells, sourceCells, layers, r.maximumCells) ||
            !charge(maxCells, sourceCells, 5 * layers + 2, max) ||
            !charge(maxConnectivity, sourceCells, 24 * layers + 8, max))
            return fail(LayeredPatchReason::OutputBudget);
        // Worst-case four boundary edges per source cell; exact sweep counts
        // are checked again by the underlying generator.
        if (!charge(out.plannedOutputBytes, nodes, sizeof(Vec3) + sizeof(mesh::Classification),
                    r.maximumOutputBytes) ||
            !charge(out.plannedOutputBytes, maxCells, sizeof(mesh::Cell), r.maximumOutputBytes) ||
            !charge(out.plannedOutputBytes, maxConnectivity, sizeof(Index), r.maximumOutputBytes))
            return fail(LayeredPatchReason::OutputBudget);
        PreparedMixedSweepRequest sweep;
        sweep.layers = static_cast<std::uint32_t>(layers);
        sweep.direction = normal;
        sweep.totalThickness = out.schedule.totalDepth;
        sweep.maximumAllowedGrowth = allowedGrowth * (r.sizeProfile ? 1 + 1e-10 : 1);
        sweep.layerFractions = out.schedule.offsets;
        for (auto &f : sweep.layerFractions)
            f /= out.schedule.totalDepth;
        sweep.layerFractions.front() = 0;
        sweep.layerFractions.back() = 1;
        sweep.volumeClassification = r.volumeClassification;
        sweep.targetClassification = r.targetClassification;
        sweep.sideClassification = r.sideClassification;
        sweep.intermediateNodeClassifications.assign(static_cast<std::size_t>(source.nodeCount()),
                                                     r.volumeClassification);
        {
            std::vector<std::array<Index, 2>> edges;
            edges.reserve(static_cast<std::size_t>(source.cellCount()) * 4);
            for (Index i = 0; i < source.cellCount(); ++i) {
                const auto *ids = source.cellNodes(i);
                const auto n = mesh::nodesPerCell(source.cell(i).type);
                for (int j = 0; j < n; ++j) {
                    std::array<Index, 2> e{ids[j], ids[(j + 1) % n]};
                    std::sort(e.begin(), e.end());
                    edges.push_back(e);
                }
            }
            std::sort(edges.begin(), edges.end());
            for (std::size_t i = 0; i < edges.size();) {
                std::size_t j = i + 1;
                while (j < edges.size() && edges[i] == edges[j])
                    ++j;
                if (j - i == 1)
                    for (auto node : edges[i])
                        sweep.intermediateNodeClassifications[static_cast<std::size_t>(node)] =
                            r.sideClassification;
                i = j;
            }
        }
        sweep.maximumNodes = r.maximumNodes;
        sweep.maximumCells = r.maximumCells;
        sweep.executionThreads = r.threads;
        sweep.fullAudit = true;
        mesh::Mesh staged;
        out.sweep = generatePreparedMixedSweep(source, sweep, staged);
        if (!out.sweep.published()) {
            if (out.sweep.reason == PreparedSweepReason::AllocationFailed)
                return fail(LayeredPatchReason::AllocationFailed);
            if (out.sweep.reason == PreparedSweepReason::BudgetExceeded ||
                out.sweep.reason == PreparedSweepReason::CountOverflow)
                return fail(LayeredPatchReason::OutputBudget);
            return fail(LayeredPatchReason::SweepRejected);
        }
        if (cancelled())
            return fail(LayeredPatchReason::Cancelled);
        out.outputBytes = static_cast<std::uint64_t>(staged.nodeCount()) *
                              (sizeof(Vec3) + sizeof(mesh::Classification)) +
                          static_cast<std::uint64_t>(staged.cellCount()) * sizeof(mesh::Cell) +
                          static_cast<std::uint64_t>(staged.connectivityCount()) * sizeof(Index);
        auto limits = r.auditLimits;
        limits.cellVolumes = {};
        const auto retainedBytes =
            retainedExcess + (layers + 1) * sizeof(double) * 2 +
            static_cast<std::uint64_t>(source.nodeCount()) * sizeof(mesh::Classification);
        limits.maximumTemporaryBytes =
            std::min(limits.maximumTemporaryBytes, r.maximumWorkingBytes - retainedBytes);
        if (!limits.cancellation)
            limits.cancellation = r.cancellation;
        out.audit = mesh::auditLinearVolumeBoundary(staged, limits);
        out.plannedWorkingBytes =
            std::max(out.plannedWorkingBytes, out.audit.temporaryBytes + retainedBytes);
        if (out.audit.cancelled || cancelled())
            return fail(LayeredPatchReason::Cancelled);
        if (out.audit.budgetExceeded)
            return fail(LayeredPatchReason::WorkingBudget);
        if (out.audit.allocationFailed)
            return fail(LayeredPatchReason::AllocationFailed);
        if (!out.audit.complete || !out.audit.valid)
            return fail(LayeredPatchReason::BoundaryAuditFailed);
        destination.swap(staged);
        out.published = true;
        out.destinationUnchanged = false;
        out.reason = LayeredPatchReason::None;
    } catch (const std::bad_alloc &) {
        out.reason = LayeredPatchReason::AllocationFailed;
    } catch (...) {
        out.reason = LayeredPatchReason::InternalError;
    }
    return out;
}
} // namespace paacs::algo
