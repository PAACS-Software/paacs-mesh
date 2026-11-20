// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/AdaptiveBoxMesher.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/mesh/CertifiedMeshStorage.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <vector>

namespace paacs::algo {
namespace {
using Clock = std::chrono::steady_clock;
constexpr auto missing = std::numeric_limits<std::uint32_t>::max();
struct Node {
    std::array<std::uint32_t, 3> p{};
    std::uint32_t width = 0, child = missing;
    std::uint8_t depth = 0;
    bool marked = false;
};
struct Key {
    std::array<std::uint32_t, 3> p{};
    friend auto operator<=>(const Key &, const Key &) = default;
};
struct Leaf {
    std::uint32_t tree = 0;
    std::uint64_t volume = 0, boundary = 0;
};
struct Patch {
    Key lower;
    std::uint32_t width = 0;
    unsigned axis = 0, side = 0;
    bool exterior = false;
};
bool finite(Vec3 p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
double seconds(Clock::time_point t) noexcept {
    return std::chrono::duration<double>(Clock::now() - t).count();
}
Vec3 position(Key k, Vec3 lo, Vec3 hi, std::uint32_t extent) noexcept {
    auto coord = [extent](double a, double b, std::uint32_t p) {
        if (!p)
            return a;
        if (p == extent)
            return b;
        return std::lerp(a, b, static_cast<double>(p) / extent);
    };
    return {coord(lo.x, hi.x, k.p[0]), coord(lo.y, hi.y, k.p[1]), coord(lo.z, hi.z, k.p[2])};
}
Key center(const Node &n) noexcept {
    return {{{n.p[0] + n.width / 2, n.p[1] + n.width / 2, n.p[2] + n.width / 2}}};
}
Key upper(const Node &n) noexcept {
    return {{{n.p[0] + n.width, n.p[1] + n.width, n.p[2] + n.width}}};
}
std::uint32_t locate(const std::vector<Node> &tree, const std::array<std::int64_t, 3> &p) noexcept {
    const auto extent = tree[0].width;
    for (auto x : p)
        if (x < 0 || x >= extent)
            return missing;
    std::uint32_t at = 0;
    while (tree[at].child != missing) {
        const auto &n = tree[at];
        unsigned octant = 0;
        for (unsigned a = 0; a < 3; ++a)
            if (p[a] >= n.p[a] + n.width / 2)
                octant |= 1U << a;
        at = n.child + octant;
    }
    return at;
}
template <class F> void patches(const std::vector<Node> &tree, const Node &n, F &&visit) {
    for (unsigned axis = 0; axis < 3; ++axis)
        for (unsigned side = 0; side < 2; ++side) {
            std::array<std::int64_t, 3> p{};
            for (unsigned a = 0; a < 3; ++a)
                p[a] = n.p[a] + n.width / 2;
            p[axis] = side ? static_cast<std::int64_t>(n.p[axis]) + n.width
                           : static_cast<std::int64_t>(n.p[axis]) - 1;
            const auto neighbor = locate(tree, p);
            const unsigned parts = neighbor != missing && tree[neighbor].depth > n.depth ? 2U : 1U;
            const auto width = n.width / parts;
            for (unsigned v = 0; v < parts; ++v)
                for (unsigned u = 0; u < parts; ++u) {
                    Key lo{n.p};
                    lo.p[axis] += side * n.width;
                    lo.p[(axis + 1) % 3] += u * width;
                    lo.p[(axis + 2) % 3] += v * width;
                    visit(Patch{lo, width, axis, side, neighbor == missing});
                }
        }
}
Key patchCenter(const Patch &p) noexcept {
    Key c = p.lower;
    c.p[(p.axis + 1) % 3] += p.width / 2;
    c.p[(p.axis + 2) % 3] += p.width / 2;
    return c;
}
template <class F> void triangles(const Patch &p, const std::vector<Key> &keys, F &&visit) {
    std::array<Key, 4> corners{p.lower, p.lower, p.lower, p.lower};
    const auto u = (p.axis + 1) % 3, v = (p.axis + 2) % 3;
    corners[1].p[u] += p.width;
    corners[2].p[u] += p.width;
    corners[2].p[v] += p.width;
    corners[3].p[v] += p.width;
    if (!p.side)
        std::swap(corners[1], corners[3]);
    const auto c = patchCenter(p);
    for (unsigned e = 0; e < 4; ++e) {
        const auto a = corners[e], b = corners[(e + 1) % 4];
        Key mid;
        for (unsigned i = 0; i < 3; ++i)
            mid.p[i] = (a.p[i] + b.p[i]) / 2;
        if (std::binary_search(keys.begin(), keys.end(), mid)) {
            visit(c, a, mid);
            visit(c, mid, b);
        } else
            visit(c, a, b);
    }
}
Index index(const std::vector<Key> &keys, Key k) {
    const auto i = std::lower_bound(keys.begin(), keys.end(), k);
    if (i == keys.end() || *i != k)
        throw std::logic_error("adaptive node missing");
    return static_cast<Index>(i - keys.begin());
}
mesh::Classification classification(Key k, std::uint32_t extent) noexcept {
    unsigned mask = 0, count = 0, face = 0;
    for (unsigned a = 0; a < 3; ++a) {
        if (k.p[a] == 0) {
            mask |= 1U << (2 * a);
            ++count;
            face = 2 * a + 1;
        } else if (k.p[a] == extent) {
            mask |= 1U << (2 * a + 1);
            ++count;
            face = 2 * a + 2;
        }
    }
    if (count == 0)
        return {3, 1};
    return {3 - static_cast<int>(count),
            count == 1 ? static_cast<int>(face) : 100 + static_cast<int>(mask)};
}
struct CellMeasure {
    double volume = 0, quality = 0;
};
CellMeasure measure(const mesh::Mesh &m, const Index *v) noexcept {
    std::array<Vec3, 4> p{};
    double scale = 0;
    for (unsigned j = 1; j < 4; ++j) {
        p[j] = m.node(v[j]) - m.node(v[0]);
        scale = std::max({scale, std::abs(p[j].x), std::abs(p[j].y), std::abs(p[j].z)});
    }
    if (!(scale > 0) || !std::isfinite(scale))
        return {};
    for (unsigned j = 1; j < 4; ++j)
        p[j] = p[j] / scale;
    const long double det =
        static_cast<long double>(p[1].x) * (static_cast<long double>(p[2].y) * p[3].z -
                                            static_cast<long double>(p[2].z) * p[3].y) -
        static_cast<long double>(p[1].y) * (static_cast<long double>(p[2].x) * p[3].z -
                                            static_cast<long double>(p[2].z) * p[3].x) +
        static_cast<long double>(p[1].z) *
            (static_cast<long double>(p[2].x) * p[3].y - static_cast<long double>(p[2].y) * p[3].x);
    double sum = 0;
    for (unsigned a = 0; a < 4; ++a)
        for (unsigned b = a + 1; b < 4; ++b)
            sum += (p[b] - p[a]).norm2();
    const double volume = static_cast<double>(det / 6 * scale * scale * scale);
    const double quality = det > 0 ? 12 * std::pow(static_cast<double>(det) / 2, 2. / 3) / sum : 0;
    return {volume, quality};
}
} // namespace

AdaptiveBoxResult meshAdaptiveBox(const AdaptiveBoxRequest &request,
                                  mesh::Mesh &destination) noexcept {
    AdaptiveBoxResult out;
    auto fail = [&](AdaptiveBoxReason r) {
        out.reason = r;
        return out;
    };
    const auto begin = Clock::now();
    const Vec3 span = request.upper - request.lower;
    const double largest = std::max({span.x, span.y, span.z});
    if (!finite(request.lower) || !finite(request.upper) || !finite(span) ||
        !(span.x > 0 && span.y > 0 && span.z > 0) || request.maximumDepth > 20 ||
        !std::isfinite(request.targetSize) || !(request.targetSize > 0) ||
        (request.field && !request.field->valid()) || request.threads < 0 ||
        !std::isfinite(request.maximumAdjacentSizeRatio) ||
        (request.maximumAdjacentSizeRatio != 0 && request.maximumAdjacentSizeRatio < 1) ||
        !std::isfinite(request.minimumQuality) || request.minimumQuality < 0 ||
        request.minimumQuality > 1)
        return fail(AdaptiveBoxReason::InvalidRequest);
    out.analyticVolume = span.x * span.y * span.z;
    if (!(out.analyticVolume > 0) || !std::isfinite(out.analyticVolume))
        return fail(AdaptiveBoxReason::NumericallyUnresolved);
    auto cancelled = [&] {
        return request.cancellation && request.cancellation->load(std::memory_order_relaxed);
    };
    if (cancelled())
        return fail(AdaptiveBoxReason::Cancelled);
    if (!request.maximumLeaves)
        return fail(AdaptiveBoxReason::LeafBudget);
    const auto queryCost = request.field ? request.field->boxQueryWorkBound() : 1;
    if (!queryCost)
        return fail(AdaptiveBoxReason::InvalidRequest);
    const std::uint32_t extent = 1U << (request.maximumDepth + 1);
    try {
        std::vector<Node> tree;
        auto working = [&](std::uint64_t bytes) {
            out.workingBytes = std::max(out.workingBytes, bytes);
            return bytes <= request.maximumWorkingBytes;
        };
        auto growTree = [&](std::size_t size) {
            if (size > missing)
                return false;
            if (size <= tree.capacity())
                return true;
            const auto capacity = std::max(size, std::max(std::size_t{1}, tree.capacity() * 2));
            if (capacity > missing || !working((capacity + tree.capacity()) * sizeof(Node)))
                return false;
            tree.reserve(capacity);
            return true;
        };
        if (!growTree(1))
            return fail(AdaptiveBoxReason::WorkingBudget);
        tree.push_back({{}, extent});
        out.leaves = 1;
        auto split = [&](std::uint32_t id) -> AdaptiveBoxReason {
            const auto old = tree[id];
            if (old.depth >= request.maximumDepth)
                return AdaptiveBoxReason::DepthLimit;
            if (request.maximumLeaves - out.leaves < 7)
                return AdaptiveBoxReason::LeafBudget;
            if (!growTree(tree.size() + 8))
                return AdaptiveBoxReason::WorkingBudget;
            const auto first = static_cast<std::uint32_t>(tree.size());
            tree[id].child = first;
            tree[id].marked = false;
            for (unsigned octant = 0; octant < 8; ++octant) {
                Node c;
                c.p = old.p;
                c.width = old.width / 2;
                c.depth = static_cast<std::uint8_t>(old.depth + 1);
                for (unsigned a = 0; a < 3; ++a)
                    if (octant & (1U << a))
                        c.p[a] += c.width;
                tree.push_back(c);
            }
            out.leaves += 7;
            return AdaptiveBoxReason::None;
        };
        for (std::uint32_t id = 0; id < tree.size(); ++id) {
            if ((id & 1023U) == 0 && cancelled())
                return fail(AdaptiveBoxReason::Cancelled);
            const auto n = tree[id];
            const auto lo = position({n.p}, request.lower, request.upper, extent);
            const auto hi = position(upper(n), request.lower, request.upper, extent);
            if (queryCost > UINT64_MAX - out.sizingWork ||
                (request.maximumSizingWork &&
                 queryCost > request.maximumSizingWork - out.sizingWork))
                return fail(AdaptiveBoxReason::SizingBudget);
            out.sizingWork += queryCost;
            const double h =
                request.field ? request.field->minimumOnBox(lo, hi) : request.targetSize;
            ++out.sizingQueries;
            if (!(h > 0) || !std::isfinite(h))
                return fail(AdaptiveBoxReason::InvalidRequest);
            if (largest * static_cast<double>(n.width) / extent > h) {
                const auto reason = split(id);
                if (reason != AdaptiveBoxReason::None)
                    return fail(reason);
            }
        }
        // Every finer leaf probes all touching face/edge/corner neighbors.
        // Marked coarse leaves are split only after the scan, preserving stable
        // node references and deterministic traversal. No spatial hash is used.
        for (;;) {
            if (cancelled())
                return fail(AdaptiveBoxReason::Cancelled);
            bool changed = false;
            const auto count = tree.size();
            for (std::uint32_t id = 0; id < count; ++id) {
                if ((id & 1023U) == 0 && cancelled())
                    return fail(AdaptiveBoxReason::Cancelled);
                const auto &n = tree[id];
                if (n.child != missing)
                    continue;
                for (int z = -1; z <= 1; ++z)
                    for (int y = -1; y <= 1; ++y)
                        for (int x = -1; x <= 1; ++x) {
                            if (!x && !y && !z)
                                continue;
                            const std::array<int, 3> side{x, y, z};
                            std::array<std::int64_t, 3> p{};
                            for (unsigned a = 0; a < 3; ++a)
                                p[a] = static_cast<std::int64_t>(n.p[a]) +
                                       (side[a] < 0   ? std::int64_t{-1}
                                        : side[a] > 0 ? static_cast<std::int64_t>(n.width)
                                                      : n.width / 2);
                            const auto at = locate(tree, p);
                            if (at != missing && tree[at].depth + 1 < n.depth) {
                                tree[at].marked = true;
                                changed = true;
                            }
                        }
            }
            if (!changed)
                break;
            ++out.balancePasses;
            for (std::uint32_t id = 0; id < count; ++id)
                if (tree[id].marked) {
                    const auto reason = split(id);
                    if (reason != AdaptiveBoxReason::None)
                        return fail(reason);
                    ++out.balanceSplits;
                }
        }
        out.treeNodes = tree.size();
        std::vector<Leaf> leaves;
        const auto treeBytes = tree.capacity() * sizeof(Node);
        if (!working(treeBytes + out.leaves * sizeof(Leaf)))
            return fail(AdaptiveBoxReason::WorkingBudget);
        leaves.reserve(static_cast<std::size_t>(out.leaves));
        std::uint64_t keyCount = 9 * out.leaves;
        out.smallestBlockWidth = largest;
        out.largestBlockWidth = 0;
        for (std::uint32_t id = 0; id < tree.size(); ++id)
            if (tree[id].child == missing) {
                leaves.push_back({id});
                out.deepestLevel =
                    std::max(out.deepestLevel, static_cast<std::uint32_t>(tree[id].depth));
                const double width = largest * static_cast<double>(tree[id].width) / extent;
                out.smallestBlockWidth = std::min(out.smallestBlockWidth, width);
                out.largestBlockWidth = std::max(out.largestBlockWidth, width);
                patches(tree, tree[id], [&](const Patch &) { ++keyCount; });
            }
        const auto leafBytes = leaves.capacity() * sizeof(Leaf);
        if (!working(treeBytes + leafBytes + keyCount * sizeof(Key)))
            return fail(AdaptiveBoxReason::WorkingBudget);
        std::vector<Key> keys;
        keys.reserve(static_cast<std::size_t>(keyCount));
        for (const auto &leaf : leaves) {
            if (cancelled())
                return fail(AdaptiveBoxReason::Cancelled);
            const auto &n = tree[leaf.tree];
            keys.push_back(center(n));
            for (unsigned octant = 0; octant < 8; ++octant) {
                Key k{n.p};
                for (unsigned a = 0; a < 3; ++a)
                    if (octant & (1U << a))
                        k.p[a] += n.width;
                keys.push_back(k);
            }
            patches(tree, n, [&](const Patch &p) { keys.push_back(patchCenter(p)); });
        }
        std::sort(keys.begin(), keys.end());
        keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
        out.nodes = keys.size();
        if (cancelled())
            return fail(AdaptiveBoxReason::Cancelled);
        for (auto &leaf : leaves) {
            leaf.volume = out.tetrahedra;
            leaf.boundary = out.boundaryTriangles;
            patches(tree, tree[leaf.tree], [&](const Patch &p) {
                triangles(p, keys, [&](Key, Key, Key) {
                    ++out.tetrahedra;
                    if (p.exterior)
                        ++out.boundaryTriangles;
                });
            });
        }
        const auto cells = out.tetrahedra + out.boundaryTriangles;
        const auto connectivity = 4 * out.tetrahedra + 3 * out.boundaryTriangles;
        out.outputBytes = out.nodes * (sizeof(Vec3) + sizeof(mesh::Classification)) +
                          cells * sizeof(mesh::Cell) + connectivity * sizeof(Index);
        if (out.nodes > request.maximumNodes || cells > request.maximumCells ||
            out.outputBytes > request.maximumOutputBytes ||
            connectivity > static_cast<std::uint64_t>(std::numeric_limits<Index>::max()))
            return fail(AdaptiveBoxReason::OutputBudget);
        out.planningSeconds = seconds(begin);
        const auto generate = Clock::now();
        mesh::detail::CertifiedMeshStorage storage(static_cast<std::size_t>(out.nodes),
                                                   static_cast<std::size_t>(cells),
                                                   static_cast<std::size_t>(connectivity));
        core::parallelForWithThreads(0, keys.size(), request.threads, [&](std::size_t i) {
            storage.setNode(i, position(keys[i], request.lower, request.upper, extent),
                            classification(keys[i], extent));
        });
        if (cancelled())
            return fail(AdaptiveBoxReason::Cancelled);
        core::parallelForWithThreads(0, leaves.size(), request.threads, [&](std::size_t i) {
            if (cancelled())
                return;
            const auto &leaf = leaves[i];
            const auto &n = tree[leaf.tree];
            const Index c = index(keys, center(n));
            auto volume = leaf.volume, boundary = leaf.boundary;
            patches(tree, n, [&](const Patch &p) {
                triangles(p, keys, [&](Key ka, Key kb, Key kc) {
                    const auto a = index(keys, ka), b = index(keys, kb), d = index(keys, kc);
                    storage.setCell(static_cast<std::size_t>(volume),
                                    static_cast<std::size_t>(volume * 4), mesh::CellType::Tet4,
                                    std::array<Index, 4>{c, a, b, d}, {3, 1});
                    ++volume;
                    if (p.exterior) {
                        storage.setCell(static_cast<std::size_t>(out.tetrahedra + boundary),
                                        static_cast<std::size_t>(4 * out.tetrahedra + 3 * boundary),
                                        mesh::CellType::Tri3, std::array<Index, 3>{a, b, d},
                                        {2, static_cast<int>(p.axis * 2 + p.side + 1)});
                        ++boundary;
                    }
                });
            });
        });
        if (cancelled())
            return fail(AdaptiveBoxReason::Cancelled);
        mesh::Mesh staged;
        storage.publishTo(staged);
        // Release planner arrays before geometric validation and the full audit.
        std::vector<Node>().swap(tree);
        std::vector<Leaf>().swap(leaves);
        std::vector<Key>().swap(keys);
        out.generationSeconds = seconds(generate);
        const auto validate = Clock::now();
        const auto volumeBytes = cells * sizeof(double);
        if (!working(volumeBytes))
            return fail(AdaptiveBoxReason::WorkingBudget);
        std::vector<double> volumes(static_cast<std::size_t>(cells), 0);
        out.minimumQuality = 1;
        out.minimumCellVolume = std::numeric_limits<double>::infinity();
        // MSVC's long double has the same precision as double. Compensate the
        // sum so volume telemetry does not drift with the number of cells.
        long double total = 0, compensation = 0;
        for (Index i = 0; i < static_cast<Index>(out.tetrahedra); ++i) {
            if ((i & 1023) == 0 && cancelled())
                return fail(AdaptiveBoxReason::Cancelled);
            const auto m = measure(staged, staged.cellNodes(i));
            if (!(m.volume > 0) || !std::isfinite(m.volume) || !(m.quality > 0) ||
                !std::isfinite(m.quality))
                return fail(AdaptiveBoxReason::InvalidCell);
            out.minimumQuality = std::min(out.minimumQuality, m.quality);
            out.minimumCellVolume = std::min(out.minimumCellVolume, m.volume);
            out.maximumCellVolume = std::max(out.maximumCellVolume, m.volume);
            const long double corrected = static_cast<long double>(m.volume) - compensation;
            const long double next = total + corrected;
            compensation = (next - total) - corrected;
            total = next;
            volumes[static_cast<std::size_t>(i)] = m.volume;
        }
        out.meshVolume = static_cast<double>(total);
        if (std::abs(total / out.analyticVolume - 1) > 1e-10L)
            return fail(AdaptiveBoxReason::NumericallyUnresolved);
        if (out.minimumQuality < request.minimumQuality)
            return fail(AdaptiveBoxReason::QualityRequirement);
        out.validationSeconds = seconds(validate);
        auto limits = request.auditLimits;
        limits.maximumTemporaryBytes =
            std::min(limits.maximumTemporaryBytes, request.maximumWorkingBytes - volumeBytes);
        limits.cellVolumes = volumes;
        if (!limits.cancellation)
            limits.cancellation = request.cancellation;
        const auto audit = Clock::now();
        out.audit = mesh::auditLinearVolumeBoundary(staged, limits);
        out.auditSeconds = seconds(audit);
        out.workingBytes = std::max(out.workingBytes, out.audit.temporaryBytes + volumeBytes);
        if (out.audit.cancelled || cancelled())
            return fail(AdaptiveBoxReason::Cancelled);
        if (out.audit.budgetExceeded)
            return fail(AdaptiveBoxReason::WorkingBudget);
        if (out.audit.allocationFailed)
            return fail(AdaptiveBoxReason::AllocationFailed);
        if (!out.audit.complete || !out.audit.valid)
            return fail(AdaptiveBoxReason::BoundaryAuditFailed);
        if (request.maximumAdjacentSizeRatio > 0 &&
            out.audit.maximumAdjacentSizeRatio > request.maximumAdjacentSizeRatio * (1 + 1e-10))
            return fail(AdaptiveBoxReason::GrowthRequirement);
        destination.swap(staged);
        out.destinationUnchanged = false;
        out.published = true;
        out.reason = AdaptiveBoxReason::None;
    } catch (const std::bad_alloc &) {
        out.reason = AdaptiveBoxReason::AllocationFailed;
    } catch (...) {
        out.reason = AdaptiveBoxReason::InternalError;
    }
    return out;
}
const char *adaptiveBoxReasonName(AdaptiveBoxReason reason) noexcept {
    switch (reason) {
#define PAACS_BOX_REASON(name)                                                                     \
    case AdaptiveBoxReason::name:                                                                  \
        return #name
        PAACS_BOX_REASON(None);
        PAACS_BOX_REASON(InvalidRequest);
        PAACS_BOX_REASON(DepthLimit);
        PAACS_BOX_REASON(LeafBudget);
        PAACS_BOX_REASON(OutputBudget);
        PAACS_BOX_REASON(WorkingBudget);
        PAACS_BOX_REASON(Cancelled);
        PAACS_BOX_REASON(NumericallyUnresolved);
        PAACS_BOX_REASON(InvalidCell);
        PAACS_BOX_REASON(QualityRequirement);
        PAACS_BOX_REASON(BoundaryAuditFailed);
        PAACS_BOX_REASON(AllocationFailed);
        PAACS_BOX_REASON(InternalError);
        PAACS_BOX_REASON(GrowthRequirement);
        PAACS_BOX_REASON(SizingBudget);
#undef PAACS_BOX_REASON
    }
    return "Unknown";
}
} // namespace paacs::algo
