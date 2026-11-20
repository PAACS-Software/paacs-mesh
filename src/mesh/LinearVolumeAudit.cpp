// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/topo/ElementTraits.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

namespace paacs::mesh {
namespace {
// The smallest node is implicit in the bucket. Remaining corners preserve
// cyclic adjacency: a crossed quad must not match a quad with the same set.
template <class NodeId> struct FaceUse {
    std::array<NodeId, 3> key;
    std::uint8_t n;
    bool boundary;
    std::int8_t orientation;
};
template <class NodeId> struct MeasuredFaceUse : FaceUse<NodeId> { double volume; };
static_assert(sizeof(FaceUse<std::uint32_t>) == 16);
static_assert(sizeof(FaceUse<std::uint64_t>) <= 32);
struct Cancelled {};
void poll(const LinearVolumeAuditLimits &limits) {
    if (limits.cancellation && limits.cancellation->load(std::memory_order_relaxed))
        throw Cancelled{};
}
const topo::ElementTraits &cellTraits(CellType t) {
    topo::ElementType e{};
    (void)topo::elementTypeOf(t, e);
    return topo::traits(e);
}
bool supported(CellType t) {
    return t == CellType::Tet4 || t == CellType::Hex8 || t == CellType::Prism6 ||
           t == CellType::Pyramid5 || t == CellType::Tri3 || t == CellType::Quad4 ||
           t == CellType::Line2;
}
// Called only after validating every connectivity row and node index.
template <class F>
void visitFaces(const Mesh &value, const LinearVolumeAuditLimits &limits, F &&visit) {
    for (Index i = 0; i < value.cellCount(); ++i) {
        if ((i & 16383) == 0)
            poll(limits);
        const auto &c = value.cell(i);
        if (c.type == CellType::Line2)
            continue;
        const auto &traits = cellTraits(c.type);
        const auto *conn = value.cellNodes(i);
        if (traits.dimension == 2) {
            std::array<Index, 4> ids{};
            for (std::uint8_t k = 0; k < traits.cornerCount; ++k)
                ids[k] = conn[k];
            visit(ids, traits.cornerCount, true, i);
        } else
            for (std::uint8_t j = 0; j < traits.faceCount; ++j) {
                const auto &f = traits.faces[j];
                std::array<Index, 4> ids{};
                for (std::uint8_t k = 0; k < f.n; ++k)
                    ids[k] = conn[f.v[k]];
                visit(ids, f.n, false, i);
            }
    }
}
void group(LinearVolumeAudit &out, std::uint64_t volumes, std::uint64_t boundaries,
           std::int64_t orientation) {
    if (volumes > 2)
        ++out.nonmanifoldFaces;
    if (boundaries > 1)
        ++out.duplicateBoundary;
    if (volumes == 1) {
        ++out.exteriorFaces;
        if (boundaries == 0)
            ++out.missingBoundary;
        // Boundary normal sense belongs to the geometry ownership provider:
        // surface adapters may orient to a chart rather than the solid.
    } else if (volumes == 2) {
        ++out.interiorFaces;
        if (orientation != 0)
            ++out.inconsistentOrientation;
        if (boundaries)
            ++out.unexpectedBoundary;
    } else if (boundaries)
        ++out.unexpectedBoundary;
}
template <class NodeId, bool Measure = false>
bool bucketAudit(const Mesh &value, const LinearVolumeAuditLimits &limits, std::uint64_t total,
                 LinearVolumeAudit &out) {
    using Face = std::conditional_t<Measure, MeasuredFaceUse<NodeId>, FaceUse<NodeId>>;
    const auto nodes = static_cast<std::uint64_t>(value.nodeCount());
    const auto maximum = std::min<std::uint64_t>(limits.maximumTemporaryBytes,
                                                 std::numeric_limits<std::size_t>::max());
    if (nodes > maximum / sizeof(std::size_t) ||
        total > (maximum - nodes * sizeof(std::size_t)) / sizeof(Face)) {
        out.budgetExceeded = true;
        return false;
    }
    out.temporaryBytes = nodes * sizeof(std::size_t) + total * sizeof(Face);
    std::vector<std::size_t> offsets(static_cast<std::size_t>(nodes), 0);
    visitFaces(value, limits, [&](const auto &ids, std::uint8_t n, bool, Index) {
        const auto least = *std::min_element(ids.begin(), ids.begin() + n);
        ++offsets[static_cast<std::size_t>(least)];
    });
    std::size_t prefix = 0;
    for (std::size_t i = 0; i < offsets.size(); ++i) {
        if ((i & 16383U) == 0)
            poll(limits);
        const auto count = offsets[i];
        out.maximumBucketSize = std::max<std::uint64_t>(out.maximumBucketSize, count);
        offsets[i] = prefix;
        prefix += count;
    }
    // Every slot is assigned before use; avoid zeroing a second large buffer.
    auto faces = std::make_unique_for_overwrite<Face[]>(static_cast<std::size_t>(total));
    visitFaces(value, limits, [&](const auto &ids, std::uint8_t n, bool boundary, Index cell) {
        unsigned first = 0;
        for (unsigned i = 1; i < n; ++i)
            if (ids[i] < ids[first])
                first = i;
        const auto next = ids[first + 1 == n ? 0 : first + 1];
        const auto prev = ids[first == 0 ? n - 1 : first - 1];
        const bool forward = next < prev;
        Face f{};
        f.n=n; f.boundary=boundary; f.orientation=static_cast<std::int8_t>(forward ? 1 : -1);
        if constexpr (Measure) f.volume=boundary ? 0 : limits.cellVolumes[static_cast<std::size_t>(cell)];
        f.key[0] = static_cast<NodeId>(forward ? next : prev);
        if (n == 3) {
            f.key[1] = static_cast<NodeId>(forward ? prev : next);
        } else {
            f.key[1] = static_cast<NodeId>(ids[first ^ 2U]);
            f.key[2] = static_cast<NodeId>(forward ? prev : next);
        }
        faces[offsets[static_cast<std::size_t>(ids[first])]++] = f;
    });
    // The scatter advanced the exclusive starts into exclusive ends. No
    // extra cursor array is needed. Empty and unused-node buckets are legal.
    std::size_t start = 0;
    for (std::size_t b = 0; b < offsets.size(); ++b) {
        if ((b & 1023U) == 0)
            poll(limits);
        const auto end = offsets[b];
        std::sort(faces.get() + start, faces.get() + end, [](const Face &a, const Face &c) {
            return a.n != c.n ? a.n < c.n : a.key < c.key;
        });
        for (std::size_t i = start; i < end;) {
            std::size_t j = i;
            std::uint64_t volumes = 0, boundaries = 0;
            std::int64_t orientation = 0;
            double minimumVolume=std::numeric_limits<double>::infinity(), maximumVolume=0;
            do {
                const auto &f = faces[j];
                if (f.boundary)
                    ++boundaries;
                else {
                    ++volumes;
                    orientation += f.orientation;
                    if constexpr (Measure) {
                        minimumVolume=std::min(minimumVolume,f.volume);
                        maximumVolume=std::max(maximumVolume,f.volume);
                    }
                }
                ++j;
            } while (j < end && faces[i].n == faces[j].n && faces[i].key == faces[j].key);
            group(out, volumes, boundaries, orientation);
            if constexpr (Measure) if (volumes==2)
                out.maximumAdjacentVolumeRatio=std::max(out.maximumAdjacentVolumeRatio,maximumVolume/minimumVolume);
            i = j;
        }
        start = end;
    }
    poll(limits);
    return true;
}
} // namespace
LinearVolumeAudit auditLinearVolumeBoundary(const Mesh &value,
                                            const LinearVolumeAuditLimits &limits) noexcept {
    LinearVolumeAudit out;
    try {
        poll(limits);
        const bool measure=!limits.cellVolumes.empty();
        if(measure && limits.cellVolumes.size()!=static_cast<std::size_t>(value.cellCount())) {
            out.invalidVolumeData=1; out.complete=true; return out;
        }
        for (Index i = 0; i < value.nodeCount(); ++i) {
            if ((i & 16383) == 0)
                poll(limits);
            const auto p = value.node(i);
            if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
                ++out.nonfiniteNodes;
        }
        std::uint64_t total = 0;
        for (Index i = 0; i < value.cellCount(); ++i) {
            if ((i & 16383) == 0)
                poll(limits);
            const auto &c = value.cell(i);
            if (!supported(c.type)) {
                ++out.unsupportedCells;
                continue;
            }
            if (c.type == CellType::Line2)
                continue;
            const auto &traits = cellTraits(c.type);
            const bool surface = traits.dimension == 2;
            const auto add = surface ? 1U : traits.faceCount;
            if (total > limits.maximumFaceUses || add > limits.maximumFaceUses - total) {
                out.budgetExceeded = true;
                return out;
            }
            total += add;
            if (surface)
                ++out.boundaryCells;
            else
                ++out.volumeCells;
            if(measure && !surface && (!(limits.cellVolumes[static_cast<std::size_t>(i)]>0) ||
                !std::isfinite(limits.cellVolumes[static_cast<std::size_t>(i)]))) ++out.invalidVolumeData;
            const auto n = static_cast<Index>(traits.nodeCount);
            if (c.first < 0 || c.first > value.connectivityCount() ||
                n > value.connectivityCount() - c.first) {
                ++out.invalidConnectivity;
                continue;
            }
            const auto *conn = value.cellNodes(i);
            bool invalid = false;
            for (Index j = 0; j < n; ++j) {
                if (conn[j] < 0 || conn[j] >= value.nodeCount())
                    invalid = true;
                for (Index k = 0; k < j; ++k)
                    if (conn[k] == conn[j])
                        invalid = true;
            }
            if (invalid)
                ++out.invalidConnectivity;
        }
        if (out.invalidConnectivity || out.unsupportedCells || out.nonfiniteNodes || out.invalidVolumeData) {
            out.complete = true;
            return out;
        }
        const bool narrow=static_cast<std::uint64_t>(value.nodeCount())<=UINT32_MAX;
        if (total && !(measure ? (narrow ? bucketAudit<std::uint32_t,true>(value,limits,total,out)
                                         : bucketAudit<std::uint64_t,true>(value,limits,total,out))
                              : (narrow ? bucketAudit<std::uint32_t>(value,limits,total,out)
                                        : bucketAudit<std::uint64_t>(value,limits,total,out))))
            return out;
        out.faceUses = total;
        out.complete = true;
        out.volumeGrowthMeasured=measure;
        if(measure) out.maximumAdjacentSizeRatio=std::cbrt(out.maximumAdjacentVolumeRatio);
        out.valid = out.volumeCells > 0 && !out.nonmanifoldFaces && !out.inconsistentOrientation &&
                    !out.missingBoundary && !out.unexpectedBoundary && !out.duplicateBoundary;
    } catch (const Cancelled &) {
        out.cancelled = true;
    } catch (const std::bad_alloc &) {
        out.allocationFailed = true;
    } catch (...) {
        out.complete = false;
        out.valid = false;
    }
    return out;
}
} // namespace paacs::mesh
