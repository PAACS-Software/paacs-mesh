// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/TetToHexMesher.hpp"
#include "PrimitiveMesherInternal.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/mesh/CertifiedMeshStorage.hpp"
#include "paacs/mesh/Hex8Jacobian.hpp"
#include <cmath>
#include <limits>
#include <new>
#include <numeric>

namespace paacs::algo {
namespace {
struct Stop {
    PrimitiveReason reason;
};
void need(bool ok, PrimitiveReason reason) {
    if (!ok)
        throw Stop{reason};
}
using Edge = std::array<Index, 2>;
using Face = std::array<Index, 3>;
Edge edge(Index a, Index b) {
    return a < b ? Edge{a, b} : Edge{b, a};
}
Face face(Index a, Index b, Index c) {
    Face f{a, b, c};
    std::sort(f.begin(), f.end());
    return f;
}
int orientation(Index a, Index b, Index c) {
    return ((a > b) + (a > c) + (b > c)) % 2 == 0 ? 1 : -1;
}
constexpr std::array<std::array<unsigned, 3>, 4> tetFaces{
    {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}}};
constexpr std::array<std::array<unsigned, 4>, 4> orders{
    {{0, 1, 2, 3}, {1, 0, 3, 2}, {2, 3, 0, 1}, {3, 2, 1, 0}}};
constexpr std::array<std::array<unsigned, 4>, 6> sub{
    {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}}};
struct Metrics {
    long double volume = 0;
    double minimumVolume = std::numeric_limits<double>::infinity(), minimumQuality = 1;
    PrimitiveQualityReport quality{};
};
} // namespace
TetToHexResult convertTet4ToHex8(const mesh::Mesh &source, const TetToHexOptions &o,
                                 mesh::Mesh &destination) noexcept {
    TetToHexResult r;
    try {
        auto cancel = [&] {
            need(!(o.cancellation && o.cancellation->load(std::memory_order_relaxed)) &&
                     !(o.auditLimits.cancellation &&
                       o.auditLimits.cancellation->load(std::memory_order_relaxed)),
                 PrimitiveReason::Cancelled);
        };
        auto work = [&](long double bytes) {
            need(bytes <= o.maximumWorkingBytes && bytes <= SIZE_MAX,
                 PrimitiveReason::WorkingBudget);
            r.workingBytes = std::max(r.workingBytes, static_cast<std::uint64_t>(bytes));
        };
        cancel();
        need(o.threads >= 1 && o.threads <= 1024 && std::isfinite(o.minimumQuality) &&
                 o.minimumQuality >= 0 && o.minimumQuality <= 1 &&
                 std::isfinite(o.warningQuality) && o.warningQuality >= 0 &&
                 o.warningQuality <= 1 &&
                 static_cast<unsigned>(o.qualityPolicy) <=
                     static_cast<unsigned>(PrimitiveQualityPolicy::Reject),
             PrimitiveReason::InvalidRequest);
        std::uint64_t nt = 0, nb = 0;
        mesh::Classification volume{};
        for (Index i = 0; i < source.cellCount(); ++i) {
            if (i % 1024 == 0)
                cancel();
            const auto c = source.cell(i);
            need((c.type == mesh::CellType::Tet4 && c.onEntity.dim == 3) ||
                     (c.type == mesh::CellType::Tri3 && c.onEntity.dim == 2),
                 PrimitiveReason::InvalidRequest);
            need(c.onEntity.tag >= 0, PrimitiveReason::InvalidRequest);
            if (c.type == mesh::CellType::Tet4) {
                if (nt == 0)
                    volume = c.onEntity;
                need(volume == c.onEntity, PrimitiveReason::InvalidRequest);
                ++nt;
            } else
                ++nb;
        }
        need(nt > 0 && nb > 0, PrimitiveReason::InvalidRequest);
        const auto nn = static_cast<std::uint64_t>(source.nodeCount());
        auto limits = o.auditLimits;
        limits.cellVolumes = {};
        if (o.cancellation)
            limits.cancellation = o.cancellation;
        limits.maximumTemporaryBytes =
            std::min(limits.maximumTemporaryBytes, o.maximumWorkingBytes);
        const auto inputAudit = mesh::auditLinearVolumeBoundary(source, limits);
        work(inputAudit.temporaryBytes);
        need(!inputAudit.cancelled, PrimitiveReason::Cancelled);
        need(!inputAudit.allocationFailed, PrimitiveReason::AllocationFailed);
        need(!inputAudit.budgetExceeded, PrimitiveReason::WorkingBudget);
        need(inputAudit.valid, PrimitiveReason::BoundaryAuditFailed);
        const long double keys =
            static_cast<long double>(nt) * (6 * sizeof(Edge) + 4 * sizeof(Face)) + nn;
        work(keys);
        std::vector<std::uint8_t> used(static_cast<std::size_t>(nn));
        std::vector<Edge> edges;
        edges.reserve(static_cast<std::size_t>(6 * nt));
        std::vector<Face> faces;
        faces.reserve(static_cast<std::size_t>(4 * nt));
        for (Index i = 0; i < source.cellCount(); ++i)
            if (source.cell(i).type == mesh::CellType::Tet4) {
                if (i % 1024 == 0)
                    cancel();
                const auto *v = source.cellNodes(i);
                const Vec3 a = source.node(v[0]), b = source.node(v[1]), c = source.node(v[2]),
                           d = source.node(v[3]);
                const double vol = (b - a).cross(c - a).dot(d - a) / 6;
                need(std::isfinite(vol) && vol > 0, PrimitiveReason::InvalidCell);
                for (unsigned j = 0; j < 4; ++j) {
                    used[static_cast<std::size_t>(v[j])] = 1;
                    for (unsigned k = j + 1; k < 4; ++k)
                        edges.push_back(edge(v[j], v[k]));
                    faces.push_back(face(v[j], v[(j + 1) % 4], v[(j + 2) % 4]));
                }
            }
        for (Index i = 0; i < source.nodeCount(); ++i) {
            const auto c = source.nodeClass(i);
            need(used[static_cast<std::size_t>(i)] && c.dim >= 0 && c.dim <= 3 && c.tag >= 0,
                 PrimitiveReason::InvalidRequest);
        }
        cancel();
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
        cancel();
        std::sort(faces.begin(), faces.end());
        faces.erase(std::unique(faces.begin(), faces.end()), faces.end());
        const auto ne = edges.size(), nf = faces.size();
        const auto edgeAt = [&](Index a, Index b) {
            const auto key = edge(a, b);
            const auto it = std::lower_bound(edges.begin(), edges.end(), key);
            need(it != edges.end() && *it == key, PrimitiveReason::InternalError);
            return static_cast<std::size_t>(it - edges.begin());
        };
        const auto faceAt = [&](Index a, Index b, Index c) {
            const auto key = face(a, b, c);
            const auto it = std::lower_bound(faces.begin(), faces.end(), key);
            need(it != faces.end() && *it == key, PrimitiveReason::InternalError);
            return static_cast<std::size_t>(it - faces.begin());
        };
        r.nodes = nn + ne + nf + nt;
        r.hexes = 4 * nt;
        r.boundaryQuads = 3 * nb;
        const auto cells = r.hexes + r.boundaryQuads, conn = 8 * r.hexes + 4 * r.boundaryQuads;
        const long double output =
            static_cast<long double>(r.nodes) * (sizeof(Vec3) + sizeof(mesh::Classification)) +
            static_cast<long double>(cells) * sizeof(mesh::Cell) +
            static_cast<long double>(conn) * sizeof(Index);
        need(r.nodes <= o.maximumNodes && r.nodes <= INT32_MAX && cells <= o.maximumCells &&
                 conn <= INT32_MAX && output <= o.maximumOutputBytes && output <= SIZE_MAX,
             PrimitiveReason::OutputBudget);
        r.outputBytes = static_cast<std::uint64_t>(output);
        const auto partsCount = (r.hexes + 1023) / 1024;
        const long double retained =
            keys + static_cast<long double>(ne + nf) * sizeof(mesh::Classification) +
            partsCount * sizeof(Metrics) + nf * sizeof(std::int8_t);
        work(retained);
        std::vector<mesh::Classification> edgeClass(ne), faceClass(nf, volume);
        // The general incidence audit checks interior volume-face winding;
        // it deliberately does not require outward exterior face winding.
        // Conversion has that stronger contract, so verify it explicitly.
        std::vector<std::int8_t> boundaryOrientation(nf);
        for (Index i = 0; i < source.cellCount(); ++i)
            if (source.cell(i).type == mesh::CellType::Tet4) {
                if (i % 1024 == 0)
                    cancel();
                const auto *v = source.cellNodes(i);
                for (const auto &f : tetFaces) {
                    const Index a = v[f[0]], b = v[f[1]], c = v[f[2]];
                    auto &value = boundaryOrientation[faceAt(a, b, c)];
                    value = static_cast<std::int8_t>(value + orientation(a, b, c));
                }
            }
        for (Index i = 0; i < source.cellCount(); ++i)
            if (source.cell(i).type == mesh::CellType::Tri3) {
                if (i % 1024 == 0)
                    cancel();
                const auto *v = source.cellNodes(i);
                const auto c = source.cell(i).onEntity;
                const auto faceIndex = faceAt(v[0], v[1], v[2]);
                need(boundaryOrientation[faceIndex] == orientation(v[0], v[1], v[2]),
                     PrimitiveReason::BoundaryAuditFailed);
                faceClass[faceIndex] = c;
                for (unsigned j = 0; j < 3; ++j) {
                    const auto a = v[j], b = v[(j + 1) % 3];
                    auto &cls = edgeClass[edgeAt(a, b)];
                    if (cls.dim < 0)
                        cls = c;
                    else if (cls.dim == 2 && cls != c) {
                        const auto na = source.nodeClass(a), nclass = source.nodeClass(b);
                        need(na == nclass && na.dim <= 1, PrimitiveReason::InvalidRequest);
                        cls = na;
                    }
                }
            }
        mesh::detail::CertifiedMeshStorage storage(r.nodes, cells, conn);
        for (Index i = 0; i < source.nodeCount(); ++i)
            storage.setNode(static_cast<std::size_t>(i), source.node(i), source.nodeClass(i));
        for (std::size_t i = 0; i < ne; ++i) {
            if (i % 1024 == 0)
                cancel();
            const auto a = source.node(edges[i][0]), b = source.node(edges[i][1]);
            storage.setNode(
                nn + i, {std::midpoint(a.x, b.x), std::midpoint(a.y, b.y), std::midpoint(a.z, b.z)},
                edgeClass[i].dim < 0 ? volume : edgeClass[i]);
        }
        for (std::size_t i = 0; i < nf; ++i) {
            if (i % 1024 == 0)
                cancel();
            const auto a = source.node(faces[i][0]), b = source.node(faces[i][1]),
                       c = source.node(faces[i][2]);
            storage.setNode(nn + ne + i, a + ((b - a) + (c - a)) / 3, faceClass[i]);
        }
        auto en = [&](Index a, Index b) { return static_cast<Index>(nn + edgeAt(a, b)); };
        auto fn = [&](Index a, Index b, Index c) {
            return static_cast<Index>(nn + ne + faceAt(a, b, c));
        };
        std::uint64_t t = 0, bound = 0;
        for (Index i = 0; i < source.cellCount(); ++i) {
            if (i % 1024 == 0)
                cancel();
            const auto *v = source.cellNodes(i);
            if (source.cell(i).type == mesh::CellType::Tet4) {
                const Index center = static_cast<Index>(nn + ne + nf + t);
                const auto p = source.node(v[0]);
                storage.setNode(static_cast<std::size_t>(center),
                                p + ((source.node(v[1]) - p) + (source.node(v[2]) - p) +
                                     (source.node(v[3]) - p)) /
                                        4,
                                volume);
                for (unsigned j = 0; j < 4; ++j) {
                    const auto &ord = orders[j];
                    const Index a = v[ord[0]], b = v[ord[1]], c = v[ord[2]], d = v[ord[3]];
                    const std::array<Index, 8> h{a,        en(a, b),    fn(a, b, c), en(a, c),
                                                 en(a, d), fn(a, b, d), center,      fn(a, c, d)};
                    storage.setCell(4 * t + j, 8 * (4 * t + j), mesh::CellType::Hex8, h, volume);
                }
                ++t;
            } else {
                for (unsigned j = 0; j < 3; ++j) {
                    const std::array<Index, 4> q{v[j], en(v[j], v[(j + 1) % 3]),
                                                 fn(v[0], v[1], v[2]), en(v[(j + 2) % 3], v[j])};
                    storage.setCell(r.hexes + bound, 8 * r.hexes + 4 * bound, mesh::CellType::Quad4,
                                    q, source.cell(i).onEntity);
                    ++bound;
                }
            }
        }
        mesh::Mesh staged;
        storage.publishTo(staged);
        std::vector<Metrics> parts(static_cast<std::size_t>(partsCount));
        core::parallelChunksWithThreads(
            0, r.hexes, 1024, o.threads, [&](std::size_t part, std::size_t lo, std::size_t hi) {
                cancel();
                auto &m = parts[part];
                for (std::size_t i = lo; i < hi; ++i) {
                    const auto *v = staged.cellNodes(static_cast<Index>(i));
                    std::array<Vec3, 8> p{};
                    bool sharp = false;
                    for (unsigned j = 0; j < 8; ++j) {
                        p[j] = staged.node(v[j]);
                        sharp = sharp || staged.nodeClass(v[j]).dim <= 1;
                    }
                    need(mesh::positiveHex8Jacobian(p), PrimitiveReason::NumericallyUnresolved);
                    long double cellVolume = 0;
                    double q = 1;
                    for (auto a : sub) {
                        const Vec3 ab = p[a[1]] - p[a[0]], ac = p[a[2]] - p[a[0]],
                                   ad = p[a[3]] - p[a[0]];
                        const double vol = ab.cross(ac).dot(ad) / 6;
                        need(vol > 0 && std::isfinite(vol), PrimitiveReason::InvalidCell);
                        const double sum =
                            ab.norm2() + ac.norm2() + ad.norm2() + (p[a[2]] - p[a[1]]).norm2() +
                            (p[a[3]] - p[a[1]]).norm2() + (p[a[3]] - p[a[2]]).norm2();
                        const double quality = 12 * std::pow(3 * vol, 2. / 3) / sum;
                        need(std::isfinite(quality) && quality > 0,
                             PrimitiveReason::NumericallyUnresolved);
                        q = std::min(q, quality);
                        cellVolume += vol;
                    }
                    m.volume += cellVolume;
                    m.minimumVolume = std::min(m.minimumVolume, static_cast<double>(cellVolume));
                    m.minimumQuality = std::min(m.minimumQuality, q);
                    primitive_detail::recordQuality(m.quality, static_cast<Index>(i), q, sharp,
                                                    o.warningQuality);
                }
            });
        long double total = 0;
        r.minimumQuality = 1;
        r.minimumCellVolume = std::numeric_limits<double>::infinity();
        for (const auto &m : parts) {
            total += m.volume;
            r.minimumQuality = std::min(r.minimumQuality, m.minimumQuality);
            r.minimumCellVolume = std::min(r.minimumCellVolume, m.minimumVolume);
            primitive_detail::mergeQuality(r.quality, m.quality);
        }
        r.volume = static_cast<double>(total);
        need(
            r.minimumQuality >= o.minimumQuality &&
                (o.qualityPolicy == PrimitiveQualityPolicy::Warn || r.quality.lowQualityCells == 0),
            PrimitiveReason::QualityRequirement);
        const long double auditBytes =
            16.L * (6 * r.hexes + r.boundaryQuads) + r.nodes * sizeof(std::size_t);
        work(retained + auditBytes);
        limits.maximumTemporaryBytes =
            std::min(limits.maximumTemporaryBytes,
                     o.maximumWorkingBytes - static_cast<std::uint64_t>(retained));
        r.audit = mesh::auditLinearVolumeBoundary(staged, limits);
        need(!r.audit.cancelled, PrimitiveReason::Cancelled);
        need(!r.audit.allocationFailed, PrimitiveReason::AllocationFailed);
        need(!r.audit.budgetExceeded, PrimitiveReason::WorkingBudget);
        need(r.audit.valid, PrimitiveReason::BoundaryAuditFailed);
        cancel();
        destination.swap(staged);
        r.published = true;
        r.destinationUnchanged = false;
        r.reason = PrimitiveReason::None;
    } catch (const Stop &e) {
        r.reason = e.reason;
    } catch (const std::bad_alloc &) {
        r.reason = PrimitiveReason::AllocationFailed;
    } catch (...) {
        r.reason = PrimitiveReason::InternalError;
    }
    return r;
}
} // namespace paacs::algo
