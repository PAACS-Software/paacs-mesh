// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/ConvexToTetMesher.hpp"
#include "paacs/mesh/CertifiedMeshStorage.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/topo/ElementTraits.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <new>
#include <vector>

namespace paacs::algo {
namespace {
using R = ConvexToTetReason;
struct Stop {
    R reason;
};
void need(bool value, R reason) {
    if (!value)
        throw Stop{reason};
}
struct FaceUse {
    std::array<Index, 4> key{};
    std::int8_t sign = 0;
    std::uint8_t count = 0;
    bool exterior = false;
};
FaceUse faceUse(std::array<Index, 4> ids, unsigned count, bool exterior) {
    need(count == 3 || count == 4, R::InvalidSource);
    const auto smallest = std::min_element(ids.begin(), ids.begin() + count);
    const auto at = static_cast<unsigned>(smallest - ids.begin());
    const std::int8_t sign = ids[(at + 1) % count] < ids[(at + count - 1) % count] ? 1 : -1;
    if (count == 3)
        ids[3] = kInvalidIndex;
    std::sort(ids.begin(), ids.end());
    return {ids, sign, static_cast<std::uint8_t>(count), exterior};
}
template <class Emit> void triangles(const std::array<Index, 4> &ids, unsigned count, Emit emit) {
    unsigned at = 0;
    for (unsigned j = 1; j < count; ++j)
        if (ids[j] < ids[at])
            at = j;
    for (unsigned j = 1; j + 1 < count; ++j)
        emit(std::array<Index, 3>{ids[at], ids[(at + j) % count], ids[(at + j + 1) % count]});
}
double magnitude(Vec3 p) {
    return std::hypot(p.x, p.y, p.z);
}
const topo::ElementTraits &traits(mesh::CellType type) {
    topo::ElementType element{};
    need(topo::elementTypeOf(type, element), R::InvalidSource);
    return topo::traits(element);
}
} // namespace

ConvexToTetResult convertConvexCellsToTet4(const mesh::Mesh &source, const ConvexToTetOptions &o,
                                           mesh::Mesh &destination) noexcept {
    ConvexToTetResult r;
    try {
        const auto cancel = [&] {
            need(!(o.cancellation && o.cancellation->load(std::memory_order_relaxed)) &&
                     !(o.auditLimits.cancellation &&
                       o.auditLimits.cancellation->load(std::memory_order_relaxed)),
                 R::Cancelled);
        };
        const auto work = [&](long double bytes) {
            need(bytes <= o.maximumWorkingBytes && bytes <= SIZE_MAX, R::WorkingBudget);
            r.workingBytes = std::max(r.workingBytes, static_cast<std::uint64_t>(bytes));
        };
        cancel();
        need(std::isfinite(o.relativePlanarityTolerance) && o.relativePlanarityTolerance > 0 &&
                 o.relativePlanarityTolerance <= 1e-6 && std::isfinite(o.minimumQuality) &&
                 o.minimumQuality >= 0 && o.minimumQuality <= 1 &&
                 std::isfinite(o.warningQuality) && o.warningQuality >= 0 && o.warningQuality <= 1,
             R::InvalidRequest);
        std::uint64_t volumeCells = 0, exterior = 0, faceUses = 0;
        mesh::Classification region{};
        for (Index i = 0; i < source.cellCount(); ++i) {
            if (i % 1024 == 0)
                cancel();
            const auto c = source.cell(i);
            need(c.onEntity.tag >= 0, R::InvalidSource);
            if (c.type == mesh::CellType::Tet4 || c.type == mesh::CellType::Prism6 ||
                c.type == mesh::CellType::Hex8) {
                need(c.onEntity.dim == 3, R::InvalidSource);
                if (volumeCells == 0)
                    region = c.onEntity;
                need(c.onEntity == region, R::InvalidSource);
                ++volumeCells;
                r.tetrahedra += c.type == mesh::CellType::Tet4     ? 1
                                : c.type == mesh::CellType::Prism6 ? 3
                                                                   : 6;
                faceUses += traits(c.type).faceCount;
            } else {
                need((c.type == mesh::CellType::Tri3 || c.type == mesh::CellType::Quad4) &&
                         c.onEntity.dim == 2,
                     R::InvalidSource);
                ++exterior;
                ++faceUses;
                r.boundaryTriangles += c.type == mesh::CellType::Tri3 ? 1 : 2;
            }
        }
        need(volumeCells > 0 && exterior > 0, R::InvalidSource);
        const auto nodes = static_cast<std::uint64_t>(source.nodeCount());
        const auto cells = r.tetrahedra + r.boundaryTriangles;
        const auto conn = 4 * r.tetrahedra + 3 * r.boundaryTriangles;
        const long double payload =
            static_cast<long double>(nodes) * (sizeof(Vec3) + sizeof(mesh::Classification)) +
            static_cast<long double>(cells) * sizeof(mesh::Cell) +
            static_cast<long double>(conn) * sizeof(Index);
        need(nodes <= o.maximumNodes && cells <= o.maximumCells && cells <= INT32_MAX &&
                 conn <= INT32_MAX && payload <= o.maximumOutputBytes && payload <= SIZE_MAX,
             R::OutputBudget);
        r.outputBytes = static_cast<std::uint64_t>(payload);
        auto limits = o.auditLimits;
        limits.cellVolumes = {};
        limits.maximumTemporaryBytes =
            std::min(limits.maximumTemporaryBytes, o.maximumWorkingBytes);
        if (o.cancellation)
            limits.cancellation = o.cancellation;
        const auto audit = [&](const mesh::Mesh &m) {
            auto a = mesh::auditLinearVolumeBoundary(m, limits);
            work(a.temporaryBytes);
            need(!a.cancelled, R::Cancelled);
            need(!a.allocationFailed, R::AllocationFailed);
            need(!a.budgetExceeded, R::WorkingBudget);
            need(a.valid, R::BoundaryAuditFailed);
            return a;
        };
        for (Index i = 0; i < source.nodeCount(); ++i) {
            if (i % 1024 == 0)
                cancel();
            const auto p = source.node(i);
            const auto c = source.nodeClass(i);
            need(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z) && c.dim >= 0 &&
                     c.dim <= 3 && c.tag >= 0,
                 R::InvalidSource);
        }
        (void)audit(source); // validates connectivity before cell-node indexing
        {
            work(static_cast<long double>(faceUses) * sizeof(FaceUse) + nodes);
            std::vector<FaceUse> uses;
            uses.reserve(static_cast<std::size_t>(faceUses));
            std::vector<std::uint8_t> used(static_cast<std::size_t>(nodes));
            for (Index i = 0; i < source.cellCount(); ++i) {
                if (i % 1024 == 0)
                    cancel();
                const auto c = source.cell(i);
                const auto *ids = source.cellNodes(i);
                const auto &t = traits(c.type);
                if (t.dimension == 2) {
                    std::array<Index, 4> f{};
                    std::copy_n(ids, mesh::nodesPerCell(c.type), f.begin());
                    uses.push_back(
                        faceUse(f, static_cast<unsigned>(mesh::nodesPerCell(c.type)), true));
                    continue;
                }
                const unsigned n = static_cast<unsigned>(mesh::nodesPerCell(c.type));
                double scale = 0, coordinate = 0;
                for (unsigned j = 0; j < n; ++j) {
                    used[static_cast<std::size_t>(ids[j])] = 1;
                    const auto p = source.node(ids[j]);
                    coordinate =
                        std::max({coordinate, std::abs(p.x), std::abs(p.y), std::abs(p.z)});
                    for (unsigned k = 0; k < j; ++k)
                        scale = std::max(scale, magnitude(p - source.node(ids[k])));
                }
                need(std::isfinite(scale) && scale > 0, R::InvalidSource);
                const double tolerance =
                    o.relativePlanarityTolerance +
                    64 * std::numeric_limits<double>::epsilon() * coordinate / scale;
                need(std::isfinite(tolerance) && tolerance < 1e-4, R::InvalidSource);
                for (unsigned j = 0; j < t.faceCount; ++j) {
                    const auto &f = t.faces[j];
                    std::array<Index, 4> face{};
                    for (unsigned k = 0; k < f.n; ++k)
                        face[k] = ids[f.v[k]];
                    uses.push_back(faceUse(face, f.n, false));
                    const auto origin = source.node(face[0]);
                    const auto a = (source.node(face[1]) - origin) / scale;
                    const auto b = (source.node(face[2]) - origin) / scale;
                    auto normal = a.cross(b);
                    const auto norm = magnitude(normal);
                    need(std::isfinite(norm) && norm > 0, R::NonConvexCell);
                    normal = normal / norm;
                    if (f.n == 4)
                        need(std::abs(((source.node(face[3]) - origin) / scale).dot(normal)) <=
                                 tolerance,
                             R::NonplanarFace);
                    for (unsigned k = 0; k < f.n; ++k) {
                        const auto p = source.node(face[k]), q = source.node(face[(k + 1) % f.n]),
                                   s = source.node(face[(k + 2) % f.n]);
                        need(((q - p) / scale).cross((s - q) / scale).dot(normal) > 0,
                             R::NonConvexCell);
                    }
                    for (unsigned k = 0; k < n; ++k)
                        if (std::find(face.begin(), face.begin() + f.n, ids[k]) ==
                            face.begin() + f.n)
                            need(((source.node(ids[k]) - origin) / scale).dot(normal) < -tolerance,
                                 R::NonConvexCell);
                }
            }
            need(std::find(used.begin(), used.end(), 0) == used.end(), R::InvalidSource);
            cancel();
            std::sort(uses.begin(), uses.end(), [](const auto &a, const auto &b) {
                if (a.count != b.count)
                    return a.count < b.count;
                return a.key < b.key;
            });
            for (std::size_t i = 0; i < uses.size();) {
                if (i % 1024 == 0)
                    cancel();
                std::size_t end = i + 1;
                while (end < uses.size() && uses[end].count == uses[i].count &&
                       uses[end].key == uses[i].key)
                    ++end;
                need(end == i + 2, R::BoundaryAuditFailed);
                const auto &a = uses[i], &b = uses[i + 1];
                if (a.exterior || b.exterior)
                    need(a.exterior != b.exterior && a.sign == b.sign, R::BoundaryAuditFailed);
                else
                    need(a.sign != b.sign, R::BoundaryAuditFailed);
                i = end;
            }
        }
        cancel();
        mesh::detail::CertifiedMeshStorage storage(nodes, cells, conn);
        for (Index i = 0; i < source.nodeCount(); ++i) {
            if (i % 1024 == 0)
                cancel();
            storage.setNode(static_cast<std::size_t>(i), source.node(i), source.nodeClass(i));
        }
        std::uint64_t tet = 0, boundary = 0;
        long double volume = 0;
        for (Index i = 0; i < source.cellCount(); ++i) {
            if (i % 1024 == 0)
                cancel();
            const auto c = source.cell(i);
            const auto *ids = source.cellNodes(i);
            const auto &t = traits(c.type);
            if (t.dimension == 2) {
                std::array<Index, 4> f{};
                const auto n = static_cast<unsigned>(mesh::nodesPerCell(c.type));
                std::copy_n(ids, n, f.begin());
                triangles(f, n, [&](const auto &tri) {
                    storage.setCell(r.tetrahedra + boundary, 4 * r.tetrahedra + 3 * boundary,
                                    mesh::CellType::Tri3, tri, c.onEntity);
                    ++boundary;
                });
            } else {
                const Index apex = *std::min_element(ids, ids + mesh::nodesPerCell(c.type));
                for (unsigned j = 0; j < t.faceCount; ++j) {
                    const auto &f = t.faces[j];
                    std::array<Index, 4> face{};
                    for (unsigned k = 0; k < f.n; ++k)
                        face[k] = ids[f.v[k]];
                    if (std::find(face.begin(), face.begin() + f.n, apex) != face.begin() + f.n)
                        continue;
                    triangles(face, f.n, [&](const auto &tri) {
                        const auto p = source.node(apex), a = source.node(tri[0]) - p,
                                   b = source.node(tri[1]) - p, d = source.node(tri[2]) - p;
                        const double scale = std::max({magnitude(a), magnitude(b), magnitude(d)});
                        const auto x = a / scale, y = b / scale, z = d / scale;
                        const double det = x.cross(y).dot(z);
                        const double v = ((det * scale) * scale) * scale / 6;
                        need(std::isfinite(v) && v > 0 &&
                                 predicates::orient3d(p, source.node(tri[0]), source.node(tri[1]),
                                                      source.node(tri[2])) ==
                                     predicates::Sign::Positive,
                             R::NonConvexCell);
                        const double sum = x.norm2() + y.norm2() + z.norm2() + (x - y).norm2() +
                                           (x - z).norm2() + (y - z).norm2();
                        const double q = 12 * std::pow(det / 2, 2. / 3) / sum;
                        need(std::isfinite(q) && q > 0, R::NonConvexCell);
                        r.minimumQuality = std::min(r.minimumQuality, q);
                        if (q < o.warningQuality)
                            ++r.lowQualityCells;
                        volume += v;
                        need(tet < r.tetrahedra, R::InternalError);
                        storage.setCell(tet, 4 * tet, mesh::CellType::Tet4,
                                        std::array<Index, 4>{apex, tri[0], tri[1], tri[2]},
                                        c.onEntity);
                        ++tet;
                    });
                }
            }
        }
        need(tet == r.tetrahedra && boundary == r.boundaryTriangles, R::InternalError);
        r.volume = static_cast<double>(volume);
        need(std::isfinite(r.volume) && r.volume > 0, R::InvalidSource);
        need(r.minimumQuality >= o.minimumQuality, R::QualityRequirement);
        mesh::Mesh staged;
        storage.publishTo(staged);
        r.audit = audit(staged);
        cancel();
        destination.swap(staged);
        r.published = true;
        r.destinationUnchanged = false;
        r.reason = R::None;
    } catch (const Stop &s) {
        r.reason = s.reason;
    } catch (const std::bad_alloc &) {
        r.reason = R::AllocationFailed;
    } catch (...) {
        r.reason = R::InternalError;
    }
    return r;
}
const char *convexToTetReasonName(ConvexToTetReason r) noexcept {
    switch (r) {
    case R::None:
        return "none";
    case R::InvalidRequest:
        return "invalid_request";
    case R::InvalidSource:
        return "invalid_source";
    case R::NonConvexCell:
        return "nonconvex_cell";
    case R::NonplanarFace:
        return "nonplanar_face";
    case R::BoundaryAuditFailed:
        return "boundary_audit_failed";
    case R::OutputBudget:
        return "output_budget";
    case R::WorkingBudget:
        return "working_budget";
    case R::QualityRequirement:
        return "quality_requirement";
    case R::Cancelled:
        return "cancelled";
    case R::AllocationFailed:
        return "allocation_failed";
    case R::InternalError:
        return "internal_error";
    }
    return "unknown";
}
} // namespace paacs::algo
