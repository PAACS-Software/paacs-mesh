// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/TransitionTemplateCatalog.hpp"

#include "paacs/mesh/CertifiedMeshAppender.hpp"
#include "paacs/predicates/Predicates.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <new>

namespace paacs::algo {
namespace {
constexpr std::array<TransitionTemplateInfo, 4> kCatalog{
    {{TransitionTemplateKind::QuadToApexPyramid5, mesh::CellType::Quad4, mesh::CellType::Pyramid5,
      1, true},
     {TransitionTemplateKind::Pyramid5ToTet4, mesh::CellType::Pyramid5, mesh::CellType::Tet4, 2,
      true},
     {TransitionTemplateKind::Wedge6ToTet4, mesh::CellType::Prism6, mesh::CellType::Tet4, 3, true},
     {TransitionTemplateKind::Hex8ToTet4, mesh::CellType::Hex8, mesh::CellType::Tet4, 6, true}}};
bool finite(const Vec3 &p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
bool validClass(mesh::Classification c, int dim) noexcept {
    return c.dim == dim && c.tag >= 0;
}
double signedTet(const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &d) noexcept {
    return (b - a).cross(c - a).dot(d - a) / 6.0;
}
bool orientTet(std::array<Index, 4> &tet, const mesh::Mesh &mesh) noexcept {
    const auto sign = predicates::orient3d(mesh.node(tet[0]), mesh.node(tet[1]), mesh.node(tet[2]),
                                           mesh.node(tet[3]));
    if (sign == predicates::Sign::Zero)
        return false;
    if (sign == predicates::Sign::Negative)
        std::swap(tet[2], tet[3]);
    return true;
}
TransitionOutcome rejected(TransitionReason reason) noexcept {
    TransitionOutcome out;
    out.status = TransitionStatus::Rejected;
    out.reason = reason;
    return out;
}
constexpr std::array<std::array<std::uint8_t, 4>, 6> kHexTets{
    {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}}};
constexpr std::array<std::array<std::uint8_t, 4>, 3> kWedgeTets{
    {{0, 1, 2, 3}, {1, 4, 2, 3}, {2, 4, 5, 3}}};
} // namespace

std::span<const TransitionTemplateInfo> transitionTemplateCatalog() noexcept {
    return kCatalog;
}

TransitionOutcome generatePyramidTransitionFan(const mesh::Mesh &quadSurface,
                                               const PyramidTransitionRequest &request,
                                               mesh::Mesh &destination) noexcept {
    if (quadSurface.nodeCount() == 0 || quadSurface.cellCount() == 0)
        return rejected(TransitionReason::EmptyInput);
    if (!finite(request.apex))
        return rejected(TransitionReason::NonfiniteGeometry);
    if (!validClass(request.volumeClassification, 3) ||
        !validClass(request.apexClassification, 3) ||
        (request.emitBoundaryFaces && !validClass(request.sideClassification, 2)))
        return rejected(TransitionReason::InvalidClassification);
    std::uint64_t quads = 0;
    for (Index ci = 0; ci < quadSurface.cellCount(); ++ci) {
        if (quadSurface.cell(ci).type != mesh::CellType::Quad4)
            return rejected(TransitionReason::UnsupportedCell);
        ++quads;
    }
    const std::uint64_t side = request.emitBoundaryFaces ? quads * 4U : 0U;
    const std::uint64_t bases = request.emitBoundaryFaces ? quads : 0U;
    const std::uint64_t cells = quads + side + bases;
    if (quads > request.maximumCells ||
        cells > static_cast<std::uint64_t>((std::numeric_limits<Index>::max)()))
        return rejected(TransitionReason::BudgetExceeded);
    TransitionOutcome out;
    try {
        mesh::Mesh local;
        mesh::detail::CertifiedMeshAppender appender(
            local, static_cast<std::size_t>(quadSurface.nodeCount()) + 1U,
            static_cast<std::size_t>(cells),
            static_cast<std::size_t>(quads * 5U + side * 3U + bases * 4U));
        for (Index n = 0; n < quadSurface.nodeCount(); ++n)
            (void)appender.addNode(quadSurface.node(n), quadSurface.nodeClass(n));
        const Index apex = appender.addNode(request.apex, request.apexClassification);
        double minimum = (std::numeric_limits<double>::infinity)();
        long double total = 0;
        for (Index ci = 0; ci < quadSurface.cellCount(); ++ci) {
            const auto &cell = quadSurface.cell(ci);
            if (!validClass(cell.onEntity, 2))
                return rejected(TransitionReason::InvalidClassification);
            const Index *src = quadSurface.cellNodes(ci);
            std::array<Index, 4> base{src[0], src[1], src[2], src[3]};
            auto sorted = base;
            std::sort(sorted.begin(), sorted.end());
            if (sorted[0] < 0 || sorted[3] >= quadSurface.nodeCount() ||
                std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
                return rejected(TransitionReason::InvalidConnectivity);
            const double a = signedTet(local.node(base[0]), local.node(base[1]),
                                       local.node(base[2]), request.apex);
            const double b = signedTet(local.node(base[0]), local.node(base[2]),
                                       local.node(base[3]), request.apex);
            if (!std::isfinite(a) || !std::isfinite(b) || !(a > 0) || !(b > 0))
                return rejected(TransitionReason::NonPositiveCell);
            minimum = std::min(minimum, a + b);
            total += a + b;
            const std::array<Index, 5> pyramid{base[0], base[1], base[2], base[3], apex};
            (void)appender.addCell(mesh::CellType::Pyramid5, pyramid, request.volumeClassification);
            if (request.emitBoundaryFaces) {
                const std::array<Index, 4> bottom{base[3], base[2], base[1], base[0]};
                (void)appender.addCell(mesh::CellType::Quad4, bottom, cell.onEntity);
                for (std::size_t e = 0; e < 4; ++e) {
                    const std::array<Index, 3> tri{base[e], base[(e + 1U) % 4U], apex};
                    (void)appender.addCell(mesh::CellType::Tri3, tri, request.sideClassification);
                }
            }
        }
        local.swap(destination);
        out.status = TransitionStatus::Published;
        out.reason = TransitionReason::None;
        out.destinationUnchanged = false;
        out.metrics.inputCells = quads;
        out.metrics.outputNodes = destination.nodeCount();
        out.metrics.outputPyramid5 = quads;
        out.metrics.boundaryTriangles = side;
        out.metrics.boundaryQuads = bases;
        out.metrics.minimumVolume = minimum;
        out.metrics.totalVolume = static_cast<double>(total);
        return out;
    } catch (const std::bad_alloc &) {
        out.reason = TransitionReason::AllocationFailed;
    } catch (...) {
        out.reason = TransitionReason::InternalError;
    }
    return out;
}

TransitionOutcome convertSupportedMixedMeshToTet4(const mesh::Mesh &source,
                                                  const MixedToTet4Request &request,
                                                  mesh::Mesh &destination) noexcept {
    if (source.nodeCount() == 0 || source.cellCount() == 0)
        return rejected(TransitionReason::EmptyInput);
    std::uint64_t tets = 0, boundary = 0;
    for (Index ci = 0; ci < source.cellCount(); ++ci) {
        switch (source.cell(ci).type) {
        case mesh::CellType::Tet4:
            ++tets;
            break;
        case mesh::CellType::Pyramid5:
            tets += 2;
            break;
        case mesh::CellType::Prism6:
            tets += 3;
            break;
        case mesh::CellType::Hex8:
            tets += 6;
            break;
        case mesh::CellType::Tri3:
            if (request.retainBoundary)
                ++boundary;
            break;
        case mesh::CellType::Quad4:
            if (request.retainBoundary) {
                if (!request.triangulateBoundaryQuads)
                    return rejected(TransitionReason::UnsupportedCell);
                boundary += 2;
            }
            break;
        default:
            return rejected(TransitionReason::UnsupportedCell);
        }
        if (tets > request.maximumTet4)
            return rejected(TransitionReason::BudgetExceeded);
    }
    if (tets + boundary > static_cast<std::uint64_t>((std::numeric_limits<Index>::max)()))
        return rejected(TransitionReason::CountOverflow);
    TransitionOutcome out;
    try {
        mesh::Mesh local;
        mesh::detail::CertifiedMeshAppender appender(
            local, static_cast<std::size_t>(source.nodeCount()),
            static_cast<std::size_t>(tets + boundary),
            static_cast<std::size_t>(tets * 4U + boundary * 3U));
        for (Index n = 0; n < source.nodeCount(); ++n) {
            if (!finite(source.node(n)))
                return rejected(TransitionReason::NonfiniteGeometry);
            (void)appender.addNode(source.node(n), source.nodeClass(n));
        }
        double minimum = (std::numeric_limits<double>::infinity)();
        long double total = 0;
        const auto emit = [&](std::array<Index, 4> tet, mesh::Classification cls) {
            if (!orientTet(tet, local))
                return false;
            const double v = signedTet(local.node(tet[0]), local.node(tet[1]), local.node(tet[2]),
                                       local.node(tet[3]));
            if (!std::isfinite(v) || !(v > 0))
                return false;
            minimum = std::min(minimum, v);
            total += v;
            (void)appender.addCell(mesh::CellType::Tet4, tet, cls);
            return true;
        };
        for (Index ci = 0; ci < source.cellCount(); ++ci) {
            const auto &c = source.cell(ci);
            const Index *n = source.cellNodes(ci);
            if (c.type == mesh::CellType::Tri3) {
                if (request.retainBoundary) {
                    const std::array<Index, 3> tri{n[0], n[1], n[2]};
                    (void)appender.addCell(mesh::CellType::Tri3, tri, c.onEntity);
                }
                continue;
            }
            if (c.type == mesh::CellType::Quad4) {
                if (request.retainBoundary) {
                    const bool d02 = std::pair{std::min(n[0], n[2]), std::max(n[0], n[2])} <
                                     std::pair{std::min(n[1], n[3]), std::max(n[1], n[3])};
                    if (d02) {
                        const std::array<Index, 3> a{n[0], n[1], n[2]}, b{n[0], n[2], n[3]};
                        (void)appender.addCell(mesh::CellType::Tri3, a, c.onEntity);
                        (void)appender.addCell(mesh::CellType::Tri3, b, c.onEntity);
                    } else {
                        const std::array<Index, 3> a{n[0], n[1], n[3]}, b{n[1], n[2], n[3]};
                        (void)appender.addCell(mesh::CellType::Tri3, a, c.onEntity);
                        (void)appender.addCell(mesh::CellType::Tri3, b, c.onEntity);
                    }
                }
                continue;
            }
            if (!validClass(c.onEntity, 3))
                return rejected(TransitionReason::InvalidClassification);
            if (c.type == mesh::CellType::Tet4) {
                if (!emit({n[0], n[1], n[2], n[3]}, c.onEntity))
                    return rejected(TransitionReason::DegenerateCell);
            } else if (c.type == mesh::CellType::Pyramid5) {
                const bool d02 = std::pair{std::min(n[0], n[2]), std::max(n[0], n[2])} <
                                 std::pair{std::min(n[1], n[3]), std::max(n[1], n[3])};
                if (d02) {
                    if (!emit({n[0], n[1], n[2], n[4]}, c.onEntity) ||
                        !emit({n[0], n[2], n[3], n[4]}, c.onEntity))
                        return rejected(TransitionReason::DegenerateCell);
                } else {
                    if (!emit({n[0], n[1], n[3], n[4]}, c.onEntity) ||
                        !emit({n[1], n[2], n[3], n[4]}, c.onEntity))
                        return rejected(TransitionReason::DegenerateCell);
                }
            } else if (c.type == mesh::CellType::Prism6) {
                for (const auto &localTet : kWedgeTets) {
                    std::array<Index, 4> tet{};
                    for (std::size_t k = 0; k < 4; ++k)
                        tet[k] = n[localTet[k]];
                    if (!emit(tet, c.onEntity))
                        return rejected(TransitionReason::DegenerateCell);
                }
            } else if (c.type == mesh::CellType::Hex8) {
                for (const auto &localTet : kHexTets) {
                    std::array<Index, 4> tet{};
                    for (std::size_t k = 0; k < 4; ++k)
                        tet[k] = n[localTet[k]];
                    if (!emit(tet, c.onEntity))
                        return rejected(TransitionReason::DegenerateCell);
                }
            }
        }
        local.swap(destination);
        out.status = TransitionStatus::Published;
        out.reason = TransitionReason::None;
        out.destinationUnchanged = false;
        out.metrics.inputCells = source.cellCount();
        out.metrics.outputNodes = destination.nodeCount();
        out.metrics.outputTet4 = tets;
        out.metrics.boundaryTriangles = boundary;
        out.metrics.minimumVolume = minimum;
        out.metrics.totalVolume = static_cast<double>(total);
        return out;
    } catch (const std::bad_alloc &) {
        out.reason = TransitionReason::AllocationFailed;
    } catch (...) {
        out.reason = TransitionReason::InternalError;
    }
    return out;
}

} // namespace paacs::algo
