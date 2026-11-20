// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "PrimitiveMesherInternal.hpp"
#include "paacs/algo/Delaunay3D.hpp"
#include "paacs/mesh/CertifiedMeshStorage.hpp"
#include <cmath>
#include <limits>
#include <new>
#include <vector>

namespace paacs::algo::primitive_detail {
namespace {
struct Stop {
    PrimitiveReason reason;
};
void need(bool ok, PrimitiveReason r) {
    if (!ok)
        throw Stop{r};
}
struct Node {
    Vec2 p;
    std::uint32_t plane;
    int rim;
};
struct Face {
    std::array<Index, 3> key, oriented;
};
constexpr std::array<std::array<unsigned, 3>, 4> faces{
    {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}}};
double jitter(int i, int j, unsigned lane) {
    std::uint64_t h = static_cast<std::uint32_t>(i) * 0x9e3779b97f4a7c15ULL +
                      static_cast<std::uint32_t>(j) * 0xbf58476d1ce4e5b9ULL + lane;
    h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
    h = (h ^ (h >> 27)) * 0x94d049bb133111ebULL;
    h ^= h >> 31;
    return (static_cast<double>(h >> 11) / 9007199254740992. - .5) * .2;
}
} // namespace
PrimitiveReason profiledCylinderSource(const PrimitiveShape &s, const PrimitiveMeshingOptions &o,
                                       PrimitiveMeshingResult &r, mesh::Mesh &result) noexcept {
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
            r.plannedWorkingBytes =
                std::max(r.plannedWorkingBytes, static_cast<std::uint64_t>(bytes));
        };
        auto output = [&](std::uint64_t n, std::uint64_t t, std::uint64_t b) {
            const long double conn = 4.L * t + 3.L * b;
            const long double bytes =
                static_cast<long double>(n) * (sizeof(Vec3) + sizeof(mesh::Classification)) +
                (static_cast<long double>(t) + b) * sizeof(mesh::Cell) + conn * sizeof(Index);
            need(n <= o.maximumNodes && n <= INT32_MAX && t + b <= o.maximumCells &&
                     conn <= INT32_MAX && bytes <= o.maximumOutputBytes && bytes <= SIZE_MAX,
                 PrimitiveReason::OutputBudget);
        };
        cancel();
        need(s.kind == PrimitiveKind::Cylinder && o.sizeProfile && o.sizeProfile->valid() &&
                 o.axialPositions.size() >= 2,
             PrimitiveReason::InvalidRequest);
        const auto positions = o.axialPositions;
        const auto planes = positions.size();
        const double roundoff = 128 * std::numeric_limits<double>::epsilon() *
                                std::max({s.radius, s.height, std::abs(s.origin.x),
                                          std::abs(s.origin.y), std::abs(s.origin.z)});
        const double tolerance = o.maximumChordError > 0 ? o.maximumChordError : .0025 * s.radius;
        need(tolerance > roundoff, PrimitiveReason::NumericallyUnresolved);
        const double angle = std::min(
            4 * std::asin(std::sqrt(std::min(.5, (tolerance - roundoff) / (2 * s.radius)))),
            std::sqrt(6 * o.maximumRelativeVolumeError));
        const double quarter =
            std::ceil(std::max(std::numbers::pi * s.radius / (2 * o.sizeProfile->minimumSize()),
                               std::numbers::pi / (2 * angle)));
        need(std::isfinite(quarter) && quarter >= 1 && quarter <= INT32_MAX / 4 &&
                 4 * quarter + 1 <= o.maximumProfileSectionNodes,
             PrimitiveReason::WorkingBudget);
        const auto rim = static_cast<std::uint32_t>(4 * quarter);
        auto section = [&](std::size_t k, auto emit) {
            const double nvalue = std::ceil(s.radius / o.sizeProfile->eval(positions[k]));
            need(std::isfinite(nvalue) && nvalue >= 1 && nvalue <= INT32_MAX / 4 &&
                     nvalue * nvalue <= o.maximumProfileSectionNodes,
                 PrimitiveReason::WorkingBudget);
            const int n = static_cast<int>(nvalue);
            const double step = s.radius / n, interior = s.radius - .65 * step;
            for (int j = -n; j <= n; ++j) {
                cancel();
                for (int i = -n; i <= n; ++i) {
                    const Vec2 p{i == 0 && j == 0 ? 0 : (i + jitter(i, j, 17)) * step,
                                 i == 0 && j == 0 ? 0 : (j + jitter(i, j, 31)) * step};
                    if (p.u * p.u + p.v * p.v < interior * interior)
                        emit(Node{p, static_cast<std::uint32_t>(k), -1});
                }
            }
            for (std::uint32_t j = 0; j < rim; ++j) {
                const double a = 2 * std::numbers::pi * static_cast<double>(j) / rim;
                emit(Node{{s.radius * std::cos(a), s.radius * std::sin(a)},
                          static_cast<std::uint32_t>(k),
                          static_cast<int>(j)});
            }
        };
        work(static_cast<long double>(planes + 1) * sizeof(std::uint64_t));
        std::vector<std::uint64_t> starts(planes + 1);
        r.minimumSectionNodes = UINT64_MAX;
        for (std::size_t k = 0; k < planes; ++k) {
            std::uint64_t n = 0;
            section(k, [&](Node) {
                ++n;
                need(n <= o.maximumProfileSectionNodes, PrimitiveReason::WorkingBudget);
            });
            starts[k + 1] = starts[k] + n;
            output(starts[k + 1], 0, 0);
            r.minimumSectionNodes = std::min(r.minimumSectionNodes, n);
            r.maximumSectionNodes = std::max(r.maximumSectionNodes, n);
        }
        const long double nodeBytes = static_cast<long double>(starts.back()) * sizeof(Node) +
                                      starts.capacity() * sizeof(std::uint64_t);
        work(nodeBytes);
        std::vector<Node> nodes(static_cast<std::size_t>(starts.back()));
        for (std::size_t k = 0; k < planes; ++k) {
            auto at = starts[k];
            section(k, [&](Node node) { nodes[at++] = node; });
            need(at == starts[k + 1], PrimitiveReason::InternalError);
        }
        std::vector<Tet> tets;
        std::size_t previousStart = 0, previousCount = 0;
        for (std::size_t k = 0; k + 1 < planes; ++k) {
            cancel();
            const auto n = starts[k + 2] - starts[k];
            need(n <= UINT64_MAX / n && n * n <= o.maximumProfileSlabWork - r.profileSlabWork,
                 PrimitiveReason::WorkingBudget);
            r.profileSlabWork += n * n;
            // Conservative O(n^2) bound for arena + scratch + slab output.
            // Keep the same admission limits when connectivity can be reused.
            const long double slabBound = 1024.L * n * n;
            work(nodeBytes + tets.capacity() * sizeof(Tet) + slabBound);
            bool reuse = o.reuseProfileSlabs && k > 0 &&
                         starts[k + 1] - starts[k] == starts[k] - starts[k - 1] &&
                         n == starts[k + 1] - starts[k - 1];
            if (reuse)
                for (std::uint64_t i = 0; i < n; ++i) {
                    const auto &a = nodes[starts[k - 1] + i], &b = nodes[starts[k] + i];
                    if (a.p.u != b.p.u || a.p.v != b.p.v || a.rim != b.rim) {
                        reuse = false;
                        break;
                    }
                }
            std::vector<Tet> slab;
            if (!reuse) {
                std::vector<Vec3> points(static_cast<std::size_t>(n));
                for (std::uint64_t i = 0; i < n; ++i) {
                    const auto a = nodes[starts[k] + i];
                    points[i] = {a.p.u / s.radius, a.p.v / s.radius,
                                 (positions[a.plane] - positions[k]) / s.radius};
                }
                Delaunay3DStats stats;
                slab = delaunay3D(points, Delaunay3DOptions{}, &stats);
                need(stats.runComplete && stats.auditPassed && stats.pointsInserted == n &&
                         stats.pointsDuplicate == 0 && stats.pointsSkipped == 0 && !slab.empty(),
                     PrimitiveReason::NumericallyUnresolved);
                r.profileDelaunayBytes =
                    std::max(r.profileDelaunayBytes, stats.arenaBytes + stats.scratchBytes);
                need(stats.arenaBytes + stats.scratchBytes + stats.outputBytes +
                             points.capacity() * sizeof(Vec3) <=
                         slabBound,
                     PrimitiveReason::WorkingBudget);
                ++r.profileDelaunaySlabs;
            } else
                ++r.profileReusedSlabs;
            const auto count = reuse ? previousCount : slab.size(), first = tets.size();
            output(nodes.size(), first + count, 0);
            if (first + count > tets.capacity()) {
                const auto capacity = std::max(first + count, 2 * tets.capacity());
                work(nodeBytes +
                     (static_cast<long double>(tets.capacity()) + capacity) * sizeof(Tet) +
                     slabBound);
                tets.reserve(capacity);
            }
            for (std::size_t j = 0; j < count; ++j) {
                auto t = reuse ? tets[previousStart + j] : slab[j];
                const auto shift = reuse ? starts[k] - starts[k - 1] : starts[k];
                for (auto &id : t.v)
                    id += static_cast<Index>(shift);
                tets.push_back(t);
            }
            previousStart = first;
            previousCount = count;
        }
        const long double retained = nodeBytes + tets.capacity() * sizeof(Tet);
        work(retained + static_cast<long double>(tets.size()) * 4 * sizeof(Face));
        std::vector<Face> uses(tets.size() * 4);
        for (std::size_t i = 0; i < tets.size(); ++i) {
            if (i % 1024 == 0)
                cancel();
            for (unsigned j = 0; j < 4; ++j) {
                auto &f = uses[4 * i + j];
                for (unsigned q = 0; q < 3; ++q)
                    f.oriented[q] = tets[i].v[faces[j][q]];
                f.key = f.oriented;
                std::sort(f.key.begin(), f.key.end());
            }
        }
        std::sort(uses.begin(), uses.end(), [](auto &a, auto &b) { return a.key < b.key; });
        auto classification = [&](const Face &f) {
            const auto a = nodes[f.key[0]], b = nodes[f.key[1]], c = nodes[f.key[2]];
            if (a.plane == 0 && b.plane == 0 && c.plane == 0)
                return 2;
            if (a.plane + 1 == planes && b.plane + 1 == planes && c.plane + 1 == planes)
                return 3;
            need(a.rim >= 0 && b.rim >= 0 && c.rim >= 0, PrimitiveReason::BoundaryAuditFailed);
            std::array<int, 3> rims{a.rim, b.rim, c.rim};
            std::sort(rims.begin(), rims.end());
            const int lo = rims.front(), hi = rims.back();
            need(lo != hi && (rims[1] == lo || rims[1] == hi) &&
                     (hi - lo == 1 || (lo == 0 && hi == static_cast<int>(rim) - 1)),
                 PrimitiveReason::BoundaryAuditFailed);
            return 1;
        };
        std::uint64_t boundaries = 0;
        for (std::size_t i = 0; i < uses.size();) {
            if (i % 1024 == 0)
                cancel();
            std::size_t j = i + 1;
            while (j < uses.size() && uses[j].key == uses[i].key)
                ++j;
            need(j - i <= 2, PrimitiveReason::BoundaryAuditFailed);
            if (j - i == 1) {
                (void)classification(uses[i]);
                ++boundaries;
            }
            i = j;
        }
        output(nodes.size(), tets.size(), boundaries);
        const auto unit = [](Vec3 p) { return p / std::hypot(p.x, p.y, p.z); };
        const Vec3 z = unit(s.axis), hint = unit(s.xDirection), x = unit(hint - z * hint.dot(z)),
                   y = z.cross(x);
        mesh::detail::CertifiedMeshStorage storage(nodes.size(), tets.size() + boundaries,
                                                   4 * tets.size() + 3 * boundaries);
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            if (i % 1024 == 0)
                cancel();
            const auto n = nodes[i];
            const int cap = n.plane == 0 ? 2 : n.plane + 1 == planes ? 3 : 0;
            const mesh::Classification tag = n.rim >= 0 && cap ? mesh::Classification{1, 99 + cap}
                                             : n.rim >= 0      ? mesh::Classification{2, 1}
                                             : cap             ? mesh::Classification{2, cap}
                                                               : mesh::Classification{3, 1};
            storage.setNode(i, s.origin + x * n.p.u + y * n.p.v + z * positions[n.plane], tag);
        }
        for (std::size_t i = 0; i < tets.size(); ++i)
            storage.setCell(i, 4 * i, mesh::CellType::Tet4, tets[i].v, {3, 1});
        std::size_t at = tets.size(), conn = 4 * tets.size();
        for (std::size_t i = 0; i < uses.size();) {
            cancel();
            std::size_t j = i + 1;
            while (j < uses.size() && uses[j].key == uses[i].key)
                ++j;
            if (j - i == 1) {
                storage.setCell(at++, conn, mesh::CellType::Tri3, uses[i].oriented,
                                {2, classification(uses[i])});
                conn += 3;
            }
            i = j;
        }
        r.nodes = nodes.size();
        r.tet4s = tets.size();
        r.boundaryTriangles = boundaries;
        r.profiledVolume = true;
        r.sweepLayers = static_cast<std::uint32_t>(planes - 1);
        r.segmentsPerQuarter = rim / 4;
        r.chordErrorBound =
            2 * s.radius * std::pow(std::sin(std::numbers::pi / (2 * rim)), 2) + roundoff;
        r.analyticVolume = std::numbers::pi * s.radius * s.radius * s.height;
        need(r.chordErrorBound <= tolerance, PrimitiveReason::GeometryTolerance);
        cancel();
        storage.publishTo(result);
        return PrimitiveReason::None;
    } catch (const Stop &e) {
        return e.reason;
    } catch (const std::bad_alloc &) {
        return PrimitiveReason::AllocationFailed;
    } catch (...) {
        return PrimitiveReason::InternalError;
    }
}
} // namespace paacs::algo::primitive_detail
