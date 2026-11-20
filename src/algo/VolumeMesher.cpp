// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Meshers.hpp"
#include "paacs/algo/Delaunay3D.hpp"
#include "paacs/algo/TetCandidate.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/core/PointGrid.hpp"
#include "paacs/mesh/Quality.hpp"
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4100)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include "paacs/topo/Workspace.hpp"
#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

#include "CheckedSizing.hpp"

#include <cmath>
#include <cstdint>
#include <map>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace paacs::algo {

using mesh::CellType;
using mesh::Classification;

namespace {

bool finite(const Vec2& value) noexcept {
    return std::isfinite(value.u) && std::isfinite(value.v);
}

std::size_t checkedAxisCount(double lower, double upper, double step) {
    if (!std::isfinite(lower) || !std::isfinite(upper) || upper < lower ||
        !detail::positiveFinite(step))
        throw detail::SizingEvaluationError{};
    const double trim = 0.4 * step;
    const double start = lower + step;
    const double stop = upper - trim;
    if (!std::isfinite(trim) || !std::isfinite(start) ||
        !std::isfinite(stop))
        throw detail::SizingEvaluationError{};
    if (upper > lower && !(start > lower))
        throw detail::SizingEvaluationError{};
    if (!(start < stop)) return 0;
    detail::requireProgress(start, step);
    const double span = stop - start;
    const double count = std::ceil(span / step);
    if (!detail::positiveFinite(span) || !detail::positiveFinite(count) ||
        count > static_cast<double>(
                    (std::numeric_limits<std::size_t>::max)()))
        throw detail::SizingEvaluationError{};
    return static_cast<std::size_t>(count);
}

std::size_t checkedCountProduct(std::size_t a, std::size_t b) {
    if (a != 0 && b > (std::numeric_limits<std::size_t>::max)() / a)
        throw detail::SizingEvaluationError{};
    return a * b;
}

void validateGridQuery(double cellSize, double radius) {
    if (!detail::positiveFinite(cellSize) ||
        !detail::positiveFinite(radius))
        throw detail::SizingEvaluationError{};
    const double radiusSquared = radius * radius;
    const double span = std::ceil(radius / cellSize);
    if (!detail::positiveFinite(radiusSquared) || !std::isfinite(span) ||
        span > static_cast<double>(
                   (std::numeric_limits<long long>::max)() / 4))
        throw detail::SizingEvaluationError{};
}

std::set<int> checkedCurvedFaceTags(const geo::Model& model) {
    std::set<int> tags;
    for (const geo::GFace* face : model.faces()) {
        if (!face) throw detail::SizingEvaluationError{};
        Vec2 lo, hi;
        face->paramBounds(lo, hi);
        if (!finite(lo) || !finite(hi) || hi.u < lo.u || hi.v < lo.v)
            throw detail::SizingEvaluationError{};
        const Vec2 samples[5] = {
            {(lo.u + hi.u) / 2, (lo.v + hi.v) / 2},
            {lo.u + 0.25 * (hi.u - lo.u),
             lo.v + 0.25 * (hi.v - lo.v)},
            {lo.u + 0.75 * (hi.u - lo.u),
             lo.v + 0.25 * (hi.v - lo.v)},
            {lo.u + 0.25 * (hi.u - lo.u),
             lo.v + 0.75 * (hi.v - lo.v)},
            {lo.u + 0.75 * (hi.u - lo.u),
             lo.v + 0.75 * (hi.v - lo.v)}};
        for (const Vec2& uv : samples) {
            if (!finite(uv)) throw detail::SizingEvaluationError{};
            const double curvature = face->curvatureAt(uv);
            if (!std::isfinite(curvature))
                throw detail::SizingEvaluationError{};
            if (curvature > 0.0) {
                tags.insert(face->tag());
                break;
            }
        }
    }
    return tags;
}

void accumulateHandoffWork(TetHandoffMetrics& target,
                           const TetHandoffMetrics& source) noexcept {
    using WorkField=std::uint64_t TetHandoffMetrics::*;
    static constexpr WorkField fields[]{
        &TetHandoffMetrics::selectedCells,
        &TetHandoffMetrics::importedNodes,
        &TetHandoffMetrics::importedCells,
        &TetHandoffMetrics::importEdgeValidationCellVisits,
        &TetHandoffMetrics::importEdgeBuildCellVisits,
        &TetHandoffMetrics::importEdgeUses,
        &TetHandoffMetrics::importUniqueEdges,
        &TetHandoffMetrics::importFaceValidationCellVisits,
        &TetHandoffMetrics::importFaceBuildCellVisits,
        &TetHandoffMetrics::importFaceUses,
        &TetHandoffMetrics::importUniqueFaces,
        &TetHandoffMetrics::importEdgeAudits,
        &TetHandoffMetrics::importFaceAudits,
        &TetHandoffMetrics::importPerCellGlobalAudits,
        &TetHandoffMetrics::importCandidateValidationCalls,
        &TetHandoffMetrics::importCandidateValidationPointRows,
        &TetHandoffMetrics::importCandidateValidationCellRows,
        &TetHandoffMetrics::importCandidateValidationFaceUses,
        &TetHandoffMetrics::importCandidateValidationFaceRows,
        &TetHandoffMetrics::importCandidateValidationBoundaryRows,
        &TetHandoffMetrics::importCandidateValidationDualNeighbourVisits,
        &TetHandoffMetrics::importCandidateValidationBoundaryEdgeUses,
        &TetHandoffMetrics::importCandidateValidationBoundaryNeighbourVisits,
        &TetHandoffMetrics::importCandidateValidationHullTests,
        &TetHandoffMetrics::importCandidateValidationSlotsInitialized,
        &TetHandoffMetrics::importSelectionRowsCopied,
        &TetHandoffMetrics::importSelectionSortComparisons,
        &TetHandoffMetrics::importSelectionDuplicateComparisons,
        &TetHandoffMetrics::importSelectionCellVisits,
        &TetHandoffMetrics::importSelectionVertexVisits,
        &TetHandoffMetrics::importSelectionPointSlotsInitialized,
        &TetHandoffMetrics::importGlobalMapRows,
        &TetHandoffMetrics::importGlobalMapComparisons,
        &TetHandoffMetrics::importGlobalMapInsertions,
        &TetHandoffMetrics::importLocalPointSlotsInitialized,
        &TetHandoffMetrics::importLocalGlobalSlotsInitialized,
        &TetHandoffMetrics::importLocalCellSlotsInitialized,
        &TetHandoffMetrics::importPersistentNodeSlotsInitialized,
        &TetHandoffMetrics::partialImportRowsCopied,
        &TetHandoffMetrics::partialImportGlobalSlotsInitialized,
        &TetHandoffMetrics::partialImportCandidateMappingVisits,
        &TetHandoffMetrics::partialImportEdgeValidationRows,
        &TetHandoffMetrics::partialImportFaceValidationRows,
        &TetHandoffMetrics::partialImportProtectedEndpointVisits,
        &TetHandoffMetrics::partialImportAuthorityPointValidations,
        &TetHandoffMetrics::partialImportDuplicateKeyComparisons,
        &TetHandoffMetrics::partialImportResolutionNodeVisits,
        &TetHandoffMetrics::partialImportEdgeResolutionLookups,
        &TetHandoffMetrics::partialImportFaceResolutionLookups,
        &TetHandoffMetrics::partialImportEdgeDataReads,
        &TetHandoffMetrics::partialImportFaceDataCopies,
        &TetHandoffMetrics::partialImportFaceIncidenceRowsCopied,
        &TetHandoffMetrics::partialImportSortComparisons};
    bool complete=source.workAccountingComplete;
    constexpr auto max=(std::numeric_limits<std::uint64_t>::max)();
    for (const WorkField field:fields) {
        auto& destination=target.*field;
        const auto value=source.*field;
        if (value>max-destination) {
            destination=max;
            complete=false;
        } else {
            destination+=value;
        }
    }
    target.workAccountingComplete=
        target.workAccountingComplete&&complete;
}

}  // namespace

void meshVolumesImpl(const geo::Model& model, const sizing::MeshSpec& spec,
                     mesh::Mesh& out, TetHandoffMetrics* handoffMetrics) {
    const bool checked = detail::usesCheckedSizing(spec);
    const bool authoritative = handoffMetrics != nullptr;
    for (const geo::GSolid* s : model.solids()) {
        std::unordered_set<int> shellFaces;
        for (int ft : s->boundaryFaceTags()) shellFaces.insert(ft);

        Vec3 lo, hi;
        s->bbox(lo, hi);
        if (checked &&
            (!detail::finite(lo) || !detail::finite(hi) || hi.x < lo.x ||
             hi.y < lo.y || hi.z < lo.z))
            throw detail::SizingEvaluationError{};
        const double hRef = checked
                                ? detail::evaluatedSize(
                                      spec, (lo + hi) * 0.5)
                                : spec.h((lo + hi) * 0.5);
        const double gridCellSize = 0.8 * hRef;
        if (checked) {
            if (!detail::positiveFinite(gridCellSize) ||
                !detail::positiveFinite(1.0 / gridCellSize))
                throw detail::SizingEvaluationError{};
        }
        core::PointGrid grid(gridCellSize);

        // 1. Seed with every existing node on the solid's shell
        //    (vertices, edge nodes and face nodes are all classified).
        std::vector<Index> meshIds;
        std::vector<Vec3>  pts;
        for (Index i = 0; i < out.nodeCount(); ++i) {
            const auto& c = out.nodeClass(i);
            const bool onShell =
                c.dim == 0 || c.dim == 1 ||             // prototype: box shell
                (c.dim == 2 && shellFaces.count(c.tag));
            if (!onShell) continue;
            if (checked && !detail::finite(out.node(i)))
                throw detail::SizingEvaluationError{};
            meshIds.push_back(i);
            pts.push_back(out.node(i));
            grid.insert(out.node(i));
        }

        // 2. Offset shield layer: one interior point behind every surface
        //    node, at ~0.75*h along the inward normal. Without this, the
        //    clearance gap between the wall point-sheet and the lattice lets
        //    Delaunay create flat wall-hugging slivers whose circumspheres
        //    fit inside the gap; the shield points break those spheres.
        //    Parallel-evaluate candidates, then commit serially in order:
        //    identical output to the serial code, on any thread count.
        {
            const size_t nShell = pts.size();
            struct ShieldCand { Vec3 q; double hp; bool ok = false; };
            std::vector<ShieldCand> sc(nShell);
            core::parallelFor(0, nShell, [&](size_t k) {
                const auto& c = out.nodeClass(meshIds[k]);
                if (c.dim != 2) return;                   // faces only
                const geo::GFace* f = model.face(c.tag);
                if (!f) return;
                const Vec3& p = pts[k];
                const Vec2 uvp = f->paramAt(p);
                if (checked && !finite(uvp))
                    throw detail::SizingEvaluationError{};
                // The sheet-gap sliver mechanism is curvature-specific;
                // shielding flat faces only perturbs an already-good layer.
                const double curvature = f->curvatureAt(uvp);
                if (checked && !std::isfinite(curvature))
                    throw detail::SizingEvaluationError{};
                if (curvature <= 1e-12) return;
                const Vec3 n = f->normalAt(uvp);
                if (checked &&
                    (!detail::finite(n) || !detail::positiveFinite(n.norm())))
                    throw detail::SizingEvaluationError{};
                // In size-gradient zones the wall spacing near p can be much
                // finer than h(p); depth must follow the LOCAL scale or the
                // shield sits too deep and re-opens the sheet-gap slivers.
                const double hp0 = checked ? detail::evaluatedSize(spec, p)
                                           : spec.h(p);
                const double firstDepth = 0.35 * hp0;
                const Vec3 insideSample = p - n * firstDepth;
                if (checked &&
                    (!detail::positiveFinite(firstDepth) ||
                     !detail::finite(insideSample)))
                    throw detail::SizingEvaluationError{};
                const double hp1 = checked
                                       ? detail::evaluatedSize(spec,
                                                               insideSample)
                                       : spec.h(insideSample);
                const double hp = std::min(hp0, hp1);
                const double jd = 0.62 +
                    0.25 * (double)((meshIds[k] * 2654435761ULL) & 0xff) / 255.0;
                const double depth = jd * hp;
                const Vec3 q = p - n * depth;
                if (checked &&
                    (!detail::positiveFinite(hp) ||
                     !detail::positiveFinite(depth) || !detail::finite(q)))
                    throw detail::SizingEvaluationError{};
                if (!s->isInside(q)) return;
                sc[k] = {q, hp, true};
            });
            for (size_t k = 0; k < nShell; ++k) {
                if (!sc[k].ok) continue;
                const double radius = 0.5477 * sc[k].hp;
                if (checked) validateGridQuery(gridCellSize, radius);
                if (grid.anyWithin(sc[k].q, radius)) continue;
                meshIds.push_back(
                    out.addNode(sc[k].q, Classification{3, s->tag()}));
                pts.push_back(sc[k].q);
                grid.insert(sc[k].q);
            }
        }

        // 3. Interior lattice driven by the sizing field, kept clear of the
        //    shell points by 0.6*h. Same parallel-evaluate / ordered-commit
        //    pattern (isInside and size-field evaluation dominate).
        {
            std::vector<Vec3> cand;
            if (checked) {
                const std::size_t nx = checkedAxisCount(lo.x, hi.x, hRef);
                const std::size_t ny = checkedAxisCount(lo.y, hi.y, hRef);
                const std::size_t nz = checkedAxisCount(lo.z, hi.z, hRef);
                const std::size_t count = checkedCountProduct(
                    checkedCountProduct(nx, ny), nz);
                if (count > cand.max_size())
                    throw detail::SizingEvaluationError{};
            }
            for (double z = lo.z + hRef; z < hi.z - 0.4 * hRef;) {
                for (double y = lo.y + hRef; y < hi.y - 0.4 * hRef;) {
                    for (double x = lo.x + hRef;
                         x < hi.x - 0.4 * hRef;) {
                        if (checked && !detail::finite(Vec3{x, y, z}))
                            throw detail::SizingEvaluationError{};
                        cand.push_back({x, y, z});
                        if (checked) detail::requireProgress(x, hRef);
                        x += hRef;
                    }
                    if (checked) detail::requireProgress(y, hRef);
                    y += hRef;
                }
                if (checked) detail::requireProgress(z, hRef);
                z += hRef;
            }

            struct LatCand { double hq = 0; bool ok = false; };
            std::vector<LatCand> lc(cand.size());
            core::parallelFor(0, cand.size(), [&](size_t k) {
                if (!s->isInside(cand[k])) return;
                lc[k] = {checked ? detail::evaluatedSize(spec, cand[k])
                                 : spec.h(cand[k]),
                         true};
            });
            for (size_t k = 0; k < cand.size(); ++k) {
                if (!lc[k].ok) continue;
                const double radius = 0.6 * lc[k].hq;
                if (checked) validateGridQuery(gridCellSize, radius);
                if (grid.anyWithin(cand[k], radius)) continue;
                meshIds.push_back(
                    out.addNode(cand[k], Classification{3, s->tag()}));
                pts.push_back(cand[k]);
                grid.insert(cand[k]);
            }
        }

        // 3. Delaunay tetrahedralization of the point set. Ties on the
        //    structured lattice resolve via the inclusive insphere tolerance
        //    (exact adaptive predicates are the robustness milestone).
        Delaunay3DStats dstats;
        std::vector<Tet> tets;
        std::unique_ptr<TetCandidate> tetCandidate;
        if (authoritative) {
            auto candidate = buildTetCandidate(pts);
            dstats = candidate.metrics.construction;
            if (!candidate.built()) {
                if (candidate.reason == TetCandidateReason::AllocationFailed)
                    throw std::bad_alloc{};
                throw std::runtime_error(
                    "meshVolumes: authoritative Tet candidate failed");
            }
            tetCandidate = std::move(candidate.candidate);
        } else {
            tets = delaunay3D(pts, &dstats);
        }
        // No silent geometric loss (product invariant #1): every generated
        // node must participate in the triangulation. Upstream deduplication
        // means duplicates or skips here are construction bugs.
        if (dstats.pointsInserted != pts.size())
            throw std::runtime_error(
                "meshVolumes: Delaunay dropped input points (inserted " +
                std::to_string(dstats.pointsInserted) + " of " +
                std::to_string(pts.size()) + ", duplicates " +
                std::to_string(dstats.pointsDuplicate) + ", skipped " +
                std::to_string(dstats.pointsSkipped) + ")");

        // 4. Keep tets whose centroid is inside; drop numerically flat ones
        //    (jitter-scale artifacts have |V| ~ amp * h^2).
        struct Cand {
            std::array<Index, 4> v{};
            bool alive = true;
            Index candidateCell = kInvalidIndex;
            Index sourceOrdinal = kInvalidIndex;
        };
        static_assert(noexcept(std::declval<std::vector<Cand>&>().swap(
            std::declval<std::vector<Cand>&>())));
        std::vector<Cand> cands;
        double minimumVolume = 0.0;
        if (checked) {
            minimumVolume = 1.0e-9 * hRef;
            if (!detail::positiveFinite(minimumVolume))
                throw detail::SizingEvaluationError{};
            minimumVolume *= hRef;
            if (!detail::positiveFinite(minimumVolume))
                throw detail::SizingEvaluationError{};
            minimumVolume *= hRef;
            if (!detail::positiveFinite(minimumVolume))
                throw detail::SizingEvaluationError{};
        }
        const std::size_t constructionCellCount =
            authoritative ? tetCandidate->cells.size() : tets.size();
        cands.reserve(constructionCellCount);
        for (std::size_t constructionCell = 0;
             constructionCell < constructionCellCount; ++constructionCell) {
            std::array<Index, 4> pointIds{};
            Index sourceOrdinal = static_cast<Index>(constructionCell);
            if (authoritative) {
                const auto& cell = tetCandidate->cells[constructionCell];
                pointIds = cell.vertices;
                sourceOrdinal = cell.sourceOrdinal;
            } else {
                pointIds = tets[constructionCell].v;
            }
            const Vec3 &a = authoritative
                                ? tetCandidate->points[size_t(pointIds[0])]
                                : pts[size_t(pointIds[0])],
                       &b = authoritative
                                ? tetCandidate->points[size_t(pointIds[1])]
                                : pts[size_t(pointIds[1])],
                       &c = authoritative
                                ? tetCandidate->points[size_t(pointIds[2])]
                                : pts[size_t(pointIds[2])],
                       &d = authoritative
                                ? tetCandidate->points[size_t(pointIds[3])]
                                : pts[size_t(pointIds[3])];
            const Vec3 cen = (a + b + c + d) * 0.25;
            if (!s->isInside(cen)) continue;
            if (std::abs(mesh::tetVolume(a, b, c, d)) <
                (checked ? minimumVolume
                         : 1e-9 * hRef * hRef * hRef))
                continue;
            std::array<Index, 4> v{};
            for (int i = 0; i < 4; ++i) {
                Index inputPoint = pointIds[static_cast<std::size_t>(i)];
                if (authoritative)
                    inputPoint = tetCandidate->pointToInput[
                        static_cast<std::size_t>(inputPoint)];
                v[static_cast<std::size_t>(i)] =
                    meshIds[static_cast<std::size_t>(inputPoint)];
            }
            if (mesh::tetVolume(a, b, c, d) < 0) std::swap(v[2], v[3]);
            cands.push_back({v, true,
                             authoritative
                                 ? static_cast<Index>(constructionCell)
                                 : kInvalidIndex,
                             sourceOrdinal});
        }
        if (authoritative)
            std::sort(cands.begin(), cands.end(),
                      [](const Cand& a, const Cand& b) {
                          return a.sourceOrdinal < b.sourceOrdinal;
                      });

        // 5. Boundary sliver peeling. On curved walls the unconstrained
        //    Delaunay hull carries near-flat tets whose four nodes all lie
        //    on the surface. Peel any sub-threshold sliver exposing a hull
        //    face, made purely of boundary nodes, AND touching a curved
        //    face: chordal caps exist only on curved boundaries, and on
        //    planar thin domains such tets ARE the domain (a 1x1x1e-3 box
        //    was previously peeled to zero tets - PM-ALG-002). Repeat,
        //    since peeling exposes new hull faces. (Superseded by
        //    constrained boundary recovery on the roadmap.)
        constexpr double kPeelQ = 0.05;
        const std::set<int> curvedTags = spec.boundaryPeeling
            ? (checked ? checkedCurvedFaceTags(model)
                       : curvedFaceTags(model))
            : std::set<int>{};
        for (int round = 0; round < 20; ++round) {
            std::map<std::array<Index, 3>, int> faceCount;
            static const int ft[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
            for (const auto& cd : cands) {
                if (!cd.alive) continue;
                for (const auto& f : ft) {
                    std::array<Index, 3> key{cd.v[size_t(f[0])],
                                             cd.v[size_t(f[1])],
                                             cd.v[size_t(f[2])]};
                    std::sort(key.begin(), key.end());
                    ++faceCount[key];
                }
            }
            int peeled = 0;
            for (auto& cd : cands) {
                if (!cd.alive) continue;
                bool allBoundary = true, touchesCurved = false;
                for (Index vv : cd.v) {
                    const auto cls = out.nodeClass(vv);
                    allBoundary = allBoundary && cls.dim < 3;
                    touchesCurved = touchesCurved ||
                        (cls.dim == 2 && curvedTags.count(cls.tag));
                }
                if (!allBoundary || !touchesCurved) continue;
                const double q = mesh::tetQuality(
                    out.node(cd.v[0]), out.node(cd.v[1]),
                    out.node(cd.v[2]), out.node(cd.v[3]));
                if (q >= kPeelQ) continue;
                bool hullFace = false;
                static const int ft2[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
                for (const auto& f : ft2) {
                    std::array<Index, 3> key{cd.v[size_t(f[0])],
                                             cd.v[size_t(f[1])],
                                             cd.v[size_t(f[2])]};
                    std::sort(key.begin(), key.end());
                    hullFace = hullFace || faceCount[key] == 1;
                }
                if (!hullFace) continue;
                cd.alive = false;
                ++peeled;
            }
            if (!peeled) break;
        }

        if (!authoritative) {
            for (const auto& cd : cands)
                if (cd.alive)
                    out.addCell(CellType::Tet4,
                                {cd.v[0], cd.v[1], cd.v[2], cd.v[3]},
                                Classification{3, s->tag()});
            continue;
        }

        std::vector<Index> selected;
        selected.reserve(cands.size());
        for (const auto& cd : cands)
            if (cd.alive) selected.push_back(cd.candidateCell);
        std::vector<Index> existingGlobal(tetCandidate->points.size(),
                                          kInvalidIndex);
        std::vector<Classification> expectedClass(
            tetCandidate->points.size());
        for (std::size_t point = 0; point < tetCandidate->points.size(); ++point) {
            const Index inputPoint = tetCandidate->pointToInput[point];
            const Index global = meshIds[static_cast<std::size_t>(inputPoint)];
            existingGlobal[point] = global;
            expectedClass[point] = out.nodeClass(global);
        }

        // Prove the accepted filter connectivity directly against the
        // canonical candidate and the explicit global mapping. This is the
        // last use of the complete filtered connectivity representation; it
        // is released before authoritative import begins.
        std::size_t selectedOrdinal = 0;
        for (const auto& cd : cands) {
            if (!cd.alive) continue;
            if (selectedOrdinal >= selected.size() ||
                selected[selectedOrdinal] != cd.candidateCell)
                throw std::logic_error(
                    "meshVolumes: authoritative selection order drift");
            const auto& cell = tetCandidate->cells[
                static_cast<std::size_t>(cd.candidateCell)];
            std::array<Index, 4> mapped{};
            for (std::size_t i = 0; i < mapped.size(); ++i) {
                const Index point = cell.vertices[i];
                if (point < 0 ||
                    static_cast<std::size_t>(point) >= existingGlobal.size())
                    throw std::logic_error(
                        "meshVolumes: candidate point map drift");
                mapped[i] = existingGlobal[static_cast<std::size_t>(point)];
            }
            if (mapped != cd.v)
                throw std::logic_error(
                    "meshVolumes: filter connectivity drift");
            ++selectedOrdinal;
        }
        if (selectedOrdinal != selected.size())
            throw std::logic_error(
                "meshVolumes: authoritative selection multiplicity drift");

        auto vectorAccounted = [](std::size_t capacity,
                                  std::size_t elementBytes) noexcept {
            const auto max = (std::numeric_limits<std::uint64_t>::max)();
            if (elementBytes != 0 && capacity > max / elementBytes)
                return max;
            return static_cast<std::uint64_t>(capacity * elementBytes);
        };
        auto addSaturating = [](std::uint64_t& target,
                                std::uint64_t value) noexcept {
            const auto max = (std::numeric_limits<std::uint64_t>::max)();
            target = value > max - target ? max : target + value;
        };
        const std::uint64_t filterConnectivityBytes =
            vectorAccounted(cands.capacity(), sizeof(Cand));
        std::uint64_t callerMappingBytes =
            vectorAccounted(selected.capacity(), sizeof(Index));
        addSaturating(callerMappingBytes,
                      vectorAccounted(existingGlobal.capacity(),
                                      sizeof(Index)));
        addSaturating(callerMappingBytes,
                      vectorAccounted(expectedClass.capacity(),
                                      sizeof(Classification)));

        {
            std::vector<Cand> released;
            released.swap(cands);
        }

        auto imported = importTetCandidate(
            *tetCandidate, selected, existingGlobal, expectedClass, out,
            Classification{3, s->tag()});
        if (!imported.imported()) {
            if (imported.reason == TetCandidateReason::AllocationFailed)
                throw std::bad_alloc{};
            throw std::runtime_error(
                "meshVolumes: authoritative Tet import failed");
        }

        auto& importedResult = *imported.result;
        importedResult.metrics.filterConnectivityStageAccountedBytes =
            filterConnectivityBytes;
        importedResult.metrics.callerMappingAccountedBytes =
            callerMappingBytes;
        importedResult.metrics.filterRepresentationReleasedBeforePublication =
            true;

        std::uint64_t publishedPayload = vectorAccounted(
            static_cast<std::size_t>(out.nodeCount()),
            sizeof(Vec3) + sizeof(Classification));
        addSaturating(publishedPayload,
                      vectorAccounted(
                          static_cast<std::size_t>(out.cellCount()),
                          sizeof(mesh::Cell)));
        std::uint64_t publishedConnectivity = 0;
        for (Index cellIndex = 0; cellIndex < out.cellCount(); ++cellIndex)
            addSaturating(
                publishedConnectivity,
                static_cast<std::uint64_t>(mesh::nodesPerCell(
                    out.cell(cellIndex).type)));
        addSaturating(
            publishedPayload,
            vectorAccounted(
                static_cast<std::size_t>(publishedConnectivity),
                sizeof(Index)));
        importedResult.metrics.publishedPayloadAccountedBytes =
            publishedPayload;
        importedResult.metrics
            .constructionCandidateConservativeAccountedBytes =
            tetCandidate->metrics.constructionCandidateConservativeBytes;
        // Import overlap: candidate + authoritative workspace/import scratch
        // + compact caller mappings + the already-published destination.
        addSaturating(importedResult.metrics.modeledOverlapPayloadBytes,
                      callerMappingBytes);
        addSaturating(importedResult.metrics.modeledOverlapPayloadBytes,
                      publishedPayload);
        // Earlier filter stage: candidate + complete filtered connectivity +
        // compact maps + published destination. The authoritative workspace
        // did not yet exist, so this is a separate stage, not a three-way sum.
        std::uint64_t filterStageOverlap = publishedPayload;
        addSaturating(filterStageOverlap,
                      importedResult.metrics.candidateAccountedBytes);
        addSaturating(filterStageOverlap, filterConnectivityBytes);
        addSaturating(filterStageOverlap, callerMappingBytes);
        importedResult.metrics.modeledOverlapPayloadBytes = std::max(
            importedResult.metrics.modeledOverlapPayloadBytes,
            filterStageOverlap);

        // Release the remaining transient complete Tet representation before
        // the destination append reserve and its no-failure publication tail.
        importedResult.metrics.candidateReleasedBeforePublication = true;
        tetCandidate.reset();
        importedResult.metrics.completeTetRepresentationsBeforePublication = 1;
        // Complete candidate/filter/authoritative connectivity never coexist:
        // filter+candidate precede import, candidate+workspace follow it.
        accumulateHandoffWork(*handoffMetrics,importedResult.metrics);
        handoffMetrics->candidateAccountedBytes = std::max(
            handoffMetrics->candidateAccountedBytes,
            importedResult.metrics.candidateAccountedBytes);
        handoffMetrics->constructionCandidateConservativeAccountedBytes =
            std::max(
            handoffMetrics->constructionCandidateConservativeAccountedBytes,
            importedResult.metrics
                .constructionCandidateConservativeAccountedBytes);
        handoffMetrics->authoritativeModeledStorageBytes = std::max(
            handoffMetrics->authoritativeModeledStorageBytes,
            importedResult.metrics.authoritativeModeledStorageBytes);
        handoffMetrics->importScratchModeledPayloadBytes = std::max(
            handoffMetrics->importScratchModeledPayloadBytes,
            importedResult.metrics.importScratchModeledPayloadBytes);
        handoffMetrics->filterConnectivityStageAccountedBytes = std::max(
            handoffMetrics->filterConnectivityStageAccountedBytes,
            importedResult.metrics.filterConnectivityStageAccountedBytes);
        handoffMetrics->callerMappingAccountedBytes = std::max(
            handoffMetrics->callerMappingAccountedBytes,
            importedResult.metrics.callerMappingAccountedBytes);
        handoffMetrics->publishedPayloadAccountedBytes = std::max(
            handoffMetrics->publishedPayloadAccountedBytes,
            importedResult.metrics.publishedPayloadAccountedBytes);
        handoffMetrics->modeledOverlapPayloadBytes = std::max(
            handoffMetrics->modeledOverlapPayloadBytes,
            importedResult.metrics.modeledOverlapPayloadBytes);
        handoffMetrics->completeTetRepresentationsAtOverlap = std::max(
            handoffMetrics->completeTetRepresentationsAtOverlap,
            importedResult.metrics.completeTetRepresentationsAtOverlap);
        handoffMetrics->completeTetRepresentationsBeforePublication = std::max(
            handoffMetrics->completeTetRepresentationsBeforePublication,
            importedResult.metrics.completeTetRepresentationsBeforePublication);
        handoffMetrics->workspaceValid =
            handoffMetrics->workspaceValid &&
            importedResult.metrics.workspaceValid;
        handoffMetrics->candidateReleasedBeforePublication =
            handoffMetrics->candidateReleasedBeforePublication &&
            importedResult.metrics.candidateReleasedBeforePublication;
        handoffMetrics->filterRepresentationReleasedBeforePublication =
            handoffMetrics->filterRepresentationReleasedBeforePublication &&
            importedResult.metrics.filterRepresentationReleasedBeforePublication;

        std::vector<Index> connectivity(4, kInvalidIndex);
        const std::size_t liveCells = selected.size();
        if (liveCells > (std::numeric_limits<std::size_t>::max)() / 4)
            throw std::length_error("meshVolumes: Tet publication overflow");
        out.reserveAppend(0, liveCells, liveCells * 4);
        for (const Index candidateCell : selected) {
            const auto localCell = importedResult.localCellByCandidate[
                static_cast<std::size_t>(candidateCell)];
            const auto cell = importedResult.workspace->cells().get(localCell);
            for (std::size_t i = 0; i < 4; ++i)
                connectivity[i] =
                    importedResult.existingGlobalByLocalNodeSlot[
                        static_cast<std::size_t>(cell.nodes[i].slot)];
            out.addCell(CellType::Tet4, connectivity,
                        Classification{3, s->tag()});
        }
    }
}

void meshVolumes(const geo::Model& model, const sizing::MeshSpec& spec,
                 mesh::Mesh& out) {
    meshVolumesImpl(model, spec, out, nullptr);
}

void meshVolumesAuthoritative(const geo::Model& model,
                              const sizing::MeshSpec& spec,
                              mesh::Mesh& out,
                              TetHandoffMetrics& metrics) {
    metrics = TetHandoffMetrics{};
    metrics.workspaceValid = true;
    metrics.candidateReleasedBeforePublication = true;
    metrics.filterRepresentationReleasedBeforePublication = true;
    metrics.workAccountingComplete = true;
    meshVolumesImpl(model, spec, out, &metrics);
    metrics.runComplete = metrics.workspaceValid &&
                           metrics.candidateReleasedBeforePublication &&
                           metrics.filterRepresentationReleasedBeforePublication &&
                           metrics.workAccountingComplete &&
                          metrics.completeTetRepresentationsBeforePublication <= 1;
}

} // namespace paacs::algo
