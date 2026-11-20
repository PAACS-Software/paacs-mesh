// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// T2-A revision 2: non-mutating, plain-value exact geometry and cavity
// planning scaffolding. This header deliberately has no workspace mutation
// access and does not include Workspace.hpp or TransactionAccess.hpp.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <queue>
#include <set>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "paacs/core/Types.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/predicates/Predicates.hpp"

namespace paacs::bnd {

using TetFacetVertexId = std::uint32_t;
inline constexpr TetFacetVertexId kInvalidTetFacetVertex =
    std::numeric_limits<TetFacetVertexId>::max();

struct TetFacetVertex {
    Vec3 position{};
    std::uint64_t sequence = 0;
};

struct TetFacetCell {
    std::array<TetFacetVertexId, 4> vertices{};
    std::uint64_t sequence = 0;
    mesh::Classification classification{3, 0};
};

struct TetFacetTriangle {
    std::array<TetFacetVertexId, 3> vertices{};
};

struct TetFacetEdgeKey {
    TetFacetVertexId lo = kInvalidTetFacetVertex;
    TetFacetVertexId hi = kInvalidTetFacetVertex;

    friend bool operator==(const TetFacetEdgeKey&,
                           const TetFacetEdgeKey&) = default;
    friend bool operator<(const TetFacetEdgeKey& a,
                          const TetFacetEdgeKey& b) {
        return std::tie(a.lo, a.hi) < std::tie(b.lo, b.hi);
    }
};

struct TetFacetFaceKey {
    std::array<TetFacetVertexId, 3> vertices{
        kInvalidTetFacetVertex, kInvalidTetFacetVertex,
        kInvalidTetFacetVertex};

    friend bool operator==(const TetFacetFaceKey&,
                           const TetFacetFaceKey&) = default;
    friend bool operator<(const TetFacetFaceKey& a,
                          const TetFacetFaceKey& b) {
        return a.vertices < b.vertices;
    }
};

struct TetFacetProtectedFace {
    TetFacetFaceKey key{};
    std::uint64_t identity = 0;
    mesh::Classification classification{2, 0};
};

struct TetFacetComplexView {
    std::span<const TetFacetVertex> vertices;
    std::span<const TetFacetCell> cells;
    TetFacetTriangle target{};
    std::span<const TetFacetEdgeKey> protectedEdges;
    std::span<const TetFacetProtectedFace> protectedFaces;
    mesh::Classification region{3, 0};
};

enum class TetFacetContact : std::uint8_t {
    Disjoint,
    ProperInteriorCut,
    DeclaredBoundaryTouch,
    UndeclaredVertexTouch,
    TangentialTouch,
    CoplanarExactFace,
    CoplanarPartialOverlap,
    InvalidInput
};

enum class TetFacetGeometryReason : std::uint8_t {
    None,
    InvalidInput,
    NonPositiveTet,
    DuplicateCell,
    NonManifoldFace,
    InvalidFaceSense,
    DisconnectedComplex,
    BoundaryNotClosed,
    InvalidVertexLink,
    ShellingFailed,
    NoProperCut,
    UnsupportedDegeneracy,
    CrossedProtectedEdge,
    CrossedProtectedFace,
    BoundaryEdgeStar,
    NonManifoldEdgeStar,
    DisconnectedEdgeStar,
    MixedRegionStar,
    BoundaryStraddlesTarget,
    NonBallClosure,
    UnsupportedNonStarCavity,
    BoundaryMismatch,
    CandidateIntersection,
    ProtectedFaceLost,
    RegionMismatch,
    OrphanNode,
    InvalidCandidate,
    AllocationFailed
};

enum class TetFacetPlanStatus : std::uint8_t { Planned, Refused, Failed };

struct TetFacetGeometryCounters {
    std::uint64_t cellsScanned = 0;
    std::uint64_t crossedEdgeStars = 0;
    std::uint64_t closureInsertions = 0;
    std::uint64_t fixedPointPasses = 0;
    std::uint64_t apexAttempts = 0;
    std::uint64_t tetTetIntersectionTests = 0;

    friend bool operator==(const TetFacetGeometryCounters&,
                           const TetFacetGeometryCounters&) = default;
};

static_assert(std::is_standard_layout_v<TetFacetGeometryCounters>);
static_assert(std::is_trivially_copyable_v<TetFacetGeometryCounters>);

// Outcomes retain partial diagnostic counters on refusal/failure. The input
// complex is const and is never mutated.

struct TetFacetCavityPlan {
    std::vector<std::size_t> seedCells;
    std::vector<std::size_t> selectedCells;
    TetFacetGeometryCounters counters{};
};

struct TetFacetPlanOutcome {
    TetFacetPlanStatus status = TetFacetPlanStatus::Failed;
    TetFacetGeometryReason reason = TetFacetGeometryReason::InvalidInput;
    TetFacetCavityPlan plan{};

    bool planned() const noexcept {
        return status == TetFacetPlanStatus::Planned;
    }
};

struct TetFacetCandidate {
    std::vector<TetFacetCell> cells;
    TetFacetFaceKey targetKey{};
    TetFacetGeometryCounters counters{};
};

struct TetFacetCandidateOutcome {
    TetFacetPlanStatus status = TetFacetPlanStatus::Failed;
    TetFacetGeometryReason reason = TetFacetGeometryReason::InvalidInput;
    TetFacetCandidate candidate{};

    bool built() const noexcept {
        return status == TetFacetPlanStatus::Planned;
    }
};

struct TetFacetAudit {
    bool valid = false;
    TetFacetGeometryReason reason = TetFacetGeometryReason::InvalidInput;
};

namespace tet_facet_detail {

using predicates::Sign;

inline int signValue(Sign s) noexcept { return static_cast<int>(s); }

inline bool opposite(Sign a, Sign b) noexcept {
    return signValue(a) * signValue(b) < 0;
}

inline TetFacetEdgeKey edgeKey(TetFacetVertexId a,
                               TetFacetVertexId b) noexcept {
    if (b < a) std::swap(a, b);
    return {a, b};
}

inline TetFacetFaceKey faceKey(std::array<TetFacetVertexId, 3> f) noexcept {
    std::sort(f.begin(), f.end());
    return {f};
}

inline TetFacetFaceKey faceKey(const TetFacetTriangle& f) noexcept {
    return faceKey(f.vertices);
}

inline bool finite(const Vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

inline constexpr std::array<std::array<int, 2>, 6> kTetEdges{{
    {{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};

inline constexpr std::array<std::array<int, 2>, 3> kTriangleEdges{{
    {{0, 1}}, {{1, 2}}, {{2, 0}}}};

// Outward for orient3d(v0,v1,v2,v3) > 0, matching ElementTraits.
inline constexpr std::array<std::array<int, 3>, 4> kTetFaces{{
    {{1, 2, 3}}, {{0, 3, 2}}, {{0, 1, 3}}, {{0, 2, 1}}}};

inline constexpr std::array<std::array<int, 3>, 1> kTriangleFaces{{
    {{0, 1, 2}}}};

inline std::array<TetFacetVertexId, 3> outwardFace(
    const TetFacetCell& c, std::size_t localFace) noexcept {
    const auto& f = kTetFaces[localFace];
    return {c.vertices[static_cast<std::size_t>(f[0])],
            c.vertices[static_cast<std::size_t>(f[1])],
            c.vertices[static_cast<std::size_t>(f[2])]};
}

inline int faceSense(const TetFacetComplexView& view,
                     const std::array<TetFacetVertexId, 3>& directed) noexcept {
    std::size_t m = 0;
    const auto order = [&](TetFacetVertexId a, TetFacetVertexId b) {
        return std::tie(view.vertices[a].sequence, a) <
               std::tie(view.vertices[b].sequence, b);
    };
    if (order(directed[1], directed[m])) m = 1;
    if (order(directed[2], directed[m])) m = 2;
    const auto next = directed[(m + 1) % 3];
    const auto prev = directed[(m + 2) % 3];
    return order(next, prev) ? +1 : -1;
}

// An exact supporting plane restricts contact to the declared common
// vertices on one face. Undeclared coplanar vertices require further tests.
template<bool PositiveTet, std::size_t N, class Charge>
inline bool separatedByTetFaceImpl(const TetFacetComplexView& view,
                                   const TetFacetCell& cell,
                                   const std::array<TetFacetVertexId, N>& other,
                                   Charge& charge) {
    for (unsigned side=0;side<4;++side) {
        const auto face=outwardFace(cell,side);
        const auto a=view.vertices[face[0]].position;
        const auto b=view.vertices[face[1]].position;
        const auto c=view.vertices[face[2]].position;
        charge();
        Sign inward = Sign::Negative;
        if constexpr (!PositiveTet)
            inward = predicates::orient3d(
                a, b, c, view.vertices[cell.vertices[side]].position);
        bool separated=inward!=predicates::Sign::Zero;
        for (auto id:other) {
            if (!separated) break;
            if (std::ranges::find(face,id)!=face.end()) continue;
            charge();
            const auto sign=predicates::orient3d(a,b,c,view.vertices[id].position);
            separated=sign!=inward&&sign!=predicates::Sign::Zero;
        }
        if (separated) return true;
    }
    return false;
}

template<std::size_t N, class Charge>
inline bool separatedByTetFace(const TetFacetComplexView& view,
                               const TetFacetCell& cell,
                               const std::array<TetFacetVertexId, N>& other,
                               Charge charge) {
    return separatedByTetFaceImpl<false>(view, cell, other, charge);
}

// The four outward face permutations of an exactly positive tet put its
// opposing vertex on the negative side. Keep the same work accounting.
template<std::size_t N, class Charge>
inline bool separatedByPositiveTetFace(
    const TetFacetComplexView& view, const TetFacetCell& cell,
    const std::array<TetFacetVertexId, N>& other, Charge charge) {
    return separatedByTetFaceImpl<true>(view, cell, other, charge);
}

inline bool canonicalDenseSequenceMapping(
    const TetFacetComplexView& view) noexcept {
    if (view.vertices.empty()) return false;
    for (std::size_t i = 0; i < view.vertices.size(); ++i) {
        if (!finite(view.vertices[i].position)) return false;
        if (i > 0 &&
            view.vertices[i - 1].sequence >= view.vertices[i].sequence)
            return false;
    }
    return true;
}

enum class Projection : std::uint8_t { XY, YZ, ZX, Invalid };

inline Vec2 project(const Vec3& p, Projection projection) noexcept {
    switch (projection) {
        case Projection::XY: return {p.x, p.y};
        case Projection::YZ: return {p.y, p.z};
        case Projection::ZX: return {p.z, p.x};
        case Projection::Invalid: break;
    }
    return {};
}

inline Projection projectionForTriangle(const Vec3& a, const Vec3& b,
                                        const Vec3& c) {
    if (predicates::orient2d({a.x, a.y}, {b.x, b.y}, {c.x, c.y}) !=
        Sign::Zero)
        return Projection::XY;
    if (predicates::orient2d({a.y, a.z}, {b.y, b.z}, {c.y, c.z}) !=
        Sign::Zero)
        return Projection::YZ;
    if (predicates::orient2d({a.z, a.x}, {b.z, b.x}, {c.z, c.x}) !=
        Sign::Zero)
        return Projection::ZX;
    return Projection::Invalid;
}

enum class PointTri2 : std::uint8_t { Outside, Boundary, Inside };

inline PointTri2 pointInTriangle2(const Vec2& p, const Vec2& a,
                                  const Vec2& b, const Vec2& c) {
    const Sign s0 = predicates::orient2d(a, b, p);
    const Sign s1 = predicates::orient2d(b, c, p);
    const Sign s2 = predicates::orient2d(c, a, p);
    const bool nonnegative = signValue(s0) >= 0 && signValue(s1) >= 0 &&
                             signValue(s2) >= 0;
    const bool nonpositive = signValue(s0) <= 0 && signValue(s1) <= 0 &&
                             signValue(s2) <= 0;
    if (!nonnegative && !nonpositive) return PointTri2::Outside;
    if (s0 == Sign::Zero || s1 == Sign::Zero || s2 == Sign::Zero)
        return PointTri2::Boundary;
    return PointTri2::Inside;
}

inline bool onSegment2(const Vec2& p, const Vec2& a, const Vec2& b) {
    if (predicates::orient2d(a, b, p) != Sign::Zero) return false;
    return p.u >= std::min(a.u, b.u) && p.u <= std::max(a.u, b.u) &&
           p.v >= std::min(a.v, b.v) && p.v <= std::max(a.v, b.v);
}

enum class Segment2Contact : std::uint8_t { None, Touch, Proper, Overlap };

inline Segment2Contact segmentContact2(const Vec2& a, const Vec2& b,
                                       const Vec2& c, const Vec2& d) {
    const Sign abC = predicates::orient2d(a, b, c);
    const Sign abD = predicates::orient2d(a, b, d);
    const Sign cdA = predicates::orient2d(c, d, a);
    const Sign cdB = predicates::orient2d(c, d, b);
    if (opposite(abC, abD) && opposite(cdA, cdB))
        return Segment2Contact::Proper;
    if (abC == Sign::Zero && abD == Sign::Zero && cdA == Sign::Zero &&
        cdB == Sign::Zero) {
        const double au0 = std::min(a.u, b.u), au1 = std::max(a.u, b.u);
        const double av0 = std::min(a.v, b.v), av1 = std::max(a.v, b.v);
        const double cu0 = std::min(c.u, d.u), cu1 = std::max(c.u, d.u);
        const double cv0 = std::min(c.v, d.v), cv1 = std::max(c.v, d.v);
        const bool uInterior = std::max(au0, cu0) < std::min(au1, cu1);
        const bool vInterior = std::max(av0, cv0) < std::min(av1, cv1);
        const bool nonPointA = a.u != b.u || a.v != b.v;
        const bool overlap = nonPointA &&
            ((a.u != b.u && uInterior) || (a.u == b.u && vInterior));
        if (overlap) return Segment2Contact::Overlap;
    }
    if ((abC == Sign::Zero && onSegment2(c, a, b)) ||
        (abD == Sign::Zero && onSegment2(d, a, b)) ||
        (cdA == Sign::Zero && onSegment2(a, c, d)) ||
        (cdB == Sign::Zero && onSegment2(b, c, d)))
        return Segment2Contact::Touch;
    return Segment2Contact::None;
}

inline bool sameNonzero(Sign a, Sign b, Sign c) noexcept {
    return a != Sign::Zero && a == b && b == c;
}

inline bool segmentCrossesOpenTriangle(
    const Vec3& p, const Vec3& q, const Vec3& a, const Vec3& b,
    const Vec3& c) {
    const Sign sp = predicates::orient3d(a, b, c, p);
    const Sign sq = predicates::orient3d(a, b, c, q);
    if (!opposite(sp, sq)) return false;
    const Sign e0 = predicates::orient3d(p, q, a, b);
    const Sign e1 = predicates::orient3d(p, q, b, c);
    const Sign e2 = predicates::orient3d(p, q, c, a);
    return sameNonzero(e0, e1, e2);
}

inline bool segmentTouchesTriangleBoundary(
    const Vec3& p, const Vec3& q, const Vec3& a, const Vec3& b,
    const Vec3& c) {
    const Sign sp = predicates::orient3d(a, b, c, p);
    const Sign sq = predicates::orient3d(a, b, c, q);
    if (!opposite(sp, sq)) return false;
    const Sign e0 = predicates::orient3d(p, q, a, b);
    const Sign e1 = predicates::orient3d(p, q, b, c);
    const Sign e2 = predicates::orient3d(p, q, c, a);
    const bool nonnegative = signValue(e0) >= 0 && signValue(e1) >= 0 &&
                             signValue(e2) >= 0;
    const bool nonpositive = signValue(e0) <= 0 && signValue(e1) <= 0 &&
                             signValue(e2) <= 0;
    return (nonnegative || nonpositive) &&
           (e0 == Sign::Zero || e1 == Sign::Zero || e2 == Sign::Zero);
}

inline bool pointStrictlyInsideTet(
    const Vec3& p, const std::array<Vec3, 4>& t) {
    for (std::size_t f = 0; f < 4; ++f) {
        const auto& lf = kTetFaces[f];
        const Sign pointSide = predicates::orient3d(
            t[static_cast<std::size_t>(lf[0])],
            t[static_cast<std::size_t>(lf[1])],
            t[static_cast<std::size_t>(lf[2])], p);
        if (pointSide != Sign::Negative) return false;
    }
    return true;
}

inline bool trianglesOpenIntersect(
    const std::array<Vec3, 3>& a, const std::array<Vec3, 3>& b) {
    for (std::size_t i = 0; i < 3; ++i)
        if (segmentCrossesOpenTriangle(
                a[i], a[(i + 1) % 3], b[0], b[1], b[2]))
            return true;
    for (std::size_t i = 0; i < 3; ++i)
        if (segmentCrossesOpenTriangle(
                b[i], b[(i + 1) % 3], a[0], a[1], a[2]))
            return true;
    return false;
}

inline bool coplanarTrianglesHaveAreaOverlap(
    const std::array<Vec3, 3>& a, const std::array<Vec3, 3>& b) {
    const Projection projection = projectionForTriangle(a[0], a[1], a[2]);
    if (projection == Projection::Invalid) return false;
    for (const Vec3& p : b)
        if (predicates::orient3d(a[0], a[1], a[2], p) != Sign::Zero)
            return false;
    std::array<Vec2, 3> ap{}, bp{};
    for (std::size_t i = 0; i < 3; ++i) {
        ap[i] = project(a[i], projection);
        bp[i] = project(b[i], projection);
    }
    const auto samePoint = [](const Vec2& lhs, const Vec2& rhs) {
        return lhs.u == rhs.u && lhs.v == rhs.v;
    };
    bool sameTriangle = true;
    for (const Vec2& p : ap) {
        bool found = false;
        for (const Vec2& q : bp) found = found || samePoint(p, q);
        sameTriangle = sameTriangle && found;
    }
    if (sameTriangle) return true;
    for (const Vec2& p : ap)
        if (pointInTriangle2(p, bp[0], bp[1], bp[2]) == PointTri2::Inside)
            return true;
    for (const Vec2& p : bp)
        if (pointInTriangle2(p, ap[0], ap[1], ap[2]) == PointTri2::Inside)
            return true;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            const auto contact = segmentContact2(
                ap[i], ap[(i + 1) % 3], bp[j], bp[(j + 1) % 3]);
            if (contact == Segment2Contact::Proper)
                return true;
        }
    return false;
}

inline bool pointOnTriangleClosed(const Vec3& point,
                                  const std::array<Vec3, 3>& triangle) {
    if (predicates::orient3d(triangle[0], triangle[1], triangle[2], point) !=
        Sign::Zero)
        return false;
    const Projection projection = projectionForTriangle(
        triangle[0], triangle[1], triangle[2]);
    if (projection == Projection::Invalid) return false;
    return pointInTriangle2(project(point, projection),
                            project(triangle[0], projection),
                            project(triangle[1], projection),
                            project(triangle[2], projection)) !=
           PointTri2::Outside;
}

inline Segment2Contact segmentContact3(const Vec3& a, const Vec3& b,
                                       const Vec3& c, const Vec3& d) {
    if (predicates::orient3d(a, b, c, d) != Sign::Zero)
        return Segment2Contact::None;
    constexpr std::array<Projection, 3> projections{
        Projection::XY, Projection::YZ, Projection::ZX};
    for (const auto projection : projections) {
        const Vec2 ap = project(a, projection);
        const Vec2 bp = project(b, projection);
        const Vec2 cp = project(c, projection);
        const Vec2 dp = project(d, projection);
        if (predicates::orient2d(ap, bp, cp) != Sign::Zero ||
            predicates::orient2d(ap, bp, dp) != Sign::Zero)
            return segmentContact2(ap, bp, cp, dp);
    }
    for (const auto projection : projections) {
        const Vec2 ap = project(a, projection);
        const Vec2 bp = project(b, projection);
        if (ap.u != bp.u || ap.v != bp.v)
            return segmentContact2(ap, bp, project(c, projection),
                                   project(d, projection));
    }
    return Segment2Contact::None;
}

struct IdentityAwareClosedContact {
    bool legalIdentityContact = false;
    bool aVertexOnBFace = false;
    bool bVertexOnAFace = false;
    bool edgeEdgeContact = false;
    bool faceAreaOverlap = false;

    bool illegal() const noexcept {
        return aVertexOnBFace || bVertexOnAFace || edgeEdgeContact ||
               faceAreaOverlap;
    }
};

template <std::size_t NA, std::size_t NEA, std::size_t NFA,
          std::size_t NB, std::size_t NEB, std::size_t NFB>
inline IdentityAwareClosedContact auditIdentityAwareClosedContact(
    const TetFacetComplexView& view,
    const std::array<TetFacetVertexId, NA>& aVertices,
    const std::array<std::array<int, 2>, NEA>& aEdges,
    const std::array<std::array<int, 3>, NFA>& aFaces,
    const std::array<TetFacetVertexId, NB>& bVertices,
    const std::array<std::array<int, 2>, NEB>& bEdges,
    const std::array<std::array<int, 3>, NFB>& bFaces) {
    IdentityAwareClosedContact audit;
    const auto faceIds = [](const auto& vertices, const auto& local) {
        return std::array<TetFacetVertexId, 3>{
            vertices[static_cast<std::size_t>(local[0])],
            vertices[static_cast<std::size_t>(local[1])],
            vertices[static_cast<std::size_t>(local[2])]};
    };
    const auto points = [&](const std::array<TetFacetVertexId, 3>& face) {
        return std::array<Vec3, 3>{view.vertices[face[0]].position,
                                   view.vertices[face[1]].position,
                                   view.vertices[face[2]].position};
    };
    const auto inFace = [](const auto& face, TetFacetVertexId vertex) {
        return std::find(face.begin(), face.end(), vertex) != face.end();
    };

    for (const auto vertex : aVertices)
        for (const auto& local : bFaces) {
            const auto face = faceIds(bVertices, local);
            if (!pointOnTriangleClosed(view.vertices[vertex].position,
                                       points(face)))
                continue;
            if (inFace(face, vertex))
                audit.legalIdentityContact = true;
            else
                audit.aVertexOnBFace = true;
        }
    for (const auto vertex : bVertices)
        for (const auto& local : aFaces) {
            const auto face = faceIds(aVertices, local);
            if (!pointOnTriangleClosed(view.vertices[vertex].position,
                                       points(face)))
                continue;
            if (inFace(face, vertex))
                audit.legalIdentityContact = true;
            else
                audit.bVertexOnAFace = true;
        }

    for (const auto& aLocal : aEdges)
        for (const auto& bLocal : bEdges) {
            const auto a0 = aVertices[static_cast<std::size_t>(aLocal[0])];
            const auto a1 = aVertices[static_cast<std::size_t>(aLocal[1])];
            const auto b0 = bVertices[static_cast<std::size_t>(bLocal[0])];
            const auto b1 = bVertices[static_cast<std::size_t>(bLocal[1])];
            const auto contact = segmentContact3(
                view.vertices[a0].position, view.vertices[a1].position,
                view.vertices[b0].position, view.vertices[b1].position);
            if (contact == Segment2Contact::None) continue;
            const auto aKey = edgeKey(a0, a1);
            const auto bKey = edgeKey(b0, b1);
            if (aKey == bKey) {
                audit.legalIdentityContact = true;
                continue;
            }
            const bool sharedVertex = a0 == b0 || a0 == b1 || a1 == b0 ||
                                      a1 == b1;
            if (contact == Segment2Contact::Touch && sharedVertex)
                audit.legalIdentityContact = true;
            else
                audit.edgeEdgeContact = true;
        }

    for (const auto& aLocal : aFaces)
        for (const auto& bLocal : bFaces) {
            const auto aFace = faceIds(aVertices, aLocal);
            const auto bFace = faceIds(bVertices, bLocal);
            if (!coplanarTrianglesHaveAreaOverlap(points(aFace),
                                                   points(bFace)))
                continue;
            if (faceKey(aFace) == faceKey(bFace))
                audit.legalIdentityContact = true;
            else
                audit.faceAreaOverlap = true;
        }
    return audit;
}

inline bool validTriangleGeometry(const TetFacetComplexView& view) {
    const auto& ids = view.target.vertices;
    if (ids[0] == ids[1] || ids[1] == ids[2] || ids[2] == ids[0])
        return false;
    for (const auto id : ids)
        if (id >= view.vertices.size() ||
            !finite(view.vertices[id].position))
            return false;
    const auto& a = view.vertices[ids[0]].position;
    const auto& b = view.vertices[ids[1]].position;
    const auto& c = view.vertices[ids[2]].position;
    return projectionForTriangle(a, b, c) != Projection::Invalid;
}

inline bool validTriangle(const TetFacetComplexView& view) {
    return canonicalDenseSequenceMapping(view) && validTriangleGeometry(view);
}

inline bool validCell(const TetFacetComplexView& view,
                      const TetFacetCell& cell) {
    auto sorted = cell.vertices;
    std::sort(sorted.begin(), sorted.end());
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
        return false;
    std::array<Vec3, 4> p{};
    for (std::size_t i = 0; i < 4; ++i) {
        if (cell.vertices[i] >= view.vertices.size()) return false;
        p[i] = view.vertices[cell.vertices[i]].position;
        if (!finite(p[i])) return false;
    }
    return predicates::orient3d(p[0], p[1], p[2], p[3]) ==
           Sign::Positive;
}

inline std::array<Vec3, 3> trianglePoints(
    const TetFacetComplexView& view, const TetFacetTriangle& triangle) {
    return {view.vertices[triangle.vertices[0]].position,
            view.vertices[triangle.vertices[1]].position,
            view.vertices[triangle.vertices[2]].position};
}

inline std::array<Vec3, 3> facePoints(
    const TetFacetComplexView& view,
    const std::array<TetFacetVertexId, 3>& face) {
    return {view.vertices[face[0]].position,
            view.vertices[face[1]].position,
            view.vertices[face[2]].position};
}

}  // namespace tet_facet_detail

inline TetFacetEdgeKey canonicalTetFacetEdge(TetFacetVertexId a,
                                             TetFacetVertexId b) noexcept {
    return tet_facet_detail::edgeKey(a, b);
}

inline TetFacetFaceKey canonicalTetFacetFace(
    std::array<TetFacetVertexId, 3> vertices) noexcept {
    return tet_facet_detail::faceKey(vertices);
}

namespace tet_facet_detail {

// Callers must have validated the finite, strictly ordered vertex table.
// The target triangle and candidate cell still receive their local checks.
inline TetFacetContact classifyTetFacetContactWithCanonicalVertices(
    const TetFacetComplexView& view, const TetFacetCell& cell) {
    if (!validTriangleGeometry(view) || !validCell(view, cell))
        return TetFacetContact::InvalidInput;

    const auto targetKey = faceKey(view.target);
    for (std::size_t f = 0; f < 4; ++f)
        if (faceKey(outwardFace(cell, f)) == targetKey)
            return TetFacetContact::CoplanarExactFace;

    const auto target = trianglePoints(view, view.target);
    std::array<Vec3, 4> tet{};
    std::array<Sign, 4> planeSigns{};
    bool anyZero = false;
    bool anyPositive = false;
    bool anyNegative = false;
    for (std::size_t i = 0; i < 4; ++i) {
        tet[i] = view.vertices[cell.vertices[i]].position;
        planeSigns[i] = predicates::orient3d(
            target[0], target[1], target[2], tet[i]);
        anyZero = anyZero || planeSigns[i] == Sign::Zero;
        anyPositive = anyPositive || planeSigns[i] == Sign::Positive;
        anyNegative = anyNegative || planeSigns[i] == Sign::Negative;
    }
    if (!anyZero && !(anyPositive && anyNegative))
        return TetFacetContact::Disjoint;

    for (const auto targetId : view.target.vertices)
        if (pointStrictlyInsideTet(view.vertices[targetId].position, tet))
            return TetFacetContact::ProperInteriorCut;

    bool tangential = false;
    for (const auto& edge : kTetEdges) {
        const Vec3& p = tet[static_cast<std::size_t>(edge[0])];
        const Vec3& q = tet[static_cast<std::size_t>(edge[1])];
        if (segmentCrossesOpenTriangle(p, q, target[0], target[1], target[2]))
            return TetFacetContact::ProperInteriorCut;
        tangential = tangential || segmentTouchesTriangleBoundary(
            p, q, target[0], target[1], target[2]);
    }
    for (std::size_t te = 0; te < 3; ++te) {
        const Vec3& p = target[te];
        const Vec3& q = target[(te + 1) % 3];
        for (std::size_t f = 0; f < 4; ++f) {
            const auto fp = facePoints(view, outwardFace(cell, f));
            if (segmentCrossesOpenTriangle(p, q, fp[0], fp[1], fp[2]))
                return TetFacetContact::ProperInteriorCut;
            tangential = tangential || segmentTouchesTriangleBoundary(
                p, q, fp[0], fp[1], fp[2]);
        }
    }

    const auto closed = auditIdentityAwareClosedContact(
        view, view.target.vertices, kTriangleEdges, kTriangleFaces,
        cell.vertices, kTetEdges, kTetFaces);
    if (closed.faceAreaOverlap)
        return TetFacetContact::CoplanarPartialOverlap;
    if (closed.bVertexOnAFace)
        return TetFacetContact::UndeclaredVertexTouch;
    if (closed.aVertexOnBFace || closed.edgeEdgeContact || tangential)
        return TetFacetContact::TangentialTouch;
    if (closed.legalIdentityContact)
        return TetFacetContact::DeclaredBoundaryTouch;
    return TetFacetContact::Disjoint;
}

}  // namespace tet_facet_detail

inline TetFacetContact classifyTetFacetContact(
    const TetFacetComplexView& view, const TetFacetCell& cell) {
    if (!tet_facet_detail::canonicalDenseSequenceMapping(view))
        return TetFacetContact::InvalidInput;
    return tet_facet_detail::classifyTetFacetContactWithCanonicalVertices(
        view, cell);
}

namespace tet_facet_detail {

struct FaceUse {
    std::size_t cell = 0;
    std::array<TetFacetVertexId, 3> directed{};
    int sense = 0;
};

using EdgeStarMap = std::map<TetFacetEdgeKey, std::vector<std::size_t>>;
using FaceUseMap = std::map<TetFacetFaceKey, std::vector<FaceUse>>;

inline std::array<TetFacetVertexId, 4> cellKey(
    const TetFacetCell& cell) noexcept {
    auto key = cell.vertices;
    std::sort(key.begin(), key.end());
    return key;
}

inline bool sameClassification(const mesh::Classification& a,
                               const mesh::Classification& b) noexcept {
    return a.dim == b.dim && a.tag == b.tag;
}

inline bool buildMaps(const TetFacetComplexView& view, EdgeStarMap& edges,
                      FaceUseMap& faces, TetFacetGeometryReason& reason) {
    std::set<std::array<TetFacetVertexId, 4>> cellKeys;
    for (std::size_t ci = 0; ci < view.cells.size(); ++ci) {
        const auto& cell = view.cells[ci];
        if (!validCell(view, cell)) {
            reason = TetFacetGeometryReason::NonPositiveTet;
            return false;
        }
        if (!cellKeys.insert(cellKey(cell)).second) {
            reason = TetFacetGeometryReason::DuplicateCell;
            return false;
        }
        for (const auto& edge : kTetEdges) {
            const auto key = edgeKey(
                cell.vertices[static_cast<std::size_t>(edge[0])],
                cell.vertices[static_cast<std::size_t>(edge[1])]);
            edges[key].push_back(ci);
        }
        for (std::size_t f = 0; f < 4; ++f) {
            const auto directed = outwardFace(cell, f);
            faces[faceKey(directed)].push_back(
                FaceUse{ci, directed, faceSense(view, directed)});
        }
    }
    for (const auto& [key, uses] : faces) {
        (void)key;
        if (uses.size() > 2) {
            reason = TetFacetGeometryReason::NonManifoldFace;
            return false;
        }
        if (uses.size() == 2 && uses[0].sense * uses[1].sense != -1) {
            reason = TetFacetGeometryReason::InvalidFaceSense;
            return false;
        }
    }
    return true;
}

inline bool protectedEdge(const TetFacetComplexView& view,
                          const TetFacetEdgeKey& key) noexcept {
    return std::find(view.protectedEdges.begin(), view.protectedEdges.end(),
                     key) != view.protectedEdges.end();
}

inline bool protectedFace(const TetFacetComplexView& view,
                          const TetFacetFaceKey& key) noexcept {
    return std::find_if(
               view.protectedFaces.begin(), view.protectedFaces.end(),
               [&](const TetFacetProtectedFace& face) {
                   return face.key == key;
               }) != view.protectedFaces.end();
}

inline std::array<TetFacetEdgeKey, 3> targetEdges(
    const TetFacetTriangle& target) noexcept {
    return {edgeKey(target.vertices[0], target.vertices[1]),
            edgeKey(target.vertices[1], target.vertices[2]),
            edgeKey(target.vertices[2], target.vertices[0])};
}

inline bool isTargetEdge(const TetFacetTriangle& target,
                         const TetFacetEdgeKey& key) noexcept {
    const auto keys = targetEdges(target);
    return std::find(keys.begin(), keys.end(), key) != keys.end();
}

inline TetFacetGeometryReason auditEdgeStar(
    const TetFacetComplexView& view, const TetFacetEdgeKey& key,
    const std::vector<std::size_t>& star) {
    if (star.size() < 3) return TetFacetGeometryReason::BoundaryEdgeStar;
    std::map<TetFacetVertexId, std::vector<TetFacetVertexId>> link;
    for (const std::size_t ci : star) {
        if (ci >= view.cells.size())
            return TetFacetGeometryReason::InvalidInput;
        const auto& cell = view.cells[ci];
        if (!sameClassification(cell.classification, view.region))
            return TetFacetGeometryReason::MixedRegionStar;
        std::array<TetFacetVertexId, 2> oppositeVertices{};
        std::size_t n = 0;
        for (const auto vertex : cell.vertices)
            if (vertex != key.lo && vertex != key.hi) {
                if (n >= oppositeVertices.size())
                    return TetFacetGeometryReason::NonManifoldEdgeStar;
                oppositeVertices[n++] = vertex;
            }
        if (n != 2 || oppositeVertices[0] == oppositeVertices[1])
            return TetFacetGeometryReason::NonManifoldEdgeStar;
        link[oppositeVertices[0]].push_back(oppositeVertices[1]);
        link[oppositeVertices[1]].push_back(oppositeVertices[0]);
    }
    bool hasDegreeOne = false;
    for (auto& [vertex, neighbours] : link) {
        (void)vertex;
        std::sort(neighbours.begin(), neighbours.end());
        if (std::adjacent_find(neighbours.begin(), neighbours.end()) !=
            neighbours.end())
            return TetFacetGeometryReason::NonManifoldEdgeStar;
        if (neighbours.size() == 1) hasDegreeOne = true;
        if (neighbours.size() != 2)
            return hasDegreeOne ? TetFacetGeometryReason::BoundaryEdgeStar
                                : TetFacetGeometryReason::NonManifoldEdgeStar;
    }
    std::set<TetFacetVertexId> visited;
    std::queue<TetFacetVertexId> pending;
    pending.push(link.begin()->first);
    visited.insert(link.begin()->first);
    while (!pending.empty()) {
        const auto vertex = pending.front();
        pending.pop();
        for (const auto neighbour : link[vertex])
            if (visited.insert(neighbour).second) pending.push(neighbour);
    }
    if (visited.size() != link.size())
        return TetFacetGeometryReason::DisconnectedEdgeStar;
    return TetFacetGeometryReason::None;
}

inline bool selectedConnected(const FaceUseMap& faces,
                              const std::vector<unsigned char>& selected,
                              std::size_t expected) {
    if (expected == 0) return false;
    std::vector<std::vector<std::size_t>> dual(selected.size());
    for (const auto& [key, uses] : faces) {
        (void)key;
        if (uses.size() != 2) continue;
        if (selected[uses[0].cell] && selected[uses[1].cell]) {
            dual[uses[0].cell].push_back(uses[1].cell);
            dual[uses[1].cell].push_back(uses[0].cell);
        }
    }
    std::size_t first = selected.size();
    for (std::size_t i = 0; i < selected.size(); ++i)
        if (selected[i]) {
            first = i;
            break;
        }
    std::vector<unsigned char> seen(selected.size(), 0);
    std::queue<std::size_t> pending;
    pending.push(first);
    seen[first] = 1;
    std::size_t count = 0;
    while (!pending.empty()) {
        const auto cell = pending.front();
        pending.pop();
        ++count;
        for (const auto neighbour : dual[cell])
            if (!seen[neighbour]) {
                seen[neighbour] = 1;
                pending.push(neighbour);
            }
    }
    return count == expected;
}

inline bool triangleStraddlesTarget(
    const TetFacetComplexView& view,
    const std::array<TetFacetVertexId, 3>& face) {
    return trianglesOpenIntersect(facePoints(view, face),
                                  trianglePoints(view, view.target));
}

inline std::vector<std::size_t> orderedCells(
    const TetFacetComplexView& view,
    const std::vector<unsigned char>& selected) {
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < selected.size(); ++i)
        if (selected[i]) out.push_back(i);
    std::sort(out.begin(), out.end(), [&](std::size_t a, std::size_t b) {
        const auto& ca = view.cells[a];
        const auto& cb = view.cells[b];
        return std::tie(ca.sequence, a) < std::tie(cb.sequence, b);
    });
    return out;
}

inline bool activeConnected(const FaceUseMap& faces,
                            const std::vector<unsigned char>& active,
                            std::size_t expected) {
    return selectedConnected(faces, active, expected);
}

inline TetFacetAudit auditBallImpl(const TetFacetComplexView& view,
                                   std::span<const std::size_t> cells) {
    TetFacetAudit audit;
    if (cells.empty() || !canonicalDenseSequenceMapping(view) ||
        view.region.dim != 3 || view.region.tag < 0) {
        audit.reason = TetFacetGeometryReason::NonBallClosure;
        return audit;
    }
    std::vector<unsigned char> selected(view.cells.size(), 0);
    std::set<std::array<TetFacetVertexId, 4>> keys;
    for (const auto ci : cells) {
        if (ci >= view.cells.size() || selected[ci] ||
            !validCell(view, view.cells[ci])) {
            audit.reason = TetFacetGeometryReason::InvalidInput;
            return audit;
        }
        if (!sameClassification(view.cells[ci].classification,
                                view.region)) {
            audit.reason = TetFacetGeometryReason::RegionMismatch;
            return audit;
        }
        selected[ci] = 1;
        if (!keys.insert(cellKey(view.cells[ci])).second) {
            audit.reason = TetFacetGeometryReason::DuplicateCell;
            return audit;
        }
    }
    FaceUseMap faces;
    for (const auto ci : cells)
        for (std::size_t f = 0; f < 4; ++f) {
            const auto directed = outwardFace(view.cells[ci], f);
            faces[faceKey(directed)].push_back(
                FaceUse{ci, directed, faceSense(view, directed)});
        }
    for (const auto& [key, uses] : faces) {
        (void)key;
        if (uses.size() > 2) {
            audit.reason = TetFacetGeometryReason::NonManifoldFace;
            return audit;
        }
        if (uses.size() == 2 && uses[0].sense * uses[1].sense != -1) {
            audit.reason = TetFacetGeometryReason::InvalidFaceSense;
            return audit;
        }
    }
    if (!activeConnected(faces, selected, cells.size())) {
        audit.reason = TetFacetGeometryReason::DisconnectedComplex;
        return audit;
    }

    std::map<TetFacetEdgeKey, std::size_t> boundaryEdges;
    std::map<TetFacetVertexId,
             std::vector<std::pair<TetFacetVertexId, TetFacetVertexId>>>
        vertexLinks;
    for (const auto& [key, uses] : faces) {
        if (uses.size() != 1) continue;
        const auto& f = key.vertices;
        for (std::size_t e = 0; e < 3; ++e)
            ++boundaryEdges[edgeKey(f[e], f[(e + 1) % 3])];
        vertexLinks[f[0]].push_back({f[1], f[2]});
        vertexLinks[f[1]].push_back({f[0], f[2]});
        vertexLinks[f[2]].push_back({f[0], f[1]});
    }
    if (boundaryEdges.empty()) {
        audit.reason = TetFacetGeometryReason::BoundaryNotClosed;
        return audit;
    }
    for (const auto& [edge, count] : boundaryEdges) {
        (void)edge;
        if (count != 2) {
            audit.reason = TetFacetGeometryReason::BoundaryNotClosed;
            return audit;
        }
    }
    for (const auto& [vertex, linkEdges] : vertexLinks) {
        (void)vertex;
        std::map<TetFacetVertexId, std::vector<TetFacetVertexId>> graph;
        for (const auto& edge : linkEdges) {
            graph[edge.first].push_back(edge.second);
            graph[edge.second].push_back(edge.first);
        }
        for (auto& [node, neighbours] : graph) {
            (void)node;
            std::sort(neighbours.begin(), neighbours.end());
            if (neighbours.size() != 2 ||
                std::adjacent_find(neighbours.begin(), neighbours.end()) !=
                    neighbours.end()) {
                audit.reason = TetFacetGeometryReason::InvalidVertexLink;
                return audit;
            }
        }
        std::set<TetFacetVertexId> seen;
        std::queue<TetFacetVertexId> pending;
        pending.push(graph.begin()->first);
        seen.insert(graph.begin()->first);
        while (!pending.empty()) {
            const auto node = pending.front();
            pending.pop();
            for (const auto neighbour : graph[node])
                if (seen.insert(neighbour).second) pending.push(neighbour);
        }
        if (seen.size() != graph.size()) {
            audit.reason = TetFacetGeometryReason::InvalidVertexLink;
            return audit;
        }
    }

    // Deterministic sufficient shelling certificate.
    auto active = selected;
    std::size_t remaining = cells.size();
    const auto order = orderedCells(view, selected);
    while (remaining > 1) {
        bool peeled = false;
        for (const auto ci : order) {
            if (!active[ci]) continue;
            std::size_t sharedFaces = 0;
            for (std::size_t f = 0; f < 4; ++f) {
                const auto it = faces.find(faceKey(outwardFace(view.cells[ci], f)));
                if (it == faces.end() || it->second.size() != 2) continue;
                const auto other = it->second[0].cell == ci
                                       ? it->second[1].cell
                                       : it->second[0].cell;
                if (active[other]) ++sharedFaces;
            }
            if (sharedFaces < 1 || sharedFaces > 3) continue;
            active[ci] = 0;
            if (activeConnected(faces, active, remaining - 1)) {
                --remaining;
                peeled = true;
                break;
            }
            active[ci] = 1;
        }
        if (!peeled) {
            audit.reason = TetFacetGeometryReason::ShellingFailed;
            return audit;
        }
    }
    audit.valid = true;
    audit.reason = TetFacetGeometryReason::None;
    return audit;
}

}  // namespace tet_facet_detail

inline TetFacetAudit auditTetFacetBall(
    const TetFacetComplexView& view,
    std::span<const std::size_t> selectedCells) noexcept {
    try {
        return tet_facet_detail::auditBallImpl(view, selectedCells);
    } catch (const std::bad_alloc&) {
        return {false, TetFacetGeometryReason::AllocationFailed};
    } catch (...) {
        return {false, TetFacetGeometryReason::InvalidInput};
    }
}

namespace tet_facet_detail {

inline TetFacetPlanOutcome planTetFacetCavityImpl(
    const TetFacetComplexView& view,
    std::span<const std::size_t> explicitProperCutSeeds = {}) noexcept {
    TetFacetPlanOutcome outcome;
    try {
        if (!validTriangle(view) || view.cells.empty() ||
            view.region.dim != 3 || view.region.tag < 0) {
            outcome.reason = TetFacetGeometryReason::InvalidInput;
            return outcome;
        }
        EdgeStarMap edgeStars;
        FaceUseMap faceUses;
        if (!buildMaps(view, edgeStars, faceUses, outcome.reason)) return outcome;

        std::vector<unsigned char> selected(view.cells.size(), 0);
        std::vector<unsigned char> seeds(view.cells.size(), 0);
        for (std::size_t ci = 0; ci < view.cells.size(); ++ci) {
            ++outcome.plan.counters.cellsScanned;
            if (!sameClassification(view.cells[ci].classification, view.region)) {
                outcome.reason = TetFacetGeometryReason::MixedRegionStar;
                outcome.status = TetFacetPlanStatus::Refused;
                return outcome;
            }
            const auto contact = classifyTetFacetContactWithCanonicalVertices(
                view, view.cells[ci]);
            if (contact == TetFacetContact::InvalidInput) {
                outcome.reason = TetFacetGeometryReason::InvalidInput;
                return outcome;
            }
            if (contact == TetFacetContact::ProperInteriorCut &&
                explicitProperCutSeeds.empty())
                seeds[ci] = 1;
            if (contact == TetFacetContact::UndeclaredVertexTouch ||
                contact == TetFacetContact::TangentialTouch ||
                contact == TetFacetContact::CoplanarPartialOverlap) {
                outcome.status = TetFacetPlanStatus::Refused;
                outcome.reason = TetFacetGeometryReason::UnsupportedDegeneracy;
                return outcome;
            }
        }
        if (!explicitProperCutSeeds.empty())
            for (const auto ci : explicitProperCutSeeds) {
                if (ci >= view.cells.size()) {
                    outcome.reason = TetFacetGeometryReason::InvalidInput;
                    return outcome;
                }
                seeds[ci] = 1;
            }
        selected = seeds;
        outcome.plan.seedCells = orderedCells(view, seeds);
        if (outcome.plan.seedCells.empty()) {
            outcome.status = TetFacetPlanStatus::Refused;
            outcome.reason = TetFacetGeometryReason::NoProperCut;
            return outcome;
        }

        const auto target = trianglePoints(view, view.target);
        for (const auto& [key, star] : edgeStars) {
            if (isTargetEdge(view.target, key)) continue;
            const Vec3& a = view.vertices[key.lo].position;
            const Vec3& b = view.vertices[key.hi].position;
            if (!segmentCrossesOpenTriangle(
                    a, b, target[0], target[1], target[2]))
                continue;
            ++outcome.plan.counters.crossedEdgeStars;
            if (protectedEdge(view, key)) {
                outcome.status = TetFacetPlanStatus::Refused;
                outcome.reason = TetFacetGeometryReason::CrossedProtectedEdge;
                return outcome;
            }
            const auto starReason = auditEdgeStar(view, key, star);
            if (starReason != TetFacetGeometryReason::None) {
                outcome.status = TetFacetPlanStatus::Refused;
                outcome.reason = starReason;
                return outcome;
            }
            for (const auto ci : star)
                if (!selected[ci]) {
                    selected[ci] = 1;
                    ++outcome.plan.counters.closureInsertions;
                }
        }

        bool changed = true;
        while (changed) {
            changed = false;
            ++outcome.plan.counters.fixedPointPasses;
            for (const auto& [key, uses] : faceUses) {
                std::size_t selectedOwners = 0;
                std::size_t selectedCell = view.cells.size();
                for (const auto& use : uses)
                    if (selected[use.cell]) {
                        ++selectedOwners;
                        selectedCell = use.cell;
                    }
                if (selectedOwners != 1) continue;
                const auto selectedUse = std::find_if(
                    uses.begin(), uses.end(), [&](const FaceUse& use) {
                        return use.cell == selectedCell;
                    });
                if (selectedUse == uses.end() ||
                    !triangleStraddlesTarget(view, selectedUse->directed))
                    continue;
                if (protectedFace(view, key)) {
                    outcome.status = TetFacetPlanStatus::Refused;
                    outcome.reason = TetFacetGeometryReason::CrossedProtectedFace;
                    return outcome;
                }
                if (uses.size() == 1) {
                    outcome.status = TetFacetPlanStatus::Refused;
                    outcome.reason =
                        TetFacetGeometryReason::BoundaryStraddlesTarget;
                    return outcome;
                }
                if (uses.size() != 2) {
                    outcome.status = TetFacetPlanStatus::Refused;
                    outcome.reason = TetFacetGeometryReason::NonManifoldFace;
                    return outcome;
                }
                if (uses[0].sense * uses[1].sense != -1) {
                    outcome.status = TetFacetPlanStatus::Refused;
                    outcome.reason = TetFacetGeometryReason::InvalidFaceSense;
                    return outcome;
                }
                const auto other = uses[0].cell == selectedCell
                                       ? uses[1].cell
                                       : uses[0].cell;
                if (!sameClassification(view.cells[other].classification,
                                        view.region)) {
                    outcome.status = TetFacetPlanStatus::Refused;
                    outcome.reason = TetFacetGeometryReason::MixedRegionStar;
                    return outcome;
                }
                if (!selected[other]) {
                    selected[other] = 1;
                    ++outcome.plan.counters.closureInsertions;
                    changed = true;
                }
            }
        }

        const auto count = static_cast<std::size_t>(
            std::count(selected.begin(), selected.end(),
                       static_cast<unsigned char>(1)));
        if (!selectedConnected(faceUses, selected, count)) {
            outcome.status = TetFacetPlanStatus::Refused;
            outcome.reason = TetFacetGeometryReason::DisconnectedComplex;
            return outcome;
        }
        outcome.plan.selectedCells = orderedCells(view, selected);
        const auto ball = auditBallImpl(view, outcome.plan.selectedCells);
        if (!ball.valid) {
            outcome.status = TetFacetPlanStatus::Refused;
            outcome.reason = TetFacetGeometryReason::NonBallClosure;
            return outcome;
        }
        outcome.status = TetFacetPlanStatus::Planned;
        outcome.reason = TetFacetGeometryReason::None;
        return outcome;
    } catch (const std::bad_alloc&) {
        outcome.status = TetFacetPlanStatus::Failed;
        outcome.reason = TetFacetGeometryReason::AllocationFailed;
        return outcome;
    } catch (...) {
        outcome.status = TetFacetPlanStatus::Failed;
        outcome.reason = TetFacetGeometryReason::InvalidInput;
        return outcome;
    }
}

}  // namespace tet_facet_detail

inline TetFacetPlanOutcome planTetFacetCavity(
    const TetFacetComplexView& view) noexcept {
    return tet_facet_detail::planTetFacetCavityImpl(view);
}

// Indexed T2-B entry.  The supplied view is the resolver-complete local
// cavity shell assembled by TetRecoveryContext; global survivor and orphan
// proofs are deliberately performed by validateTetFacetCandidateIndexed()
// from separate indexed query/use-count views.  Keeping this as an additive
// entry prevents the accepted T2-A public call from being fed a queried
// subset or from acquiring different semantics.
inline TetFacetPlanOutcome planTetFacetCavityIndexed(
    const TetFacetComplexView& resolverCompleteLocalView) noexcept {
    return tet_facet_detail::planTetFacetCavityImpl(
        resolverCompleteLocalView);
}

#if defined(PAACS_TOPO_TESTING)
inline TetFacetPlanOutcome planTetFacetCavityWithSeedsForTesting(
    const TetFacetComplexView& view,
    std::span<const std::size_t> explicitProperCutSeeds) noexcept {
    return tet_facet_detail::planTetFacetCavityImpl(
        view, explicitProperCutSeeds);
}
#endif

namespace tet_facet_detail {

inline bool tetImproperlyIntersects(const TetFacetComplexView& view,
                                    const TetFacetCell& a,
                                    const TetFacetCell& b) {
    if (cellKey(a) == cellKey(b)) return true;
    std::array<Vec3, 4> ap{}, bp{};
    for (std::size_t i = 0; i < 4; ++i) {
        ap[i] = view.vertices[a.vertices[i]].position;
        bp[i] = view.vertices[b.vertices[i]].position;
    }
    const auto closed = auditIdentityAwareClosedContact(
        view, a.vertices, kTetEdges, kTetFaces, b.vertices, kTetEdges,
        kTetFaces);
    if (closed.illegal()) return true;
    for (const Vec3& p : ap)
        if (pointStrictlyInsideTet(p, bp)) return true;
    for (const Vec3& p : bp)
        if (pointStrictlyInsideTet(p, ap)) return true;
    for (const auto& edge : kTetEdges)
        for (std::size_t f = 0; f < 4; ++f) {
            const auto face = outwardFace(b, f);
            const auto fp = facePoints(view, face);
            if (segmentCrossesOpenTriangle(
                    ap[static_cast<std::size_t>(edge[0])],
                    ap[static_cast<std::size_t>(edge[1])], fp[0], fp[1],
                    fp[2]))
                return true;
        }
    for (const auto& edge : kTetEdges)
        for (std::size_t f = 0; f < 4; ++f) {
            const auto face = outwardFace(a, f);
            const auto fp = facePoints(view, face);
            if (segmentCrossesOpenTriangle(
                    bp[static_cast<std::size_t>(edge[0])],
                    bp[static_cast<std::size_t>(edge[1])], fp[0], fp[1],
                    fp[2]))
                return true;
        }
    for (std::size_t af = 0; af < 4; ++af)
        for (std::size_t bf = 0; bf < 4; ++bf) {
            const auto aFace = outwardFace(a, af);
            const auto bFace = outwardFace(b, bf);
            if (faceKey(aFace) == faceKey(bFace)) {
                if (faceSense(view, aFace) * faceSense(view, bFace) != -1)
                    return true;
            }
        }
    return false;
}

inline FaceUseMap selectedFaceMap(const TetFacetComplexView& view,
                                  std::span<const std::size_t> cells) {
    FaceUseMap faces;
    for (const auto ci : cells)
        for (std::size_t f = 0; f < 4; ++f) {
            const auto directed = outwardFace(view.cells[ci], f);
            faces[faceKey(directed)].push_back(
                FaceUse{ci, directed, faceSense(view, directed)});
        }
    return faces;
}

inline FaceUseMap candidateFaceMap(const TetFacetComplexView& view,
                                   std::span<const TetFacetCell> cells) {
    FaceUseMap faces;
    for (std::size_t ci = 0; ci < cells.size(); ++ci)
        for (std::size_t f = 0; f < 4; ++f) {
            const auto directed = outwardFace(cells[ci], f);
            faces[faceKey(directed)].push_back(
                FaceUse{ci, directed, faceSense(view, directed)});
        }
    return faces;
}

inline std::map<TetFacetFaceKey, int> boundarySenses(
    const FaceUseMap& faces) {
    std::map<TetFacetFaceKey, int> boundary;
    for (const auto& [key, uses] : faces)
        if (uses.size() == 1) boundary.emplace(key, uses[0].sense);
    return boundary;
}

inline int facePlaneSide(const TetFacetComplexView& view,
                         const TetFacetFaceKey& face) {
    const auto target = trianglePoints(view, view.target);
    int side = 0;
    for (const auto id : face.vertices) {
        const auto sign = predicates::orient3d(
            target[0], target[1], target[2], view.vertices[id].position);
        const int value = signValue(sign);
        if (value == 0) continue;
        if (side != 0 && side != value) return 2;
        side = value;
    }
    return side;
}

inline std::optional<TetFacetCell> makePositiveCell(
    const TetFacetComplexView& view, TetFacetVertexId apex,
    const std::array<TetFacetVertexId, 3>& face,
    std::uint64_t sequence) {
    TetFacetCell cell;
    cell.vertices = {apex, face[0], face[1], face[2]};
    cell.sequence = sequence;
    cell.classification = view.region;
    const auto& a = view.vertices[cell.vertices[0]].position;
    const auto& b = view.vertices[cell.vertices[1]].position;
    const auto& c = view.vertices[cell.vertices[2]].position;
    const auto& d = view.vertices[cell.vertices[3]].position;
    const auto orientation = predicates::orient3d(a, b, c, d);
    if (orientation == Sign::Zero) return std::nullopt;
    if (orientation == Sign::Negative)
        std::swap(cell.vertices[1], cell.vertices[2]);
    return cell;
}

inline bool boundaryMatches(
    const FaceUseMap& candidateFaces,
    const std::map<TetFacetFaceKey, int>& expected,
    const TetFacetFaceKey& targetKey) {
    const auto actual = boundarySenses(candidateFaces);
    if (actual.size() != expected.size()) return false;
    for (const auto& [key, sense] : expected) {
        const auto it = actual.find(key);
        if (it == actual.end()) return false;
        if (!(key == targetKey) && it->second != sense) return false;
    }
    return true;
}

inline std::vector<TetFacetVertexId> cavityVertices(
    const TetFacetComplexView& view,
    std::span<const std::size_t> selectedCells) {
    std::vector<TetFacetVertexId> vertices;
    for (const auto ci : selectedCells)
        for (const auto vertex : view.cells[ci].vertices)
            vertices.push_back(vertex);
    std::sort(vertices.begin(), vertices.end());
    vertices.erase(std::unique(vertices.begin(), vertices.end()),
                   vertices.end());
    std::sort(vertices.begin(), vertices.end(),
              [&](TetFacetVertexId a, TetFacetVertexId b) {
                  return std::tie(view.vertices[a].sequence, a) <
                         std::tie(view.vertices[b].sequence, b);
              });
    return vertices;
}

inline TetFacetAudit validateCandidateImpl(
    const TetFacetComplexView& source,
    std::span<const std::size_t> removed,
    std::span<const TetFacetCell> candidate,
    TetFacetGeometryCounters& counters) {
    TetFacetAudit audit;
    if (candidate.empty()) {
        audit.reason = TetFacetGeometryReason::InvalidInput;
        return audit;
    }
    TetFacetComplexView candidateView{
        source.vertices, candidate, source.target, source.protectedEdges,
        source.protectedFaces, source.region};
    std::vector<std::size_t> all(candidate.size());
    for (std::size_t i = 0; i < all.size(); ++i) all[i] = i;
    audit = auditBallImpl(candidateView, all);
    if (!audit.valid) return audit;

    const auto sourceBoundary =
        boundarySenses(selectedFaceMap(source, removed));
    const auto candidateFaces = candidateFaceMap(candidateView, candidate);
    const auto candidateBoundary = boundarySenses(candidateFaces);
    if (sourceBoundary != candidateBoundary) {
        audit.valid = false;
        audit.reason = TetFacetGeometryReason::BoundaryMismatch;
        return audit;
    }
    const auto targetKey = faceKey(source.target);
    const auto targetIt = candidateFaces.find(targetKey);
    if (targetIt == candidateFaces.end() || targetIt->second.size() != 2 ||
        targetIt->second[0].sense * targetIt->second[1].sense != -1) {
        audit.valid = false;
        audit.reason = TetFacetGeometryReason::BoundaryMismatch;
        return audit;
    }

    std::vector<unsigned char> removedMask(source.cells.size(), 0);
    std::vector<unsigned char> removedNodeUse(source.vertices.size(), 0);
    std::vector<unsigned char> finalNodeUse(source.vertices.size(), 0);
    for (const auto ci : removed) {
        if (ci >= source.cells.size() || removedMask[ci]) {
            audit.valid = false;
            audit.reason = TetFacetGeometryReason::InvalidInput;
            return audit;
        }
        removedMask[ci] = 1;
        for (const auto vertex : source.cells[ci].vertices)
            removedNodeUse[vertex] = 1;
    }
    for (std::size_t ci = 0; ci < source.cells.size(); ++ci)
        if (!removedMask[ci])
            for (const auto vertex : source.cells[ci].vertices)
                finalNodeUse[vertex] = 1;
    for (const auto& cell : candidate)
        for (const auto vertex : cell.vertices)
            finalNodeUse[vertex] = 1;
    for (std::size_t vertex = 0; vertex < removedNodeUse.size(); ++vertex)
        if (removedNodeUse[vertex] && !finalNodeUse[vertex]) {
            audit.valid = false;
            audit.reason = TetFacetGeometryReason::OrphanNode;
            return audit;
        }

    for (const auto& protectedFaceSpec : source.protectedFaces) {
        bool affected = false;
        const auto sourceFaces = selectedFaceMap(source, removed);
        affected = sourceFaces.find(protectedFaceSpec.key) != sourceFaces.end();
        if (affected && candidateFaces.find(protectedFaceSpec.key) ==
                            candidateFaces.end()) {
            audit.valid = false;
            audit.reason = TetFacetGeometryReason::ProtectedFaceLost;
            return audit;
        }
    }

    for (std::size_t i = 0; i < candidate.size(); ++i)
        for (std::size_t j = i + 1; j < candidate.size(); ++j) {
            ++counters.tetTetIntersectionTests;
            if (tetImproperlyIntersects(candidateView, candidate[i],
                                        candidate[j])) {
                audit.valid = false;
                audit.reason = TetFacetGeometryReason::CandidateIntersection;
                return audit;
            }
        }
    for (const auto& created : candidate)
        for (std::size_t ci = 0; ci < source.cells.size(); ++ci) {
            if (removedMask[ci]) continue;
            ++counters.tetTetIntersectionTests;
            if (tetImproperlyIntersects(source, created, source.cells[ci])) {
                audit.valid = false;
                audit.reason = TetFacetGeometryReason::CandidateIntersection;
                return audit;
            }
        }
    audit.valid = true;
    audit.reason = TetFacetGeometryReason::None;
    return audit;
}

}  // namespace tet_facet_detail

inline TetFacetAudit validateTetFacetCandidate(
    const TetFacetComplexView& source,
    std::span<const std::size_t> removedCells,
    std::span<const TetFacetCell> candidateCells,
    TetFacetGeometryCounters& counters) noexcept {
    try {
        return tet_facet_detail::validateCandidateImpl(
            source, removedCells, candidateCells, counters);
    } catch (const std::bad_alloc&) {
        return {false, TetFacetGeometryReason::AllocationFailed};
    } catch (...) {
        return {false, TetFacetGeometryReason::InvalidInput};
    }
}

// Additive indexed validation entry.  Its `completeProofView` is assembled
// from the selected cavity, complete node stars and every survivor returned
// by the closed candidate-AABB queries.  Untouched-node use counts are
// audited by TetRecoveryContext before this call.
inline TetFacetAudit validateTetFacetCandidateIndexed(
    const TetFacetComplexView& completeProofView,
    std::span<const std::size_t> removedCells,
    std::span<const TetFacetCell> candidateCells,
    TetFacetGeometryCounters& counters) noexcept {
    try {
        return tet_facet_detail::validateCandidateImpl(
            completeProofView, removedCells, candidateCells, counters);
    } catch (const std::bad_alloc&) {
        return {false, TetFacetGeometryReason::AllocationFailed};
    } catch (...) {
        return {false, TetFacetGeometryReason::InvalidInput};
    }
}

inline TetFacetCandidateOutcome reconstructTetFacetCandidate(
    const TetFacetComplexView& view,
    const TetFacetCavityPlan& plan) noexcept {
    using namespace tet_facet_detail;
    TetFacetCandidateOutcome outcome;
    outcome.candidate.counters = plan.counters;
    try {
        if (plan.selectedCells.empty()) {
            outcome.reason = TetFacetGeometryReason::InvalidInput;
            return outcome;
        }
        const auto sourceFaces = selectedFaceMap(view, plan.selectedCells);
        const auto sourceBoundary = boundarySenses(sourceFaces);
        const auto targetKey = faceKey(view.target);
        outcome.candidate.targetKey = targetKey;
        if (sourceBoundary.find(targetKey) != sourceBoundary.end()) {
            outcome.status = TetFacetPlanStatus::Refused;
            outcome.reason = TetFacetGeometryReason::UnsupportedDegeneracy;
            return outcome;
        }

        std::array<std::map<TetFacetFaceKey, int>, 2> sideBoundary;
        for (const auto& [key, sense] : sourceBoundary) {
            const int side = facePlaneSide(view, key);
            if (side == 2 || side == 0) {
                outcome.status = TetFacetPlanStatus::Refused;
                outcome.reason = TetFacetGeometryReason::UnsupportedDegeneracy;
                return outcome;
            }
            sideBoundary[side > 0 ? 1U : 0U].emplace(key, sense);
        }
        sideBoundary[0].emplace(targetKey, 0);
        sideBoundary[1].emplace(targetKey, 0);

        const auto cavityNodes = cavityVertices(view, plan.selectedCells);
        struct SideTrial {
            TetFacetVertexId apex = kInvalidTetFacetVertex;
            std::vector<TetFacetCell> cells;
        };
        std::array<std::vector<SideTrial>, 2> trials;
        const auto targetPointsValue = trianglePoints(view, view.target);
        for (std::size_t sideIndex = 0; sideIndex < 2; ++sideIndex) {
            const int desiredSide = sideIndex == 0 ? -1 : +1;
            for (const auto apex : cavityNodes) {
                const auto sign = predicates::orient3d(
                    targetPointsValue[0], targetPointsValue[1],
                    targetPointsValue[2], view.vertices[apex].position);
                if (signValue(sign) != desiredSide) continue;
                ++outcome.candidate.counters.apexAttempts;
                std::vector<TetFacetCell> trial;
                std::set<std::array<TetFacetVertexId, 4>> keys;
                bool valid = true;
                for (const auto& [key, sense] : sideBoundary[sideIndex]) {
                    (void)sense;
                    if (std::find(key.vertices.begin(), key.vertices.end(),
                                  apex) != key.vertices.end())
                        continue;
                    auto cell = makePositiveCell(
                        view, apex, key.vertices, 1 + trial.size());
                    if (!cell || !keys.insert(cellKey(*cell)).second) {
                        valid = false;
                        break;
                    }
                    trial.push_back(*cell);
                }
                if (!valid || trial.empty()) continue;
                TetFacetComplexView trialView{
                    view.vertices, trial, view.target, view.protectedEdges,
                    view.protectedFaces, view.region};
                std::vector<std::size_t> trialIds(trial.size());
                for (std::size_t i = 0; i < trialIds.size(); ++i)
                    trialIds[i] = i;
                const auto ball = auditBallImpl(trialView, trialIds);
                if (!ball.valid) continue;
                const auto trialFaces = candidateFaceMap(trialView, trial);
                if (!boundaryMatches(trialFaces, sideBoundary[sideIndex],
                                     targetKey))
                    continue;
                trials[sideIndex].push_back(
                    SideTrial{apex, std::move(trial)});
            }
            if (trials[sideIndex].empty()) {
                outcome.status = TetFacetPlanStatus::Refused;
                outcome.reason =
                    TetFacetGeometryReason::UnsupportedNonStarCavity;
                return outcome;
            }
        }

        bool sawProtectedFaceLoss = false;
        for (const auto& negative : trials[0])
            for (const auto& positive : trials[1]) {
                (void)negative.apex;
                (void)positive.apex;
                std::vector<TetFacetCell> combined;
                combined.reserve(negative.cells.size() +
                                 positive.cells.size());
                combined.insert(combined.end(), negative.cells.begin(),
                                negative.cells.end());
                combined.insert(combined.end(), positive.cells.begin(),
                                positive.cells.end());
                for (std::size_t i = 0; i < combined.size(); ++i)
                    combined[i].sequence = 1 + i;
                const auto audit = validateCandidateImpl(
                    view, plan.selectedCells, combined,
                    outcome.candidate.counters);
                if (audit.valid) {
                    outcome.candidate.cells = std::move(combined);
                    outcome.status = TetFacetPlanStatus::Planned;
                    outcome.reason = TetFacetGeometryReason::None;
                    return outcome;
                }
                sawProtectedFaceLoss =
                    sawProtectedFaceLoss ||
                    audit.reason == TetFacetGeometryReason::ProtectedFaceLost;
            }
        if (sawProtectedFaceLoss) {
            outcome.status = TetFacetPlanStatus::Refused;
            outcome.reason = TetFacetGeometryReason::ProtectedFaceLost;
            return outcome;
        }
        outcome.status = TetFacetPlanStatus::Failed;
        outcome.reason = TetFacetGeometryReason::InvalidCandidate;
        return outcome;
    } catch (const std::bad_alloc&) {
        outcome.status = TetFacetPlanStatus::Failed;
        outcome.reason = TetFacetGeometryReason::AllocationFailed;
        return outcome;
    } catch (...) {
        outcome.status = TetFacetPlanStatus::Failed;
        outcome.reason = TetFacetGeometryReason::InvalidInput;
        return outcome;
    }
}

}  // namespace paacs::bnd
