// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// ============================================================================
// BND-2A - constrained Steiner-enabled segment recovery (2D, Tri3).
// Design authority: docs/PMBND002_DESIGN.md rev 3.
//
// This FOUNDATION slice of Increment 3 provides the parametric (BND-layer)
// artifacts that the single-Steiner recovery route composes, WITHOUT any
// workspace mutation of its own:
//   * the GeometryProvider abstraction (identity + revision + pointAt), §4;
//   * the portable EXACT symbolic parameter model AuthoritativeParam
//     (Stored | Affine), §5.2-5.4, evaluated through the accepted
//     predicates::orient2dAffine / orderOnSegmentAffine kernel (no rounded
//     mapping coordinate ever enters a predicate or ordering decision);
//   * canonical undirected endpoint identity + portable dyadic ordering
//     (uint64, no __int128), §5.2-5.3;
//   * the SteinerRecoveryCandidate value types (SteinerNodeId, SteinerNode,
//     TriNodeRef, TriNodeRefTriangle), §1.1;
//   * the §9 outcome-reason taxonomy for the Steiner route;
//   * lowering a validated candidate into exactly one topo::SteinerCavitySpec
//     (staged authoritative param -> derived Vec3 via pointAt; localIds
//     assigned from the structural SteinerNodeId in deterministic order),
//     §1.2 - the ONLY bridge into paacs::topo. No ChartId / AuthoritativeParam
//     / orient2d over a parametric domain crosses into topology.
//
// The recovery ROUTE, its six geometric gates over this model, the batched
// ASCR commit and the allocation-free shadow-metadata publication build on
// these types (subsequent slice). Mutation happens exclusively through the
// accepted topo::SteinerCavityTransaction primitive.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "paacs/core/Types.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/topo/SteinerCavity.hpp"
#include "paacs/bnd/SegmentRecovery.hpp"   // ConstraintVertexId, ChartId,
                                           // ParamPoint, OwnershipStamp,
                                           // ConstraintSegment, ConstraintSet,
                                           // preflightConstraintSet
#include "paacs/bnd/SegmentRecoverySession.hpp"   // RecoverySession::open
                                                 // (full A1+A2+PSLG+index)
#include "paacs/bnd/SegmentRecoveryCommit.hpp"   // routeEntity + recoverSegments
                                                 // (mark-legality routing reuse)

namespace paacs::bnd {

// ---------------------------------------------------------------------------
// §4 - geometry provider (mapping authority). Identity+revision stamped;
// pointAt is used ONLY to derive the physical Vec3 of a mapping coordinate,
// NEVER for a topology predicate or ordering decision.
// ---------------------------------------------------------------------------
struct GeometryProviderStamp {
    std::uint64_t identity = 0;
    std::uint64_t revision = 0;
    bool operator==(const GeometryProviderStamp& o) const {
        return identity == o.identity && revision == o.revision;
    }
};
// §4 NARROW invalid-input signal. A provider that can recognize an
// out-of-domain parameter throws THIS specific type; the route classifies it
// as Rejected(ProviderDomainRejected) (input, not a tier-1 failure). Any
// OTHER exception is Failed(GeometryProviderError); std::bad_alloc is
// Failed(CandidateAllocFailed). This is the entire, closed taxonomy — there
// is no implicit "provider contract identifies invalid input" back-channel.
struct GeometryProviderDomainError : std::runtime_error {
    explicit GeometryProviderDomainError(const std::string& what)
        : std::runtime_error(what) {}
};
class GeometryProvider {
public:
    virtual ~GeometryProvider() = default;
    virtual GeometryProviderStamp stamp() const = 0;
    // May throw GeometryProviderDomainError (-> Rejected(ProviderDomainRejected)),
    // std::bad_alloc (-> Failed(CandidateAllocFailed)), or any other exception
    // (-> Failed(GeometryProviderError)) per the §4 taxonomy; the route owns
    // the classification. Must return a finite Vec3 for a finite in-domain
    // param, else the route classifies Rejected(NonfiniteMappedPosition).
    virtual Vec3 pointAt(ChartId chart, const Vec2& uv) const = 0;
};

// ---------------------------------------------------------------------------
// §5.2 - authoritative symbolic parameter. Stored = an existing node's
// arbitrary stored uv (need not be dyadic). Affine = a staged node at
// t = num/2^level from the CANONICAL undirected endpoints (lo,hi) = (min,max)
// by ConstraintVertexId, so reversing input (b,a) with t' = 1-t yields the
// SAME reduced tuple (invariant 15). num reduced (odd); the level is capped at
// kMaxSteinerLevel == predicates::kMaxAffineLevel (52) because EVERY one of the
// six gates evaluates through the exact orient2dAffine kernel, which is
// qualified only through level 52 (its integer weights must stay exactly
// representable as double). The portable uint64 dyadic ordering/collinearity
// would fit higher levels, but the model does NOT permit them.
// ---------------------------------------------------------------------------
struct StoredParam {
    Vec2 uv{};
};
struct AffineParam {
    ConstraintVertexId lo = 0, hi = 0;   // canonical: lo < hi
    std::uint64_t num = 0;               // reduced (odd), 0 < num < 2^level
    std::uint32_t level = 0;             // 1 <= level <= kMaxSteinerLevel
};
using AuthoritativeParam = std::variant<StoredParam, AffineParam>;

// The recovery level limit is the exact-kernel limit (blocker 2): every one
// of the six gates evaluates through predicates::orient2dAffine, which is
// exactly qualified ONLY through kMaxAffineLevel (its integer weights must
// stay exactly representable as double). Levels above it are NOT supported.
inline constexpr std::uint32_t kMaxSteinerLevel = predicates::kMaxAffineLevel;

inline void rejAffine(const char* why) {
    throw std::invalid_argument(std::string("bnd::AffineParam: ") + why);
}
// Centralized DYADIC-level check used BEFORE any shift or narrowing. `1u<<L`
// is only ever formed after L is confirmed <= kMaxSteinerLevel (<=52), so no
// shift can be undefined and no uint32->uint8 narrowing can alias.
inline void requireValidLevel(std::uint32_t level) {
    if (level < 1 || level > kMaxSteinerLevel)
        rejAffine("level must be in [1, kMaxSteinerLevel]");
}
inline void requireValidNumerator(std::uint64_t num, std::uint32_t level) {
    requireValidLevel(level);                       // guards the shift below
    if (num == 0 || num >= (std::uint64_t(1) << level))
        rejAffine("numerator must satisfy 0 < num < 2^level");
}
// Fully validate a canonical (already-reduced) AffineParam BEFORE it maps,
// compares, or enters a predicate. CANONICAL endpoint order (lo < hi, so the
// undirected identity is unique); valid level; interior reduced-odd numerator.
inline void validateAffineParam(const AffineParam& ap) {
    if (!(ap.lo < ap.hi))
        rejAffine("endpoints must be canonical (lo < hi)");
    requireValidNumerator(ap.num, ap.level);        // level + range (guarded)
    if ((ap.num & 1u) == 0) rejAffine("numerator must be reduced (odd)");
}
inline void requireFiniteVec2(const Vec2& p, const char* what) {
    if (!(std::isfinite(p.u) && std::isfinite(p.v)))
        throw std::invalid_argument(
            std::string("bnd::SteinerRecovery: non-finite ") + what);
}

// Canonicalize a staged split of the UNDIRECTED segment {a,b} at parameter
// t = num/2^level measured from `a`. VALIDATES its raw input first (distinct
// ids, valid level, 0 < num < 2^level), then keys by (min,max) and reduces to
// odd; reversing (a,b) yields the identical tuple. Throws std::invalid_argument
// on any malformed input (before any shift).
inline AffineParam makeAffine(ConstraintVertexId a, ConstraintVertexId b,
                              std::uint64_t num, std::uint32_t level) {
    if (a == b) rejAffine("endpoints must be distinct");
    requireValidNumerator(num, level);              // level + range (guarded)
    AffineParam ap;
    if (a < b) {
        ap.lo = a; ap.hi = b; ap.num = num;
    } else {
        ap.lo = b; ap.hi = a;
        ap.num = (std::uint64_t(1) << level) - num;   // t' = 1 - t
    }
    ap.level = level;
    // reduce to odd numerator (canonical): strip common factors of two
    while (ap.level > 1 && (ap.num & 1u) == 0) {
        ap.num >>= 1;
        --ap.level;
    }
    return ap;
}

// §5.3 portable dyadic compare of t1=num1/2^L1 vs t2=num2/2^L2 WITHOUT
// __int128: align to Lmax; aligned numerators num_i * 2^(Lmax-Li) < 2^Lmax
// <= 2^kMaxSteinerLevel fit uint64. VALIDATES both tokens BEFORE any shift.
// Returns -1/0/+1.
inline int compareDyadic(std::uint64_t num1, std::uint32_t l1,
                         std::uint64_t num2, std::uint32_t l2) {
    requireValidNumerator(num1, l1);
    requireValidNumerator(num2, l2);
    const std::uint32_t lmax = l1 > l2 ? l1 : l2;
    const std::uint64_t a = num1 << (lmax - l1);
    const std::uint64_t b = num2 << (lmax - l2);
    return a < b ? -1 : (a > b ? 1 : 0);
}

// ---------------------------------------------------------------------------
// BND-2A Increment 4 - CANONICAL, REVERSAL-SYMMETRIC dyadic subdivision.
// Staged Steiner positions are ALWAYS canonical dyadic parameters (never the
// non-dyadic crossed-edge intersections). For a subdivision demand of D staged
// nodes on the undirected segment {a,b}, choose the coarsest level
// L = ceil(log2(D+1)) and place D numerators SYMMETRIC about the midpoint
// 2^(L-1) at denominator 2^L: the centre (if D is odd) plus symmetric pairs
// (c-k, c+k) outward. makeAffine() then keys each to canonical (min,max)
// identity and reduces to an odd numerator, so:
//   * reversing the input endpoints (b,a) maps every numerator num -> 2^L-num,
//     which is exactly the same symmetric set -> identical reduced tokens
//     (reversal symmetry, invariant 15);
//   * the result is independent of input order and allocator history.
// D==1 -> {1/2}; D==2 -> {1/4, 3/4}; D==3 -> {1/4, 1/2, 3/4}; D==4 ->
// {1/4, 3/8, 5/8, 3/4}; etc. Returns params ASCENDING by t. Throws
// std::invalid_argument if the required level exceeds kMaxSteinerLevel.
// The dyadic level for a subdivision demand D (smallest L with 2^L >= D+1).
// Throws std::invalid_argument if it would exceed the exact-kernel limit.
inline std::uint32_t dyadicLevelFor(std::uint32_t D) {
    std::uint32_t L = 1;
    while ((std::uint64_t(1) << L) < static_cast<std::uint64_t>(D) + 1) {
        ++L;
        if (L > kMaxSteinerLevel)
            throw std::invalid_argument(
                "bnd::dyadicLevelFor: subdivision demand exceeds the "
                "exact-kernel level limit");
    }
    return L;
}
// Visit the D canonical dyadic numerators (at level L) for demand D WITHOUT
// allocating - used by the budget/counting pass so the staged demand (after
// reuse coincidences are removed) can be enforced BEFORE any subdivision vector
// / provider call / workspace mutation (blocker 2, budget-before-allocation).
template <class Fn>
inline void forEachCanonicalDyadicNumerator(std::uint32_t D, Fn&& fn) {
    if (D == 0) return;
    const std::uint32_t L = dyadicLevelFor(D);
    const std::uint64_t c = (std::uint64_t(1) << L) / 2;
    std::uint32_t n = 0;
    if (D % 2 == 1) { fn(c, L); ++n; }
    for (std::uint64_t k = 1; n < D; ++k) {
        fn(c - k, L); ++n;
        if (n < D) { fn(c + k, L); ++n; }
    }
}

inline std::vector<AffineParam> canonicalDyadicParams(ConstraintVertexId a,
                                                      ConstraintVertexId b,
                                                      std::uint32_t D) {
    std::vector<AffineParam> out;
    if (a == b)
        rejAffine("endpoints must be distinct");   // even when D == 0
    if (D == 0) return out;
    const std::uint32_t L = dyadicLevelFor(D);
    const std::uint64_t den = std::uint64_t(1) << L;
    const std::uint64_t c = den / 2;
    std::vector<std::uint64_t> nums;
    nums.reserve(D);
    if (D % 2 == 1) nums.push_back(c);
    for (std::uint64_t k = 1; nums.size() < D; ++k) {
        nums.push_back(c - k);
        if (nums.size() < D) nums.push_back(c + k);
    }
    std::sort(nums.begin(), nums.end());
    for (const std::uint64_t num : nums)
        out.push_back(makeAffine(a, b, num, L));   // canonical + reduced
    // ascending by parameter (portable dyadic compare on the canonical tokens)
    std::sort(out.begin(), out.end(),
              [](const AffineParam& x, const AffineParam& y) {
                  return compareDyadic(x.num, x.level, y.num, y.level) < 0;
              });
    return out;
}

// Exact minimum-feature test: are two canonical dyadic parameters separated by
// at least 2^-minLevel along the segment's [0,1] parameter? Portable uint64
// (no __int128, no rounded mapping coordinate): align both numerators and the
// threshold to Lmax = max(l1,l2,minLevel) (<= kMaxSteinerLevel so shifts fit).
inline bool dyadicSeparated(const AffineParam& p, const AffineParam& q,
                            std::uint32_t minLevel) {
    validateAffineParam(p);                 // both tokens valid BEFORE shifting
    validateAffineParam(q);
    if (!(p.lo == q.lo && p.hi == q.hi))
        throw std::invalid_argument(
            "bnd::dyadicSeparated: parameters from different undirected "
            "segments (endpoint identities differ)");
    if (minLevel < 1 || minLevel > kMaxSteinerLevel)
        throw std::invalid_argument(
            "bnd::dyadicSeparated: minLevel out of [1, kMaxSteinerLevel]");
    const std::uint32_t lmax =
        std::max(std::max(p.level, q.level), minLevel);
    const std::uint64_t ap = p.num << (lmax - p.level);
    const std::uint64_t aq = q.num << (lmax - q.level);
    const std::uint64_t diff = ap > aq ? ap - aq : aq - ap;
    return diff >= (std::uint64_t(1) << (lmax - minLevel));
}

// ---------------------------------------------------------------------------
// Increment 4 - staged-aware GEOMETRIC ear clipping over an INDEX ring into a
// parallel AUTHORITATIVE-coordinate array AP of predicates::AffinePoint2 (an
// existing node is its exact stored param via storedPoint(); a staged node is
// its exact AffineParam). EVERY topology decision - orientation, apex ranking,
// point-in-triangle, on-edge, diagonal visibility - is made with the EXACT
// affine kernel (orient2dAffine / compareOnAxisAffine); a staged node's rounded
// mapping coordinate NEVER controls topology (blocker 1). Faithful port of
// BND-1's earClipOriented (13.9): convex apex in the polygon's OWN winding, no
// vertex inside/on the diagonal, diagonal crosses no nonadjacent edge; exact
// lexicographic-apex clipping order (allocator-independent); then normalize each
// triangle to the chart orientation `want`. nullopt on degenerate/non-simple.
// ---------------------------------------------------------------------------
namespace scdetail {
using predicates::AffinePoint2;
using predicates::Sign;

// Exact lexicographic order of two affine points: <0 / 0 / >0.
inline int lexCmp(const AffinePoint2& a, const AffinePoint2& b) {
    const Sign su = predicates::compareOnAxisAffine(a, b, 0);
    if (su != Sign::Zero) return su == Sign::Negative ? -1 : 1;
    const Sign sv = predicates::compareOnAxisAffine(a, b, 1);
    return sv == Sign::Zero ? 0 : (sv == Sign::Negative ? -1 : 1);
}
// p on the CLOSED segment [a,b], exactly (collinear + between). Uses an axis on
// which a != b (a,b are distinct ring vertices).
inline bool onClosedSeg(const AffinePoint2& a, const AffinePoint2& b,
                        const AffinePoint2& p) {
    if (predicates::orient2dAffine(a, b, p) != Sign::Zero) return false;
    int ax = predicates::compareOnAxisAffine(a, b, 0) != Sign::Zero ? 0 : 1;
    const Sign pa = predicates::compareOnAxisAffine(p, a, ax);
    const Sign pb = predicates::compareOnAxisAffine(p, b, ax);
    if (pa == Sign::Zero || pb == Sign::Zero) return true;   // p == a or p == b
    return pa != pb;                                         // strictly between
}
// segments (a,b) and (c,d) cross at an interior point of BOTH (strict).
inline bool properCross(const AffinePoint2& a, const AffinePoint2& b,
                        const AffinePoint2& c, const AffinePoint2& d) {
    const Sign d1 = predicates::orient2dAffine(a, b, c);
    const Sign d2 = predicates::orient2dAffine(a, b, d);
    const Sign d3 = predicates::orient2dAffine(c, d, a);
    const Sign d4 = predicates::orient2dAffine(c, d, b);
    return d1 != d2 && d1 != Sign::Zero && d2 != Sign::Zero && d3 != d4 &&
           d3 != Sign::Zero && d4 != Sign::Zero;
}
inline bool earValid(const std::vector<AffinePoint2>& AP,
                     const std::vector<std::size_t>& ring, std::size_t i,
                     Sign want) {
    const std::size_t n = ring.size();
    const std::size_t ip = ring[(i + n - 1) % n], ic = ring[i],
                      in = ring[(i + 1) % n];
    const AffinePoint2 &pp = AP[ip], &pc = AP[ic], &pn = AP[in];
    if (predicates::orient2dAffine(pp, pc, pn) != want) return false;   // (a)
    for (std::size_t j = 0; j < n; ++j) {
        const std::size_t rj = ring[j];
        if (rj == ip || rj == ic || rj == in) continue;
        const AffinePoint2& pq = AP[rj];
        const Sign s1 = predicates::orient2dAffine(pp, pc, pq);
        const Sign s2 = predicates::orient2dAffine(pc, pn, pq);
        const Sign s3 = predicates::orient2dAffine(pn, pp, pq);
        if (s1 != Sign::Zero && s1 == s2 && s2 == s3) return false;    // (b)
        if (onClosedSeg(pp, pn, pq)) return false;                     // (c)
    }
    for (std::size_t e = 0; e < n; ++e) {                              // (d)
        const std::size_t r = ring[e], s = ring[(e + 1) % n];
        if (r == ip || r == in || s == ip || s == in) continue;
        if (properCross(pp, pn, AP[r], AP[s])) return false;
    }
    return true;
}
inline Sign winding(const std::vector<AffinePoint2>& AP,
                    const std::vector<std::size_t>& ring) {
    const std::size_t n = ring.size();
    std::size_t lo = 0;
    for (std::size_t i = 1; i < n; ++i)
        if (lexCmp(AP[ring[i]], AP[ring[lo]]) < 0) lo = i;
    return predicates::orient2dAffine(AP[ring[(lo + n - 1) % n]], AP[ring[lo]],
                                      AP[ring[(lo + 1) % n]]);
}
inline std::optional<std::vector<std::array<std::size_t, 3>>> earClipOriented(
    const std::vector<AffinePoint2>& AP, std::vector<std::size_t> ring,
    Sign want) {
    if (ring.size() < 3) return std::nullopt;
    const Sign wind = winding(AP, ring);
    if (wind == Sign::Zero) return std::nullopt;
    std::vector<std::array<std::size_t, 3>> out;
    while (ring.size() > 3) {
        std::optional<std::size_t> best;
        for (std::size_t i = 0; i < ring.size(); ++i) {
            if (!earValid(AP, ring, i, wind)) continue;
            if (!best || lexCmp(AP[ring[i]], AP[ring[*best]]) < 0) best = i;
        }
        if (!best) return std::nullopt;
        const std::size_t i = *best, n = ring.size();
        out.push_back({ring[(i + n - 1) % n], ring[i], ring[(i + 1) % n]});
        ring.erase(ring.begin() + static_cast<std::ptrdiff_t>(i));
    }
    if (ring.size() == 3) out.push_back({ring[0], ring[1], ring[2]});
    for (auto& t : out) {
        const Sign s = predicates::orient2dAffine(AP[t[0]], AP[t[1]], AP[t[2]]);
        if (s == Sign::Zero) return std::nullopt;
        if (s != want) std::swap(t[1], t[2]);
    }
    return out;
}

}  // namespace scdetail

// Plain-value GENERAL multi/reuse geometry validator.  This is the same
// validator used by the public route.  Vertex identity is an index into
// `vertices`, so the contract is independently testable without workspace
// handles or staged-id lookup.
enum class GeneralSteinerGeometryIssue {
    None,
    InvalidVertexReference,
    DegenerateTriangle,
    InconsistentOrientation,
    ForeignVertexInside,
    UnintendedEdgeTJunction,
    NonadjacentEdgeCrossing,
    StagedPointOffSegment,
    StagedPointBeyondSegment
};

struct GeneralSteinerGeometryVertex {
    predicates::AffinePoint2 point{};
    bool staged = false;
};

struct GeneralSteinerGeometry {
    std::vector<GeneralSteinerGeometryVertex> vertices;
    std::vector<std::array<std::size_t, 3>> triangles;
};

inline GeneralSteinerGeometryIssue validateGeneralSteinerGeometry(
    const Vec2& uvA, const Vec2& uvB,
    const predicates::AffinePoint2& pA,
    const predicates::AffinePoint2& pB,
    const GeneralSteinerGeometry& g) {
    using predicates::Sign;
    Sign referenceOrientation = Sign::Zero;
    for (const auto& tri : g.triangles) {
        for (const std::size_t v : tri)
            if (v >= g.vertices.size())
                return GeneralSteinerGeometryIssue::InvalidVertexReference;
        const Sign s = predicates::orient2dAffine(
            g.vertices[tri[0]].point, g.vertices[tri[1]].point,
            g.vertices[tri[2]].point);
        if (s == Sign::Zero)
            return GeneralSteinerGeometryIssue::DegenerateTriangle;
        if (referenceOrientation == Sign::Zero)
            referenceOrientation = s;
        else if (s != referenceOrientation)
            return GeneralSteinerGeometryIssue::InconsistentOrientation;
    }

    for (const auto& tri : g.triangles) {
        const auto& a = g.vertices[tri[0]].point;
        const auto& b = g.vertices[tri[1]].point;
        const auto& c = g.vertices[tri[2]].point;
        const Sign w = predicates::orient2dAffine(a, b, c);
        for (std::size_t vi = 0; vi < g.vertices.size(); ++vi) {
            if (vi == tri[0] || vi == tri[1] || vi == tri[2]) continue;
            const auto& q = g.vertices[vi].point;
            const Sign s0 = predicates::orient2dAffine(a, b, q);
            const Sign s1 = predicates::orient2dAffine(b, c, q);
            const Sign s2 = predicates::orient2dAffine(c, a, q);
            if (s0 == w && s1 == w && s2 == w)
                return GeneralSteinerGeometryIssue::ForeignVertexInside;
            if (scdetail::onClosedSeg(a, b, q) ||
                scdetail::onClosedSeg(b, c, q) ||
                scdetail::onClosedSeg(c, a, q))
                return GeneralSteinerGeometryIssue::UnintendedEdgeTJunction;
        }
    }

    for (std::size_t i = 0; i < g.triangles.size(); ++i)
        for (std::size_t j = i + 1; j < g.triangles.size(); ++j)
            for (std::size_t a = 0; a < 3; ++a)
                for (std::size_t b = 0; b < 3; ++b) {
                    const std::size_t i0 = g.triangles[i][a];
                    const std::size_t i1 = g.triangles[i][(a + 1) % 3];
                    const std::size_t j0 = g.triangles[j][b];
                    const std::size_t j1 = g.triangles[j][(b + 1) % 3];
                    if (i0 == j0 || i0 == j1 || i1 == j0 || i1 == j1)
                        continue;
                    if (scdetail::properCross(
                            g.vertices[i0].point, g.vertices[i1].point,
                            g.vertices[j0].point, g.vertices[j1].point))
                        return GeneralSteinerGeometryIssue::
                            NonadjacentEdgeCrossing;
                }

    for (const auto& vertex : g.vertices) {
        if (!vertex.staged) continue;
        if (predicates::orient2dAffine(pA, pB, vertex.point) != Sign::Zero)
            return GeneralSteinerGeometryIssue::StagedPointOffSegment;
        if (predicates::orderOnSegmentAffine(uvA, uvB, pA, vertex.point) !=
                Sign::Negative ||
            predicates::orderOnSegmentAffine(uvA, uvB, vertex.point, pB) !=
                Sign::Negative)
            return GeneralSteinerGeometryIssue::StagedPointBeyondSegment;
    }
    return GeneralSteinerGeometryIssue::None;
}

#if defined(PAACS_TOPO_TESTING)
namespace scdetail {
inline thread_local const GeneralSteinerGeometry*
    generalSteinerGeometryTestOverride = nullptr;
}
// Gated seam used only to substitute plain geometry immediately before the
// production validator call.  The public route still maps the real validator
// result through Failed(InvalidSteinerGeometry).
class GeneralSteinerGeometryTestScope {
public:
    explicit GeneralSteinerGeometryTestScope(
        const GeneralSteinerGeometry& replacement) noexcept
        : previous_(scdetail::generalSteinerGeometryTestOverride) {
        scdetail::generalSteinerGeometryTestOverride = &replacement;
    }
    ~GeneralSteinerGeometryTestScope() {
        scdetail::generalSteinerGeometryTestOverride = previous_;
    }
    GeneralSteinerGeometryTestScope(const GeneralSteinerGeometryTestScope&) =
        delete;
    GeneralSteinerGeometryTestScope& operator=(
        const GeneralSteinerGeometryTestScope&) = delete;

private:
    const GeneralSteinerGeometry* previous_ = nullptr;
};
#endif

// ---------------------------------------------------------------------------
// §1.1 - candidate value types (parametric). Corners are Existing (a current
// NodeHandle) or Staged (a SteinerNodeId). Lowering resolves Staged ids to
// topo localIds in deterministic (ascending SteinerNodeId) order.
// ---------------------------------------------------------------------------
using SteinerNodeId = std::uint64_t;   // structural, allocator-history free

struct SteinerNode {
    SteinerNodeId id = 0;
    AuthoritativeParam param;
    ChartId chart = 0;
    Vec3 physical{};                    // derived via provider.pointAt (§4)
    mesh::Classification cls{};         // dim-1 segment classification (§6.2)
};

struct TriNodeRef {
    bool staged = false;
    topo::NodeHandle existing{};        // when !staged
    SteinerNodeId stagedId = 0;         // when staged
    static TriNodeRef Existing(const topo::NodeHandle& h) {
        return TriNodeRef{false, h, 0};
    }
    static TriNodeRef Staged(SteinerNodeId id) {
        return TriNodeRef{true, {}, id};
    }
};
struct TriNodeRefTriangle {
    std::array<TriNodeRef, 3> nodes{};
    mesh::Classification cls{};
};
struct SteinerEdgeReq {
    TriNodeRef a, b;
    mesh::Classification cls{};
};

struct SteinerRecoveryCandidate {
    OwnershipStamp stamp;                       // workspace + session identity
    GeometryProviderStamp providerStamp;        // §4 geometry identity
    ConstraintEntityId entity = 0;
    ChartId chart = 0;
    // Ordered start..end chain of REFERENCES (§1.1): a staged node cannot be
    // a NodeHandle before commit, so the chain is TriNodeRef. The all-current
    // NodeHandle chain is produced only AFTER ASCR returns its
    // localId->NodeHandle mapping.
    std::vector<TriNodeRef> resultChain;
    std::vector<SteinerNode> staged;            // interior additions
    std::vector<topo::CellHandle> removed;
    std::vector<TriNodeRefTriangle> added;
    std::vector<SteinerEdgeReq> requirements;
    mesh::Classification cellClass{};           // owning FACE class (§6.2)
};

// ---------------------------------------------------------------------------
// §9 outcome-reason taxonomy for the Steiner route. Reuses the five public
// Outcomes; these are the Steiner-specific reasons (a superset of BND-1's,
// kept as a separate enum so BND-1 stays byte-identical).
// ---------------------------------------------------------------------------
enum class SteinerReject {
    None,
    StaleGeometryProvider,       // §4: identity/revision mismatch
    NonfiniteMappedPosition,     // §4: provider returned non-finite Vec3
    ProviderDomainRejected,      // §4: provider signalled out-of-domain input
    NonInjectiveStagedId,        // §1.2
    UnknownStagedRef,            // §1.2: a ref names an undeclared staged id
    StaleCandidateStamp,         // candidate stamp != current workspace stamp
    ChartMismatch,
    InvalidSplitPolicy,
    InvalidClassification,       // §6.2
    PslgRejected,                // whole-set preflight
};
enum class SteinerRefuse {
    None,
    SteinerDisabled,
    SplitBudgetExceeded,         // §5.5
    MinimumFeatureReached,       // §5.5
    UnsupportedChart,
    UnsupportedSeam,
    CrossedRequiredEdge,
    PartialChain,                // BND-1 13.3-0 case 3: partial existing chain
    StripNotRecoverable,         // the Steiner strip walk hit a degeneracy
                                 // outside the qualified scope (tangential
                                 // touch, walk guard, boundary reached, no
                                 // unique initial cell, collinear overlap) -
                                 // NOT a single-cavity limitation (multi-cell,
                                 // multi-Steiner and existing-node reuse are
                                 // all supported).
};
enum class SteinerFail {
    None,
    InvalidSteinerGeometry,      // six-gate reject                  (tier 1)
    CavityRejected,              // ASCR combinatorial reject        (tier 1)
    CandidateAllocFailed,        // off-workspace alloc / provider bad_alloc (t1)
    GeometryProviderError,       // provider other-throw             (tier 1)
    StagedNodeAllocFailed,       // node capacity reservation fail   (tier 2)
    CavityThrew,                 // cell/edge/incidence reservation fail (t2)
    IndexRebuildFailed,          // post-commit prefix attribution to the
                                 // NEXT entity (committed prefix preserved)
};

// ---------------------------------------------------------------------------
// Symbolic evaluation helpers - the ONLY place the exact predicate kernel is
// touched. A resolver maps ConstraintVertexId -> stored Vec2 (its param).
// ---------------------------------------------------------------------------
template <class Resolver>
inline predicates::AffinePoint2 toAffinePoint2(const AuthoritativeParam& ap,
                                               Resolver&& resolve) {
    if (const auto* s = std::get_if<StoredParam>(&ap))
        return predicates::storedPoint(s->uv);
    const auto& a = std::get<AffineParam>(ap);
    validateAffineParam(a);   // BEFORE the uint32->uint8 narrowing below, so a
                              // level like 257 can never alias level 1
    return predicates::AffinePoint2{resolve(a.lo), resolve(a.hi), a.num,
                                    static_cast<std::uint8_t>(a.level)};
}

// §5.4 deterministic mapping Vec2 (for pointAt ONLY, never a predicate).
// VALIDATES the token and rejects non-finite Stored values, resolved
// endpoints, and mapping results BEFORE the shift / before the caller hands
// the value to the provider - the same discipline as the predicate path.
template <class Resolver>
inline Vec2 toMappingVec2(const AuthoritativeParam& ap, Resolver&& resolve) {
    if (const auto* s = std::get_if<StoredParam>(&ap)) {
        requireFiniteVec2(s->uv, "stored mapping param");
        return s->uv;
    }
    const auto& a = std::get<AffineParam>(ap);
    validateAffineParam(a);                       // level + canonical + reduced
    const Vec2 lo = resolve(a.lo), hi = resolve(a.hi);
    requireFiniteVec2(lo, "resolved lo endpoint");
    requireFiniteVec2(hi, "resolved hi endpoint");
    const double t = double(a.num) / double(std::uint64_t(1) << a.level);
    const Vec2 m{lo.u + t * (hi.u - lo.u), lo.v + t * (hi.v - lo.v)};
    requireFiniteVec2(m, "mapped position");
    return m;
}

// Exact ordering of two authoritative params along a directed segment
// loSeg -> hiSeg (Sign of t_p - t_q on the dominant axis). Wraps the accepted
// predicates::orderOnSegmentAffine over the resolved corners.
template <class Resolver>
inline predicates::Sign orderOnSegment(const Vec2& loSeg, const Vec2& hiSeg,
                                       const AuthoritativeParam& p,
                                       const AuthoritativeParam& q,
                                       Resolver&& resolve) {
    return predicates::orderOnSegmentAffine(
        loSeg, hiSeg, toAffinePoint2(p, resolve), toAffinePoint2(q, resolve));
}

// Exact collinearity of an authoritative param with the segment supporting
// line through (loSeg, hiSeg): Sign::Zero iff the point lies on the line.
template <class Resolver>
inline predicates::Sign collinearWithSegment(const Vec2& loSeg,
                                             const Vec2& hiSeg,
                                             const AuthoritativeParam& p,
                                             Resolver&& resolve) {
    return predicates::orient2dAffine(predicates::storedPoint(loSeg),
                                      predicates::storedPoint(hiSeg),
                                      toAffinePoint2(p, resolve));
}

// ---------------------------------------------------------------------------
// §1.2 lowering: a VALIDATED candidate -> exactly one topo::SteinerCavitySpec.
// Staged authoritative params are mapped to their derived physical Vec3
// (already computed in SteinerNode::physical by the route via the provider);
// localIds are assigned from the structural SteinerNodeId in ASCENDING order
// (deterministic, allocator-history independent). Returns nullopt on a
// non-injective staged id (a lowering precondition the route also checks).
// This function performs NO geometry-kernel call and NO workspace mutation.
// ---------------------------------------------------------------------------
struct LoweredSpec {
    topo::SteinerCavitySpec spec;
    // stagedId -> localId, ascending; the route maps the ASCR result's
    // localId->NodeHandle back to SteinerNodeId through this.
    std::vector<std::pair<SteinerNodeId, std::uint32_t>> idToLocal;
    // The ordered chain PRE-LOWERED to Existing|Staged(local) refs, so that
    // after the ASCR commit supplies localId->NodeHandle the conversion to a
    // current-handle chain is bounded and non-throwing (no late resolution of
    // an unknown staged id could fail post-commit).
    std::vector<topo::SteinerNodeRef> loweredChain;
};

inline std::optional<LoweredSpec> lowerCandidate(
    const SteinerRecoveryCandidate& c, const topo::WorkspaceStamp& wsStamp) {
    // (0) stamp gate (blocker 4): a caller-supplied fresh workspace stamp must
    //     NOT launder a candidate built against an older workspace/revision.
    //     Every pool identity AND the revision must match.
    if (!(c.stamp.nodePool == wsStamp.nodePool &&
          c.stamp.cellPool == wsStamp.cellPool &&
          c.stamp.edgePool == wsStamp.edgePool &&
          c.stamp.facePool == wsStamp.facePool &&
          c.stamp.revision == wsStamp.revision))
        return std::nullopt;

    // (0b) local-id width guard (blocker 4): localIds are uint32; refuse a
    //      batch that could not be indexed without an overflowing cast.
    if (c.staged.size() >
        static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
        return std::nullopt;

    LoweredSpec out;
    out.spec.stamp = wsStamp;
    out.spec.removed = c.removed;

    // deterministic localId assignment: ascending SteinerNodeId
    std::vector<SteinerNodeId> ids;
    ids.reserve(c.staged.size());
    for (const auto& s : c.staged) ids.push_back(s.id);
    std::sort(ids.begin(), ids.end());
    for (std::size_t i = 1; i < ids.size(); ++i)
        if (ids[i - 1] == ids[i]) return std::nullopt;   // non-injective
    // EQUALITY-checked resolution (blocker 3): an undeclared id must NOT
    // silently alias a neighbour via lower_bound's insertion point.
    auto localOf = [&](SteinerNodeId id) -> std::optional<std::uint32_t> {
        const auto it = std::lower_bound(ids.begin(), ids.end(), id);
        if (it == ids.end() || *it != id) return std::nullopt;
        return static_cast<std::uint32_t>(it - ids.begin());
    };
    for (SteinerNodeId id : ids)
        out.idToLocal.emplace_back(id, *localOf(id));

    // staged NodeData (Vec3 physical + classification), in localId order
    out.spec.steinerNodes.resize(c.staged.size());
    for (const auto& s : c.staged) {
        const auto l = localOf(s.id);
        if (!l) return std::nullopt;
        topo::SteinerNodeSpec ns;
        ns.localId = *l;
        ns.coord = s.physical;
        ns.cls = s.cls;
        out.spec.steinerNodes[*l] = ns;
    }
    auto lowerRef =
        [&](const TriNodeRef& r) -> std::optional<topo::SteinerNodeRef> {
        if (!r.staged) return topo::SteinerNodeRef::Existing(r.existing);
        const auto l = localOf(r.stagedId);
        if (!l) return std::nullopt;                     // unknown staged id
        return topo::SteinerNodeRef::Staged(*l);
    };
    out.spec.added.reserve(c.added.size());
    for (const auto& t : c.added) {
        topo::SteinerCavityCell cell;
        for (int i = 0; i < 3; ++i) {
            const auto lr = lowerRef(t.nodes[std::size_t(i)]);
            if (!lr) return std::nullopt;
            cell.nodes[std::size_t(i)] = *lr;
        }
        cell.cls = t.cls;
        out.spec.added.push_back(cell);
    }
    out.spec.requirements.reserve(c.requirements.size());
    for (const auto& r : c.requirements) {
        const auto la = lowerRef(r.a);
        const auto lb = lowerRef(r.b);
        if (!la || !lb) return std::nullopt;
        topo::SteinerEdgeRequirement req;
        req.a = *la;
        req.b = *lb;
        req.cls = r.cls;
        out.spec.requirements.push_back(req);
    }
    // Pre-lower the ordered chain (blocker 3): an unknown staged id used ONLY
    // in resultChain must fail lowering NOW, never post-commit.
    out.loweredChain.reserve(c.resultChain.size());
    for (const auto& r : c.resultChain) {
        const auto lr = lowerRef(r);
        if (!lr) return std::nullopt;
        out.loweredChain.push_back(*lr);
    }
    return out;
}

// ---------------------------------------------------------------------------
// Increment 3 - single-Steiner recovery route (BND layer). Bounded scope: one
// interior Steiner node in a single QUAD cavity (two triangles whose shared
// diagonal the segment crosses) - the canonical single-Steiner scenario
// (design §10 "First implementation may support one Steiner point per
// constraint"). All six geometric gates run in BND over the exact symbolic
// model BEFORE ASCR preparation; provider calls + exception classification are
// entirely in BND; the metadata shadow / chain capacity / result capacity are
// reserved BEFORE commit; after ASCR begins committing only bounded,
// non-throwing handle substitution + publication happen (no map insertion,
// vector growth, provider call, predicate, or allocation post-commit).
// ---------------------------------------------------------------------------
// Blocker 5: audit the persistent authoritative-parameter store. Every record
// must reference a current, owned node with a matching ordinary param entry and
// chart; stored values finite; affine tokens canonical/reduced with existing
// endpoint identities on the same chart; and the authoritative token must AGREE
// with the stored mapping coordinate. A single bad record fails the whole set.
inline bool auditAuthParams(const topo::TopologyWorkspace& ws,
                            const ConstraintSet& cs) {
    for (const auto& [h, rec] : cs.authParams) {
        if (h.pool != ws.nodes().poolId() || !ws.nodes().isCurrent(h))
            return false;
        const ParamPoint* pp = cs.paramOf(h);
        if (!pp) return false;                 // must have an ordinary entry
        if (pp->chart != rec.chart) return false;
        if (!(std::isfinite(pp->p.u) && std::isfinite(pp->p.v))) return false;
        if (!rec.affine) {
            if (!(std::isfinite(rec.storedUv.u) &&
                  std::isfinite(rec.storedUv.v)))
                return false;
            if (!(rec.storedUv.u == pp->p.u && rec.storedUv.v == pp->p.v))
                return false;
            continue;
        }
        // canonical reduced affine token
        if (!(rec.lo < rec.hi)) return false;
        if (rec.level < 1 || rec.level > kMaxSteinerLevel) return false;
        if (rec.num == 0 || rec.num >= (std::uint64_t(1) << rec.level))
            return false;
        if ((rec.num & 1u) == 0) return false;
        const auto ila = cs.vertices.find(rec.lo);
        const auto ilb = cs.vertices.find(rec.hi);
        if (ila == cs.vertices.end() || ilb == cs.vertices.end()) return false;
        const ParamPoint* pl = cs.paramOf(ila->second);
        const ParamPoint* ph = cs.paramOf(ilb->second);
        if (!pl || !ph) return false;
        if (pl->chart != rec.chart || ph->chart != rec.chart) return false;
        const double t = double(rec.num) / double(std::uint64_t(1) << rec.level);
        const Vec2 m{pl->p.u + t * (ph->p.u - pl->p.u),
                     pl->p.v + t * (ph->p.v - pl->p.v)};
        if (!(m.u == pp->p.u && m.v == pp->p.v)) return false;   // agreement
    }
    return true;
}

// Blocker 6: a standalone PLAIN-VALUE six-gate preview validator, independent
// of the workspace/topology. It takes resolved parameter coordinates for the
// quad corners A,P,B,Q, the Steiner point as an exact AffinePoint2, the M
// triangles as exact AffinePoint2 triples, and the chain shape. One
// discriminating adversary per gate exists in tests (bnd2.steiner_gate_*).
enum class SteinerGate {
    Ok,
    Degenerate,        // gate 1: an M triangle has zero area
    Orientation,       // gate 2: M triangles disagree in orientation
    Interiority,       // gate 3: S not exactly on A-B strictly between
    BoundaryNonSimple, // gate 4: quad A-P-B-Q diagonals do not properly cross
    OverlapTJunction,  // gate 5: P,Q not strictly on opposite sides of A-B
    ChainNonConform,   // gate 6: chain is not exactly Existing,Staged,Existing
};
inline SteinerGate steinerSixGate(
    const Vec2& uvA, const Vec2& uvB, const Vec2& pP, const Vec2& pQ,
    const predicates::AffinePoint2& sAP,
    const std::vector<std::array<predicates::AffinePoint2, 3>>& mtris,
    std::size_t chainSize, bool chainFrontStaged, bool chainMidStaged,
    bool chainBackStaged) {
    using Sign = predicates::Sign;
    // gate 1 nondegeneracy + gate 2 consistent orientation
    Sign ref0 = Sign::Zero;
    for (const auto& t : mtris) {
        const Sign s = predicates::orient2dAffine(t[0], t[1], t[2]);
        if (s == Sign::Zero) return SteinerGate::Degenerate;
        if (ref0 == Sign::Zero) ref0 = s;
        else if (s != ref0) return SteinerGate::Orientation;
    }
    // gate 3 interiority: S exactly collinear with A-B and strictly between
    {
        const auto pa = predicates::storedPoint(uvA);
        const auto pb = predicates::storedPoint(uvB);
        if (predicates::orient2dAffine(pa, pb, sAP) != Sign::Zero)
            return SteinerGate::Interiority;
        if (predicates::orderOnSegmentAffine(uvA, uvB, pa, sAP) !=
                Sign::Negative ||
            predicates::orderOnSegmentAffine(uvA, uvB, sAP, pb) !=
                Sign::Negative)
            return SteinerGate::Interiority;
    }
    // gate 5 (T-junction / overlap) is the PREREQUISITE and is checked FIRST so
    // it is independently reachable: P and Q must be strictly on OPPOSITE sides
    // of the line A-B. Same side, or a boundary node ON A-B, is an overlap /
    // unintended T-junction.
    {
        using predicates::orient2d;
        const Sign sp = orient2d(uvA, uvB, pP), sq = orient2d(uvA, uvB, pQ);
        if (sp == Sign::Zero || sq == Sign::Zero || sp == sq)
            return SteinerGate::OverlapTJunction;
    }
    // gate 4 (boundary simplicity): given P,Q straddle the line A-B, the finite
    // diagonals A-B and P-Q must actually cross - i.e. A and B must also
    // straddle the line P-Q. Opposite sides of A-B but no finite crossing is a
    // non-simple boundary.
    {
        using predicates::orient2d;
        const Sign e1 = orient2d(pP, pQ, uvA), e2 = orient2d(pP, pQ, uvB);
        if (e1 == Sign::Zero || e2 == Sign::Zero || e1 == e2)
            return SteinerGate::BoundaryNonSimple;
    }
    // gate 6 chain conformity: exactly Existing, Staged, Existing
    if (chainSize != 3 || chainFrontStaged || chainBackStaged ||
        !chainMidStaged)
        return SteinerGate::ChainNonConform;
    return SteinerGate::Ok;
}

struct SteinerPolicy {
    bool steinerEnabled = true;
    std::uint32_t splitBudget = 1;                 // max staged nodes / segment
    std::uint32_t minFeatureLevel = kMaxSteinerLevel;
};
enum class SteinerStatus { AlreadyPresent, Recovered, Refused, Failed, Rejected };
struct SteinerResult {
    SteinerStatus status = SteinerStatus::Failed;
    SteinerReject reject = SteinerReject::None;
    SteinerRefuse refuse = SteinerRefuse::None;
    SteinerFail fail = SteinerFail::None;
    std::vector<topo::NodeHandle> recoveredChain;   // all-current, post-commit
    // (handle, authoritative param) for each staged node the route CREATED and
    // ALREADY published into the ConstraintSet (cs.params + cs.authParams via the
    // atomic copy/swap on commit). This vector is a report of what was created,
    // NOT a to-do for the caller - the metadata is already persisted.
    std::vector<std::pair<topo::NodeHandle, AuthoritativeParam>> createdParams;
    std::uint64_t revision = 0;
};

// ---------------------------------------------------------------------------
// Increment 4 - the STEINER strip walk (blocker 2, existing-node reuse). It is
// an adaptation of BND-1's FROZEN transverse walk (SegmentRecoveryCommit.hpp
// 13.8) that, instead of REFUSING at an interior on-segment mesh vertex,
// records that vertex as an existing-node reuse WAYPOINT and continues from it.
// EVERY other degeneracy and the crossed-required-edge policy are preserved
// verbatim (CrossedRequiredEdge, TangentialTouch, WalkGuard, BoundaryReached,
// NoUniqueInitialCell, CollinearOverlapWithEdge). The result is the ordered
// waypoint chain A, X1.., B and, per LEG (between consecutive waypoints), the
// crossed strip, both boundary chains, and the crossed-edge demand. The reused
// node needs NO new identity: it is a current NodeHandle with a StoredParam
// (parameter totality already satisfied by the A1 preflight) and is referenced
// as TriNodeRef::Existing; staged nodes stay anchored to the ORIGINAL segment
// endpoints.
// ---------------------------------------------------------------------------
namespace swalk {
struct StripLeg {
    std::vector<topo::CellHandle> strip;
    std::vector<topo::NodeHandle> left, right;   // wi .. wj inclusive, per side
    std::size_t crossed = 0;                      // crossed interior edges
};
struct StripWalkResult {
    bool ok = false;
    RefusedReason why = RefusedReason::None;
    std::vector<topo::NodeHandle> waypoints;      // A, reused interiors..., B
    std::vector<StripLeg> legs;
};

inline StripWalkResult steinerStripWalk(const topo::TopologyWorkspace& ws,
                                        const ConstraintSet& cs,
                                        const TriAdjacencyIndex& idx,
                                        const topo::NodeHandle& A,
                                        const topo::NodeHandle& B) {
    StripWalkResult R;
    R.waypoints.push_back(A);
    const Vec2& pAseg = detail::pOf(cs, A);
    const Vec2& pBseg = detail::pOf(cs, B);
    const std::size_t guard = ws.cells().activeBySequence().size() + 1;
    topo::NodeHandle legStart = A;
    std::size_t legGuard = 0;
    while (true) {
        if (++legGuard > guard + 1) { R.why = RefusedReason::WalkGuard; return R; }
        topo::CellHandle cur{};
        if (!detail::initialCell(ws, cs, idx, legStart, B, cur, R.why)) return R;
        StripLeg leg;
        leg.left.push_back(legStart);
        leg.right.push_back(legStart);
        topo::NodeHandle fromU{}, fromV{};
        bool haveEdge = false, reachedB = false, hitInterior = false;
        topo::NodeHandle stopX{};
        while (leg.strip.size() <= guard) {
            leg.strip.push_back(cur);
            const auto cd = ws.cells().get(cur);
            bool hasB = false;
            for (int i = 0; i < 3; ++i)
                if (cd.nodes[std::size_t(i)] == B) hasB = true;
            if (hasB) { reachedB = true; break; }
            int crossedLocal = -1;
            topo::NodeHandle cu{}, cv{};
            int crossings = 0;
            for (int i = 0; i < 3 && !hitInterior; ++i) {
                const auto& x = cd.nodes[std::size_t(i)];
                const auto& y = cd.nodes[std::size_t((i + 1) % 3)];
                if (haveEdge && ((x == fromU && y == fromV) ||
                                 (x == fromV && y == fromU)))
                    continue;                            // the entry edge
                const Vec2& px = detail::pOf(cs, x);
                const Vec2& py = detail::pOf(cs, y);
                const auto d1 = predicates::orient2d(pAseg, pBseg, px);
                const auto d2 = predicates::orient2d(pAseg, pBseg, py);
                const auto d3 = predicates::orient2d(px, py, pAseg);
                const auto d4 = predicates::orient2d(px, py, pBseg);
                if (d1 == predicates::Sign::Zero &&
                    d2 == predicates::Sign::Zero) {
                    // BOTH endpoints of this edge are on the segment line: the
                    // edge overlaps the segment. This is a collinear overlap
                    // (BND-1 CollinearOverlapWithEdge), NOT node reuse - handle
                    // it BEFORE the single-endpoint waypoint cases.
                    R.why = RefusedReason::CollinearOverlapWithEdge;
                    return R;
                }
                if (d1 == predicates::Sign::Zero ||
                    d2 == predicates::Sign::Zero) {
                    // exactly ONE endpoint lies on the segment line: a
                    // non-endpoint one is an existing-node reuse waypoint (where
                    // BND-1 would refuse with InteriorNodeHit).
                    if (d1 == predicates::Sign::Zero &&
                        !(x == legStart || x == B)) {
                        stopX = x; hitInterior = true;
                    } else if (d2 == predicates::Sign::Zero &&
                               !(y == legStart || y == B)) {
                        stopX = y; hitInterior = true;
                    }
                    continue;
                }
                const bool s1 = d1 != d2;
                const bool s2 = d3 != d4 && d3 != predicates::Sign::Zero &&
                                d4 != predicates::Sign::Zero;
                if (s1 && s2) { crossedLocal = i; cu = x; cv = y; ++crossings; }
            }
            if (hitInterior) break;                      // terminate leg at X
            if (crossings == 0) { R.why = RefusedReason::TangentialTouch; return R; }
            if (crossings > 1) { R.why = RefusedReason::WalkGuard; return R; }
            {
                const auto sideU =
                    predicates::orient2d(pAseg, pBseg, detail::pOf(cs, cu));
                if (sideU == predicates::Sign::Positive) {
                    if (!(leg.left.back() == cu)) leg.left.push_back(cu);
                    if (!(leg.right.back() == cv)) leg.right.push_back(cv);
                } else {
                    if (!(leg.right.back() == cu)) leg.right.push_back(cu);
                    if (!(leg.left.back() == cv)) leg.left.push_back(cv);
                }
            }
            {
                const auto eh = ws.edges().find(ws.nodes(), cu, cv);
                if (eh) {
                    if (ws.edges().get(*eh).requiredConstraint) {
                        R.why = RefusedReason::CrossedRequiredEdge;
                        return R;
                    }
                    ++leg.crossed;
                }
            }
            const auto nb = idx.neighbors.find(cur);
            if (nb == idx.neighbors.end()) { R.why = RefusedReason::BoundaryReached; return R; }
            const auto next = nb->second[std::size_t(crossedLocal)];
            if (!ws.cells().isCurrent(next)) { R.why = RefusedReason::BoundaryReached; return R; }
            fromU = cu; fromV = cv; haveEdge = true; cur = next;
        }
        if (!reachedB && !hitInterior) { R.why = RefusedReason::WalkGuard; return R; }
        const topo::NodeHandle stop = reachedB ? B : stopX;
        // exact safety: a reuse waypoint must be strictly AFTER the current leg
        // start (not merely somewhere between the global A and B) AND strictly
        // before B, so the walk makes monotone forward progress every leg.
        if (!reachedB) {
            const Vec2& pLeg = detail::pOf(cs, legStart);
            const auto sp = predicates::storedPoint(detail::pOf(cs, stop));
            if (predicates::orderOnSegmentAffine(
                    pAseg, pBseg, predicates::storedPoint(pLeg), sp) !=
                    predicates::Sign::Negative ||
                predicates::orderOnSegmentAffine(
                    pAseg, pBseg, sp, predicates::storedPoint(pBseg)) !=
                    predicates::Sign::Negative) {
                R.why = RefusedReason::TangentialTouch;
                return R;
            }
        }
        leg.left.push_back(stop);
        leg.right.push_back(stop);
        R.legs.push_back(std::move(leg));
        R.waypoints.push_back(stop);
        if (reachedB) { R.ok = true; return R; }
        legStart = stopX;
    }
}
}  // namespace swalk

// Per-entity recovery using a caller-OWNED, already-validated session. The
// caller (the single-segment opener OR the set wrapper) has already run the
// whole-set preflight (A1 mapping + authParams audit + RecoverySession::open)
// and passes the validated session in; this helper NEVER opens a session, so
// the set wrapper can own exactly one and rebuild it only after a real commit.
inline SteinerResult recoverSteinerEntityValidated(
    topo::TopologyWorkspace& ws, ConstraintSet& cs,
    const ConstraintSegment& seg, const GeometryProvider& provider,
    const GeometryProviderStamp& expectedProvider, const SteinerPolicy& policy,
    const RecoverySession& sess) {
    SteinerResult r;
    auto reject = [&](SteinerReject x) {
        r.status = SteinerStatus::Rejected; r.reject = x; return r;
    };
    auto refuse = [&](SteinerRefuse x) {
        r.status = SteinerStatus::Refused; r.refuse = x; return r;
    };
    auto fail = [&](SteinerFail x) {
        r.status = SteinerStatus::Failed; r.fail = x; return r;
    };

    // (3) the supplied segment must belong to the validated set (exact,
    //     including classification).
    {
        bool inSet = false;
        for (const auto& s : cs.segments)
            if (s.id == seg.id && s.a == seg.a && s.b == seg.b &&
                s.chart == seg.chart && s.cls.dim == seg.cls.dim &&
                s.cls.tag == seg.cls.tag)
                inSet = true;
        if (!inSet) return reject(SteinerReject::PslgRejected);
    }
    if (seg.cls.dim != 1 || seg.cls.tag < 0)
        return reject(SteinerReject::InvalidClassification);

    // (2b) split-policy validity is checked BEFORE any provider call / mutation
    //      so an out-of-range policy is a structured Rejected(InvalidSplitPolicy)
    //      with the provider PROVABLY untouched (never an escaped
    //      std::invalid_argument). splitBudget has no invalid value (0 simply
    //      refuses a create route); only the dyadic level has a valid domain.
    if (policy.minFeatureLevel < 1 || policy.minFeatureLevel > kMaxSteinerLevel)
        return reject(SteinerReject::InvalidSplitPolicy);

    // endpoints + params.
    const topo::NodeHandle A = cs.vertices.at(seg.a);
    const topo::NodeHandle B = cs.vertices.at(seg.b);
    const ParamPoint* ppA = cs.paramOf(A);
    const ParamPoint* ppB = cs.paramOf(B);
    if (!ppA || !ppB) return reject(SteinerReject::PslgRejected);
    if (ppA->chart != seg.chart || ppB->chart != seg.chart)
        return reject(SteinerReject::ChartMismatch);
    const Vec2 uvA = ppA->p, uvB = ppB->p;

    // NOTE (blocker 3): the geometry-provider identity gate (§4) is DEFERRED to
    // the Steiner create path, AFTER demand/budget/min-feature, so that an
    // over-budget or minimum-feature refusal touches NO provider call (not even
    // stamp()). The mark-only / already-present routes create no nodes and do
    // not consult the provider at all.

    // (5) route the entity through the NORMATIVE routing (13.3-0), reusing
    //     BND-1: AlreadyPresent recognizes the recovered A-S-B chain;
    //     MarkOnly/PartialChain reuse BND-1 mark-only recovery with EXACT reason
    //     preservation; Transverse enters the Steiner create path.
    std::vector<topo::NodeHandle> chain;
    const Routing route =
        routeEntity(ws, sess.input(), sess.index(), seg, chain);
    if (route == Routing::AlreadyPresent) {
        r.status = SteinerStatus::AlreadyPresent;
        r.recoveredChain = chain;
        r.revision = ws.topologyRevision();
        return r;
    }
    if (route == Routing::MarkOnly || route == Routing::PartialChain) {
        ConstraintSet one;
        one.vertices = cs.vertices;
        one.params = cs.params;
        one.authParams = cs.authParams;
        one.stamp = cs.stamp;
        one.segments = {seg};
        const RecoveryReport rep = recoverSegments(ws, one);
        if (rep.rejected || rep.processed.empty())
            return reject(SteinerReject::PslgRejected);
        const EntityResult er = rep.processed.front();
        switch (er.outcome) {
            case Outcome::AlreadyPresent:
                r.status = SteinerStatus::AlreadyPresent;
                r.recoveredChain = chain; break;
            case Outcome::Recovered:
                r.status = SteinerStatus::Recovered;
                r.recoveredChain = chain; break;
            case Outcome::Refused:
                // Preserve BND-1 refusal semantics distinctly (blocker 5).
                switch (er.refused) {
                    case RefusedReason::ConflictingEdgeMark:
                    case RefusedReason::MixedCellClassification:
                        // §6.2: a conflicting/illegal mark is a CLASSIFICATION
                        // reject, never a crossed-edge refusal.
                        r.status = SteinerStatus::Rejected;
                        r.reject = SteinerReject::InvalidClassification;
                        break;
                    case RefusedReason::PartialChain:
                        r.status = SteinerStatus::Refused;
                        r.refuse = SteinerRefuse::PartialChain;
                        break;
                    case RefusedReason::CrossedRequiredEdge:
                        r.status = SteinerStatus::Refused;
                        r.refuse = SteinerRefuse::CrossedRequiredEdge;
                        break;
                    default:
                        r.status = SteinerStatus::Refused;
                        r.refuse = SteinerRefuse::StripNotRecoverable;
                        break;
                }
                break;
            case Outcome::Failed:
                // Preserve BND-1 failure tiers distinctly (blocker 5).
                r.status = SteinerStatus::Failed;
                switch (er.failed) {
                    case FailedReason::CandidateAllocFailed:
                        r.fail = SteinerFail::CandidateAllocFailed; break;
                    case FailedReason::CavityThrew:
                        r.fail = SteinerFail::CavityThrew; break;
                    case FailedReason::IndexRebuildFailed:
                        r.fail = SteinerFail::IndexRebuildFailed; break;
                    default:   // CavityRejected / InvalidCandidateGeometry / etc.
                        r.fail = SteinerFail::CavityRejected; break;
                }
                break;
            case Outcome::Rejected:
            default:
                r.status = SteinerStatus::Rejected;
                r.reject = SteinerReject::InvalidClassification; break;
        }
        if (r.status == SteinerStatus::Recovered ||
            r.status == SteinerStatus::AlreadyPresent) {
            r.revision = ws.topologyRevision();
            cs.stamp = OwnershipStamp::of(ws);
        }
        return r;
    }
    // route == Routing::Transverse: the general multi-Steiner strip CREATE
    // path. The Steiner strip walk splits the segment at every EXISTING on-
    // segment mesh node the walk passes through (existing-node reuse), yielding
    // one LEG per span between consecutive waypoints; each leg is a transverse
    // strip. All legs' cavities + staged nodes commit in ONE atomic ASCR call
    // with ONE revision bump. The single-strip (no reuse) case is the one-leg
    // specialization; K==1 on a single leg reproduces the accepted single-
    // Steiner 4-triangle fan and chain a-S-b exactly.
    //
    // REUSE model: a reused node needs NO new identity - it is a current
    // NodeHandle referenced as TriNodeRef::Existing with its StoredParam (an
    // ANONYMOUS mesh node not declared in cs.vertices is reusable, since A1
    // parameter totality already requires it to have a StoredParam; it need not
    // be a constraint vertex). Staged identity stays anchored to the ORIGINAL
    // segment endpoints (a single global canonical dyadic subdivision on A-B),
    // and a canonical token that COINCIDES with a reused waypoint is served by
    // that node (dropped, no staged creation). Each reused waypoint is made a
    // strictly INTERIOR existing vertex of the cavity (§7b) so the region is one
    // simple disk. The staged/reused minimum-feature is an exact staged/Stored
    // dyadic-boundary test; a reused node's own classification must be
    // compatible with the segment (checked below, §6.2).
    if (!policy.steinerEnabled) return refuse(SteinerRefuse::SteinerDisabled);

    // resolvers (checked; never fabricate a coordinate).
    bool resolveOk = true;
    auto paramOfHandle = [&](const topo::NodeHandle& h) -> Vec2 {
        const ParamPoint* p = cs.paramOf(h);
        if (!p) { resolveOk = false; return Vec2{}; }
        return p->p;
    };
    auto resolve = [&](ConstraintVertexId id) -> Vec2 {
        const auto it = cs.vertices.find(id);
        const ParamPoint* p =
            it != cs.vertices.end() ? cs.paramOf(it->second) : nullptr;
        if (!p) { resolveOk = false; return Vec2{}; }
        return p->p;
    };

    const auto wantOpt = sess.orientationOf(seg.chart);
    if (!wantOpt) return reject(SteinerReject::ChartMismatch);
    const predicates::Sign want = *wantOpt;

    const auto pAsp = predicates::storedPoint(uvA);
    const auto pBsp = predicates::storedPoint(uvB);

    // (7) STEINER STRIP WALK (blocker 2): the ordered waypoint chain A, reused
    //     interior nodes..., B, and per LEG the crossed strip + both boundary
    //     chains + the crossed-edge demand. An existing on-segment mesh node is
    //     REUSED as a chain waypoint (TriNodeRef::Existing + its StoredParam);
    //     it needs no new identity. Every non-reuse degeneracy + the
    //     crossed-required-edge policy are preserved from the frozen BND-1 walk.
    const auto walkRes =
        swalk::steinerStripWalk(ws, sess.input(), sess.index(), A, B);
    if (!walkRes.ok) {
        if (walkRes.why == RefusedReason::CrossedRequiredEdge)
            return refuse(SteinerRefuse::CrossedRequiredEdge);
        return refuse(SteinerRefuse::StripNotRecoverable);
    }
    struct Waypt { topo::NodeHandle h; Vec2 uv; };
    std::vector<Waypt> way;
    way.reserve(walkRes.waypoints.size());
    for (const auto& h : walkRes.waypoints) way.push_back({h, paramOfHandle(h)});
    if (!resolveOk) return reject(SteinerReject::PslgRejected);

    struct Leg {
        std::vector<topo::CellHandle> strip;
        std::vector<topo::NodeHandle> left, right;
        std::uint32_t demand = 0;
        std::vector<AffineParam> params;    // this leg's staged params, A->B
    };
    std::vector<Leg> legs;
    legs.reserve(walkRes.legs.size());
    bool haveClass = false;
    mesh::Classification cellClass{};
    std::size_t Ktotal = 0;
    // PHASE A1 - demand + classification only. NO subdivision allocation and NO
    // provider call until the STAGED budget + min-feature pass. Overflow-safe.
    for (const auto& wl : walkRes.legs) {
        for (const auto& ch : wl.strip) {   // uniform owning-face class (§6.2)
            const auto cc = ws.cells().get(ch).cls;
            if (!haveClass) { cellClass = cc; haveClass = true; }
            else if (cc.dim != cellClass.dim || cc.tag != cellClass.tag)
                return reject(SteinerReject::InvalidClassification);
        }
        const std::size_t Ki = wl.crossed;
        if (Ki == 0)   // a direct-edge leg = partial existing chain (13.3-0 c3)
            return refuse(SteinerRefuse::PartialChain);
        if (Ki > static_cast<std::size_t>(
                     std::numeric_limits<std::uint32_t>::max()))
            return refuse(SteinerRefuse::SplitBudgetExceeded);
        if (Ktotal > static_cast<std::size_t>(
                         std::numeric_limits<std::uint32_t>::max()) - Ki)
            return refuse(SteinerRefuse::SplitBudgetExceeded);
        Ktotal += Ki;
        Leg leg;
        leg.strip = wl.strip;
        leg.left = wl.left;
        leg.right = wl.right;
        leg.demand = static_cast<std::uint32_t>(Ki);
        legs.push_back(std::move(leg));
    }
    if (Ktotal == 0) return refuse(SteinerRefuse::StripNotRecoverable);

    // (blocker 1, §6.2) every INTERIOR reused waypoint must ALREADY carry a
    // EXACT classification of the required segment - it is NOT silently
    // reclassified. No vertex-to-curve ownership relation exists here, so an
    // arbitrary dimension-0 node cannot prove membership in this curve.
    // Endpoints are outside this loop and remain unaffected.
    for (std::size_t wi = 1; wi + 1 < way.size(); ++wi) {
        const auto vc = ws.nodes().get(way[wi].h).cls;
        const bool compatible =
            vc.dim == seg.cls.dim && vc.tag == seg.cls.tag;
        if (!compatible) return reject(SteinerReject::InvalidClassification);
    }

    // exact test: does a canonical dyadic token COINCIDE with a reused waypoint?
    // (same parameter on A-B => same point, since both are on the segment.)
    auto coincidesWaypoint = [&](const AffineParam& t) {
        for (std::size_t wi = 1; wi + 1 < way.size(); ++wi)
            if (predicates::orderOnSegmentAffine(
                    uvA, uvB, toAffinePoint2(AuthoritativeParam{t}, resolve),
                    predicates::storedPoint(way[wi].uv)) ==
                predicates::Sign::Zero)
                return true;
        return false;
    };

    // (blocker 2) STAGED-ONLY budget. The subdivision demand is Ktotal canonical
    // dyadic tokens, but a token that coincides EXACTLY with a reused waypoint is
    // SERVED by that existing node (no staged node is created there). Count the
    // reuse coincidences with a NON-ALLOCATING token pass, derive the staged
    // demand, and enforce the budget BEFORE the subdivision vector, any provider
    // call, or any workspace mutation (budget-before-allocation).
    std::uint32_t stagedDemand = 0;
    try {
        forEachCanonicalDyadicNumerator(
            static_cast<std::uint32_t>(Ktotal),
            [&](std::uint64_t num, std::uint32_t L) {
                if (!coincidesWaypoint(makeAffine(seg.a, seg.b, num, L)))
                    ++stagedDemand;
            });
    } catch (const std::invalid_argument&) {
        return refuse(SteinerRefuse::SplitBudgetExceeded);   // level > kMax
    }
    if (!resolveOk) return reject(SteinerReject::PslgRejected);
    if (stagedDemand > policy.splitBudget)
        return refuse(SteinerRefuse::SplitBudgetExceeded);

    // PHASE A2 - materialize the canonical dyadic subdivision (anchored to the
    // ORIGINAL segment endpoints), DROP tokens coincident with a reused waypoint
    // (served by that node), partition the remaining staged tokens into legs by
    // EXACT ordering, then enforce the exact minimum-feature policy. Still NO
    // provider. Staged identity never depends on reused-node positions.
    std::vector<AffineParam> gp;
    try {
        gp = canonicalDyadicParams(seg.a, seg.b,
                                   static_cast<std::uint32_t>(Ktotal));
    } catch (const std::invalid_argument&) {
        return refuse(SteinerRefuse::SplitBudgetExceeded);   // level > kMax
    }
    std::sort(gp.begin(), gp.end(), [&](const AffineParam& x,
                                        const AffineParam& y) {
        return predicates::orderOnSegmentAffine(
                   uvA, uvB, toAffinePoint2(AuthoritativeParam{x}, resolve),
                   toAffinePoint2(AuthoritativeParam{y}, resolve)) ==
               predicates::Sign::Negative;
    });
    // partition the globally-ordered staged params into their legs (each is
    // STRICTLY between the leg's two waypoints along A->B); a token coincident
    // with a waypoint is DROPPED (the reused node serves that split).
    {
        std::size_t li = 0;
        for (const auto& g : gp) {
            const auto gap = toAffinePoint2(AuthoritativeParam{g}, resolve);
            while (li + 1 < way.size() &&
                   predicates::orderOnSegmentAffine(
                       uvA, uvB, gap,
                       predicates::storedPoint(way[li + 1].uv)) !=
                       predicates::Sign::Negative)
                ++li;                                // advance past waypoint li+1
            if (coincidesWaypoint(g)) continue;      // reused: no staged node
            if (li >= legs.size())
                return fail(SteinerFail::InvalidSteinerGeometry);
            legs[li].params.push_back(g);
        }
        if (!resolveOk) return reject(SteinerReject::PslgRejected);
    }
    // exact minimum feature: staged/endpoint (t=0,1), staged/staged (dyadic),
    // and staged/arbitrary-Stored reused waypoint (exact dyadic-boundary test).
    auto farFromEnds = [&](const AffineParam& p) {
        const std::uint32_t L = std::max(p.level, policy.minFeatureLevel);
        const std::uint64_t num = p.num << (L - p.level);
        const std::uint64_t den = std::uint64_t(1) << L;
        const std::uint64_t thr = std::uint64_t(1) << (L - policy.minFeatureLevel);
        return num >= thr && (den - num) >= thr;
    };
    // is staged param p at least 2^-minLevel (in A-B parameter) from Stored X?
    auto farFromStored = [&](const AffineParam& p, const Vec2& xuv) {
        const std::uint32_t L = std::max(p.level, policy.minFeatureLevel);
        const std::uint64_t base = p.num << (L - p.level);   // t = base/2^L
        const std::uint64_t step = std::uint64_t(1) << (L - policy.minFeatureLevel);
        const std::uint64_t den = std::uint64_t(1) << L;
        const auto Xsp = predicates::storedPoint(xuv);
        bool ok = false;
        if (base > step) {   // lower dyadic boundary is strictly inside (0,1)
            const AffineParam lo = makeAffine(seg.a, seg.b, base - step, L);
            if (predicates::orderOnSegmentAffine(
                    uvA, uvB, Xsp,
                    toAffinePoint2(AuthoritativeParam{lo}, resolve)) !=
                predicates::Sign::Positive)
                ok = true;   // X <= t-2^-m
        }
        if (!ok && base + step < den) {   // upper boundary inside (0,1)
            const AffineParam hi = makeAffine(seg.a, seg.b, base + step, L);
            if (predicates::orderOnSegmentAffine(
                    uvA, uvB, Xsp,
                    toAffinePoint2(AuthoritativeParam{hi}, resolve)) !=
                predicates::Sign::Negative)
                ok = true;   // X >= t+2^-m
        }
        return ok;
    };
    for (std::size_t li = 0; li < legs.size(); ++li) {
        const Leg& leg = legs[li];
        for (std::size_t j = 0; j < leg.params.size(); ++j) {
            if (!farFromEnds(leg.params[j]))         // staged / global endpoint
                return refuse(SteinerRefuse::MinimumFeatureReached);
            if (j + 1 < leg.params.size()) {         // staged / staged
                if (!dyadicSeparated(leg.params[j], leg.params[j + 1],
                                     policy.minFeatureLevel))
                    return refuse(SteinerRefuse::MinimumFeatureReached);
            } else if (li + 1 < way.size() &&
                       !(way[li + 1].h == B)) {       // last staged / reused wj
                if (!farFromStored(leg.params[j], way[li + 1].uv))
                    return refuse(SteinerRefuse::MinimumFeatureReached);
            }
            if (j == 0 && li > 0) {                   // first staged / reused wi
                if (!farFromStored(leg.params[0], way[li].uv))
                    return refuse(SteinerRefuse::MinimumFeatureReached);
            }
        }
    }

    // (§4) geometry-provider identity gate - reached ONLY after budget +
    // min-feature pass, so a refusal above touches no provider call.
    GeometryProviderStamp provStamp;
    try {
        provStamp = provider.stamp();
    } catch (const GeometryProviderDomainError&) {
        return reject(SteinerReject::ProviderDomainRejected);
    } catch (const std::bad_alloc&) {
        return fail(SteinerFail::CandidateAllocFailed);
    } catch (...) {
        return fail(SteinerFail::GeometryProviderError);
    }
    if (!(provStamp == expectedProvider))
        return reject(SteinerReject::StaleGeometryProvider);

    // Canonical, DIRECTION-INDEPENDENT SteinerNodeId assignment: order ALL
    // staged params by their reversal-symmetric affine identity (lo, hi, dyadic
    // value). lowerCandidate creates nodes in ascending SteinerNodeId order, so
    // committed node SLOTS - hence the semantic hash - are identical whether the
    // segment is recovered a->b or b->a. Chain order follows leg traversal.
    std::vector<std::vector<SteinerNodeId>> legStagedId(legs.size());
    {
        struct Ent { std::size_t li, j; };
        std::vector<Ent> flat;
        for (std::size_t li = 0; li < legs.size(); ++li)
            for (std::size_t j = 0; j < legs[li].params.size(); ++j)
                flat.push_back({li, j});
        std::sort(flat.begin(), flat.end(), [&](const Ent& x, const Ent& y) {
            const AffineParam& a = legs[x.li].params[x.j];
            const AffineParam& b = legs[y.li].params[y.j];
            if (a.lo != b.lo) return a.lo < b.lo;
            if (a.hi != b.hi) return a.hi < b.hi;
            return compareDyadic(a.num, a.level, b.num, b.level) < 0;
        });
        for (std::size_t li = 0; li < legs.size(); ++li)
            legStagedId[li].resize(legs[li].params.size());
        SteinerNodeId id = 1;
        for (const auto& e : flat) legStagedId[e.li][e.j] = id++;
    }

    // ----- PHASE B: provider mapping + per-leg retriangulation. Staged ids are
    // the canonical ones just assigned; the merged chain follows leg traversal
    // order (allocator-history independent, direction-canonical). -----
    SteinerRecoveryCandidate c;
    c.stamp = OwnershipStamp::of(ws);
    c.providerStamp = expectedProvider;
    c.entity = seg.id;
    c.chart = seg.chart;
    c.cellClass = cellClass;

    std::vector<TriNodeRef> chainRefs;          // full merged A..B chain refs
    chainRefs.push_back(TriNodeRef::Existing(A));
    for (std::size_t li = 0; li < legs.size(); ++li) {
        Leg& leg = legs[li];
        const Waypt& wj = way[li + 1];
        // staged nodes for this leg, in traversal (wi->wj) order; each carries
        // its canonical direction-independent id.
        std::vector<SteinerNode> legStaged;
        legStaged.reserve(leg.params.size());
        for (std::size_t j = 0; j < leg.params.size(); ++j) {
            const AffineParam& ap = leg.params[j];
            SteinerNode s;
            s.id = legStagedId[li][j];
            s.param = ap;
            s.chart = seg.chart;
            s.cls = seg.cls;
            Vec2 mv;
            try {
                mv = toMappingVec2(s.param, resolve);
            } catch (const std::invalid_argument&) {
                return reject(SteinerReject::NonfiniteMappedPosition);
            }
            try {
                s.physical = provider.pointAt(seg.chart, mv);
            } catch (const GeometryProviderDomainError&) {
                return reject(SteinerReject::ProviderDomainRejected);
            } catch (const std::bad_alloc&) {
                return fail(SteinerFail::CandidateAllocFailed);
            } catch (...) {
                return fail(SteinerFail::GeometryProviderError);
            }
            if (!(std::isfinite(s.physical.x) && std::isfinite(s.physical.y) &&
                  std::isfinite(s.physical.z)))
                return reject(SteinerReject::NonfiniteMappedPosition);
            legStaged.push_back(s);
        }
        if (!resolveOk) return reject(SteinerReject::PslgRejected);

        // unified node table + parallel EXACT affine coords for the ear-clip:
        // an existing node is its stored param; a staged node is its exact
        // AffineParam. No rounded mapping coordinate enters the ear-clip.
        std::vector<TriNodeRef> nodes;
        std::vector<predicates::AffinePoint2> AP;
        std::map<topo::NodeHandle, std::size_t> existIdx;
        std::map<SteinerNodeId, std::size_t> stagedIdx;
        auto idxExisting = [&](const topo::NodeHandle& h) {
            auto it = existIdx.find(h);
            if (it != existIdx.end()) return it->second;
            const std::size_t k = nodes.size();
            nodes.push_back(TriNodeRef::Existing(h));
            AP.push_back(predicates::storedPoint(paramOfHandle(h)));
            existIdx.emplace(h, k);
            return k;
        };
        auto idxStaged = [&](std::size_t si) {
            const SteinerNodeId id = legStaged[si].id;
            auto it = stagedIdx.find(id);
            if (it != stagedIdx.end()) return it->second;
            const std::size_t k = nodes.size();
            nodes.push_back(TriNodeRef::Staged(id));
            AP.push_back(toAffinePoint2(legStaged[si].param, resolve));
            stagedIdx.emplace(id, k);
            return k;
        };
        auto buildRing = [&](const std::vector<topo::NodeHandle>& side) {
            std::vector<std::size_t> ring;
            for (const auto& h : side) ring.push_back(idxExisting(h));  // wi..wj
            for (std::size_t s = legStaged.size(); s-- > 0;)
                ring.push_back(idxStaged(s));                          // rev
            return ring;
        };
        const std::vector<std::size_t> leftRing = buildRing(leg.left);
        const std::vector<std::size_t> rightRing = buildRing(leg.right);
        if (!resolveOk) return reject(SteinerReject::PslgRejected);
        auto leftTris = scdetail::earClipOriented(AP, leftRing, want);
        auto rightTris = scdetail::earClipOriented(AP, rightRing, want);
        if (!leftTris || !rightTris)
            return fail(SteinerFail::InvalidSteinerGeometry);
        auto emitTri = [&](const std::array<std::size_t, 3>& t) {
            TriNodeRefTriangle tri;
            tri.nodes = {nodes[t[0]], nodes[t[1]], nodes[t[2]]};
            tri.cls = cellClass;
            c.added.push_back(tri);
        };
        for (const auto& t : *leftTris) emitTri(t);
        for (const auto& t : *rightTris) emitTri(t);

        for (const auto& ch : leg.strip) c.removed.push_back(ch);
        for (const auto& s : legStaged) {
            chainRefs.push_back(TriNodeRef::Staged(s.id));
            c.staged.push_back(s);
        }
        chainRefs.push_back(TriNodeRef::Existing(wj.h));
    }

    // (7b) existing-node reuse (blocker 2): the segment passes THROUGH each
    //      interior waypoint, so leaving it on the cavity boundary would pinch
    //      the cavity (its boundary would visit the node twice - not a
    //      topological disk, which the frozen ASCR requires). Pull each reused
    //      node's remaining (non-crossed) incident cells into the cavity and
    //      re-add them UNCHANGED (all-existing triangles, orientation-normalized
    //      to `want`), so the reused node becomes a strictly INTERIOR existing
    //      vertex and the cavity boundary is one simple cycle. The frozen ASCR
    //      disk/node-link audit is rank-based and admits interior existing nodes
    //      (Cavity.hpp interiorNodes), so this remains ONE atomic commit.
    if (way.size() > 2) {
        std::set<topo::CellHandle> removedSet(c.removed.begin(), c.removed.end());
        for (std::size_t wi = 1; wi + 1 < way.size(); ++wi) {
            const auto sit = sess.index().star.find(way[wi].h);
            if (sit == sess.index().star.end())
                return fail(SteinerFail::InvalidSteinerGeometry);
            for (const auto& ch : sit->second) {
                if (!removedSet.insert(ch).second) continue;   // already in R
                const auto cd = ws.cells().get(ch);
                if (cd.cls.dim != cellClass.dim || cd.cls.tag != cellClass.tag)
                    return reject(SteinerReject::InvalidClassification);
                c.removed.push_back(ch);
                std::array<topo::NodeHandle, 3> n3{cd.nodes[0], cd.nodes[1],
                                                   cd.nodes[2]};
                const auto s = predicates::orient2dAffine(
                    predicates::storedPoint(paramOfHandle(n3[0])),
                    predicates::storedPoint(paramOfHandle(n3[1])),
                    predicates::storedPoint(paramOfHandle(n3[2])));
                if (s == predicates::Sign::Zero)
                    return fail(SteinerFail::InvalidSteinerGeometry);
                if (s != want) std::swap(n3[1], n3[2]);
                TriNodeRefTriangle tri;
                tri.nodes = {TriNodeRef::Existing(n3[0]),
                             TriNodeRef::Existing(n3[1]),
                             TriNodeRef::Existing(n3[2])};
                tri.cls = cellClass;
                c.added.push_back(tri);
            }
        }
        if (!resolveOk) return reject(SteinerReject::PslgRejected);
    }

    // (8) requirements = every consecutive edge of the full merged chain; the
    //     result chain is the ordered A .. B references.
    for (std::size_t i = 0; i + 1 < chainRefs.size(); ++i)
        c.requirements.push_back({chainRefs[i], chainRefs[i + 1], seg.cls});
    c.resultChain = chainRefs;

    // (9) GENERAL geometric candidate validator (blocker 4) - INDEPENDENT of
    //     the ear-clip and of the ASCR combinatorial audit, over authoritative
    //     AFFINE coordinates. It proves the added triangles form a valid planar
    //     triangulation (a geometric embedding): (a) each triangle nondegenerate
    //     + consistently oriented; (b) NO vertex strictly interior to any other
    //     triangle (overlap); (c) NO vertex strictly interior to any triangle
    //     EDGE it is not an endpoint of (unintended T-junction); (d) NO proper
    //     crossing between any two non-adjacent triangle edges; and (e) every
    //     staged node exactly on A-B and strictly interior. (a)+(b)+(c)+(d)
    //     together imply boundary simplicity and diagonal visibility.
    {
        auto sameRef = [](const TriNodeRef& a, const TriNodeRef& b) {
            if (a.staged != b.staged) return false;
            return a.staged ? a.stagedId == b.stagedId : a.existing == b.existing;
        };
        auto pointOf = [&](const TriNodeRef& ref,
                           bool& found) -> predicates::AffinePoint2 {
            if (!ref.staged)
                return predicates::storedPoint(paramOfHandle(ref.existing));
            for (const auto& s : c.staged)
                if (s.id == ref.stagedId)
                    return toAffinePoint2(s.param, resolve);
            found = false;
            return predicates::storedPoint(Vec2{});
        };

        GeneralSteinerGeometry geometry;
        std::vector<TriNodeRef> refs;
        geometry.triangles.reserve(c.added.size());
        for (const auto& t : c.added) {
            std::array<std::size_t, 3> tri{};
            for (std::size_t corner = 0; corner < 3; ++corner) {
                std::size_t index = refs.size();
                for (std::size_t i = 0; i < refs.size(); ++i)
                    if (sameRef(refs[i], t.nodes[corner])) {
                        index = i;
                        break;
                    }
                if (index == refs.size()) {
                    bool found = true;
                    const auto point = pointOf(t.nodes[corner], found);
                    if (!found || !resolveOk)
                        return fail(SteinerFail::InvalidSteinerGeometry);
                    refs.push_back(t.nodes[corner]);
                    geometry.vertices.push_back(
                        GeneralSteinerGeometryVertex{point,
                                                     t.nodes[corner].staged});
                }
                tri[corner] = index;
            }
            geometry.triangles.push_back(tri);
        }

        const GeneralSteinerGeometry* geometryToValidate = &geometry;
#if defined(PAACS_TOPO_TESTING)
        if (scdetail::generalSteinerGeometryTestOverride != nullptr)
            geometryToValidate =
                scdetail::generalSteinerGeometryTestOverride;
#endif
        if (validateGeneralSteinerGeometry(
                uvA, uvB, pAsp, pBsp, *geometryToValidate) !=
            GeneralSteinerGeometryIssue::None)
            return fail(SteinerFail::InvalidSteinerGeometry);
    }

    // (10) lower to exactly one SteinerCavitySpec.
    const topo::WorkspaceStamp wsStamp = topo::WorkspaceStamp::of(ws);
    auto low = lowerCandidate(c, wsStamp);
    if (!low) return fail(SteinerFail::CavityRejected);
    r.recoveredChain.reserve(low->loweredChain.size());
    r.createdParams.reserve(c.staged.size());

    // (11) ONE atomic ASCR commit (typed tiers) + copy/swap metadata
    //      publication for ALL staged nodes.
    topo::SteinerCavityCommit out;
    std::unordered_map<topo::NodeHandle, ParamPoint> nextParams;
    std::unordered_map<topo::NodeHandle, AuthoritativeParamRecord> nextAuth;
    {
        topo::SteinerCavityTransaction tx = ws.beginSteinerCavityTransaction();
        try {
            (void)tx.stage(low->spec);
        } catch (const std::invalid_argument&) {
            return fail(SteinerFail::CavityRejected);
        } catch (const std::bad_alloc&) {
            return fail(SteinerFail::CandidateAllocFailed);
        }
        try {
            tx.prepare();
        } catch (const topo::StagedNodeReservationError&) {
            tx.abandon();
            return fail(SteinerFail::StagedNodeAllocFailed);
        } catch (const topo::CavityPreparationError&) {
            tx.abandon();
            return fail(SteinerFail::CavityThrew);
        } catch (const std::bad_alloc&) {
            tx.abandon();
            return fail(SteinerFail::CavityThrew);
        }
        try {
            nextParams = cs.params;
            nextAuth = cs.authParams;
            const auto predicted = tx.predictedStagedNodes();
            for (const auto& sn : c.staged) {
                std::uint32_t local = 0;
                for (const auto& kv : low->idToLocal)
                    if (kv.first == sn.id) local = kv.second;
                const topo::NodeHandle h = predicted[local];
                const Vec2 muv = toMappingVec2(sn.param, resolve);
                if (!nextParams.emplace(h, ParamPoint{muv, sn.chart}).second)
                    throw std::logic_error("shadow param collision");
                AuthoritativeParamRecord rec;
                if (const auto* ap = std::get_if<AffineParam>(&sn.param)) {
                    rec.affine = true;
                    rec.lo = ap->lo; rec.hi = ap->hi;
                    rec.num = ap->num; rec.level = ap->level;
                } else {
                    rec.affine = false;
                    rec.storedUv = std::get<StoredParam>(sn.param).uv;
                }
                rec.chart = sn.chart;
                if (!nextAuth.emplace(h, rec).second)
                    throw std::logic_error("shadow auth collision");
            }
            out = tx.commit();
        } catch (...) {
            tx.abandon();
            return fail(SteinerFail::CavityThrew);
        }
    }
    cs.params.swap(nextParams);
    cs.authParams.swap(nextAuth);
    for (const auto& ref : low->loweredChain)
        r.recoveredChain.push_back(ref.isStaged ? out.nodeOf(ref.localId)
                                                : ref.existing);
    for (const auto& s : c.staged) {
        std::uint32_t local = 0;
        for (const auto& kv : low->idToLocal)
            if (kv.first == s.id) local = kv.second;
        r.createdParams.emplace_back(out.nodeOf(local), s.param);
    }
    cs.stamp = OwnershipStamp::of(ws);
    r.revision = out.revision;
    r.status = SteinerStatus::Recovered;
    return r;
}

// Single-segment opener: runs the whole-set preflight (fine-grained A1 mapping
// for specific reasons, authoritative-store audit, full RecoverySession::open)
// ONCE, then delegates to the validated-session entity helper.
inline SteinerResult recoverSteinerSegmentImpl(
    topo::TopologyWorkspace& ws, ConstraintSet& cs,
    const ConstraintSegment& seg, const GeometryProvider& provider,
    const GeometryProviderStamp& expectedProvider,
    const SteinerPolicy& policy) {
    SteinerResult r;
    auto reject = [&](SteinerReject x) {
        r.status = SteinerStatus::Rejected; r.reject = x; return r;
    };
    auto fail = [&](SteinerFail x) {
        r.status = SteinerStatus::Failed; r.fail = x; return r;
    };
    {
        const auto issues = preflightConstraintSet(ws, cs);
        auto has = [&](InvalidInputReason why) {
            for (const auto& i : issues)
                if (i.reason == why) return true;
            return false;
        };
        if (has(InvalidInputReason::StaleStamp))
            return reject(SteinerReject::StaleCandidateStamp);
        if (has(InvalidInputReason::InvalidClassification))
            return reject(SteinerReject::InvalidClassification);
        if (has(InvalidInputReason::ChartMismatch))
            return reject(SteinerReject::ChartMismatch);
        if (!issues.empty())
            return reject(SteinerReject::PslgRejected);
    }
    if (!auditAuthParams(ws, cs))
        return reject(SteinerReject::PslgRejected);
    SessionStatus ss{};
    auto sessOpt = RecoverySession::open(ws, cs, ss);
    if (!sessOpt) {
        switch (ss) {
            case SessionStatus::InputStampMismatch:
                return reject(SteinerReject::StaleCandidateStamp);
            case SessionStatus::IndexBuildFailed:
                return fail(SteinerFail::IndexRebuildFailed);
            default:
                return reject(SteinerReject::PslgRejected);
        }
    }
    return recoverSteinerEntityValidated(ws, cs, seg, provider,
                                         expectedProvider, policy, *sessOpt);
}

// Public entry: LEAK-PROOF over ANY stdlib. The off-workspace candidate build
// can allocate; a std::bad_alloc escaping there is classified
// Failed(CandidateAllocFailed) (§9, tier 1) - nothing is committed because the
// whole allocating candidate build precedes the ASCR commit, and the
// post-commit publication is allocation-free (all buffers pre-reserved). ASCR
// stage/prepare/commit exceptions are already classified inside the impl.
inline SteinerResult recoverSteinerSegment(
    topo::TopologyWorkspace& ws, ConstraintSet& cs,
    const ConstraintSegment& seg, const GeometryProvider& provider,
    const GeometryProviderStamp& expectedProvider,
    const SteinerPolicy& policy = {}) {
    try {
        return recoverSteinerSegmentImpl(ws, cs, seg, provider,
                                         expectedProvider, policy);
    } catch (const std::bad_alloc&) {
        SteinerResult r;
        r.status = SteinerStatus::Failed;
        r.fail = SteinerFail::CandidateAllocFailed;
        return r;
    }
}

// ---------------------------------------------------------------------------
// Session/set-level wrapper (§9 completeness). Processes every segment of the
// ConstraintSet in ascending entity id with a COMPLETE report: each input
// entity appears exactly once in `processed` XOR `unattempted`. Prefix rule:
// CONTINUE after AlreadyPresent/Recovered; STOP at the first Refused / Failed /
// Rejected, listing all later entities as unattempted. Leak-proof (each
// per-entity call is). The report is reserved up front so it cannot itself
// fail to record an entity.
// ---------------------------------------------------------------------------
struct SteinerEntityResult {
    ConstraintEntityId entity = 0;
    SteinerResult result;
};
struct SteinerSegmentsReport {
    std::vector<SteinerEntityResult> processed;   // ascending entity id
    std::vector<ConstraintEntityId> unattempted;  // prefix-stopped remainder
};
inline SteinerSegmentsReport recoverSteinerSegmentsImpl(
    topo::TopologyWorkspace& ws, ConstraintSet& cs,
    const GeometryProvider& provider,
    const GeometryProviderStamp& expectedProvider,
    const SteinerPolicy& policy) {
    SteinerSegmentsReport rep;
    rep.processed.reserve(cs.segments.size());
    rep.unattempted.reserve(cs.segments.size());

    std::vector<ConstraintSegment> segs = cs.segments;
    std::sort(segs.begin(), segs.end(),
              [](const ConstraintSegment& a, const ConstraintSegment& b) {
                  return a.id < b.id;
              });

    // Itemise EVERY entity as Rejected(reason), ascending, unattempted empty
    // (§9 whole-set preflight vector) - distinct from a per-entity prefix stop.
    auto fillAllRejected = [&](SteinerReject rj) {
        rep.processed.clear();
        rep.unattempted.clear();
        for (const auto& s : segs) {
            SteinerResult r;
            r.status = SteinerStatus::Rejected;
            r.reject = rj;
            rep.processed.push_back({s.id, r});
        }
    };

    // (A) ONE whole-set preflight BEFORE the loop and BEFORE any provider call
    //     (blockers 1+2). Fine-grained A1 mapping, then authoritative-store
    //     audit, then the full session (A2 + PSLG + index).
    {
        const auto issues = preflightConstraintSet(ws, cs);
        auto has = [&](InvalidInputReason why) {
            for (const auto& i : issues)
                if (i.reason == why) return true;
            return false;
        };
        if (has(InvalidInputReason::StaleStamp)) {
            fillAllRejected(SteinerReject::StaleCandidateStamp);
            return rep;
        }
        if (has(InvalidInputReason::InvalidClassification)) {
            fillAllRejected(SteinerReject::InvalidClassification);
            return rep;
        }
        if (has(InvalidInputReason::ChartMismatch)) {
            fillAllRejected(SteinerReject::ChartMismatch);
            return rep;
        }
        if (!issues.empty()) {
            fillAllRejected(SteinerReject::PslgRejected);
            return rep;
        }
    }
    if (!auditAuthParams(ws, cs)) {
        fillAllRejected(SteinerReject::PslgRejected);
        return rep;
    }
    // The set wrapper OWNS exactly one validated session and passes it into
    // each entity; it is rebuilt only after a genuine commit (never inside an
    // entity, never on AlreadyPresent).
    std::optional<RecoverySession> session;
    {
        SessionStatus ss{};
        session = RecoverySession::open(ws, cs, ss);
        if (!session) {
            fillAllRejected(ss == SessionStatus::InputStampMismatch
                                ? SteinerReject::StaleCandidateStamp
                                : SteinerReject::PslgRejected);
            return rep;
        }
    }

    // (B) per-entity loop with prefix stopping + 13.2a index-rebuild
    //     attribution. Each entity uses the OWNED validated session (no re-open
    //     inside). The 13.2a rebuild happens ONLY after a real commit and only
    //     when another entity remains; it is guarded so a rebuild allocation
    //     failure becomes the NEXT entity's Failed(IndexRebuildFailed) rather
    //     than an escaped exception. AlreadyPresent changes nothing and never
    //     triggers a rebuild (fix: no false attribution).
    bool stop = false, indexFailedForNext = false;
    for (std::size_t i = 0; i < segs.size(); ++i) {
        const ConstraintSegment& seg = segs[i];
        if (stop) {
            rep.unattempted.push_back(seg.id);
            continue;
        }
        if (indexFailedForNext) {
            SteinerResult r;
            r.status = SteinerStatus::Failed;
            r.fail = SteinerFail::IndexRebuildFailed;
            rep.processed.push_back({seg.id, std::move(r)});
            stop = true;   // committed prefix preserved; later unattempted
            continue;
        }
        SteinerResult r = recoverSteinerEntityValidated(
            ws, cs, seg, provider, expectedProvider, policy, *session);
        const bool committed = r.status == SteinerStatus::Recovered;
        const bool ok = committed || r.status == SteinerStatus::AlreadyPresent;
        rep.processed.push_back({seg.id, std::move(r)});
        if (!ok) {
            stop = true;   // Refused/Failed/Rejected -> prefix STOP
            continue;
        }
        // Rebuild the session ONLY after a genuine commit and only if a later
        // entity will use it. AlreadyPresent retains the current session.
        if (committed && i + 1 < segs.size()) {
            try {
                SessionStatus ss2{};
                auto s2 = RecoverySession::open(ws, cs, ss2);
                if (!s2) indexFailedForNext = true;
                else session = std::move(s2);
            } catch (const std::bad_alloc&) {
                indexFailedForNext = true;
            }
        }
    }
    return rep;
}

// Public set-level entry: LEAK-PROOF. Every allocatable setup step (report
// reservation, segment copy/sort, whole-set preflight) is inside the backstop;
// a std::bad_alloc with no prefix recorded becomes a COMPLETE structured
// whole-set rejection (retried once), falling back to a bare empty report only
// under true exhaustion.
inline SteinerSegmentsReport recoverSteinerSegments(
    topo::TopologyWorkspace& ws, ConstraintSet& cs,
    const GeometryProvider& provider,
    const GeometryProviderStamp& expectedProvider,
    const SteinerPolicy& policy = {}) {
    try {
        return recoverSteinerSegmentsImpl(ws, cs, provider, expectedProvider,
                                          policy);
    } catch (const std::bad_alloc&) {
        SteinerSegmentsReport rep;
        try {
            std::vector<ConstraintEntityId> ids;
            ids.reserve(cs.segments.size());
            for (const auto& s : cs.segments) ids.push_back(s.id);
            std::sort(ids.begin(), ids.end());
            rep.processed.reserve(ids.size());
            for (const auto id : ids) {
                SteinerResult r;
                r.status = SteinerStatus::Failed;
                r.fail = SteinerFail::CandidateAllocFailed;
                rep.processed.push_back({id, r});
            }
        } catch (...) {
            rep = SteinerSegmentsReport{};   // bare fallback (true exhaustion)
        }
        return rep;
    }
}

}  // namespace paacs::bnd
