// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// paacs::predicates - the ONLY sanctioned source of geometric sign
// decisions (M3). Every topological choice (orientation, circumsphere /
// circumcircle membership) must route through these four predicates; ad hoc
// epsilon comparisons on determinants are forbidden in new code.
//
// Implementation: filtered exact evaluation. A fast double-precision path
// with a conservative forward error bound answers the overwhelming
// majority of calls; uncertain cases fall back to exact arithmetic using
// Shewchuk-style floating-point expansions (two-sum / FMA two-product /
// zero-eliminating grow, sum and scale). The exact path returns the true
// sign, including exact Zero for genuinely degenerate configurations
// (cospherical, coplanar, duplicate) - callers define deterministic Zero
// policies and must document them.
//
// Conventions (locked by unit tests):
//   orient2d(a,b,c)      > 0  iff  c lies to the LEFT of the directed line
//                              a->b (a,b,c counter-clockwise).
//   incircle(a,b,c,d)    > 0  iff  d lies strictly INSIDE the circle
//                              through a,b,c, with (a,b,c) CCW.
//   orient3d(a,b,c,d)    > 0  iff  (b-a) x (c-a) . (d-a) > 0
//                              (d on the positive side of plane abc).
//   insphere(a,b,c,d,e)  > 0  iff  e lies strictly INSIDE the sphere
//                              through a,b,c,d, with orient3d(a,b,c,d) > 0.
//
// Numerical envelope (part of the public contract): results are exact for
// the double values PASSED IN, provided intermediate products of
// coordinate differences stay within the NORMAL double range. Differences
// whose products underflow to subnormal zero (external differential
// testing observed this for perturbations at the 5e-324 minimum
// subnormal) may report Zero for a mathematically nonzero determinant.
// Callers must pre-scale into a sane frame; the Delaunay kernels do this
// via their normalized working frames. Note the corollary: a caller that
// normalizes gets answers exact RELATIVE TO ITS ROUNDED NORMALIZED
// COORDINATES, not to the pre-transform inputs. A formally specified
// power-of-two WorkingFrame is the M3-B refinement.
//
#include <cstdint>
#include <type_traits>

#include "paacs/core/Types.hpp"

namespace paacs::predicates {

enum class Sign : int { Negative = -1, Zero = 0, Positive = 1 };

Sign orient2d(const Vec2& a, const Vec2& b, const Vec2& c);
Sign incircle(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& d);
Sign orient3d(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d);
// Metric-only estimate of the exact binary64-coordinate determinant, computed
// with expansions to retain cancellation tails. Divide by six for tet volume.
// Topological decisions must still use orient3d(), never this rounded value.
// Same normal-product numerical envelope as orient3d; nonfinite inputs give NaN.
// Does not update predicate counters and does not allocate.
double orient3dDeterminantEstimate(const Vec3& a, const Vec3& b,
                                  const Vec3& c, const Vec3& d) noexcept;
Sign insphere(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d,
              const Vec3& e);

// Delaunay-private instrumentation-free entry points. They execute the exact
// same filter and exact fallback as orient3d()/insphere(), but omit the global
// process statistics increments from the performance-critical path. Tests and
// diagnostic tools retain the instrumented public functions above.
Sign orient3dFiltered(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d);
Sign insphereFiltered(const Vec3& a, const Vec3& b, const Vec3& c,
                      const Vec3& d, const Vec3& e);

// Exact checked power predicates used by bounded interior Tet refinement.
// They evaluate the original binary64 values directly. Invalid geometry has a
// structured status and always carries Sign::Zero; allocation failure from the
// exact arbitrary-width path propagates as std::bad_alloc.
enum class PowerPredicateStatus : std::uint8_t {
    Valid,
    NonFiniteInput,
    DegenerateEdge,
    DegenerateTriangle
};

struct PowerPredicateResult {
    PowerPredicateStatus status = PowerPredicateStatus::Valid;
    Sign sign = Sign::Zero;
};

PowerPredicateResult edgeDiametralPower(
    const Vec3& a, const Vec3& b, const Vec3& p);
PowerPredicateResult triangleCircumballPower(
    const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& p);

static_assert(std::is_standard_layout_v<PowerPredicateResult>);
static_assert(std::is_trivially_copyable_v<PowerPredicateResult>);
static_assert(std::is_nothrow_copy_constructible_v<PowerPredicateResult>);
static_assert(std::is_nothrow_move_constructible_v<PowerPredicateResult>);

// ---------------------------------------------------------------------------
// Exact AFFINE-point predicates (BND-2A / T1-A support; additive, generic -
// no downstream identifiers appear here).
//
// An AffinePoint2 denotes the exact point  lo + t*(hi - lo)  with the DYADIC
// parameter  t = numerator / 2^level  measured from lo to hi. A STORED point
// (an ordinary Vec2) is represented as lo == hi, numerator == 0, level == 0
// (t = 0, so the point is lo). A GENUINE affine point has 0 < numerator <
// 2^level with numerator reduced (odd).
//
// level is bounded by kMaxAffineLevel = 52 so that numerator, 2^level, and
// 2^level - numerator are all EXACTLY representable as doubles. orient2dAffine
// returns the EXACT sign of orient2d over the authoritative affine points: it
// expands each argument into its (one or two) stored corners and combines the
// up-to-eight exact stored-point orient2d determinants with exact
// integer-weighted expansion arithmetic (Shewchuk scale/sum) - it NEVER uses
// the rounded mapping coordinate. Zero/one/two/three affine arguments are all
// supported. Same numerical envelope as the stored predicates above:
// coordinates and their weighted products must stay in the normal double
// range.
// Valid representations are EXACTLY: (stored) level==0, numerator==0,
// lo==hi, finite coordinates; or (genuine affine) 1<=level<=52,
// 0<numerator<2^level, ODD numerator (canonical reduced form), finite
// endpoints. Both entry points below VALIDATE every argument first and
// throw std::invalid_argument on any other representation - an invalid token
// is never silently read as an exact point.
struct AffinePoint2 {
    Vec2          lo{};
    Vec2          hi{};
    std::uint64_t numerator = 0;
    std::uint8_t  level = 0;
};
inline constexpr std::uint8_t kMaxAffineLevel = 52;

// Convenience: an AffinePoint2 for a plain stored point.
inline AffinePoint2 storedPoint(const Vec2& p) { return AffinePoint2{p, p, 0, 0}; }

Sign orient2dAffine(const AffinePoint2& a, const AffinePoint2& b,
                    const AffinePoint2& c);

// Exact 1-D ordering along the directed segment loSeg -> hiSeg: the sign of
// (t_p - t_q), where t is each point's parameter on the segment's DOMINANT
// axis. Supports Stored and Affine in any combination (both directions).
// Points are assumed to lie on the segment's supporting line (the BND chain
// invariant); ordering is exact via expansion arithmetic on the dominant
// coordinate. Generic - takes only Vec2 / AffinePoint2.
Sign orderOnSegmentAffine(const Vec2& loSeg, const Vec2& hiSeg,
                          const AffinePoint2& p, const AffinePoint2& q);

// Exact sign of (value(p) - value(q)) on a single axis (axis 0 = u, 1 = v),
// evaluated over the AUTHORITATIVE affine coordinates with no rounding. Enables
// exact lexicographic ordering of affine points and exact closed-segment
// containment (collinear + between) without a Vec2 segment - used by the
// Steiner ear-clip / candidate validator so that staged (affine) node positions
// never enter topology decisions through their rounded mapping coordinates.
Sign compareOnAxisAffine(const AffinePoint2& p, const AffinePoint2& q, int axis);

// Instrumentation: how often the fast filter decided vs. the exact
// fallback ran (per process, monotonic; relaxed atomics).
struct Stats {
    std::uint64_t orient2dFast = 0, orient2dExact = 0;
    std::uint64_t incircleFast = 0, incircleExact = 0;
    std::uint64_t orient3dFast = 0, orient3dExact = 0;
    std::uint64_t insphereFast = 0, insphereExact = 0;
};
Stats stats();
void resetStats();

} // namespace paacs::predicates
