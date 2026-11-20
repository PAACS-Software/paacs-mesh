// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/predicates/Predicates.hpp"

#include <atomic>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace paacs::predicates {

namespace {

// A deliberately small private exact dyadic accumulator. Every finite double
// is decoded as integer*2^exponent; additions align integer magnitudes and
// products use base-2^32 schoolbook multiplication. This path is independent
// of coordinate translation/normalization and therefore remains exact across
// the complete binary64 exponent range. Vector growth is intentionally
// fallible: std::bad_alloc is part of the checked-predicate contract.
struct BigUnsigned {
    std::vector<std::uint32_t> words;

    bool zero() const noexcept { return words.empty(); }

    void trim() noexcept {
        while (!words.empty() && words.back() == 0u) words.pop_back();
    }
};

int compare(const BigUnsigned& a, const BigUnsigned& b) noexcept {
    if (a.words.size() != b.words.size())
        return a.words.size() < b.words.size() ? -1 : 1;
    for (std::size_t i = a.words.size(); i != 0; --i) {
        const auto av = a.words[i - 1], bv = b.words[i - 1];
        if (av != bv) return av < bv ? -1 : 1;
    }
    return 0;
}

BigUnsigned addMagnitude(const BigUnsigned& a, const BigUnsigned& b) {
    BigUnsigned out;
    const std::size_t n = (std::max)(a.words.size(), b.words.size());
    out.words.resize(n + 1u, 0u);
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const std::uint64_t av = i < a.words.size() ? a.words[i] : 0u;
        const std::uint64_t bv = i < b.words.size() ? b.words[i] : 0u;
        const std::uint64_t sum = av + bv + carry;
        out.words[i] = static_cast<std::uint32_t>(sum);
        carry = sum >> 32u;
    }
    out.words[n] = static_cast<std::uint32_t>(carry);
    out.trim();
    return out;
}

// Precondition a >= b.
BigUnsigned subtractMagnitude(const BigUnsigned& a, const BigUnsigned& b) {
    BigUnsigned out;
    out.words.resize(a.words.size(), 0u);
    std::uint64_t borrow = 0;
    for (std::size_t i = 0; i < a.words.size(); ++i) {
        const std::uint64_t av = a.words[i];
        const std::uint64_t bv =
            (i < b.words.size() ? b.words[i] : 0u) + borrow;
        out.words[i] = static_cast<std::uint32_t>(av - bv);
        borrow = av < bv ? 1u : 0u;
    }
    out.trim();
    return out;
}

BigUnsigned shiftLeft(const BigUnsigned& a, unsigned bits) {
    if (a.zero() || bits == 0u) return a;
    const std::size_t wordShift = bits / 32u;
    const unsigned bitShift = bits % 32u;
    BigUnsigned out;
    out.words.assign(a.words.size() + wordShift + (bitShift ? 1u : 0u), 0u);
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < a.words.size(); ++i) {
        const std::uint64_t value =
            (static_cast<std::uint64_t>(a.words[i]) << bitShift) | carry;
        out.words[i + wordShift] = static_cast<std::uint32_t>(value);
        carry = value >> 32u;
    }
    if (bitShift) out.words[a.words.size() + wordShift] =
        static_cast<std::uint32_t>(carry);
    out.trim();
    return out;
}

BigUnsigned multiplyMagnitude(const BigUnsigned& a, const BigUnsigned& b) {
    BigUnsigned out;
    if (a.zero() || b.zero()) return out;
    out.words.assign(a.words.size() + b.words.size(), 0u);
    for (std::size_t i = 0; i < a.words.size(); ++i) {
        std::uint64_t carry = 0;
        for (std::size_t j = 0; j < b.words.size(); ++j) {
            const std::size_t k = i + j;
            const std::uint64_t product =
                static_cast<std::uint64_t>(a.words[i]) * b.words[j];
            const std::uint64_t sum = product + out.words[k] + carry;
            out.words[k] = static_cast<std::uint32_t>(sum);
            carry = sum >> 32u;
        }
        std::size_t k = i + b.words.size();
        while (carry != 0u) {
            const std::uint64_t sum =
                static_cast<std::uint64_t>(out.words[k]) + carry;
            out.words[k] = static_cast<std::uint32_t>(sum);
            carry = sum >> 32u;
            ++k;
        }
    }
    out.trim();
    return out;
}

struct ExactDyadic {
    BigUnsigned magnitude;
    int exponent = 0;
    bool negative = false;

    bool zero() const noexcept { return magnitude.zero(); }
};

ExactDyadic exactDouble(double value) {
    ExactDyadic out;
    const std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
    const std::uint64_t fraction = bits & ((std::uint64_t{1} << 52u) - 1u);
    const unsigned biased = static_cast<unsigned>((bits >> 52u) & 0x7ffu);
    std::uint64_t significand = fraction;
    if (biased != 0u) significand |= std::uint64_t{1} << 52u;
    if (significand == 0u) return out;
    out.magnitude.words.push_back(static_cast<std::uint32_t>(significand));
    out.magnitude.words.push_back(static_cast<std::uint32_t>(significand >> 32u));
    out.magnitude.trim();
    out.exponent = biased == 0u ? -1074 : static_cast<int>(biased) - 1023 - 52;
    out.negative = (bits >> 63u) != 0u;
    return out;
}

ExactDyadic negate(ExactDyadic value) noexcept {
    if (!value.zero()) value.negative = !value.negative;
    return value;
}

ExactDyadic addExact(const ExactDyadic& a, const ExactDyadic& b) {
    if (a.zero()) return b;
    if (b.zero()) return a;
    const int exp = (std::min)(a.exponent, b.exponent);
    const auto am = shiftLeft(a.magnitude,
        static_cast<unsigned>(a.exponent - exp));
    const auto bm = shiftLeft(b.magnitude,
        static_cast<unsigned>(b.exponent - exp));
    ExactDyadic out;
    out.exponent = exp;
    if (a.negative == b.negative) {
        out.magnitude = addMagnitude(am, bm);
        out.negative = a.negative;
    } else {
        const int order = compare(am, bm);
        if (order == 0) return out;
        if (order > 0) {
            out.magnitude = subtractMagnitude(am, bm);
            out.negative = a.negative;
        } else {
            out.magnitude = subtractMagnitude(bm, am);
            out.negative = b.negative;
        }
    }
    return out;
}

ExactDyadic subtractExact(const ExactDyadic& a, const ExactDyadic& b) {
    return addExact(a, negate(b));
}

ExactDyadic multiplyExact(const ExactDyadic& a, const ExactDyadic& b) {
    ExactDyadic out;
    if (a.zero() || b.zero()) return out;
    out.magnitude = multiplyMagnitude(a.magnitude, b.magnitude);
    out.exponent = a.exponent + b.exponent;
    out.negative = a.negative != b.negative;
    return out;
}

Sign exactSign(const ExactDyadic& value) noexcept {
    if (value.zero()) return Sign::Zero;
    return value.negative ? Sign::Negative : Sign::Positive;
}

using ExactVec3 = std::array<ExactDyadic, 3>;

ExactVec3 exactVector(const Vec3& value) {
    return {exactDouble(value.x), exactDouble(value.y), exactDouble(value.z)};
}

ExactVec3 subtractExact(const ExactVec3& a, const ExactVec3& b) {
    return {subtractExact(a[0], b[0]), subtractExact(a[1], b[1]),
            subtractExact(a[2], b[2])};
}

ExactVec3 crossExact(const ExactVec3& a, const ExactVec3& b) {
    return {
        subtractExact(multiplyExact(a[1], b[2]), multiplyExact(a[2], b[1])),
        subtractExact(multiplyExact(a[2], b[0]), multiplyExact(a[0], b[2])),
        subtractExact(multiplyExact(a[0], b[1]), multiplyExact(a[1], b[0]))};
}

ExactDyadic dotExact(const ExactVec3& a, const ExactVec3& b) {
    return addExact(addExact(multiplyExact(a[0], b[0]),
                             multiplyExact(a[1], b[1])),
                    multiplyExact(a[2], b[2]));
}

bool finite(const Vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

// ---------------------------------------------------------------------------
// Error-free transformations. twoSum is Knuth's branch-free version (valid
// for any rounding-to-nearest doubles); twoProd uses FMA, which yields the
// exact residual of a product in one instruction.
// ---------------------------------------------------------------------------
inline void twoSum(double a, double b, double& x, double& y) {
    x = a + b;
    const double bv = x - a;
    const double av = x - bv;
    y = (a - av) + (b - bv);
}
inline void twoDiff(double a, double b, double& x, double& y) {
    twoSum(a, -b, x, y);   // negation of a double is exact
}
inline void twoProd(double a, double b, double& x, double& y) {
    x = a * b;
    y = std::fma(a, b, -x);
}
// Precondition: |a| >= |b| (guaranteed at its call sites below).
inline void fastTwoSum(double a, double b, double& x, double& y) {
    x = a + b;
    y = b - (x - a);
}

// ---------------------------------------------------------------------------
// Expansions: arrays of doubles, nonoverlapping, increasing magnitude,
// zero components eliminated. The sign of a (zero-eliminated) expansion is
// the sign of its last component. All routines follow Shewchuk's
// GROW-EXPANSION / EXPANSION-SUM / SCALE-EXPANSION with zero elimination.
// ---------------------------------------------------------------------------

// h = e + b (single double). h may not alias e.
int growZE(int elen, const double* e, double b, double* h) {
    double Q = b;
    int hlen = 0;
    for (int i = 0; i < elen; ++i) {
        double Qn, hh;
        twoSum(Q, e[i], Qn, hh);
        Q = Qn;
        if (hh != 0.0) h[hlen++] = hh;
    }
    if (Q != 0.0 || hlen == 0) h[hlen++] = Q;
    return hlen;
}

// h = e + f, via repeated grow (Shewchuk's EXPANSION-SUM; O(m*n) but our
// operands are small and zero elimination keeps them smaller in practice).
// h may not alias e or f.
int sumZE(int elen, const double* e, int flen, const double* f, double* h,
          double* scratch) {
    if (elen == 0) { std::memcpy(h, f, size_t(flen) * 8); return flen; }
    int len = elen;
    const double* src = e;
    double* bufs[2] = {h, scratch};
    int which = 0;
    for (int j = 0; j < flen; ++j) {
        double* dst = bufs[which];
        len = growZE(len, src, f[j], dst);
        src = dst;
        which ^= 1;
    }
    if (src != h) std::memcpy(h, src, size_t(len) * 8);
    return len;
}

// h = e * b. h may not alias e.
int scaleZE(int elen, const double* e, double b, double* h) {
    if (elen == 0 || b == 0.0) { h[0] = 0.0; return 1; }
    double Q, hh;
    int hlen = 0;
    twoProd(e[0], b, Q, hh);
    if (hh != 0.0) h[hlen++] = hh;
    for (int i = 1; i < elen; ++i) {
        double T, t, Qn;
        twoProd(e[i], b, T, t);
        twoSum(Q, t, Qn, hh);
        if (hh != 0.0) h[hlen++] = hh;
        fastTwoSum(T, Qn, Q, hh);
        if (hh != 0.0) h[hlen++] = hh;
    }
    if (Q != 0.0 || hlen == 0) h[hlen++] = Q;
    return hlen;
}

Sign signOf(int len, const double* e) {
    // Zero-eliminated expansions keep at most one trailing zero (the case
    // where the whole expansion is zero).
    const double top = e[len - 1];
    if (top > 0.0) return Sign::Positive;
    if (top < 0.0) return Sign::Negative;
    return Sign::Zero;
}

// Fixed-capacity scratch for the small exact paths (orient2d/3d, incircle);
// insphere's larger products use thread_local vectors.
struct Buf192 { double d[192]; };

// (u0,u1) * (v0,v1) for 2-expansions -> up to 8 components.
int mul22(const double* u, const double* v, double* h, double* s1,
          double* s2) {
    const int l1 = scaleZE(2, u, v[0], s1);
    const int l2 = scaleZE(2, u, v[1], s2);
    double scratch[16];
    return sumZE(l1, s1, l2, s2, h, scratch);
}

// det2 of 2-expansions: u*w - v*z  -> up to 16 components.
int det2exp(const double* u, const double* w, const double* v,
            const double* z, double* h) {
    double p1[8], p2[8], s1[8], s2[8], scratch[24];
    const int l1 = mul22(u, w, p1, s1, s2);
    int l2 = mul22(v, z, p2, s1, s2);
    for (int i = 0; i < l2; ++i) p2[i] = -p2[i];
    return sumZE(l1, p1, l2, p2, h, scratch);
}

// det3 of rows of 2-expansions r0,r1,r2 (each row = 3 two-expansions laid
// out as row[0..1], row[2..3], row[4..5]) -> up to 192 components.
int det3exp(const double* r0, const double* r1, const double* r2,
            double* out, double* scratchA, double* scratchB) {
    double m0[16], m1[16], m2[16];
    const int lm0 = det2exp(r1 + 2, r2 + 4, r1 + 4, r2 + 2, m0);
    const int lm1 = det2exp(r1 + 0, r2 + 4, r1 + 4, r2 + 0, m1);
    const int lm2 = det2exp(r1 + 0, r2 + 2, r1 + 2, r2 + 0, m2);

    // t0 = r0x * m0 ; t1 = -r0y * m1 ; t2 = r0z * m2 ; out = t0 + t1 + t2
    double t[3][64];
    int tl[3];
    for (int k = 0; k < 3; ++k) {
        const double* m = k == 0 ? m0 : k == 1 ? m1 : m2;
        const int lm = k == 0 ? lm0 : k == 1 ? lm1 : lm2;
        double a[32], b[32];
        const int la = scaleZE(lm, m, r0[2 * k + 0], a);
        const int lb = scaleZE(lm, m, r0[2 * k + 1], b);
        double sc[64];
        tl[k] = sumZE(la, a, lb, b, t[k], sc);
        if (k == 1)
            for (int i = 0; i < tl[k]; ++i) t[k][i] = -t[k][i];
    }
    const int l01 = sumZE(tl[0], t[0], tl[1], t[1], scratchA, scratchB);
    return sumZE(l01, scratchA, tl[2], t[2], out, scratchB);
}

// ---------------------------------------------------------------------------
// Statistics.
// ---------------------------------------------------------------------------
std::atomic<std::uint64_t> gO2F{0}, gO2E{0}, gICF{0}, gICE{0}, gO3F{0},
    gO3E{0}, gISF{0}, gISE{0};
std::atomic<bool> gDelaunayZeroSeen{false};
inline void bump(std::atomic<std::uint64_t>& c) {
    c.fetch_add(1, std::memory_order_relaxed);
}

// Conservative static filter constants: Shewchuk's first-stage error bounds
// are ~(4..17)*eps ~ 1e-15 x permanent; we use margins >10x looser, trading
// a few unnecessary exact evaluations for an unarguably safe filter.
constexpr double kO2Bound = 1e-14;
constexpr double kO3Bound = 1e-13;
constexpr double kICBound = 1e-13;
constexpr double kISBound = 1e-12;

} // namespace

PowerPredicateResult edgeDiametralPower(
    const Vec3& a, const Vec3& b, const Vec3& p) {
    if (!finite(a) || !finite(b) || !finite(p))
        return {PowerPredicateStatus::NonFiniteInput, Sign::Zero};

    const ExactVec3 ea = exactVector(a);
    const ExactVec3 eb = exactVector(b);
    const ExactVec3 ep = exactVector(p);
    const ExactVec3 edge = subtractExact(eb, ea);
    if (edge[0].zero() && edge[1].zero() && edge[2].zero())
        return {PowerPredicateStatus::DegenerateEdge, Sign::Zero};
    const ExactDyadic value =
        dotExact(subtractExact(ep, ea), subtractExact(ep, eb));
    return {PowerPredicateStatus::Valid, exactSign(value)};
}

PowerPredicateResult triangleCircumballPower(
    const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& p) {
    if (!finite(a) || !finite(b) || !finite(c) || !finite(p))
        return {PowerPredicateStatus::NonFiniteInput, Sign::Zero};

    const ExactVec3 ea = exactVector(a);
    const ExactVec3 u = subtractExact(exactVector(b), ea);
    const ExactVec3 v = subtractExact(exactVector(c), ea);
    const ExactVec3 w = subtractExact(exactVector(p), ea);
    const ExactVec3 n = crossExact(u, v);
    if (n[0].zero() && n[1].zero() && n[2].zero())
        return {PowerPredicateStatus::DegenerateTriangle, Sign::Zero};

    const ExactDyadic uu = dotExact(u, u);
    const ExactDyadic vv = dotExact(v, v);
    const ExactDyadic nn = dotExact(n, n);
    const ExactDyadic ww = dotExact(w, w);
    const ExactVec3 vCrossN = crossExact(v, n);
    const ExactVec3 nCrossU = crossExact(n, u);
    ExactVec3 centreNumerator;
    for (std::size_t axis = 0; axis < centreNumerator.size(); ++axis) {
        centreNumerator[axis] = addExact(
            multiplyExact(uu, vCrossN[axis]),
            multiplyExact(vv, nCrossU[axis]));
    }
    const ExactDyadic value = subtractExact(
        multiplyExact(nn, ww), dotExact(w, centreNumerator));
    return {PowerPredicateStatus::Valid, exactSign(value)};
}

// ---------------------------------------------------------------------------
Sign orient2d(const Vec2& a, const Vec2& b, const Vec2& c) {
    const double l = (b.u - a.u) * (c.v - a.v);
    const double r = (b.v - a.v) * (c.u - a.u);
    const double det = l - r;
    const double perm = std::abs(l) + std::abs(r);
    if (std::abs(det) > kO2Bound * perm) {
        bump(gO2F);
        return det > 0 ? Sign::Positive : Sign::Negative;
    }
    bump(gO2E);
    double bu[2], bv[2], cu[2], cv[2];
    twoDiff(b.u, a.u, bu[1], bu[0]);
    twoDiff(b.v, a.v, bv[1], bv[0]);
    twoDiff(c.u, a.u, cu[1], cu[0]);
    twoDiff(c.v, a.v, cv[1], cv[0]);
    double h[16];
    const int hl = det2exp(bu, cv, bv, cu, h);
    return signOf(hl, h);
}

// ---------------------------------------------------------------------------
static Sign orient3dImpl(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d, bool countStats) {
    const double adx = a.x - d.x, ady = a.y - d.y, adz = a.z - d.z;
    const double bdx = b.x - d.x, bdy = b.y - d.y, bdz = b.z - d.z;
    const double cdx = c.x - d.x, cdy = c.y - d.y, cdz = c.z - d.z;

    const double bdxcdy = bdx * cdy, cdxbdy = cdx * bdy;
    const double cdxady = cdx * ady, adxcdy = adx * cdy;
    const double adxbdy = adx * bdy, bdxady = bdx * ady;

    // Rows (v - d) give the NEGATIVE of the declared triple-product
    // convention; negate so fast and exact paths agree.
    const double det = -(adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) +
                         cdz * (adxbdy - bdxady));
    const double perm = std::abs(adz) * (std::abs(bdxcdy) + std::abs(cdxbdy)) +
                        std::abs(bdz) * (std::abs(cdxady) + std::abs(adxcdy)) +
                        std::abs(cdz) * (std::abs(adxbdy) + std::abs(bdxady));
    if (std::abs(det) > kO3Bound * perm) {
        if (countStats) bump(gO3F);
        return det > 0 ? Sign::Positive : Sign::Negative;
    }
    if (countStats) bump(gO3E);
    // Exact: rows are (v - a) so the convention matches
    // (b-a) x (c-a) . (d-a); det3exp computes r0.(r1 x r2) via cofactors of
    // the first row, i.e. det [r0; r1; r2]. det[b-a; c-a; d-a] equals the
    // triple product convention above.
    double r0[6], r1[6], r2[6];
    twoDiff(b.x, a.x, r0[1], r0[0]);
    twoDiff(b.y, a.y, r0[3], r0[2]);
    twoDiff(b.z, a.z, r0[5], r0[4]);
    twoDiff(c.x, a.x, r1[1], r1[0]);
    twoDiff(c.y, a.y, r1[3], r1[2]);
    twoDiff(c.z, a.z, r1[5], r1[4]);
    twoDiff(d.x, a.x, r2[1], r2[0]);
    twoDiff(d.y, a.y, r2[3], r2[2]);
    twoDiff(d.z, a.z, r2[5], r2[4]);
    static thread_local Buf192 out, sa, sb;
    const int hl = det3exp(r0, r1, r2, out.d, sa.d, sb.d);
    return signOf(hl, out.d);
}

Sign orient3d(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    return orient3dImpl(a, b, c, d, true);
}
Sign orient3dFiltered(const Vec3& a, const Vec3& b, const Vec3& c,
                      const Vec3& d) {
    const Sign sign = orient3dImpl(a, b, c, d, false);
    if (sign == Sign::Zero)
        gDelaunayZeroSeen.store(true, std::memory_order_relaxed);
    return sign;
}

double orient3dDeterminantEstimate(const Vec3& a, const Vec3& b,
                                  const Vec3& c, const Vec3& d) noexcept {
    for (const auto p : std::array{a, b, c, d})
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
            return std::numeric_limits<double>::quiet_NaN();
    double r0[6], r1[6], r2[6];
    twoDiff(b.x, a.x, r0[1], r0[0]);
    twoDiff(b.y, a.y, r0[3], r0[2]);
    twoDiff(b.z, a.z, r0[5], r0[4]);
    twoDiff(c.x, a.x, r1[1], r1[0]);
    twoDiff(c.y, a.y, r1[3], r1[2]);
    twoDiff(c.z, a.z, r1[5], r1[4]);
    twoDiff(d.x, a.x, r2[1], r2[0]);
    twoDiff(d.y, a.y, r2[3], r2[2]);
    twoDiff(d.z, a.z, r2[5], r2[4]);
    static thread_local Buf192 out, sa, sb;
    const int count = det3exp(r0, r1, r2, out.d, sa.d, sb.d);
    double estimate = 0;
    for (int i = 0; i < count; ++i) estimate += out.d[i];
    return estimate;
}

// ---------------------------------------------------------------------------
Sign incircle(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& d) {
    const double adx = a.u - d.u, ady = a.v - d.v;
    const double bdx = b.u - d.u, bdy = b.v - d.v;
    const double cdx = c.u - d.u, cdy = c.v - d.v;
    const double al = adx * adx + ady * ady;
    const double bl = bdx * bdx + bdy * bdy;
    const double cl = cdx * cdx + cdy * cdy;
    const double m0 = bdx * cdy - cdx * bdy;
    const double m1 = adx * cdy - cdx * ady;
    const double m2 = adx * bdy - bdx * ady;
    const double det = al * m0 - bl * m1 + cl * m2;
    const double perm = al * (std::abs(bdx * cdy) + std::abs(cdx * bdy)) +
                        bl * (std::abs(adx * cdy) + std::abs(cdx * ady)) +
                        cl * (std::abs(adx * bdy) + std::abs(bdx * ady));
    if (std::abs(det) > kICBound * perm) {
        bump(gICF);
        return det > 0 ? Sign::Positive : Sign::Negative;
    }
    bump(gICE);
    // Exact 3x3 with lifted first column, expanded along it.
    double ax[2], ay[2], bx[2], by[2], cx[2], cy[2];
    twoDiff(a.u, d.u, ax[1], ax[0]);
    twoDiff(a.v, d.v, ay[1], ay[0]);
    twoDiff(b.u, d.u, bx[1], bx[0]);
    twoDiff(b.v, d.v, by[1], by[0]);
    twoDiff(c.u, d.u, cx[1], cx[0]);
    twoDiff(c.v, d.v, cy[1], cy[0]);

    auto lift = [&](const double* x, const double* y, double* L,
                    double* s1, double* s2) {
        double xx[8], yy[8], scr[24];
        const int lx = mul22(x, x, xx, s1, s2);
        const int ly = mul22(y, y, yy, s1, s2);
        return sumZE(lx, xx, ly, yy, L, scr);
    };
    double s1[8], s2[8];
    double La[24], Lb[24], Lc[24];
    const int lla = lift(ax, ay, La, s1, s2);
    const int llb = lift(bx, by, Lb, s1, s2);
    const int llc = lift(cx, cy, Lc, s1, s2);

    double M0[16], M1[16], M2[16];
    const int lm0 = det2exp(bx, cy, cx, by, M0);
    const int lm1 = det2exp(ax, cy, cx, ay, M1);
    const int lm2 = det2exp(ax, by, bx, ay, M2);

    static thread_local std::vector<double> t0, t1, t2, acc, scr, icPiece,
        icTmp, icScr;
    auto mulLM = [&](int ll, const double* L, int lm, const double* M,
                     std::vector<double>& out) {
        out.assign(size_t(2 * ll * lm) + 8, 0.0);
        icPiece.assign(size_t(2 * lm) + 4, 0.0);
        icTmp.assign(out.size(), 0.0);
        icScr.assign(out.size(), 0.0);
        auto& piece = icPiece; auto& tmp = icTmp; auto& scr2 = icScr;
        int len = 0;
        for (int i = 0; i < ll; ++i) {
            const int lp = scaleZE(lm, M, L[i], piece.data());
            len = sumZE(len, out.data(), lp, piece.data(), tmp.data(),
                        scr2.data());
            std::memcpy(out.data(), tmp.data(), size_t(len) * 8);
        }
        return len;
    };
    const int lt0 = mulLM(lla, La, lm0, M0, t0);
    int lt1 = mulLM(llb, Lb, lm1, M1, t1);
    for (int i = 0; i < lt1; ++i) t1[size_t(i)] = -t1[size_t(i)];
    const int lt2 = mulLM(llc, Lc, lm2, M2, t2);

    acc.assign(size_t(lt0 + lt1 + lt2) + 8, 0.0);
    scr.assign(acc.size(), 0.0);
    int la = sumZE(lt0, t0.data(), lt1, t1.data(), acc.data(), scr.data());
    std::vector<double> fin(size_t(la + lt2) + 8);
    la = sumZE(la, acc.data(), lt2, t2.data(), fin.data(), scr.data());
    return signOf(la, fin.data());
}

// ---------------------------------------------------------------------------
static Sign insphereImpl(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d,
                         const Vec3& e, bool countStats) {
    const double aex = a.x - e.x, aey = a.y - e.y, aez = a.z - e.z;
    const double bex = b.x - e.x, bey = b.y - e.y, bez = b.z - e.z;
    const double cex = c.x - e.x, cey = c.y - e.y, cez = c.z - e.z;
    const double dex = d.x - e.x, dey = d.y - e.y, dez = d.z - e.z;

    auto det3 = [](double x0, double y0, double z0, double x1, double y1,
                   double z1, double x2, double y2, double z2) {
        return x0 * (y1 * z2 - z1 * y2) - y0 * (x1 * z2 - z1 * x2) +
               z0 * (x1 * y2 - y1 * x2);
    };
    auto perm3 = [](double x0, double y0, double z0, double x1, double y1,
                    double z1, double x2, double y2, double z2) {
        return std::abs(x0) * (std::abs(y1 * z2) + std::abs(z1 * y2)) +
               std::abs(y0) * (std::abs(x1 * z2) + std::abs(z1 * x2)) +
               std::abs(z0) * (std::abs(x1 * y2) + std::abs(y1 * x2));
    };
    const double La = aex * aex + aey * aey + aez * aez;
    const double Lb = bex * bex + bey * bey + bez * bez;
    const double Lc = cex * cex + cey * cey + cez * cez;
    const double Ld = dex * dex + dey * dey + dez * dez;

    const double Ma = det3(bex, bey, bez, cex, cey, cez, dex, dey, dez);
    const double Mb = det3(aex, aey, aez, cex, cey, cez, dex, dey, dez);
    const double Mc = det3(aex, aey, aez, bex, bey, bez, dex, dey, dez);
    const double Md = det3(aex, aey, aez, bex, bey, bez, cex, cey, cez);

    const double det = -La * Ma + Lb * Mb - Lc * Mc + Ld * Md;
    const double perm =
        La * perm3(bex, bey, bez, cex, cey, cez, dex, dey, dez) +
        Lb * perm3(aex, aey, aez, cex, cey, cez, dex, dey, dez) +
        Lc * perm3(aex, aey, aez, bex, bey, bez, dex, dey, dez) +
        Ld * perm3(aex, aey, aez, bex, bey, bez, cex, cey, cez);
    if (std::abs(det) > kISBound * perm) {
        if (countStats) bump(gISF);
        // For orient3d(a,b,c,d) > 0 (triple-product convention) this lifted
        // determinant is NEGATIVE when e is inside (locked by unit test).
        return det < 0 ? Sign::Positive : Sign::Negative;
    }
    if (countStats) bump(gISE);

    // Exact: rows (v - e) as 2-expansions.
    double ra[6], rb[6], rc[6], rd[6];
    auto row = [](const Vec3& v, const Vec3& e0, double* r) {
        twoDiff(v.x, e0.x, r[1], r[0]);
        twoDiff(v.y, e0.y, r[3], r[2]);
        twoDiff(v.z, e0.z, r[5], r[4]);
    };
    row(a, e, ra);
    row(b, e, rb);
    row(c, e, rc);
    row(d, e, rd);

    auto lift3 = [](const double* r, double* L) {
        double s1[8], s2[8], xx[8], yy[8], zz[8], t[16], scr[32];
        const int lx = mul22(r + 0, r + 0, xx, s1, s2);
        const int ly = mul22(r + 2, r + 2, yy, s1, s2);
        const int lz = mul22(r + 4, r + 4, zz, s1, s2);
        const int lt = sumZE(lx, xx, ly, yy, t, scr);
        return sumZE(lt, t, lz, zz, L, scr);
    };
    double LA[24], LB[24], LC[24], LD[24];
    const int la = lift3(ra, LA), lb = lift3(rb, LB), lc = lift3(rc, LC),
              ld = lift3(rd, LD);

    static thread_local Buf192 mA, mB, mC, mD, sa, sb;
    const int lma = det3exp(rb, rc, rd, mA.d, sa.d, sb.d);
    const int lmb = det3exp(ra, rc, rd, mB.d, sa.d, sb.d);
    const int lmc = det3exp(ra, rb, rd, mC.d, sa.d, sb.d);
    const int lmd = det3exp(ra, rb, rc, mD.d, sa.d, sb.d);

    static thread_local std::vector<double> term, accum, tmp, scr, piece,
        t2, s2;
    auto addTerm = [&](int lL, const double* L, int lM, const double* M,
                       double sign, int accLen) {
        // term = sign * L * M, accumulated into accum.
        term.assign(size_t(2 * lL * lM) + 8, 0.0);
        piece.assign(size_t(2 * lM) + 4, 0.0);
        t2.assign(term.size(), 0.0);
        s2.assign(term.size(), 0.0);
        int len = 0;
        for (int i = 0; i < lL; ++i) {
            const int lp = scaleZE(lM, M, sign * L[i], piece.data());
            len = sumZE(len, term.data(), lp, piece.data(), t2.data(),
                        s2.data());
            std::memcpy(term.data(), t2.data(), size_t(len) * 8);
        }
        tmp.assign(size_t(accLen + len) + 8, 0.0);
        scr.assign(tmp.size(), 0.0);
        const int newLen = sumZE(accLen, accum.data(), len, term.data(),
                                 tmp.data(), scr.data());
        accum.assign(tmp.begin(), tmp.begin() + newLen);
        return newLen;
    };
    accum.assign(1, 0.0);
    int accLen = 1;
    accLen = addTerm(la, LA, lma, mA.d, -1.0, accLen);
    accLen = addTerm(lb, LB, lmb, mB.d, +1.0, accLen);
    accLen = addTerm(lc, LC, lmc, mC.d, -1.0, accLen);
    accLen = addTerm(ld, LD, lmd, mD.d, +1.0, accLen);

    const Sign s = signOf(accLen, accum.data());
    if (s == Sign::Positive) return Sign::Negative;   // same flip as fast path
    if (s == Sign::Negative) return Sign::Positive;
    return Sign::Zero;
}

Sign insphere(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d,
              const Vec3& e) {
    return insphereImpl(a, b, c, d, e, true);
}
Sign insphereFiltered(const Vec3& a, const Vec3& b, const Vec3& c,
                      const Vec3& d, const Vec3& e) {
    const Sign sign = insphereImpl(a, b, c, d, e, false);
    if (sign == Sign::Zero)
        gDelaunayZeroSeen.store(true, std::memory_order_relaxed);
    return sign;
}

namespace detail {
void resetDelaunayZeroSeen() noexcept {
    gDelaunayZeroSeen.store(false, std::memory_order_relaxed);
}
bool delaunayZeroSeen() noexcept {
    return gDelaunayZeroSeen.load(std::memory_order_relaxed);
}
} // namespace detail

// ---------------------------------------------------------------------------
// Exact AFFINE-point predicates (BND-2A / T1-A support).
//
// An AffinePoint2 is  lo + t*(hi-lo),  t = numerator/2^level. orient2d is
// AFFINE (degree <= 1) in each of its three point arguments (translation
// invariant; the self-product terms cancel), so it is MULTIAFFINE in the
// three points. Multiaffine interpolation gives the exact identity
//
//   2^(lA+lB+lC) * orient2d(A,B,C)
//     = sum over corner choices a in {Alo,Ahi}, b in {Blo,Bhi}, c in {Clo,Chi}
//         WA(a) * WB(b) * WC(c) * orient2d(a,b,c)
//
// with integer weights WX(lo)=2^lX-numX, WX(hi)=numX (each in [0,2^52], so
// EXACT as doubles). 2^(lA+lB+lC) > 0, so the sign of orient2d(A,B,C) equals
// the sign of that integer-weighted sum of stored-point orient2d
// determinants - computed here with exact Shewchuk scale/sum expansions
// (never the rounded mapping coordinate).
namespace {

// Validate an AffinePoint2 BEFORE any field enters ldexp / weight / expansion
// arithmetic. Valid representations are EXACTLY: (stored) level==0,
// numerator==0, lo==hi, finite; or (genuine affine) 1<=level<=52,
// 0<numerator<2^level, ODD numerator (canonical reduced form), finite
// endpoints. Anything else fails loudly rather than being misread as an
// exact point.
inline void requireFinite(const Vec2& p, const char* what) {
    if (!std::isfinite(p.u) || !std::isfinite(p.v))
        throw std::invalid_argument(what);
}
void validateAffine(const AffinePoint2& p) {
    requireFinite(p.lo, "AffinePoint2: non-finite lo coordinate");
    if (p.level == 0) {
        if (p.numerator != 0 || p.lo.u != p.hi.u || p.lo.v != p.hi.v)
            throw std::invalid_argument(
                "AffinePoint2: stored form requires level==0, numerator==0, "
                "lo==hi");
        return;
    }
    if (p.level > kMaxAffineLevel)
        throw std::invalid_argument(
            "AffinePoint2: level exceeds kMaxAffineLevel (52)");
    requireFinite(p.hi, "AffinePoint2: non-finite hi coordinate");
    const std::uint64_t s = std::uint64_t(1) << p.level;   // safe: level<=52
    if (p.numerator == 0 || p.numerator >= s)
        throw std::invalid_argument(
            "AffinePoint2: numerator must satisfy 0 < numerator < 2^level");
    if ((p.numerator & 1ULL) == 0ULL)
        throw std::invalid_argument(
            "AffinePoint2: numerator must be odd (canonical reduced form)");
}

// Exact orient2d expansion of three STORED points -> h (needs >=16), len out.
int orient2dStoredExp(const Vec2& a, const Vec2& b, const Vec2& c, double* h) {
    double bu[2], bv[2], cu[2], cv[2];
    twoDiff(b.u, a.u, bu[1], bu[0]);
    twoDiff(b.v, a.v, bv[1], bv[0]);
    twoDiff(c.u, a.u, cu[1], cu[0]);
    twoDiff(c.v, a.v, cv[1], cv[0]);
    return det2exp(bu, cv, bv, cu, h);
}

// Active stored corners of an AffinePoint2 with their integer weights
// (weight = 2^level*cornerWeight). Stored point (num==0,level==0) -> one
// corner (lo) weight 1. Returns corner count (1 or 2).
int affineCorners(const AffinePoint2& X, Vec2 pts[2], double w[2]) {
    const double s = std::ldexp(1.0, X.level);       // 2^level exact (level<=52)
    const double nhi = static_cast<double>(X.numerator);
    const double nlo = s - nhi;                       // exact integer
    int k = 0;
    if (nlo != 0.0) { pts[k] = X.lo; w[k] = nlo; ++k; }
    if (nhi != 0.0) { pts[k] = X.hi; w[k] = nhi; ++k; }
    if (k == 0) { pts[0] = X.lo; w[0] = 0.0; k = 1; }  // guard (t==0 & s==0)
    return k;
}

// value(P) on axis d, scaled by 2^level so it is an exact expansion:
//   2^level * (P.lo[d]) + numerator * (P.hi[d] - P.lo[d]).
int axisValueScaled(const AffinePoint2& P, int axis, double* h) {
    const double lo = axis == 0 ? P.lo.u : P.lo.v;
    const double hi = axis == 0 ? P.hi.u : P.hi.v;
    const double s = std::ldexp(1.0, P.level);        // exact power of two
    const double t1 = s * lo;                          // exact (power-of-two)
    double diff[2];
    twoDiff(hi, lo, diff[1], diff[0]);
    double scaled[8];
    const int ls = scaleZE(2, diff, static_cast<double>(P.numerator), scaled);
    return growZE(ls, scaled, t1, h);                  // h needs >= ls+1
}

inline Sign flipSign(Sign s) {
    return s == Sign::Positive   ? Sign::Negative
           : s == Sign::Negative ? Sign::Positive
                                 : Sign::Zero;
}

}  // namespace

Sign orient2dAffine(const AffinePoint2& A, const AffinePoint2& B,
                    const AffinePoint2& C) {
    validateAffine(A);
    validateAffine(B);
    validateAffine(C);
    Vec2 ap[2], bp[2], cp[2];
    double aw[2], bw[2], cw[2];
    const int na = affineCorners(A, ap, aw);
    const int nb = affineCorners(B, bp, bw);
    const int nc = affineCorners(C, cp, cw);
    static thread_local std::vector<double> acc, nxt, scr;
    acc.assign(1, 0.0);
    int accLen = 1;
    for (int i = 0; i < na; ++i)
        for (int j = 0; j < nb; ++j)
            for (int k = 0; k < nc; ++k) {
                double h[16];
                const int hl = orient2dStoredExp(ap[i], bp[j], cp[k], h);
                double s1[64], s2[192], s3[512];
                const int l1 = scaleZE(hl, h, aw[i], s1);
                const int l2 = scaleZE(l1, s1, bw[j], s2);
                const int l3 = scaleZE(l2, s2, cw[k], s3);
                const std::size_t need = static_cast<std::size_t>(accLen + l3);
                if (nxt.size() < need) nxt.resize(need);
                if (scr.size() < need) scr.resize(need);
                accLen = sumZE(accLen, acc.data(), l3, s3, nxt.data(),
                               scr.data());
                acc.swap(nxt);
                if (acc.size() < static_cast<std::size_t>(accLen))
                    acc.resize(static_cast<std::size_t>(accLen));
            }
    bump(gO2E);
    return signOf(accLen, acc.data());
}

Sign orderOnSegmentAffine(const Vec2& loSeg, const Vec2& hiSeg,
                          const AffinePoint2& p, const AffinePoint2& q) {
    requireFinite(loSeg, "orderOnSegmentAffine: non-finite segment endpoint");
    requireFinite(hiSeg, "orderOnSegmentAffine: non-finite segment endpoint");
    if (loSeg.u == hiSeg.u && loSeg.v == hiSeg.v)
        throw std::invalid_argument(
            "orderOnSegmentAffine: zero-length segment");
    validateAffine(p);
    validateAffine(q);
    const double du = hiSeg.u - loSeg.u, dv = hiSeg.v - loSeg.v;
    const int axis = std::abs(du) >= std::abs(dv) ? 0 : 1;
    double hp[8], hq[8];
    const int lp = axisValueScaled(p, axis, hp);   // = value(p) * 2^p.level
    const int lq = axisValueScaled(q, axis, hq);   // = value(q) * 2^q.level
    // sign(value_p - value_q) = sign(hp*2^q.level - hq*2^p.level)
    double a[24], b[24];
    const int la = scaleZE(lp, hp, std::ldexp(1.0, q.level), a);
    int lb = scaleZE(lq, hq, std::ldexp(1.0, p.level), b);
    for (int i = 0; i < lb; ++i) b[i] = -b[i];
    double sum[64], scratch[64];
    const int ls = sumZE(la, a, lb, b, sum, scratch);
    Sign s = signOf(ls, sum);
    // fold in the segment direction on the dominant axis so the result is
    // sign(t_p - t_q) measured lo->hi, not the raw coordinate order.
    const double dir = axis == 0 ? du : dv;
    if (dir < 0.0) s = flipSign(s);
    return s;
}

Sign compareOnAxisAffine(const AffinePoint2& p, const AffinePoint2& q,
                         int axis) {
    if (axis != 0 && axis != 1)
        throw std::invalid_argument("compareOnAxisAffine: axis must be 0 or 1");
    validateAffine(p);
    validateAffine(q);
    double hp[8], hq[8];
    const int lp = axisValueScaled(p, axis, hp);   // = value(p) * 2^p.level
    const int lq = axisValueScaled(q, axis, hq);   // = value(q) * 2^q.level
    // sign(value_p - value_q) = sign(hp*2^q.level - hq*2^p.level)
    double a[24], b[24];
    const int la = scaleZE(lp, hp, std::ldexp(1.0, q.level), a);
    int lb = scaleZE(lq, hq, std::ldexp(1.0, p.level), b);
    for (int i = 0; i < lb; ++i) b[i] = -b[i];
    double sum[64], scratch[64];
    const int ls = sumZE(la, a, lb, b, sum, scratch);
    return signOf(ls, sum);
}

// ---------------------------------------------------------------------------
Stats stats() {
    Stats s;
    s.orient2dFast = gO2F.load();
    s.orient2dExact = gO2E.load();
    s.incircleFast = gICF.load();
    s.incircleExact = gICE.load();
    s.orient3dFast = gO3F.load();
    s.orient3dExact = gO3E.load();
    s.insphereFast = gISF.load();
    s.insphereExact = gISE.load();
    return s;
}
void resetStats() {
    for (auto* c : {&gO2F, &gO2E, &gICF, &gICE, &gO3F, &gO3E, &gISF, &gISE})
        c->store(0);
}

} // namespace paacs::predicates
