// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Tet10JacobianCertificate.hpp"

#include "paacs/predicates/Predicates.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>

static_assert(std::numeric_limits<double>::is_iec559,
              "Tet10 certification requires IEC-559 binary64");

namespace paacs::algo {
namespace {

using I = Tet10Interval;
constexpr double inf = std::numeric_limits<double>::infinity();

double down(double x) noexcept {
    return std::nextafter(x, -inf);
}
double up(double x) noexcept {
    return std::nextafter(x, inf);
}
I exact(double x) noexcept {
    return {x, x};
}
I add(I a, I b) noexcept {
    return {down(a.lo + b.lo), up(a.hi + b.hi)};
}
I neg(I a) noexcept {
    return {-a.hi, -a.lo};
}
I sub(I a, I b) noexcept {
    return add(a, neg(b));
}
I mul(I a, I b) noexcept {
    const std::array<double, 4> v{a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi};
    return {down(*std::min_element(v.begin(), v.end())), up(*std::max_element(v.begin(), v.end()))};
}
I scale(I a, double s) noexcept {
    return s < 0.0 ? mul(a, {s, s}) : mul(a, {s, s});
}
bool finite(I a) noexcept {
    return std::isfinite(a.lo) && std::isfinite(a.hi) && a.lo <= a.hi;
}

struct Exp {
    unsigned r, s, t;
};
constexpr std::array<Exp, 20> exps{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {2, 0, 0},
                                    {1, 1, 0}, {1, 0, 1}, {0, 2, 0}, {0, 1, 1}, {0, 0, 2},
                                    {3, 0, 0}, {2, 1, 0}, {2, 0, 1}, {1, 2, 0}, {1, 1, 1},
                                    {1, 0, 2}, {0, 3, 0}, {0, 2, 1}, {0, 1, 2}, {0, 0, 3}}};

int indexOf(unsigned r, unsigned s, unsigned t) noexcept {
    for (std::size_t i = 0; i < exps.size(); ++i)
        if (exps[i].r == r && exps[i].s == s && exps[i].t == t)
            return static_cast<int>(i);
    return -1;
}

struct Poly {
    std::array<I, 20> c{};
};
Poly constant(double x) noexcept {
    Poly p;
    p.c[0] = exact(x);
    return p;
}
Poly variable(unsigned axis) noexcept {
    Poly p;
    p.c[axis == 0 ? 1 : (axis == 1 ? 2 : 3)] = exact(1.0);
    return p;
}
Poly padd(const Poly &a, const Poly &b) noexcept {
    Poly p;
    for (std::size_t i = 0; i < 20; ++i)
        p.c[i] = add(a.c[i], b.c[i]);
    return p;
}
Poly pscale(const Poly &a, double s) noexcept {
    Poly p;
    for (std::size_t i = 0; i < 20; ++i)
        p.c[i] = scale(a.c[i], s);
    return p;
}
Poly pmul(const Poly &a, const Poly &b) noexcept {
    Poly p;
    for (std::size_t i = 0; i < 20; ++i)
        for (std::size_t j = 0; j < 20; ++j) {
            const unsigned r = exps[i].r + exps[j].r;
            const unsigned s = exps[i].s + exps[j].s;
            const unsigned t = exps[i].t + exps[j].t;
            if (r + s + t > 3)
                continue;
            const int k = indexOf(r, s, t);
            p.c[static_cast<std::size_t>(k)] =
                add(p.c[static_cast<std::size_t>(k)], mul(a.c[i], b.c[j]));
        }
    return p;
}
Poly deriv(const Poly &a, unsigned axis) noexcept {
    Poly p;
    for (std::size_t i = 0; i < 20; ++i) {
        const unsigned n = axis == 0 ? exps[i].r : (axis == 1 ? exps[i].s : exps[i].t);
        if (n == 0)
            continue;
        unsigned r = exps[i].r, s = exps[i].s, t = exps[i].t;
        if (axis == 0)
            --r;
        else if (axis == 1)
            --s;
        else
            --t;
        const int k = indexOf(r, s, t);
        p.c[static_cast<std::size_t>(k)] =
            add(p.c[static_cast<std::size_t>(k)], scale(a.c[i], double(n)));
    }
    return p;
}

unsigned choose(unsigned n, unsigned k) noexcept {
    if (k > n)
        return 0;
    if (k == 0 || k == n)
        return 1;
    if (n == 2)
        return 2;
    if (n == 3)
        return k == 1 || k == 2 ? 3 : 1;
    return 1;
}
unsigned factorial(unsigned n) noexcept {
    return n < 2 ? 1 : (n == 2 ? 2 : 6);
}
unsigned tetraDenominator(const Exp &e) noexcept {
    const unsigned degree = e.r + e.s + e.t;
    return 6 / (factorial(e.r) * factorial(e.s) * factorial(e.t) * factorial(3 - degree));
}

constexpr std::array<std::array<unsigned, 4>, 20> alpha{
    {{{3, 0, 0, 0}}, {{2, 1, 0, 0}}, {{2, 0, 1, 0}}, {{2, 0, 0, 1}}, {{1, 2, 0, 0}},
     {{1, 1, 1, 0}}, {{1, 1, 0, 1}}, {{1, 0, 2, 0}}, {{1, 0, 1, 1}}, {{1, 0, 0, 2}},
     {{0, 3, 0, 0}}, {{0, 2, 1, 0}}, {{0, 2, 0, 1}}, {{0, 1, 2, 0}}, {{0, 1, 1, 1}},
     {{0, 1, 0, 2}}, {{0, 0, 3, 0}}, {{0, 0, 2, 1}}, {{0, 0, 1, 2}}, {{0, 0, 0, 3}}}};

std::array<I, 20> powerToBernstein(const Poly &p) noexcept {
    std::array<I, 20> out{};
    for (std::size_t a = 0; a < 20; ++a) {
        I sum{};
        for (std::size_t bi = 0; bi < 20; ++bi) {
            const auto e = exps[bi];
            if (e.r > alpha[a][1] || e.s > alpha[a][2] || e.t > alpha[a][3])
                continue;
            const unsigned num =
                choose(alpha[a][1], e.r) * choose(alpha[a][2], e.s) * choose(alpha[a][3], e.t);
            const unsigned den = tetraDenominator(e);
            sum = add(sum, scale(p.c[bi], double(num) / double(den)));
        }
        out[a] = sum;
    }
    return out;
}

int alphaIndex(const std::array<unsigned, 4> &a, unsigned degree) noexcept {
    if (degree == 3) {
        for (std::size_t i = 0; i < 20; ++i)
            if (alpha[i] == a)
                return static_cast<int>(i);
    }
    return -1;
}

using Bary = std::array<I, 4>;
constexpr std::array<std::array<double, 4>, 10> baryVertices{{{{1, 0, 0, 0}},
                                                              {{0, 1, 0, 0}},
                                                              {{0, 0, 1, 0}},
                                                              {{0, 0, 0, 1}},
                                                              {{.5, .5, 0, 0}},
                                                              {{0, .5, .5, 0}},
                                                              {{.5, 0, .5, 0}},
                                                              {{.5, 0, 0, .5}},
                                                              {{0, .5, 0, .5}},
                                                              {{0, 0, .5, .5}}}};
constexpr std::array<std::array<unsigned, 4>, 8> childVertices{{{{0, 4, 6, 7}},
                                                                {{1, 4, 8, 5}},
                                                                {{2, 6, 5, 9}},
                                                                {{3, 7, 9, 8}},
                                                                {{4, 9, 6, 7}},
                                                                {{4, 9, 7, 8}},
                                                                {{4, 9, 8, 5}},
                                                                {{4, 9, 5, 6}}}};

I blossom(const std::array<I, 20> &coeff, const std::array<Bary, 3> &args) noexcept {
    // Sparse fixed-size maps for degrees 3,2,1,0. This deliberately mirrors
    // the lexicographic rational oracle used by the focused tests.
    std::array<I, 20> cur = coeff;
    unsigned degree = 3;
    for (unsigned step = 0; step < 3; ++step) {
        std::array<I, 20> next{};
        std::size_t out = 0;
        for (unsigned a0 = degree - 1;; --a0) {
            for (unsigned a1 = degree - 1 - a0;; --a1) {
                for (unsigned a2 = degree - 1 - a0 - a1;; --a2) {
                    const unsigned a3 = degree - 1 - a0 - a1 - a2;
                    const std::array<unsigned, 4> beta{{a0, a1, a2, a3}};
                    I v{};
                    for (unsigned q = 0; q < 4; ++q) {
                        auto parent = beta;
                        ++parent[q];
                        int pi = -1;
                        if (degree == 3)
                            pi = alphaIndex(parent, 3);
                        else {
                            // Earlier reductions are stored in the same
                            // descending lexicographic order.
                            std::size_t idx = 0;
                            for (unsigned b0 = degree;; --b0) {
                                for (unsigned b1 = degree - b0;; --b1) {
                                    for (unsigned b2 = degree - b0 - b1;; --b2) {
                                        const unsigned b3 = degree - b0 - b1 - b2;
                                        if (parent == std::array<unsigned, 4>{{b0, b1, b2, b3}})
                                            pi = static_cast<int>(idx);
                                        ++idx;
                                        if (b2 == 0)
                                            break;
                                    }
                                    if (b1 == 0)
                                        break;
                                }
                                if (b0 == 0)
                                    break;
                            }
                        }
                        v = add(v, mul(args[step][q], cur[static_cast<std::size_t>(pi)]));
                    }
                    next[out++] = v;
                    if (a2 == 0)
                        break;
                }
                if (a1 == 0)
                    break;
            }
            if (a0 == 0)
                break;
        }
        cur = next;
        --degree;
    }
    return cur[0];
}

bool restrictChild(const std::array<I, 20> &parent, unsigned child,
                   std::array<I, 20> &result) noexcept {
    if (child >= 8)
        return false;
    for (std::size_t ai = 0; ai < 20; ++ai) {
        std::array<Bary, 3> args{};
        unsigned at = 0;
        for (unsigned v = 0; v < 4; ++v)
            for (unsigned n = 0; n < alpha[ai][v]; ++n) {
                for (unsigned k = 0; k < 4; ++k)
                    args[at][k] = exact(baryVertices[childVertices[child][v]][k]);
                ++at;
            }
        result[ai] = blossom(parent, args);
    }
    return true;
}

bool normalize(const std::array<Vec3, 10> &in, std::array<std::array<I, 3>, 10> &out,
               double &scaleOut) noexcept {
    Vec3 lo = in[0], hi = in[0];
    for (const auto &p : in) {
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
            return false;
        lo.x = std::min(lo.x, p.x);
        lo.y = std::min(lo.y, p.y);
        lo.z = std::min(lo.z, p.z);
        hi.x = std::max(hi.x, p.x);
        hi.y = std::max(hi.y, p.y);
        hi.z = std::max(hi.z, p.z);
    }
    const Vec3 c{std::midpoint(lo.x, hi.x), std::midpoint(lo.y, hi.y), std::midpoint(lo.z, hi.z)};
    I m{0, 0};
    for (const auto &p : in)
        for (unsigned axis = 0; axis < 3; ++axis) {
            const double x = axis == 0 ? p.x : (axis == 1 ? p.y : p.z);
            const double cc = axis == 0 ? c.x : (axis == 1 ? c.y : c.z);
            const I d = sub(exact(x), exact(cc));
            m.hi = std::max(m.hi, std::max(std::abs(d.lo), std::abs(d.hi)));
        }
    if (!(m.hi > 0.0) || !std::isfinite(m.hi))
        return false;
    const int k = std::ilogb(m.hi) + 1;
    // Record the finite normalized upper magnitude.  Never materialize
    // 2^-k: that value can overflow even when scalbn(value,-k) is valid.
    scaleOut = std::scalbn(m.hi, -k);
    if (!std::isfinite(scaleOut) || scaleOut == 0.0)
        return false;
    for (std::size_t n = 0; n < 10; ++n)
        for (unsigned axis = 0; axis < 3; ++axis) {
            const double x = axis == 0 ? in[n].x : (axis == 1 ? in[n].y : in[n].z);
            const double cc = axis == 0 ? c.x : (axis == 1 ? c.y : c.z);
            I d = sub(exact(x), exact(cc));
            I q{down(std::scalbn(d.lo, -k)), up(std::scalbn(d.hi, -k))};
            if (!finite(q))
                return false;
            out[n][axis] = q;
        }
    return true;
}

struct Pending {
    std::array<I, 20> det{}, shift{};
    std::uint64_t path = 0;
    std::uint32_t depth = 0;
};

} // namespace

Tet10CertificateResult certifyTet10Jacobian(const std::array<Vec3, 10> &nodes,
                                            const Tet10CertificatePolicy &policy) noexcept {
    Tet10CertificateResult out;
    try {
        if (!std::isfinite(policy.minimumJacobianRatio) || policy.minimumJacobianRatio < 0.0 ||
            policy.minimumJacobianRatio > 1.0 || policy.maxDepth > 20)
            return out;
        std::array<std::array<I, 3>, 10> p{};
        if (!normalize(nodes, p, out.normalizationScale)) {
            out.status = Tet10CertificateStatus::RangeUnrepresentable;
            return out;
        }
        // Evaluate the topology-affecting sign from the original binary64
        // corners.  The frozen predicate has a documented normal-product
        // envelope, so exact Zero at an extreme exponent receives one
        // uniform power-of-two retry.  This does not translate or round the
        // input values and cannot alter a nonzero determinant's sign.
        auto cornerSign = predicates::orient3d(nodes[0], nodes[1], nodes[2], nodes[3]);
        bool cornerScaleUnrepresentable = false;
        if (cornerSign == predicates::Sign::Zero) {
            double magnitude = 0.0;
            for (unsigned n = 0; n < 4; ++n) {
                magnitude = std::max(magnitude, std::abs(nodes[n].x));
                magnitude = std::max(magnitude, std::abs(nodes[n].y));
                magnitude = std::max(magnitude, std::abs(nodes[n].z));
            }
            if (magnitude > 0.0 && std::isfinite(magnitude)) {
                const int exponent = -std::ilogb(magnitude);
                std::array<Vec3, 4> scaled{};
                bool valid = true;
                for (unsigned n = 0; n < 4; ++n) {
                    scaled[n] = {std::scalbn(nodes[n].x, exponent),
                                 std::scalbn(nodes[n].y, exponent),
                                 std::scalbn(nodes[n].z, exponent)};
                    valid = valid && std::isfinite(scaled[n].x) && std::isfinite(scaled[n].y) &&
                            std::isfinite(scaled[n].z);
                    for (unsigned axis = 0; axis < 3; ++axis) {
                        const double original =
                            axis == 0 ? nodes[n].x : (axis == 1 ? nodes[n].y : nodes[n].z);
                        const double working =
                            axis == 0 ? scaled[n].x : (axis == 1 ? scaled[n].y : scaled[n].z);
                        // A nonzero normal result proves the power-of-two
                        // scaling did not round through underflow.
                        if (original != 0.0 && !std::isnormal(working))
                            valid = false;
                    }
                }
                for (unsigned a = 0; a < 4; ++a)
                    for (unsigned b = a + 1; b < 4; ++b)
                        for (unsigned axis = 0; axis < 3; ++axis) {
                            const double oa =
                                axis == 0 ? nodes[a].x : (axis == 1 ? nodes[a].y : nodes[a].z);
                            const double ob =
                                axis == 0 ? nodes[b].x : (axis == 1 ? nodes[b].y : nodes[b].z);
                            const double sa =
                                axis == 0 ? scaled[a].x : (axis == 1 ? scaled[a].y : scaled[a].z);
                            const double sb =
                                axis == 0 ? scaled[b].x : (axis == 1 ? scaled[b].y : scaled[b].z);
                            if (oa != ob && sa == sb)
                                valid = false;
                        }
                if (valid)
                    cornerSign = predicates::orient3d(scaled[0], scaled[1], scaled[2], scaled[3]);
                else
                    cornerScaleUnrepresentable = true;
            }
        }
        if (cornerSign != predicates::Sign::Positive) {
            if (cornerScaleUnrepresentable)
                out.status = Tet10CertificateStatus::RangeUnrepresentable;
            return out;
        }
        // One root-polynomial work unit is one fixed interval-array operation
        // (shape construction, differentiation, convolution, determinant, or
        // Bernstein conversion).  The evaluator executes this constant amount
        // before any coefficient can be classified.
        constexpr std::uint64_t rootPolynomialWork = 2048;
        auto consumeWork = [&](std::uint64_t units) noexcept {
            if (units > policy.maxWorkUnits || out.metrics.workUnits > policy.maxWorkUnits - units)
                return false;
            out.metrics.workUnits += units;
            return true;
        };
        if (!consumeWork(rootPolynomialWork)) {
            out.status = Tet10CertificateStatus::IncompleteWork;
            return out;
        }
        const Poly r = variable(0), s = variable(1), t = variable(2);
        Poly l0 = padd(constant(1.0), pscale(padd(padd(r, s), t), -1.0));
        const std::array<Poly, 4> l{{l0, r, s, t}};
        std::array<Poly, 10> shape{};
        for (unsigned i = 0; i < 4; ++i)
            shape[i] = padd(pscale(pmul(l[i], l[i]), 2.0), pscale(l[i], -1.0));
        constexpr std::array<std::array<unsigned, 2>, 6> edges{
            {{{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};
        for (unsigned e = 0; e < 6; ++e)
            shape[4 + e] = pscale(pmul(l[edges[e][0]], l[edges[e][1]]), 4.0);
        std::array<Poly, 3> x{};
        for (unsigned axis = 0; axis < 3; ++axis)
            for (unsigned n = 0; n < 10; ++n) {
                Poly q = shape[n];
                for (auto &value : q.c)
                    value = mul(value, p[n][axis]);
                x[axis] = padd(x[axis], q);
            }
        Poly j[3][3];
        for (unsigned axis = 0; axis < 3; ++axis)
            for (unsigned d = 0; d < 3; ++d)
                j[axis][d] = deriv(x[axis], d);
        Poly det = padd(
            padd(pmul(j[0][0], padd(pmul(j[1][1], j[2][2]), pscale(pmul(j[1][2], j[2][1]), -1.0))),
                 pscale(pmul(j[0][1],
                             padd(pmul(j[1][0], j[2][2]), pscale(pmul(j[1][2], j[2][0]), -1.0))),
                        -1.0)),
            pmul(j[0][2], padd(pmul(j[1][0], j[2][1]), pscale(pmul(j[1][1], j[2][0]), -1.0))));
        Poly corner = constant(0.0);
        const std::array<I, 3> a{
            {sub(p[1][0], p[0][0]), sub(p[1][1], p[0][1]), sub(p[1][2], p[0][2])}};
        const std::array<I, 3> b{
            {sub(p[2][0], p[0][0]), sub(p[2][1], p[0][1]), sub(p[2][2], p[0][2])}};
        const std::array<I, 3> c{
            {sub(p[3][0], p[0][0]), sub(p[3][1], p[0][1]), sub(p[3][2], p[0][2])}};
        I dc = add(sub(mul(a[0], sub(mul(b[1], c[2]), mul(b[2], c[1]))),
                       mul(a[1], sub(mul(b[0], c[2]), mul(b[2], c[0])))),
                   mul(a[2], sub(mul(b[0], c[1]), mul(b[1], c[0]))));
        if (dc.lo <= 0.0) {
            out.status = dc.hi <= 0.0 ? Tet10CertificateStatus::InvalidInput
                                      : Tet10CertificateStatus::RangeUnrepresentable;
            return out;
        }
        corner.c[0] = dc;
        Poly shifted = padd(det, pscale(corner, -policy.minimumJacobianRatio));
        out.rootDeterminant = powerToBernstein(det);
        out.rootShifted = powerToBernstein(shifted);
        // maxDepth is at most 20, so DFS high-water is 1+7*20=141. Fixed
        // storage makes certificate execution allocation-free after entry.
        std::array<Pending, 141> stack{};
        std::size_t used = 1;
        stack[0] = {out.rootDeterminant, out.rootShifted, 0, 0};
        while (used != 0) {
            Pending cur = stack[--used];
            if (out.metrics.visitedSubcells >= policy.maxVisitedSubcells) {
                out.status = Tet10CertificateStatus::IncompleteSubcells;
                return out;
            }
            ++out.metrics.visitedSubcells;
            out.metrics.maximumDepth = std::max(out.metrics.maximumDepth, cur.depth);
            bool positive = true, allShiftNonpositive = true, allDetNonpositive = true;
            for (unsigned i = 0; i < 20; ++i) {
                if (out.metrics.coefficientTests >= policy.maxCoefficientTests) {
                    out.status = Tet10CertificateStatus::IncompleteCoefficients;
                    return out;
                }
                if (!consumeWork(1)) {
                    out.status = Tet10CertificateStatus::IncompleteWork;
                    return out;
                }
                ++out.metrics.coefficientTests;
                out.conservativeMinimumLowerBound =
                    std::min(out.conservativeMinimumLowerBound, cur.shift[i].lo);
                positive = positive && cur.shift[i].lo > 0.0;
                allShiftNonpositive = allShiftNonpositive && cur.shift[i].hi <= 0.0;
                allDetNonpositive = allDetNonpositive && cur.det[i].hi <= 0.0;
                if ((allShiftNonpositive || allDetNonpositive) && i == 19) {
                    out.witness = {cur.path, cur.depth, i, cur.shift[i], cur.det[i]};
                }
            }
            if (positive)
                continue;
            if (allDetNonpositive) {
                out.status = Tet10CertificateStatus::Nonpositive;
                return out;
            }
            if (allShiftNonpositive) {
                out.status = Tet10CertificateStatus::PositiveBelowThreshold;
                return out;
            }
            if (cur.depth >= policy.maxDepth) {
                out.status = Tet10CertificateStatus::IncompleteDepth;
                return out;
            }
            for (int child = 7; child >= 0; --child) {
                // One child restriction produces 20 paired determinant/
                // shifted coefficients.  Charge it before touching the DFS
                // stack so a hard budget is never crossed.
                if (!consumeWork(20)) {
                    out.status = Tet10CertificateStatus::IncompleteWork;
                    return out;
                }
                Pending next;
                next.depth = cur.depth + 1;
                next.path = cur.path | (std::uint64_t(child) << (3u * cur.depth));
                restrictChild(cur.det, static_cast<unsigned>(child), next.det);
                restrictChild(cur.shift, static_cast<unsigned>(child), next.shift);
                if (used >= stack.size()) {
                    out.status = Tet10CertificateStatus::InvalidInput;
                    return out;
                }
                stack[used++] = next;
            }
        }
        out.status = Tet10CertificateStatus::Certified;
        out.complete = true;
        return out;
    } catch (const std::bad_alloc &) {
        out.status = Tet10CertificateStatus::AllocationFailed;
        return out;
    } catch (...) {
        out.status = Tet10CertificateStatus::InvalidInput;
        return out;
    }
}

#if defined(PAACS_TOPO_TESTING)
bool restrictTet10BernsteinChildForTesting(const std::array<Tet10Interval, 20> &parent,
                                           std::uint32_t child,
                                           std::array<Tet10Interval, 20> &result) noexcept {
    return restrictChild(parent, child, result);
}
#endif

} // namespace paacs::algo
