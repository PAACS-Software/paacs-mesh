// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/mesh/Quality.hpp"
#include "paacs/predicates/Predicates.hpp"

#include <algorithm>
#include "paacs/core/Parallel.hpp"
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <sstream>

namespace paacs::mesh {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kInvSqrt2 = 0.70710678118654752440;  // regular-tet corner J
} // namespace

double tetVolume(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    const Vec3 x=b-a, y=c-a, z=d-a;
    const double determinant=x.cross(y).dot(z);
    const double permanent=(std::abs(x.y*y.z)+std::abs(x.z*y.y))*std::abs(z.x)+
        (std::abs(x.z*y.x)+std::abs(x.x*y.z))*std::abs(z.y)+
        (std::abs(x.x*y.y)+std::abs(x.y*y.x))*std::abs(z.z);
    // A nearly flat cell can lose its sign when edge differences and products
    // cancel. Measure from the original doubles in that uncommon case.
    if (std::isfinite(permanent) &&
        std::abs(determinant)<=32*std::numeric_limits<double>::epsilon()*permanent)
        return predicates::orient3dDeterminantEstimate(a,b,c,d)/6.;
    return determinant/6.;
}

double tetScaledJacobian(const Vec3& a, const Vec3& b, const Vec3& c,
                         const Vec3& d) {
    const Vec3* v[4] = {&a, &b, &c, &d};
    // Corner k Jacobian: det of the three edges leaving corner k, ordered so
    // a positively-oriented tet gives positive determinants at every corner.
    static const int e[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
    double worst = 1e300;
    for (int k = 0; k < 4; ++k) {
        const Vec3 e0 = *v[e[k][0]] - *v[k];
        const Vec3 e1 = *v[e[k][1]] - *v[k];
        const Vec3 e2 = *v[e[k][2]] - *v[k];
        const double den = e0.norm() * e1.norm() * e2.norm();
        if (den < 1e-300) return 0.0;
        worst = std::min(worst, e0.cross(e1).dot(e2) / den);
    }
    return worst / kInvSqrt2;  // normalize: regular tet -> 1
}

double tetSkewness(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    const double V = tetVolume(a, b, c, d);
    if (V <= 0) return 1.0;

    // Circumradius R via |num| / (12 V) (same identity as tetQuality).
    const Vec3 ba = b - a, ca = c - a, da = d - a;
    const Vec3 num = ba.norm2() * ca.cross(da) +
                     ca.norm2() * da.cross(ba) +
                     da.norm2() * ba.cross(ca);
    const double R = num.norm() / (12.0 * V);
    if (R <= 0) return 1.0;

    // Regular tet inscribed in sphere of radius R: edge a_eq = R*sqrt(8/3),
    // V_eq = a_eq^3 / (6*sqrt(2)).
    const double aEq = R * std::sqrt(8.0 / 3.0);
    const double vEq = aEq * aEq * aEq / (6.0 * std::sqrt(2.0));
    return std::clamp((vEq - V) / vEq, 0.0, 1.0);
}

void tetDihedralAngles(const Vec3& a, const Vec3& b, const Vec3& c,
                       const Vec3& d, double& minDeg, double& maxDeg) {
    const Vec3* v[4] = {&a, &b, &c, &d};
    // Outward face normals for a positively-oriented tet; the face opposite
    // vertex k, with outward orientation.
    auto n = [&](int i, int j, int k) {
        return (*v[j] - *v[i]).cross(*v[k] - *v[i]).normalized();
    };
    const Vec3 N[4] = {n(1, 3, 2),   // opposite 0
                       n(0, 2, 3),   // opposite 1
                       n(0, 3, 1),   // opposite 2
                       n(0, 1, 2)};  // opposite 3
    minDeg = 360; maxDeg = 0;
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j) {
            // Dihedral along the shared edge = pi - angle(N_i, N_j).
            const double cosang = std::clamp(N[i].dot(N[j]), -1.0, 1.0);
            const double dih = 180.0 - std::acos(cosang) * 180.0 / kPi;
            minDeg = std::min(minDeg, dih);
            maxDeg = std::max(maxDeg, dih);
        }
}

double tetQuality(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    const double V = std::abs(tetVolume(a, b, c, d));
    if (V < 1e-300) return 0.0;

    auto faceArea = [](const Vec3& p, const Vec3& q, const Vec3& r) {
        return 0.5 * (q - p).cross(r - p).norm();
    };
    const double A = faceArea(a, b, c) + faceArea(a, b, d) +
                     faceArea(a, c, d) + faceArea(b, c, d);
    const double rIn = 3.0 * V / A;

    const Vec3 ba = b - a, ca = c - a, da = d - a;
    const Vec3 num = ba.norm2() * ca.cross(da) +
                     ca.norm2() * da.cross(ba) +
                     da.norm2() * ba.cross(ca);
    const double R = num.norm() / (12.0 * V);
    if (R <= 0) return 0.0;
    return std::clamp(3.0 * rIn / R, 0.0, 1.0);
}

QualityReport analyzeTets(const Mesh& m) {
    // Deterministic parallel reduce: fixed 512-tet chunks accumulated
    // sequentially, merged in chunk order => bit-identical on any thread
    // count. Metrics per element are independent, so this scales linearly.
    struct Acc {
        MetricStats rr, sj, sk, t10j;
        double dmin = 180, dmax = 0, vol = 0, t10vol = 0;
        Index nTets = 0, nInv = 0, nT10 = 0, nT10Inv = 0, nNF = 0;
        std::array<Index, 10> hist{};
    };
    auto accMetric = [](MetricStats& s, double x) {
        s.min = std::min(s.min, x);
        s.max = std::max(s.max, x);
        s.mean += x;                       // running sum; divided at the end
        ++s.n;
    };

    const Index nc = m.cellCount();
    Acc total = core::parallelReduce(
        std::size_t{0}, static_cast<std::size_t>(nc), 512, Acc{},
        [&](Acc& a, std::size_t i) {
            const auto& c = m.cell(static_cast<Index>(i));
            auto finite = [](const Vec3& p) {
                return std::isfinite(p.x) && std::isfinite(p.y) &&
                       std::isfinite(p.z);
            };
            if (c.type == CellType::Tet10) {
                const Index* v = m.cellNodes(static_cast<Index>(i));
                Vec3 p[10];
                bool ok = true;
                for (int k = 0; k < 10; ++k) {
                    p[k] = m.node(v[k]);
                    ok = ok && finite(p[k]);
                }
                if (!ok) { ++a.nNF; return; }
                const double jr = tet10MinJacobianRatio(p);
                accMetric(a.t10j, jr);
                if (jr <= 0.0) ++a.nT10Inv;
                a.t10vol += std::abs(tetVolume(p[0], p[1], p[2], p[3]));
                ++a.nT10;
                return;
            }
            if (c.type != CellType::Tet4) return;
            const Index* v = m.cellNodes(static_cast<Index>(i));
            const Vec3 &p0 = m.node(v[0]), &p1 = m.node(v[1]),
                       &p2 = m.node(v[2]), &p3 = m.node(v[3]);
            if (!finite(p0) || !finite(p1) || !finite(p2) || !finite(p3)) {
                ++a.nNF;
                return;
            }

            const double vol = tetVolume(p0, p1, p2, p3);
            if (predicates::orient3d(p0,p1,p2,p3) != predicates::Sign::Positive)
                ++a.nInv;
            a.vol += std::abs(vol);
            ++a.nTets;

            const double rr = tetQuality(p0, p1, p2, p3);
            accMetric(a.rr, rr);
            accMetric(a.sj, tetScaledJacobian(p0, p1, p2, p3));
            accMetric(a.sk, tetSkewness(p0, p1, p2, p3));

            double dmin, dmax;
            tetDihedralAngles(p0, p1, p2, p3, dmin, dmax);
            a.dmin = std::min(a.dmin, dmin);
            a.dmax = std::max(a.dmax, dmax);
            ++a.hist[static_cast<size_t>(
                std::clamp(static_cast<int>(rr * 10.0), 0, 9))];
        },
        [](Acc& out, const Acc& a) {
            auto mergeMetric = [](MetricStats& s, const MetricStats& t) {
                s.min = std::min(s.min, t.min);
                s.max = std::max(s.max, t.max);
                s.mean += t.mean;
                s.n += t.n;
            };
            mergeMetric(out.rr, a.rr);
            mergeMetric(out.sj, a.sj);
            mergeMetric(out.sk, a.sk);
            mergeMetric(out.t10j, a.t10j);
            out.dmin = std::min(out.dmin, a.dmin);
            out.dmax = std::max(out.dmax, a.dmax);
            out.vol += a.vol;
            out.t10vol += a.t10vol;
            out.nTets += a.nTets;
            out.nInv += a.nInv;
            out.nT10 += a.nT10;
            out.nT10Inv += a.nT10Inv;
            out.nNF += a.nNF;
            for (size_t k = 0; k < 10; ++k) out.hist[k] += a.hist[k];
        });

    QualityReport r;
    r.radiusRatio = total.rr;
    r.scaledJacobian = total.sj;
    r.skewness = total.sk;
    r.minDihedralDeg = total.dmin;
    r.maxDihedralDeg = total.dmax;
    r.totalVolume = total.vol;
    r.nTets = total.nTets;
    r.nInverted = total.nInv;
    r.radiusRatioHist = total.hist;
    r.nTet10 = total.nT10;
    r.nTet10Invalid = total.nT10Inv;
    r.tet10Jacobian = total.t10j;
    r.tet10CornerVolume = total.t10vol;
    r.nNonFinite = total.nNF;
    for (MetricStats* s : {&r.radiusRatio, &r.scaledJacobian, &r.skewness,
                           &r.tet10Jacobian})
        if (s->n) s->mean /= static_cast<double>(s->n);
    return r;
}

std::string QualityReport::toString() const {
    // Dynamic formatting (the previous fixed buffer could truncate and, via
    // size_t underflow of the remaining-space argument, crash - PM-QUAL-001).
    // Empty and Tet10-only populations are explicit states, never sentinel
    // numbers.
    std::ostringstream os;
    os.setf(std::ios::fixed);
    auto metric = [&](const char* name, const MetricStats& s, int prec = 3) {
        os << "  " << name;
        if (!s.n) { os << " (no elements)\n"; return; }
        os << std::setprecision(prec)
           << " min " << s.min << "  mean " << s.mean
           << "  max " << s.max << "\n";
    };
    if (nTets == 0 && nTet10 == 0) {
        if (nNonFinite > 0) {
            os << "mesh contains no analyzable tetrahedra ("
               << nNonFinite << " with non-finite coordinates)\n";
            return os.str();
        }
        return "mesh contains no tetrahedra (nothing analyzed)\n";
    }
    if (nNonFinite > 0)
        os << "WARNING: " << nNonFinite
           << " element(s) with non-finite coordinates excluded\n";

    if (nTets > 0) {
        os << std::setprecision(6)
           << "tet4 " << nTets << " (inverted " << nInverted
           << "), volume " << totalVolume << "\n";
        metric("radius-ratio   ", radiusRatio);
        metric("scaled Jacobian", scaledJacobian);
        metric("skewness       ", skewness);
        os << std::setprecision(1)
           << "  dihedral angles [" << minDihedralDeg << ", "
           << maxDihedralDeg << "] deg\n  rr hist |";
        for (Index h : radiusRatioHist) os << ' ' << h;
        os << " |\n";
    }
    if (nTet10 > 0) {
        os << std::setprecision(6)
           << "tet10 " << nTet10 << " (invalid " << nTet10Invalid
           << "), corner volume " << tet10CornerVolume << "\n";
        metric("sampled J ratio", tet10Jacobian);
    }
    return os.str();
}


QualityStats tetQualityStats(const Mesh& m) {
    auto r = analyzeTets(m);
    QualityStats s;
    s.min = r.radiusRatio.n ? r.radiusRatio.min : 0;
    s.max = r.radiusRatio.n ? r.radiusRatio.max : 0;
    s.mean = r.radiusRatio.mean;
    s.n = r.radiusRatio.n;
    return s;
}

double totalTetVolume(const Mesh& m) { return analyzeTets(m).totalVolume; }


double tet10MinJacobianRatio(const Vec3 p[10]) {
    // Barycentric gradients (constant): L0 = 1-x-y-z, L1 = x, L2 = y, L3 = z.
    static const double gL[4][3] = {
        {-1, -1, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    // Edge node -> its two corner indices (VTK/MSH tet10 ordering).
    static const int ed[6][2] = {{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};

    const double v6 = 6.0 * tetVolume(p[0], p[1], p[2], p[3]);
    if (std::abs(v6) < 1e-300) return 0.0;

    auto detAt = [&](double L0, double L1, double L2, double L3) {
        const double L[4] = {L0, L1, L2, L3};
        double J[3][3] = {};
        auto addGrad = [&](const Vec3& x, const double g[3]) {
            J[0][0] += x.x * g[0]; J[0][1] += x.x * g[1]; J[0][2] += x.x * g[2];
            J[1][0] += x.y * g[0]; J[1][1] += x.y * g[1]; J[1][2] += x.y * g[2];
            J[2][0] += x.z * g[0]; J[2][1] += x.z * g[1]; J[2][2] += x.z * g[2];
        };
        for (int i = 0; i < 4; ++i) {
            const double s = 4.0 * L[i] - 1.0;
            const double g[3] = {s * gL[i][0], s * gL[i][1], s * gL[i][2]};
            addGrad(p[i], g);
        }
        for (int e = 0; e < 6; ++e) {
            const int a = ed[e][0], b = ed[e][1];
            const double g[3] = {
                4.0 * (L[a] * gL[b][0] + L[b] * gL[a][0]),
                4.0 * (L[a] * gL[b][1] + L[b] * gL[a][1]),
                4.0 * (L[a] * gL[b][2] + L[b] * gL[a][2])};
            addGrad(p[4 + e], g);
        }
        return J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
               J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
               J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    };

    // 15 sample points: 4 corners, 6 edge midpoints, 4 face centroids,
    // 1 centroid.
    static const double S[15][4] = {
        {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1},
        {.5,.5,0,0}, {0,.5,.5,0}, {.5,0,.5,0},
        {.5,0,0,.5}, {0,.5,0,.5}, {0,0,.5,.5},
        {1./3,1./3,1./3,0}, {1./3,1./3,0,1./3},
        {1./3,0,1./3,1./3}, {0,1./3,1./3,1./3},
        {.25,.25,.25,.25}};
    double mn = 1e300;
    for (const auto& s : S)
        mn = std::min(mn, detAt(s[0], s[1], s[2], s[3]) / v6);
    return mn;
}

Index countInvalidTet10(const Mesh& m, double threshold) {
    Index n = 0;
    for (Index ci = 0; ci < m.cellCount(); ++ci) {
        if (m.cell(ci).type != CellType::Tet10) continue;
        const Index* v = m.cellNodes(ci);
        Vec3 p[10];
        for (int k = 0; k < 10; ++k) p[k] = m.node(v[k]);
        if (tet10MinJacobianRatio(p) <= threshold) ++n;
    }
    return n;
}

} // namespace paacs::mesh
