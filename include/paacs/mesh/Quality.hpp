// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include <array>
#include <string>
#include "paacs/mesh/Mesh.hpp"

namespace paacs::mesh {

// Signed volume of tetrahedron (a,b,c,d); positive if d on the +side of abc.
double tetVolume(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d);

// Radius-ratio quality in [0,1]: 3*r_in / r_circ. 1 = regular tet, 0 = sliver.
double tetQuality(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d);

// Minimum corner scaled Jacobian, normalized so the regular tet scores 1.
// Range: 1 = regular, -> 0 degenerate, < 0 inverted. (Verdict-compatible up
// to the sqrt(2) normalization, which we fold in so "1 is perfect".)
double tetScaledJacobian(const Vec3& a, const Vec3& b, const Vec3& c,
                         const Vec3& d);

// Equilateral-volume skewness (ANSYS convention):
//   skew = (V_eq - V) / V_eq,  V_eq = volume of the regular tet sharing the
// circumsphere. 0 = perfect, -> 1 degenerate.
double tetSkewness(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d);

// Min/max dihedral angles in degrees (regular tet: 70.53 on all six edges).
void tetDihedralAngles(const Vec3& a, const Vec3& b, const Vec3& c,
                       const Vec3& d, double& minDeg, double& maxDeg);

struct MetricStats { double min = 1e300, max = -1e300, mean = 0; Index n = 0; };

struct QualityReport {
    MetricStats radiusRatio;     // higher better, [0,1]
    MetricStats scaledJacobian;  // higher better, 1 = regular
    MetricStats skewness;        // lower better, 0 = regular

    // Quadratic population (Tet10), reported separately: sampled min
    // Jacobian ratio distribution, invalid count, corner-tet volume.
    Index nTet10 = 0, nTet10Invalid = 0;
    MetricStats tet10Jacobian;
    double tet10CornerVolume = 0;

    // Elements with any non-finite node coordinate. They are excluded from
    // every statistic above (a single NaN would otherwise poison min/mean/
    // max silently) and reported explicitly.
    Index nNonFinite = 0;
    double minDihedralDeg = 180, maxDihedralDeg = 0;
    double totalVolume = 0;
    Index  nTets = 0, nInverted = 0;

    std::array<Index, 10> radiusRatioHist{};  // bins [0,0.1) ... [0.9,1.0]
    std::string toString() const;
};

QualityReport analyzeTets(const Mesh& m);

// Back-compat helpers.
struct QualityStats { double min = 1, max = 0, mean = 0; Index n = 0; };
QualityStats tetQualityStats(const Mesh& m);
double totalTetVolume(const Mesh& m);


// Minimum sampled Jacobian ratio of a curved 10-node tet (VTK/MSH node
// ordering: corners 0-3, then edge midsides 01,12,20,03,13,23).
// Returns min over sample points of det J(xi) / det J_linear, where
// J_linear comes from the corner tet: 1.0 for a straight element, <= 0 for
// an invalid (locally inverted) curved element. Sampling (corners, edge
// midpoints, face centroids, centroid) is standard practice but not a
// certificate; Bezier-coefficient bounds are the roadmap upgrade.
double tet10MinJacobianRatio(const Vec3 p[10]);

// Count Tet10 cells whose sampled Jacobian ratio falls at or below the
// threshold (0.0 = locally inverted somewhere).
Index countInvalidTet10(const Mesh& m, double threshold = 0.0);

} // namespace paacs::mesh
