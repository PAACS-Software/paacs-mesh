// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Bounded interval-Bernstein certificate for one quadratic Tet10 mapping.
// Connectivity order is corners 0,1,2,3 followed by edges
// 01,12,20,03,13,23.  The result certifies the complete reference tetrahedron;
// sampled Jacobians are never an acceptance criterion.

#include "paacs/core/Types.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace paacs::algo {

struct Tet10Interval {
    double lo = 0.0;
    double hi = 0.0;
};

enum class Tet10CertificateStatus : std::uint8_t {
    Certified,
    PositiveBelowThreshold,
    Nonpositive,
    IncompleteDepth,
    IncompleteSubcells,
    IncompleteCoefficients,
    IncompleteWork,
    InvalidInput,
    RangeUnrepresentable,
    AllocationFailed
};

struct Tet10CertificatePolicy {
    double minimumJacobianRatio = 0.0;
    std::uint32_t maxDepth = 0;
    std::uint64_t maxVisitedSubcells = 1;
    std::uint64_t maxCoefficientTests = 20;
    std::uint64_t maxWorkUnits = 4096;
};

struct Tet10CertificateWitness {
    std::uint64_t pathBits = 0;
    std::uint32_t depth = 0;
    std::uint32_t coefficient = 0;
    Tet10Interval shifted{};
    Tet10Interval determinant{};
};

struct Tet10CertificateMetrics {
    std::uint64_t visitedSubcells = 0;
    std::uint64_t coefficientTests = 0;
    std::uint64_t workUnits = 0;
    std::uint32_t maximumDepth = 0;
};

struct Tet10CertificateResult {
    Tet10CertificateStatus status = Tet10CertificateStatus::InvalidInput;
    Tet10CertificateWitness witness{};
    Tet10CertificateMetrics metrics{};
    std::array<Tet10Interval,20> rootDeterminant{};
    std::array<Tet10Interval,20> rootShifted{};
    double normalizationScale = 0.0;
    double conservativeMinimumLowerBound =
        (std::numeric_limits<double>::infinity)();
    bool complete = false;
};

Tet10CertificateResult certifyTet10Jacobian(
    const std::array<Vec3,10>& nodes,
    const Tet10CertificatePolicy& policy) noexcept;

#if defined(PAACS_TOPO_TESTING)
// Plain-value child restriction seam. child is 0..7 and uses the exact same
// fixed lexicographic simplex-blossom recursion as production.
bool restrictTet10BernsteinChildForTesting(
    const std::array<Tet10Interval,20>& parent, std::uint32_t child,
    std::array<Tet10Interval,20>& result) noexcept;
#endif

static_assert(std::is_standard_layout_v<Tet10Interval>);
static_assert(std::is_trivially_copyable_v<Tet10Interval>);
static_assert(std::is_standard_layout_v<Tet10CertificateMetrics>);
static_assert(std::is_trivially_copyable_v<Tet10CertificateMetrics>);
static_assert(std::is_standard_layout_v<Tet10CertificateResult>);
static_assert(std::is_trivially_copyable_v<Tet10CertificateResult>);

} // namespace paacs::algo
