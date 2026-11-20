// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/core/Types.hpp"
#include <array>
namespace paacs::algo {
// Sufficient interval test on the whole reference tetrahedron. The Jacobian
// matrix is affine, so each entry lies in the hull of its four vertex values.
// A false result is unresolved, NOT evidence of an invalid element. Particularly
// inexpensive for straight-sided elevation (including midpoint rounding).
// Tet10 order: corners, then edges 01,12,20,03,13,23. Strict IEEE FP required.
bool positiveTet10JacobianFast(const std::array<Vec3, 10> &nodes) noexcept;
} // namespace paacs::algo
