// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/core/Types.hpp"
#include <array>
namespace paacs::mesh {
// Sufficient positivity check for the COMPLETE trilinear [0,1]^3 map in
// canonical Hex8 order. No allocations or subdivisions. Coordinates are
// translated and scaled before evaluating its 27 tensor Bernstein values.
// Uses strict IEEE round-to-nearest arithmetic, as required by PAACS; a
// 2^18*epsilon guard in normalized coordinates covers roundoff conservatively.
// false includes valid but numerically unresolved/highly distorted cells;
// it does not mean the cell is necessarily inverted. It is not an exact
// arithmetic certificate, nor a check for overlap with OTHER elements.
bool positiveHex8Jacobian(const std::array<Vec3, 8> &points) noexcept;
} // namespace paacs::mesh
