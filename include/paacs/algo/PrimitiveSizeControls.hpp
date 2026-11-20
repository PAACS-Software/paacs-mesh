// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/algo/PrimitiveMesher.hpp"
namespace paacs::algo {
// Cylinder/Tube analytic cap tags: 2 = start, 3 = end. Unsupported geometry
// or tags return an invalid control which compilation refuses. No silent
// projection of a curved wall or imported CAD entity onto the sweep axis.
sizing::ProfileSizeControl primitiveFaceSizeControl(const PrimitiveShape &, int faceTag,
                                                    double size, double growth = 1.3) noexcept;
} // namespace paacs::algo
