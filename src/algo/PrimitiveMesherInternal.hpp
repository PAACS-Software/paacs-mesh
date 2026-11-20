// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/algo/PrimitiveMesher.hpp"
#include <algorithm>
namespace paacs::algo::primitive_detail {
inline void recordQuality(PrimitiveQualityReport &r, Index cell, double q, bool sharp,
                          double threshold) {
    ++r.histogram[std::min(9u, static_cast<unsigned>(std::max(0., q) * 10))];
    if (r.worstCell == kInvalidIndex || q < r.worstQuality) {
        r.worstCell = cell;
        r.worstQuality = q;
    }
    if (q < threshold) {
        ++r.lowQualityCells;
        r.sharpFeatureLowQualityCells += sharp;
        if (r.exampleCount < r.examples.size())
            r.examples[r.exampleCount++] = {cell, q, sharp};
        else
            r.examplesTruncated = true;
    }
}
inline void mergeQuality(PrimitiveQualityReport &a, const PrimitiveQualityReport &b) {
    for (unsigned i = 0; i < a.histogram.size(); ++i)
        a.histogram[i] += b.histogram[i];
    if (b.worstCell != kInvalidIndex &&
        (a.worstCell == kInvalidIndex || b.worstQuality < a.worstQuality)) {
        a.worstCell = b.worstCell;
        a.worstQuality = b.worstQuality;
    }
    a.lowQualityCells += b.lowQualityCells;
    a.sharpFeatureLowQualityCells += b.sharpFeatureLowQualityCells;
    for (unsigned i = 0; i < b.exampleCount && a.exampleCount < a.examples.size(); ++i)
        a.examples[a.exampleCount++] = b.examples[i];
    a.examplesTruncated = a.lowQualityCells > a.exampleCount;
}
PrimitiveMeshingResult tetrahedra(const PrimitiveShape &, const PrimitiveMeshingOptions &,
                                  mesh::Mesh &) noexcept;
PrimitiveMeshingResult profiled(const PrimitiveShape &, const PrimitiveMeshingOptions &,
                                mesh::Mesh &);
PrimitiveReason profiledCylinderSource(const PrimitiveShape &, const PrimitiveMeshingOptions &,
                                       PrimitiveMeshingResult &, mesh::Mesh &) noexcept;
} // namespace paacs::algo::primitive_detail
