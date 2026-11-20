// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/sizing/SizeField.hpp"
#include <cstdint>

namespace paacs::sizing {
// A conservative lower bound over the WHOLE closed box, not corner sampling.
// Implementations must return a finite positive value for valid boxes, or NaN
// when the query cannot be supported. Immutable implementations are shareable.
class BoundedSizeField : public SizeField {
  public:
    virtual bool valid() const noexcept = 0;
    virtual double minimumOnBox(const Vec3 &lower, const Vec3 &upper) const noexcept = 0;
    // Algorithmic query-work units, not CPU instructions or elapsed time.
    virtual std::uint64_t boxQueryWorkBound() const noexcept { return 1; }
};
} // namespace paacs::sizing
