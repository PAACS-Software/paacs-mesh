// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <cstdint>
#include <span>
#include <vector>

namespace paacs::algo {

// Structured grading is compiled once before node generation.  The hot loops
// consume only the resulting normalized boundary coordinates; no logarithms,
// sizing-field calls or associative lookups occur per node/cell.
enum class StructuredAxisGrowthMode : std::uint8_t {
    Uniform,
    // Ratio of the final cell width to the first cell width.  Adjacent cells
    // follow one constant geometric ratio derived from this total expansion.
    EndToStartCellWidthRatio
};

struct StructuredAxisGrowth {
    StructuredAxisGrowthMode mode = StructuredAxisGrowthMode::Uniform;
    double endToStartCellWidthRatio = 1.0;
    // Symmetric bound: both q and 1/q must be <= this value.
    double maximumAdjacentCellWidthRatio = 1.3;
    // Optional caller-owned 0...1 boundaries, exactly cells+1 entries.
    // Overrides geometric grading; validated before copying. Retain until
    // the meshing call returns. Used by reusable capped size profiles.
    std::span<const double> boundaryFractions{};
};

struct StructuredBlockGrowth {
    StructuredAxisGrowth u{};
    StructuredAxisGrowth v{};
    StructuredAxisGrowth w{};
};

enum class StructuredGrowthReason : std::uint8_t {
    None,
    InvalidCellCount,
    NonfiniteRatio,
    NonPositiveRatio,
    InvalidAdjacentLimit,
    AdjacentLimitExceeded,
    AllocationFailed,
    NumericFailure,
    InvalidCoordinates
};

struct StructuredAxisDistribution {
    bool valid = false;
    StructuredGrowthReason reason = StructuredGrowthReason::NumericFailure;
    std::uint32_t cells = 0;
    // Normalized cell-boundary coordinates.  Size is cells+1, begins exactly
    // at 0 and ends exactly at 1.
    std::vector<double> boundaries;
    double firstCellFraction = 0.0;
    double lastCellFraction = 0.0;
    double minimumCellFraction = 0.0;
    double maximumCellFraction = 0.0;
    double maximumAdjacentRatio = 0.0;

    double coordinate(std::uint32_t ordinal) const noexcept {
        return boundaries[static_cast<std::size_t>(ordinal)];
    }
};

StructuredAxisDistribution
buildStructuredAxisDistribution(std::uint32_t cells, const StructuredAxisGrowth &growth) noexcept;

} // namespace paacs::algo
