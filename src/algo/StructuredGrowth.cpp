// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/StructuredGrowth.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>

namespace paacs::algo {

StructuredAxisDistribution
buildStructuredAxisDistribution(std::uint32_t cells, const StructuredAxisGrowth &growth) noexcept {
    StructuredAxisDistribution out;
    out.cells = cells;
    if (cells == 0) {
        out.reason = StructuredGrowthReason::InvalidCellCount;
        return out;
    }
    if (!std::isfinite(growth.endToStartCellWidthRatio)) {
        out.reason = StructuredGrowthReason::NonfiniteRatio;
        return out;
    }
    if (!(growth.endToStartCellWidthRatio > 0.0)) {
        out.reason = StructuredGrowthReason::NonPositiveRatio;
        return out;
    }
    if (!std::isfinite(growth.maximumAdjacentCellWidthRatio) ||
        growth.maximumAdjacentCellWidthRatio < 1.0) {
        out.reason = StructuredGrowthReason::InvalidAdjacentLimit;
        return out;
    }

    try {
        if (!growth.boundaryFractions.empty()) {
            const auto f = growth.boundaryFractions;
            if (f.size() != static_cast<std::size_t>(cells) + 1 || f.front() != 0 ||
                f.back() != 1) {
                out.reason = StructuredGrowthReason::InvalidCoordinates;
                return out;
            }
            double low = 1, high = 0, adjacent = 1, previous = 0;
            for (std::size_t i = 1; i < f.size(); ++i) {
                const double width = f[i] - f[i - 1];
                if (!std::isfinite(f[i]) || !(width > 0)) {
                    out.reason = StructuredGrowthReason::InvalidCoordinates;
                    return out;
                }
                low = std::min(low, width);
                high = std::max(high, width);
                if (i > 1)
                    adjacent = std::max({adjacent, width / previous, previous / width});
                previous = width;
            }
            if (adjacent > growth.maximumAdjacentCellWidthRatio * (1 + 1e-10)) {
                out.reason = StructuredGrowthReason::AdjacentLimitExceeded;
                return out;
            }
            out.boundaries.assign(f.begin(), f.end());
            out.firstCellFraction = f[1];
            out.lastCellFraction = 1 - f[f.size() - 2];
            out.minimumCellFraction = low;
            out.maximumCellFraction = high;
            out.maximumAdjacentRatio = adjacent;
            out.valid = true;
            out.reason = StructuredGrowthReason::None;
            return out;
        }
        out.boundaries.resize(static_cast<std::size_t>(cells) + 1U);
        out.boundaries.front() = 0.0;
        out.boundaries.back() = 1.0;

        if (cells == 1 || growth.mode == StructuredAxisGrowthMode::Uniform ||
            growth.endToStartCellWidthRatio == 1.0) {
            const double inv = 1.0 / static_cast<double>(cells);
            for (std::uint32_t i = 1; i < cells; ++i)
                out.boundaries[static_cast<std::size_t>(i)] = static_cast<double>(i) * inv;
            out.firstCellFraction = inv;
            out.lastCellFraction = inv;
            out.minimumCellFraction = inv;
            out.maximumCellFraction = inv;
            out.maximumAdjacentRatio = 1.0;
            out.valid = true;
            out.reason = StructuredGrowthReason::None;
            return out;
        }

        const long double totalRatio = static_cast<long double>(growth.endToStartCellWidthRatio);
        const long double logQ = std::log(totalRatio) / static_cast<long double>(cells - 1U);
        const long double q = std::exp(logQ);
        const long double reciprocal = 1.0L / q;
        const long double adjacent = std::max(q, reciprocal);
        if (!std::isfinite(static_cast<double>(q)) || !(q > 0.0L)) {
            out.reason = StructuredGrowthReason::NumericFailure;
            return out;
        }
        if (adjacent > static_cast<long double>(growth.maximumAdjacentCellWidthRatio) *
                           (1.0L + 32.0L * std::numeric_limits<long double>::epsilon())) {
            out.reason = StructuredGrowthReason::AdjacentLimitExceeded;
            return out;
        }

        long double sum = 0.0L;
        if (std::abs(q - 1.0L) <= 64.0L * std::numeric_limits<long double>::epsilon()) {
            sum = static_cast<long double>(cells);
        } else {
            const long double qN = std::exp(logQ * static_cast<long double>(cells));
            sum = (qN - 1.0L) / (q - 1.0L);
        }
        if (!std::isfinite(static_cast<double>(sum)) || !(sum > 0.0L)) {
            out.reason = StructuredGrowthReason::NumericFailure;
            return out;
        }

        long double width = 1.0L / sum;
        long double position = 0.0L;
        double minWidth = (std::numeric_limits<double>::infinity)();
        double maxWidth = 0.0;
        for (std::uint32_t i = 0; i < cells; ++i) {
            const double current = static_cast<double>(width);
            if (!std::isfinite(current) || !(current > 0.0)) {
                out.reason = StructuredGrowthReason::NumericFailure;
                return out;
            }
            if (i == 0)
                out.firstCellFraction = current;
            if (i + 1U == cells)
                out.lastCellFraction = current;
            minWidth = std::min(minWidth, current);
            maxWidth = std::max(maxWidth, current);
            position += width;
            if (i + 1U < cells) {
                const double next = static_cast<double>(position);
                if (!std::isfinite(next) || !(next > out.boundaries[static_cast<std::size_t>(i)])) {
                    out.reason = StructuredGrowthReason::NumericFailure;
                    return out;
                }
                out.boundaries[static_cast<std::size_t>(i + 1U)] = next;
            }
            width *= q;
        }
        // Keep the public block boundary exact and absorb only accumulated
        // round-off in the last cell.
        out.boundaries.back() = 1.0;
        out.lastCellFraction = 1.0 - out.boundaries[static_cast<std::size_t>(cells - 1U)];
        minWidth = std::min(minWidth, out.lastCellFraction);
        maxWidth = std::max(maxWidth, out.lastCellFraction);
        if (!(out.lastCellFraction > 0.0) || !std::isfinite(minWidth) || !std::isfinite(maxWidth)) {
            out.reason = StructuredGrowthReason::NumericFailure;
            return out;
        }
        out.minimumCellFraction = minWidth;
        out.maximumCellFraction = maxWidth;
        out.maximumAdjacentRatio = static_cast<double>(adjacent);
        out.valid = true;
        out.reason = StructuredGrowthReason::None;
        return out;
    } catch (const std::bad_alloc &) {
        out.reason = StructuredGrowthReason::AllocationFailed;
    } catch (...) {
        out.reason = StructuredGrowthReason::NumericFailure;
    }
    return out;
}

} // namespace paacs::algo
