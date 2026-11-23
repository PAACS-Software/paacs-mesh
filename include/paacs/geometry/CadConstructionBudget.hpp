// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include <cstdint>
#include <limits>

namespace paacs::geo {
enum class CadResourcePolicy : std::uint8_t { Bounded, Normal };
// Normal requests keep accounting and overflow guards without an operational
// logical-work ceiling. Explicit caller limits still take precedence.
inline constexpr std::uint64_t cadNormalWorkAllowance =
    (std::numeric_limits<std::uint64_t>::max)()/4;
enum class CadConstructionEffort : std::uint8_t { Legacy, Automatic, Thorough };

// Stage allowances remain subordinate to the enclosing node, cell, work,
// modeled-memory, output and cancellation limits. Zero selects the policy's
// automatic allowance; positive values are explicit stage ceilings.
struct CadConstructionBudget {
    CadConstructionEffort effort = CadConstructionEffort::Legacy;
    std::uint64_t surfacePointsPerFace = 0, surfaceChordPointsPerFace = 0;
    std::uint64_t interiorCandidates = 0, interiorSites = 0;
    CadResourcePolicy resourcePolicy = CadResourcePolicy::Bounded;
};

// Queue termination describes construction effort, not a certificate that
// every physical edge meets its requested size or that the CAD is valid.
struct CadSurfaceSpacingReport {
    int faceTag = -1;
    std::uint64_t pointLimit = 0, chordPointLimit = 0, chartNodeLimit = 0;
    std::uint64_t addedPoints = 0, chordAddedPoints = 0;
    std::uint64_t pendingEntries = 0, blockedCandidates = 0;
    const char* termination = "not_attempted";
    bool queueExhausted = false, budgetLimited = false;
};
} // namespace paacs::geo
