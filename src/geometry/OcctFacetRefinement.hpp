// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private CAD-evaluated facet candidates. Native PLC recovery can select these
// after boundary barriers block completion; final volume admission remains a caller gate.
#include "OcctStarMesher.hpp"
#include "paacs/geometry/OcctModel.hpp"
#include <span>

namespace paacs::geo::detail {
enum class OcctFacetRefinementReason : std::uint8_t {
    None, InvalidRequest, UnsupportedGeometry, BudgetExceeded, Cancelled,
    InvalidBoundary, NativeQueryFailed, AllocationFailed, InternalError
};
struct OcctFacetBoundary { int bodyTag=-1; mesh::Mesh* boundary=nullptr; };
struct OcctFacetSelection { int bodyTag=-1; Index triangle=-1; };
struct OcctFacetRefinementOptions {
    std::span<const OcctFacetSelection> selections{};
    double maximumProjectionDistance=0;
    std::uint64_t maximumAddedNodes=1024;
    OcctStarOptions limits{};
};
struct OcctFacetRefinementOutcome {
    bool published=false;
    const char* diagnostic="internal-error";
    OcctFacetRefinementReason reason=OcctFacetRefinementReason::InternalError;
    std::uint64_t workUsed=0, addedNodes=0, addedTriangles=0, nativeQueries=0;
    double maximumDisplacement=0;
    int failedFaceTag=-1;
    Index failedTriangle=-1;
};
// For each selected original triangle, projects its canonical centroid onto its
// actual trimmed CAD face and replaces it with three triangles. The projected
// point's planar foot must lie strictly inside the source triangle. All old
// nodes/classes and all unaffected directed triangles are preserved. Shared
// facets require bitwise matching source triangles in every actual coowner;
// one point is evaluated once and reused across owners. Exact closure and
// surface embedding are certified before/after; all publication is atomic.
// The caller must audit CAD facet error, material containment and volume.
// This does not claim curved-patch conformity against continuous CAD, or repair
// missing trims, source cracks, mismatched shared triangulations or bad pcurves.
OcctFacetRefinementOutcome refineOcctFacets(const OcctModel&,
    std::span<const OcctFacetBoundary>,const OcctFacetRefinementOptions&) noexcept;
} // namespace paacs::geo::detail
