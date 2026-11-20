// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private, opt-in preparation for PLC recovery. Nonlinear wire candidates are
// explicitly admitted and still require exact surface and CAD facet gates.
#include "OcctStarMesher.hpp"
#include "paacs/geometry/OcctModel.hpp"
#include <span>

namespace paacs::geo::detail {
enum class OcctWireRefinementReason : std::uint8_t {
    None, InvalidRequest, UnsupportedGeometry, BudgetExceeded, Cancelled,
    InvalidBoundary, NativeQueryFailed, AllocationFailed, InternalError
};
struct OcctWireBoundary {
    int bodyTag = -1;
    mesh::Mesh* boundary = nullptr;
};
struct OcctWireRefinementOptions {
    int edgeTag = -1;
    unsigned segmentsPerExistingSegment = 2;
    bool allowNonlinearWireCandidates = false;
    bool allowClosedWireCandidates = false;
    std::uint64_t maximumAddedNodes = 4096;
    OcctStarOptions limits{};
};
struct OcctWireRefinementOutcome {
    bool published = false;
    const char* diagnostic = "internal-error";
    OcctWireRefinementReason reason = OcctWireRefinementReason::InternalError;
    std::uint64_t workUsed = 0, addedNodes = 0, addedTriangles = 0;
    std::uint64_t nativeUseChecks = 0;
};
// Atomically subdivides a selected actual CAD wire in every supplied
// incident body. Original coordinates/classes remain unchanged. All coowners
// of every affected face must be present, and all incident triangles use one
// common chain. Exact closed/embedded surface certificates are required before
// and after. Owning pcurves are checked at source tolerances. CAD facet error,
// whole volume validity and persistent cross-call identities remain caller
// gates. Work limits count explicit visits, not OCCT allocator/predicate work.
// Closed wires require explicit admission and a complete, uniquely ordered
// native parameter cycle. Trimmed periodic arcs retain their finite interval.
// A fan split on
// a nonlinear curve is a candidate, not a general trimmed-surface mesher.
OcctWireRefinementOutcome refineOcctWireRefinement(const OcctModel&,
    std::span<const OcctWireBoundary>, const OcctWireRefinementOptions&) noexcept;
} // namespace paacs::geo::detail
