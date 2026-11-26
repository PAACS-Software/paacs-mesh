// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "OcctBoundaryInternal.hpp"
#include "paacs/geometry/OcctMesher.hpp"

namespace paacs::geo::detail {
struct OcctSourceFallbackDecision {
    bool permitted=false;
    const char* diagnostic="disabled";
};

// Provider substitution is not a relaxation of the source certificate or a
// way to evade the chart's explicit allowances. A boundary reason alone is
// insufficient: historical chart wrappers also use it for resource refusals.
inline OcctSourceFallbackDecision nativeSourceFallbackDecision(
    const NativeCadMeshOptions& options,std::size_t solids,bool spatialSizing,
    bool advancingInteriorFront,OcctBoundaryReason reason,
    const OcctBoundaryConstructionWitness& witness,bool constructionBudgetLimited,
    bool attemptAvailable) noexcept {
    if(!options.allowNativeSourceFallback)return {};
    if(options.geometry.cancellation&&options.geometry.cancellation->load(std::memory_order_relaxed))
        return {false,"cancelled"};
    if(solids!=1||!spatialSizing||!advancingInteriorFront||!options.useNativeSurfaceFallback||
       options.elements!=algo::PrimitiveElements::Tet4)
        return {false,"requires_single_body_live_constructor"};
    if(!options.faceSizes.empty()||!options.sizeBoxes.empty()||options.geometry.spatialField||
       !options.geometry.bodyControls.empty()||options.geometry.maximumElementSizeRatio>0||
       options.geometry.warningElementSizeRatio>0)
        return {false,"explicit_sizing_controls_not_supported"};
    if(options.construction.surfacePointsPerFace||options.construction.surfaceChordPointsPerFace)
        return {false,"explicit_surface_caps_not_supported"};
    if(constructionBudgetLimited)return {false,"source_construction_budget"};
    if(reason!=OcctBoundaryReason::NonconformingEdge&&reason!=OcctBoundaryReason::OpenBoundary&&
       reason!=OcctBoundaryReason::MissingTriangulation)
        return {false,"source_reason_not_eligible"};
    if(!witness.geometricRefusal)return {false,"no_geometric_refusal_witness"};
    if(!attemptAvailable)return {false,"shared_source_attempt_limit"};
    return {true,"named_geometric_source_refusal"};
}
} // namespace paacs::geo::detail
