// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Last native fallback: immutable PLC front construction and bounded recovery.
#include "paacs/geometry/OcctMesher.hpp"
#include "OcctNativeValidation.hpp"
namespace paacs::geo::detail {
enum class OcctFrontSurface { Standard, Spaced, Regular, RegularQuality, Curvature, SpatialSizing };
struct OcctFrontOutcome {
    bool published=false;
    bool originalFacetCentersCertified=false;
    NativeCadMeshReason reason=NativeCadMeshReason::MeshingFailed;
    const char* diagnostic="not-started";
    const char* surfaceDiagnostic=nullptr;
    std::uint64_t workUsed=0,surfaceWork=0,nativeQueries=0,cavitiesClosed=0;
    std::uint64_t nativeRayAdmissions=0;
    std::uint64_t removedInteriorNodes=0,removedIncidentCells=0,removedCentroidCells=0;
    std::uint64_t randomRetries=0,randomCandidateProposals=0,localFrontSeeds=0;
    std::uint64_t mergedComponents=0,sourceBlockedCandidates=0,otherComponentBlockedCandidates=0,sourceFacetNodes=0;
    NativeCadRegionRecoveryReport regionRecovery;
    NativeCadConstraintRecoveryReport constraintRecovery;
    std::uint64_t liveMaterialInitialRejects=0,liveMaterialRemainingRejects=0,liveMaterialPatches=0;
    std::uint64_t localDelaunaySupportNodes=0,localDelaunayRefinementNodes=0;
    unsigned localDelaunayQualityEdits=0;
    const char* interiorFrontDiagnostic="not_attempted";
    std::uint64_t interiorFrontNodes=0,interiorFrontCandidates=0;
    std::uint64_t interiorFrontWork=0,interiorFrontQueries=0;
    std::uint64_t interiorFrontRemainingOversized=0;
    std::uint64_t interiorFrontActivations=0,interiorVoidActivations=0,interiorFrontAdvances=0;
    double interiorFrontSeconds=0;
    bool localDelaunayUsedGlobalFallback=false;
    double surfaceSeconds=0,volumeSeconds=0,meshVolume=0,chord=0;
    double maximumNodeError=0,maximumFacetError=0;
    std::uint64_t surfaceTriangles=0,rebuiltSurfaceCharts=0;
    OcctNativeValidationOutcome nativeValidation;
    OcctNativeValidationOutcome chordSampleValidation;
};
OcctFrontOutcome buildOcctFrontVolume(const OcctModel&,int solidTag,
    const NativeCadMeshOptions&,std::uint64_t maximumWork,mesh::Mesh&,
    std::uint64_t earlyConstraintWork=0,OcctFrontSurface surface=OcctFrontSurface::Standard,
    bool liveMaterialTrial=false) noexcept;
}
