// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Complete bounded cavities in an existing positive Tet4 mesh of an immutable
// closed triangular PLC. This is polyhedral recovery, not CAD certification.
#include "paacs/mesh/Mesh.hpp"
#include <cstdint>
#include <stop_token>
#include <functional>
#include <vector>
namespace paacs::bnd {
enum class TetPlcCavityReason : std::uint8_t {None,InvalidRequest,InvalidSource,InvalidPartialMesh,Unavailable,BudgetExceeded,Cancelled,AllocationFailed,InternalError};
struct TetPlcPatchSearchOptions {
    // Zero disables the local, backtracking region refill. Existing vertices
    // are reused; source triangles and the outer patch boundary stay fixed.
    // Candidate count bounds retained admissible cells; rejected geometry has
    // a separate cache of at most this size and still spends maximumWork.
    unsigned maximumStates=0,maximumCandidateTets=8192,maximumVertices=64;
    unsigned maximumShells=2;
};
struct TetPlcCavityOptions {
    mesh::Classification region{3,0};
    // maximumCells counts Tet4 cells; the unchanged Tri3 source is additional.
    std::uint64_t maximumNodes=200000,maximumCells=200000,maximumSourceTriangles=20000;
    std::uint64_t maximumWork=100000000,maximumWorkingBytes=512ULL*1024*1024;
    unsigned maximumExpansionCells=128,maximumPatchCells=512,maximumComponentFaces=512;
    unsigned maximumCandidatePoints=32,maximumPatches=32;
    // Optional bounded interior proposals; zero preserves the established
    // deterministic candidate set. The seed never perturbs source coordinates.
    unsigned maximumRandomCandidatePoints=0;
    std::uint64_t randomSeed=1;
    // Additional unfilled front components that one visibility patch may
    // absorb. Original source triangles remain hard barriers. Zero disables.
    unsigned maximumMergedComponents=0;
    TetPlcPatchSearchOptions patchSearch;
    // A caller retrying an already exhausted fan may go straight to the region
    // search. Enable visibility expansion, region search, or local advance.
    bool useVisibilityExpansion=true;
    double warningQuality=.25;
    std::stop_token cancellation;
    // Optional caller-owned material predicate for introduced interior points
    // and every Tet4 centroid. It only adds refusals to exact PLC validation.
    // Throwing callbacks refuse without publication. The caller accounts for
    // predicate work and keeps captured geometry alive through this call.
    std::function<bool(Vec3)> interiorAdmission;
    // Optional cell-only policy. Empty retains interiorAdmission for centroids.
    // This never changes admission of introduced interior nodes or exact PLC,
    // orientation, conformity and contact checks.
    std::function<bool(Vec3)> centroidAdmission;
    // Bounded local advance before broader cavity searches. Zero disables it.
    unsigned maximumLocalFrontTets=0;
};
struct TetPlcCavityOutcome {
    bool published=false,destinationUnchanged=true;
    TetPlcCavityReason reason=TetPlcCavityReason::InternalError;
    const char* diagnostic="not-started";
    std::uint64_t workUsed=0,workingPayloadEstimate=0,inputCells=0,outputCells=0;
    std::uint64_t initialFrontFaces=0,remainingFrontFaces=0,cavitiesClosed=0,removedCells=0,addedCells=0;
    std::uint64_t candidatePointsTested=0,randomCandidatePointsProposed=0;
    std::uint64_t patchesCommitted=0,mergedComponents=0,sourceBlockedCandidates=0,otherComponentBlockedCandidates=0;
    std::uint64_t regionSearchAttempts=0,regionSearchStates=0,regionSearchBacktracks=0;
    std::uint64_t regionSearchCandidates=0,regionSearchPairCacheHits=0,regionSearchPatches=0;
    std::uint64_t regionSearchBoundHits=0;
    std::uint64_t regionVertexBoundHits=0,regionStateBoundHits=0;
    std::uint64_t largestRegionVertices=0,largestRegionFaces=0;
    // Elapsed wall time, including refused work; these are not CPU counters.
    double preparationSeconds=0,visibilitySeconds=0,regionSearchSeconds=0;
    double publicationSeconds=0,validationSeconds=0;
    double sourceVolume=0,meshVolume=0,minimumQuality=0;
    bool qualityWarning=false;
    // On an unavailable recovery, indices of original source triangles still
    // exposed to an unfilled component. Diagnostic only; no partial mesh.
    std::vector<Index> unresolvedSourceTriangles;
    // Original-source barriers hit by visibility expansion, highest hit count
    // first. Includes covered source facets, which can still block a repair.
    std::vector<Index> blockedSourceTriangles;
    std::uint64_t localFrontTets=0,localFrontClosed=0;
    double localFrontSeconds=0;
};
// boundary contains Tri3 cells only. partial contains Tet4 cells only and uses
// exactly the same original node identities/classes/coordinates as boundary
// in its node prefix; appended nodes belong to options.region. An empty partial
// is allowed for suitable star-shaped PLCs. Coordinate-only welding is never
// performed. Source face/node classifications and the source node prefix are
// retained; unused appended points are omitted from the completed output.
//
// Source/partial/destination aliases are refused. Both inputs and destination
// remain unchanged on every refusal. Success returns a complete, positive,
// nonoverlapping PLC mesh with the exact original oriented boundary. There is
// no quality rejection threshold. Unsupported topology or exhausted allowances
// refuse instead of publishing a partially recovered mesh.
//
// Work counts modelled scans and contact tests, not CPU instructions. Working
// storage admission is conservative payload accounting, not a process-RSS cap.
// CAD projection, body identity beyond the supplied region, geometric tolerance
// and CAD material classification remain the adapter's responsibility.
TetPlcCavityOutcome recoverTetPlcCavities(const mesh::Mesh& boundary,
    const mesh::Mesh& partial,const TetPlcCavityOptions&,mesh::Mesh& destination) noexcept;
} // namespace paacs::bnd
