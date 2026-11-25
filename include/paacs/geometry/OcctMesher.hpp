// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/geometry/OcctModel.hpp"
#include "paacs/geometry/CadConstructionBudget.hpp"
#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/algo/PrimitiveMesher.hpp"
#include <array>
#include <cstddef>
#include <optional>
#include <string>

namespace paacs::geo {
enum class NativeCadMeshReason : std::uint8_t {
    None, InvalidRequest, UnsupportedBody, UnresolvedBodyInterface,
    UnsupportedControl, RecognitionFailed, MeshingFailed, GeometryError,
    BudgetExceeded, Cancelled, AllocationFailed, InternalError,
    QualityRequirementNotMet, SizeRequirementNotMet
};
struct CadFaceSizeControl { int faceTag = -1; double size = 0, growth = 1.3; };
struct CadBoxSizeControl {
    Vec3 lower{}, upper{};
    double size = 0, growth = 1.3;
};
enum class NativeCadAssemblyPolicy : std::uint8_t { ConformingDomain, IndependentParts };
enum class NativeCadRecoveryMode : std::uint8_t { Standard, Advanced };
enum class NativeCadInterferenceState : std::uint8_t { NoMaterialOverlap, MaterialOverlap, Unresolved };
struct NativeCadInterference {
    int firstSolidTag=-1,secondSolidTag=-1;
    NativeCadInterferenceState state=NativeCadInterferenceState::Unresolved;
    double nativeOverlapVolume=0;
};
struct NativeCadPartitionAttempt {
    const char* policy = "none";
    const char* diagnostic = "not_attempted";
    int failedPieceTag = -1;
    std::uint64_t pieces = 0, interfaces = 0, work = 0;
    double seconds = 0, surfaceDeflection = 0;
};
struct NativeCadStageEvent {
    const char* stage;
    const char* diagnostic;
    int entityTag;
    bool beginning;
};
using NativeCadStageObserver = void (*)(const NativeCadStageEvent&) noexcept;
struct NativeCadMeshOptions {
    algo::GeometryMeshingOptions geometry;
    // IndependentParts preserves every original solid and its distinct nodes.
    // Parts may overlap; this certifies each part, not a global material domain.
    // Interference reports accompany the result. Default conformity is unchanged.
    NativeCadAssemblyPolicy assemblyPolicy=NativeCadAssemblyPolicy::ConformingDomain;
    // Automatic interior placement applies to compatible single-solid Tet4
    // requests. Other element families retain their established routes.
    algo::PrimitiveElements elements = algo::PrimitiveElements::Tet4;
    double maximumChordError = 0, maximumRelativeVolumeError = .01;
    // Explicit approximation policy for the exact source-conforming Tet4 PLC
    // fallback only. Interior nodes remain IN/ON. A centroid classified OUT
    // must be no farther from the original solid than the body's chord bound.
    // Four additional interior samples per tetrahedron must satisfy the same
    // native distance bound. Unknown states still refuse. Closure, contacts, native face/curve
    // identity, facet error and volume gates are unchanged. Default is strict.
    bool allowCentroidsWithinChordTolerance = false;
    double minimumSize = 0;
    std::vector<CadFaceSizeControl> faceSizes;
    // Select the native size-field route for one Tet4 solid. A box also
    // selects this route. Existing analytic cap controls remain available.
    bool useSpatialSizing = false;
    std::vector<CadBoxSizeControl> sizeBoxes;
    // Native construction and validation share this logical work admission.
    // Match the geometry meshing default so an unset request admits both phases.
    std::uint64_t maximumValidationWork = 50'000'000;
    std::uint64_t maximumPrimitiveWorkingBytes = 512ULL * 1024 * 1024;
    bool optimizeProfileInterior = false;
    std::uint64_t maximumQualityNodeTrials = 100'000;
    // Point/edge contacts may retain separate body nodes. Other contact
    // policies are explicit; ConformingDomain refuses overlapping material.
    bool separateLowerDimensionalContacts = false;
    bool useNativeSurfaceFallback = true;
    // Controls native Boolean partition attempts. Explicit native-kernel
    // recovery can still fill an existing CAD surface when cuts are disabled.
    bool useNativePartitionFallback = true;
    // Optional small-piece recovery after direct/swept routes. Longitudinal
    // native/Cartesian cuts share exact interfaces; each piece uses a bounded Tet4 front.
    // Requires native kernel admission and obeys useNativePartitionFallback.
    bool useLongitudinalRecovery = false;
    std::uint64_t maximumNativePartitionPieces = 512;
    // Opt-in separate full-face contacts for native axis-aligned boxes at the
    // eight-corner Tet4 resolution. Exact opposite triangles must cover the
    // complete native contact area. Body nodes and identities stay distinct.
    bool separateMatchedBoxFaceContacts = false;
    // Native validation workers use independent owned geometry/classifiers.
    // Meshing thread count remains a separate control; one retains serial behavior.
    unsigned maximumNativeValidationThreads = 1;
    // Optional per-call, one-surface projector reuse. Every point is solved
    // and all trim/branch checks still run; native solver state is thread-local.
    bool reuseNativeQueries = false;
    // Skip expensive facet distance checks for privately discarded direct-star
    // candidates. Every accepted candidate still completes all facet checks.
    bool deferDiscardedNativeFacetChecks = false;
    // Explicit separate contact surfaces for general native assemblies. The
    // fallback splits one copied assembly and requires exact opposite shared
    // triangles while retaining each original body's nodes and identity.
    bool separateMatchedNativeFaceContacts = false;
    // Additional bounded native admission and recovery routes, tried after
    // established routes. Fan proposals and recovered Tet4 centroids require
    // original-body inclusion. Recovery may change artificial exterior-face
    // diagonals; native curve segments and shared partition triangles stay fixed.
    bool useNativeKernelAdmission = false;
    // Bounded interior-candidate retries after deterministic PLC recovery fails.
    // Source coordinates stay fixed; zero retries disables diversified search.
    unsigned maximumRecoveryRandomRetries=2;
    std::uint64_t recoverySeed=1;
    // Optional existing-node edge-cavity fallback using one staged point.
    // Zero preserves the previous fixed-PLC route; maximum 64 proposals.
    unsigned maximumEdgeRecoverySteinerCandidates=0;
    // Final bounded region refill after the established PLC repairs refuse.
    // States are bounded per local region/shell; zero disables this route.
    unsigned maximumRecoveryPatchStates=2048;
    // Advanced tries to recover source edges and triangles in the full
    // Delaunay complex before selecting material for the audited Tet4 fallback.
    NativeCadRecoveryMode recoveryMode=NativeCadRecoveryMode::Standard;
    // Try bounded Tet4 recovery on a fixed partition when the coarse split
    // leaves a piece with a handle. Original CAD validation still applies.
    bool useLocalPartitionRecovery = false;
    // Experimental curve and Tet4 recovery, route scheduling and exact query reuse.
    bool useAdaptiveTetRecovery = false;
    // Try a discrete fill followed by local original-CAD cavity repair.
    bool useLiveMaterialRecovery = false;
    // Try bounded, boundary-preserving Tet4 quality improvement on native pieces
    // and after spatial sizing. The latter has an independent postfill budget;
    // it does not reduce the insertion budget or change the prescribed surface.
    bool improveNativeTetQuality = false;
    // Synchronous diagnostics on the calling thread. Stage names and
    // diagnostics are string literals; a null observer disables tracing.
    NativeCadStageObserver stageObserver = nullptr;
    // Try size-graded native CAD face meshing before the single-body Tet4 fallback.
    bool improveNativeSurfaceSpacing = false;
    // Use a graded advancing front for native Tet4 surface spacing.
    bool useRegularSurfaceFront = false;
    // An unset override enables this bounded Tet4 route in Advanced recovery.
    std::optional<bool> improveRegularMeshQuality;
    // Opt into curvature-sized native curves and faces for a bounded
    // single-solid Tet4 recovery trial. Zero keeps the default route order.
    unsigned curvaturePointsPerTurn = 0;
    // Experimental positive-only original CAD ray witnesses in the serial
    // strict live Tet4 route. Other routes use ordinary validation.
    bool useTrimmedRayWitnesses = false;
    // Try local interior-point insertion on a recovered Tet4 mesh in the
    // strict live material recovery trial.
    bool useLocalDelaunayRefinement = false;
    // Unset selects initial wall/void construction automatically for compatible
    // single-solid Tet4 requests. True forces the existing native size-field
    // route and its control restrictions; false retains the legacy dispatch.
    // Optional sizing controls can trigger a separate adaptive stage.
    std::optional<bool> useAdvancingInteriorFront;
    [[nodiscard]] bool advancingInteriorFrontEnabled(std::size_t solidCount) const noexcept {
        if(useAdvancingInteriorFront.has_value())return *useAdvancingInteriorFront;
        // Legacy analytic cap controls can be coarser than the global ceiling;
        // the native size-field route deliberately accepts only finer controls.
        for(const auto& control:faceSizes)if(control.size>geometry.targetSize)return false;
        return solidCount==1&&elements==algo::PrimitiveElements::Tet4&&useNativeSurfaceFallback&&
            (assemblyPolicy==NativeCadAssemblyPolicy::ConformingDomain||
             assemblyPolicy==NativeCadAssemblyPolicy::IndependentParts)&&!separateMatchedNativeFaceContacts&&
            geometry.strategy==algo::GeometryStrategy::Auto&&geometry.bodyControls.empty()&&
            !geometry.spatialField&&geometry.maximumElementSizeRatio==0&&geometry.warningElementSizeRatio==0;
    }
    [[nodiscard]] bool regularMeshQualityEnabled() const noexcept {
        return improveRegularMeshQuality.value_or(recoveryMode==NativeCadRecoveryMode::Advanced);
    }
    CadConstructionBudget construction;
    // Finalize numerical slivers after point insertion without moving nodes.
    // Unset enables the bounded stage for Normal admission only.
    std::optional<bool> initialTopologyQuality;
    // Experimental provider substitution after a named geometric chart/wire
    // refusal. Retains the live constructor and all Native validity checks.
    // Explicit local surface controls/caps are not bypassed.
    // A source refusal remains on this pipeline rather than retrying a legacy
    // whole-method constructor. This toggle does not force initial-front ON.
    bool allowNativeSourceFallback=false;
    // Experimental order within the live size-field constructor: recover the
    // boundary before placing interior sites. Does not force that constructor.
    // A failed source is not replaced by legacy whole-method retry in AUTO.
    bool recoverBoundaryBeforeInteriorSites=false;
};

// Existing aggregate options retain their bounded defaults. This preset opts
// into progressive construction; explicit work fields assigned afterwards
// remain authoritative, including a caller's deliberately small limit.
inline NativeCadMeshOptions normalNativeCadMeshOptions() {
    NativeCadMeshOptions options;
    options.construction.resourcePolicy=CadResourcePolicy::Normal;
    options.geometry.maximumMeshingWork=cadNormalWorkAllowance;
    options.maximumValidationWork=cadNormalWorkAllowance;
    return options;
}
struct NativeCadConstraintRecoveryReport {
    const char* diagnostic="not_attempted";
    std::uint64_t candidateCells=0,missingEdgesBefore=0,missingFacesBefore=0;
    std::uint64_t missingEdgesAfter=0,missingFacesAfter=0,patches=0,seedCells=0,work=0;
    double seconds=0;
    std::uint64_t edgeWalks=0,facetWalks=0,walkFallbacks=0;
    std::uint64_t steinerAttempts=0,steinerWork=0,steinerNodes=0,retainedSteinerNodes=0;
    std::uint64_t piercingPatches=0;
    bool initialTopologyUsed=false;
    const char* initialTopologyDiagnostic="not_available";
    std::uint64_t initialTopologyCells=0,initialTopologyImportWork=0,initialTopologyPeakModeledBytes=0;
};
struct NativeCadRegionRecoveryReport {
    std::uint64_t attempts=0,states=0,backtracks=0,candidates=0,pairCacheHits=0,patches=0,boundHits=0;
    std::uint64_t vertexBoundHits=0,stateBoundHits=0,largestVertices=0,largestFaces=0;
    double searchSeconds=0,preparationSeconds=0,visibilitySeconds=0,publicationSeconds=0,validationSeconds=0;
    std::uint64_t localFrontTets=0,localFrontClosed=0;
    double localFrontSeconds=0;
};
struct NativeCadInitialFrontReport {
    const char* diagnostic="not_attempted";
    std::uint64_t nodes=0,candidates=0,work=0,queries=0;
    std::uint64_t wallNodes=0,voidNodes=0,voidCandidates=0;
    // Queue exhaustion is not a coverage certificate: refused candidates can
    // remain. The unresolved ratio also includes historical refusals.
    bool sourceSweepComplete=false,voidQueueExhausted=false;
    double maximumUnresolvedVoidRatio=0;
    std::uint64_t cones=0,reusedApices=0,sourceActivations=0,propagatedActivations=0;
    // Propagation ancestry is not a physical distance from the CAD boundary.
    std::uint64_t propagatedSites=0,maximumDepth=0,remainingFaces=0;
    std::uint64_t delaunayEvaluations=0,delaunayCommits=0,fullCertifications=0;
    // A zero cavityLimit selects the complete admitted numerical complex.
    std::uint64_t cavityUnavailable=0,cavityLimit=0,maximumPlannedCavity=0,maximumAcceptedCavity=0;
    std::uint64_t storageRejections=0,compactions=0,retiredCellsReclaimed=0,compactionWork=0;
    std::uint64_t auxiliaryStoredCells=0,auxiliaryLiveCells=0,auxiliaryStoredPeak=0,auxiliaryStorageCapacity=0;
    // Modeled auxiliary payload during compaction, not peak process RSS.
    std::uint64_t compactionPeakModeledBytes=0,auxiliaryWorkingBytes=0;
    double seconds=0;
    const char* budgetReason="none";
    std::uint64_t candidateLimit=0,siteLimit=0,queryLimit=0,workLimit=0;
    std::uint64_t reservedBulkCandidates=0,reservedBulkSites=0;
    std::uint64_t sourceCandidates=0,sourceDeferredFaces=0;
    bool sourcePhaseDeferred=false;
    std::uint64_t auxiliaryImportedCells=0,auxiliaryImportWork=0,auxiliaryImportModeledBytes=0;
    std::uint64_t auxiliaryInitialCapacity=0,auxiliaryHardCapacity=0;
    std::uint64_t auxiliaryGrowthAttempts=0,auxiliaryGrowths=0,auxiliaryGrowthWork=0;
    std::uint64_t maximumCombinedModeledBytes=0,maximumHeapCapacity=0,resourceRetries=0;
    const char* auxiliaryResourceDiagnostic="not_attempted";
    std::uint64_t geometryRefusalCacheHits=0,cachedGeometryRefusals=0;
};
struct NativeCadInitialTopologyQualityReport {
    bool attempted=false,completed=false,changed=false,cancelled=false,budgetLimited=false;
    bool summaryAvailable=false,candidateListTruncated=false;
    const char* diagnostic="not_attempted";
    std::uint64_t work=0,maximumModeledBytes=0,candidates=0,proposals=0,certifications=0;
    std::uint64_t stagedEdits=0,edits=0,flips23=0,flips32=0,flips44=0;
    std::uint64_t liveBefore=0,liveAfter=0,severeBefore=0,severeAfter=0;
    std::uint64_t protectedSourceFaces=0,protectedSourceEdges=0;
    double minimumMeanRatioBefore=1,minimumMeanRatioAfter=1;
    double minimumScaledJacobianBefore=1,minimumScaledJacobianAfter=1;
    double seconds=0;
    double searchBelowMeanRatio=1e-6;
    std::uint64_t workLimit=0;
};
struct NativeCadSourceConstructionWitness {
    const char* stage="not_reported";
    const char* detail="not_reported";
    int edgeTag=-1,faceTag=-1;
    std::optional<std::uint64_t> useId;
    bool geometricRefusal=false;
};
struct NativeCadSourceShellAttempt {
    double chordErrorLimit=0,seconds=0;
    double signedVolume=0,cadVolume=0,relativeVolumeError=0;
    std::uint64_t work=0,nodes=0,triangles=0;
    bool published=false,sourceModelUnmodified=false,volumeCriterionMet=false;
    bool constructionBudgetLimited=false,constructionStalled=false;
    const char* surfaceDiagnostic="not_attempted";
    const char* provider="not_reported";
    // Configuration evidence, not achieved spacing or a point-count certificate.
    bool globalSizingApplied=false;
    unsigned curvaturePointsPerTurn=0;
    double angularDeflection=0;
    NativeCadSourceConstructionWitness witness;
};
struct NativeCadSourceShellReport {
    // Source geometry retries are bounded independently of volume recovery.
    // Queue exhaustion is not a requested element-size certificate.
    std::array<NativeCadSourceShellAttempt,3> attempts{};
    unsigned count=0;
    double requestedGlobalSize=0,requestedChordError=0,maximumRelativeVolumeError=0;
    bool accepted=false;
    const char* diagnostic="not_attempted";
    bool nativeFallbackRequested=false,nativeFallbackAttempted=false;
    const char* nativeFallbackDiagnostic="disabled";
};
// Single-body AUTO may retry legacy after a bounded source-stage refusal,
// before initial point construction. Uncontrolled Normal multi-part AUTO
// can also retain a failed front before one guarded legacy retry. This is
// the discarded trial, separate from the published route's initial front.
struct NativeCadInitialFrontFallbackReport {
    bool attempted=false;
    NativeCadMeshReason reason=NativeCadMeshReason::None;
    std::uint64_t work=0;
    double seconds=0;
    const char* surfaceDiagnostic="not_attempted";
    const char* volumeDiagnostic="not_attempted";
    NativeCadSourceShellReport sourceShell;
};
struct NativeCadRouteAttemptReport {
    const char* method="not_attempted";
    NativeCadMeshReason reason=NativeCadMeshReason::None;
    bool published=false;
    std::uint64_t work=0;
    double seconds=0;
};
struct NativeCadRouteResolutionReport {
    const char* reason="not_reported";
    const char* candidate="none";
    bool preflightAttempted=false,wholeBodyVerified=false;
    std::uint64_t preflightWork=0;
    double preflightSeconds=0;
    std::optional<double> fallbackChordLimit;
    std::array<NativeCadRouteAttemptReport,2> attempts{};
    unsigned count=0;
};
// Cumulative material-admission batches within one live volume attempt. The
// completion/source flags describe the last batch, not publication of a mesh.
struct NativeCadMaterialClassificationReport {
    const char* diagnostic="not_attempted";
    bool completed=false,sourceModelUnmodified=false;
    // A native query is counted when started; interruption need not yield a
    // decision. These attempts never stand in for the completed batch mask.
    std::uint64_t attempts=0,requestedPoints=0,attemptedPoints=0;
    std::uint64_t nativeQueries=0,strictReuses=0,rayAdmissions=0,work=0;
    double seconds=0;
    // Only explicitly profiled strict serial batches contribute these fields.
    // Zero profiledBatches means unavailable, not zero classification cost.
    std::uint64_t profiledBatches=0,rayWitnessAttempts=0;
    double classifierSetupSeconds=0,nativeQuerySeconds=0,rayWitnessSeconds=0;
};
// Experimental recovery ordering evidence. Source-PLC inclusion never admits
// an original-CAD rejection. Counts accumulate across private route attempts;
// an interrupted attempt may contain only partial diagnostic measurements.
struct NativeCadSourceMismatchPriorityReport {
    const char* diagnostic="not_attempted";
    const char* fallbackReason="not_attempted";
    std::uint64_t attempts=0,completedAttempts=0,requestedPoints=0;
    std::uint64_t insideSourceOutsideCad=0,outsideSource=0,ambiguousSource=0;
    std::uint64_t prioritizedFacets=0,insideWitnessesNearMissingFacets=0;
    std::uint64_t work=0,maximumWorkingPayloadEstimate=0;
    std::uint64_t sourceGeneration=0,sourceMeshHash=0;
    bool priorityApplied=false;
    double seconds=0;
};
struct NativeCadBodyMeshReport {
    int solidTag = -1;
    const char* route = "unsupported";
    NativeCadMeshReason reason = NativeCadMeshReason::UnsupportedBody;
    double cadVolume = 0, meshVolume = 0, relativeVolumeError = 0;
    double maximumBoundaryNodeError = 0, maximumBoundarySampleError = 0;
    int startFaceTag = -1, endFaceTag = -1;
    algo::PrimitiveMeshingResult primitive;
    algo::GeometryMeshingOutcome planar;
    bool optimizationAttempted = false;
    algo::MixedMeshOptimizationOutcome optimization;
    const char* nativeSurfaceReason = "not_attempted";
    const char* nativeVolumeReason = "not_attempted";
    double surfaceDeflection = 0, surfaceChordErrorLimit = 0;
    NativeCadSourceShellReport sourceShell;
    NativeCadMaterialClassificationReport materialClassification;
    NativeCadSourceMismatchPriorityReport sourceMismatchPriority;
    // An absent value means this route has not supplied construction evidence.
    // True records use of a volume size field, not verification of every cell.
    std::optional<bool> volumeSizingApplied;
    const char* volumeSizingDiagnostic="not_reported";
    std::uint64_t surfaceTriangles = 0, surfaceAttempts = 0, rebuiltSurfaceCharts = 0;
    std::uint64_t nativeSurfaceWork = 0, nativeVolumeWork = 0;
    double nativeSurfaceSeconds = 0, nativeVolumeSeconds = 0;
    double nativePartitionSeconds = 0, nativeAssemblySeconds = 0;
    std::uint64_t exactSurfacePairs = 0, kernelTrials = 0;
    const char* partitionReason = "not_attempted";
    const char* partitionPolicy = "none";
    std::uint64_t partitionAttempts = 0, partitionPieces = 0, partitionInterfaces = 0;
    std::uint64_t partitionWork = 0;
    double partitionVolumeChange = 0;
    double partitionVolumeAdmission = 1e-6;
    bool sampledPartitionTrimWitness = false;
    std::uint64_t partitionTrimSamples = 0;
    double maximumPartitionTrimDrift = 0, maximumPartitionTrimToleranceRatio = 0;
    std::uint64_t surfaceDiagonalRepairs = 0;
    std::vector<NativeCadPartitionAttempt> partitionAttemptReports;
    bool boundaryEmbeddingCertified = false;
    const char* contactCylinderReason = "not_attempted";
    int contactSourceEdgeTag = -1, contactTargetEdgeTag = -1;
    std::uint64_t contactRingNodes = 0, contactCylinderWork = 0;
    std::uint64_t boxContactVerticesEvaluated = 0, boxContactWork = 0;
    double boxContactMaximumEndpointCorrection = 0;
    unsigned nativeCentroidThreads = 0;
    std::uint64_t nativeCentroidsChecked = 0, nativeCentroidReferenceQueries = 0;
    std::uint64_t nativeCentroidExactReuses = 0;
    std::uint64_t nativeCentroidGeometryPairs = 0;
    bool nativeCentroidSourceGuard = false;
    double centroidOutsideDistanceLimit = 0, maximumOutsideCentroidDistance = 0;
    std::uint64_t centroidsWithinChordTolerance = 0, nativeCentroidDistanceQueries = 0;
    std::uint64_t nativeChordSamplesChecked = 0, chordSamplesWithinTolerance = 0;
    double maximumOutsideChordSampleDistance = 0;
    std::uint64_t nativeCandidateQueries = 0;
    std::uint64_t nativeCandidateRayAdmissions = 0, nativeFinalRayAdmissions = 0;
    std::uint64_t recoveryRemovedInteriorNodes=0,recoveryRemovedIncidentCells=0,recoveryRemovedCentroidCells=0;
    std::uint64_t recoveryRandomRetries=0,recoveryRandomCandidateProposals=0;
    std::uint64_t recoverySteinerNodes=0;
    std::uint64_t recoveryMergedComponents=0,recoverySourceBlockedCandidates=0,recoveryOtherComponentBlockedCandidates=0,recoverySourceFacetNodes=0;
    std::uint64_t liveMaterialInitialRejects=0,liveMaterialRemainingRejects=0,liveMaterialPatches=0;
    std::uint64_t localDelaunaySupportNodes=0,localDelaunayRefinementNodes=0;
    unsigned localDelaunayQualityEdits=0;
    NativeCadInitialFrontReport initialFront;
    NativeCadInitialTopologyQualityReport initialTopologyQuality;
    NativeCadInitialFrontFallbackReport initialFrontFallback;
    const char* interiorFrontDiagnostic="not_attempted";
    std::uint64_t interiorFrontNodes=0,interiorFrontCandidates=0;
    std::uint64_t interiorFrontWork=0,interiorFrontQueries=0;
    std::uint64_t interiorFrontRemainingOversized=0;
    std::uint64_t interiorFrontActivations=0,interiorVoidActivations=0,interiorFrontAdvances=0;
    double interiorFrontSeconds=0;
    bool localDelaunayUsedGlobalFallback=false;
    NativeCadRegionRecoveryReport recoveryRegions;
    NativeCadConstraintRecoveryReport constraintRecovery;
    std::vector<CadSurfaceSpacingReport> surfaceSpacing;
    bool planarAttempted=false;
    std::optional<double> effectiveMaximumChordError;
    const char* chordErrorPolicy="not_reported";
    NativeCadRouteResolutionReport routeResolution;
};
struct NativeCadMeshOutcome {
    bool published = false, destinationUnchanged = true;
    NativeCadMeshReason reason = NativeCadMeshReason::InternalError;
    NativeCadAssemblyPolicy assemblyPolicy=NativeCadAssemblyPolicy::ConformingDomain;
    std::uint64_t validatedParts=0;
    std::vector<NativeCadInterference> interferences;
    std::uint64_t validationWork = 0;
    double preparationSeconds = 0, generationSeconds = 0, validationSeconds = 0;
    double auditSeconds = 0, qualitySeconds = 0, totalSeconds = 0;
    std::vector<NativeCadBodyMeshReport> bodies;
    NativeCadInitialFrontFallbackReport initialFrontFallback;
    mesh::LinearVolumeAudit audit;
    algo::MixedMeshQualityReport quality;
    bool qualityWarning = false;
    std::uint64_t separateContactPairs = 0, contactTrianglePairs = 0;
    std::uint64_t separateBoxFaceContactPairs = 0, matchedBoxFaceTriangles = 0;
    std::uint64_t nativeProjectionQueries = 0, nativeProjectionInitializations = 0;
    std::uint64_t nativeTrimmedProjectionQueries = 0, nativeTrimmedProjectionMemoHits = 0;
    std::uint64_t nativeProjectionReuses = 0;
    const char* nativeContactDiagnostic = "not_attempted";
    bool nativeContactSourceGuard = false;
    bool nativeContactSourceAuditAvailable = false;
    std::uint64_t nativeContactPieces = 0, nativeContactPlanes = 0;
    std::uint64_t matchedNativeContactTriangles = 0;
    std::uint64_t nativeContactMaximumWorkingPayloadEstimate = 0;
    NativeCadRouteResolutionReport routeResolution;
};
// Native BRep remains the geometric authority. Existing planar/analytic routes
// run first. The optional surface fallback retains native trim/edge identities
// and accepts a Tet4 fan only after exact embedding, closed vertex-manifold
// topology and a strictly certified star kernel. A bounded partition fallback
// splits a copied native body and joins only native shared topology with exact
// matching, oppositely wound cut triangles. Original CAD identities survive.
// Each assembled Tet4 centroid is checked inside its original body. The
// explicit PLC chord-band option additionally measures OUT centroid distance
// against the original solid; it never changes interior-node admission. These
// native geometry samples and
// volume error are measured; they are not a Hausdorff certificate. Unsupported
// controls and unresolved body interfaces refuse transactionally. The opt-in
// coarse box full-face route evaluates original CAD vertex identities and
// requires unique exact opposite contact triangles covering native area;
// it retains distinct body nodes rather than bonding a material interface.
// OCCT recognition/projection allocations are not a process RSS limit.
// geometry controls are applied per body. Analytic routes refuse spatial fields,
// strict edge-size audits and body strategy overrides. Their work admission is
// a node+cell construction floor plus bounded profile work, not a full kernel
// operation counter. Native fallback source/contact/proposal operations share
// the validation work admission; OCCT internals remain outside it. Source-cloud
// limits and recovery flags concern planar routes. Fallback global size controls
// the surface request; volume edge sizes are not certified by the star route.
NativeCadMeshOutcome meshNativeCad(const OcctModel&, const NativeCadMeshOptions&,
                                   mesh::Mesh& destination) noexcept;
const char* nativeCadMeshReasonName(NativeCadMeshReason) noexcept;
} // namespace paacs::geo
