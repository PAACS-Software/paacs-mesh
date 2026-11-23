// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/OcctMesher.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "NumberArgument.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <stdexcept>

namespace {
struct JsonFinite { double value; };
std::ostream& operator<<(std::ostream& stream,JsonFinite number) {
    if(std::isfinite(number.value))stream<<number.value;
    else stream<<"null";
    return stream;
}
template<class T>T number(const char* text) {
    T value{};
    if (!paacs::cli::parseNumber(text, value)) throw std::invalid_argument("invalid number");
    return value;
}
void traceNativeStage(const paacs::geo::NativeCadStageEvent& event) noexcept {
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
    std::fprintf(stderr, "PAACS_STAGE %.9f %s %d %s %s\n", seconds,
        event.beginning ? "begin" : "end", event.entityTag, event.stage, event.diagnostic);
    std::fflush(stderr);
}
void writeInitialFrontFallback(const paacs::geo::NativeCadInitialFrontFallbackReport& report){
    std::cout<<"{\"attempted\":"<<(report.attempted?"true":"false")
        <<",\"reason\":\""<<paacs::geo::nativeCadMeshReasonName(report.reason)<<"\""
        <<",\"work\":"<<report.work<<",\"seconds\":"<<report.seconds
        <<",\"surface_diagnostic\":\""<<report.surfaceDiagnostic<<"\"}";
}
}
int runNativeCadCli(int argc,char** argv) {
    using namespace paacs;
    geo::NativeCadImportOptions io;geo::NativeCadImportOutcome cad;
    geo::NativeCadMeshOptions options;geo::NativeCadMeshOutcome result;mesh::Mesh mesh;
    bool report=false,archivePublished=false;double sizeFactor=0,importSeconds=0,growth=1.3;
    const auto started=std::chrono::steady_clock::now();
    const auto receipt=[&](const char* status,const char* reason) {
        if(!report)return;
        const auto& construction=options.construction;
        const auto cap=[](std::uint64_t value){return value?std::to_string(value):std::string("null");};
        const bool interiorFrontAttempted=std::any_of(result.bodies.begin(),result.bodies.end(),[](const auto& body){
            return body.initialFront.diagnostic&&std::string_view(body.initialFront.diagnostic)!="not_attempted";
        });
        double volumeReference=cad.cadVolume;
        bool auditedReference=result.published&&result.bodies.size()==cad.solids;
        double bodyVolumeSum=0;
        if(auditedReference)for(const auto& body:result.bodies){
            if(!(body.cadVolume>0)||!std::isfinite(body.cadVolume)){auditedReference=false;break;}
            bodyVolumeSum+=body.cadVolume;
        }
        if(auditedReference&&std::isfinite(bodyVolumeSum))volumeReference=bodyVolumeSum;
        else auditedReference=false;
        std::cout<<std::setprecision(17)<<"{\"schema\":\"paacs-native-cad-report-v1\",\"status\":\""<<status
            <<"\",\"reason\":\""<<reason<<"\",\"mesh_published\":"<<(result.published?"true":"false")
            <<",\"archive_published\":"<<(archivePublished?"true":"false")<<",\"solids\":"<<cad.solids
            <<",\"mesh_scope\":\""<<(options.assemblyPolicy==geo::NativeCadAssemblyPolicy::IndependentParts?"independent_parts":"conforming_domain")<<"\""
            <<",\"interbody_nonoverlap_certified\":"<<(result.published&&options.assemblyPolicy==geo::NativeCadAssemblyPolicy::ConformingDomain?"true":"false")
            <<",\"validated_parts\":"<<result.validatedParts
            <<",\"construction\":{\"effort\":\""<<(construction.effort==geo::CadConstructionEffort::Legacy?"legacy":
                construction.effort==geo::CadConstructionEffort::Automatic?"auto":"thorough")
            <<"\",\"max_surface_points\":"<<cap(construction.surfacePointsPerFace)
            <<",\"max_surface_chord_points\":"<<cap(construction.surfaceChordPointsPerFace)
            <<",\"max_interior_candidates\":"<<cap(construction.interiorCandidates)
            <<",\"max_interior_sites\":"<<cap(construction.interiorSites)<<"}"
            <<",\"faces\":"<<cad.faces<<",\"curved_faces\":"<<cad.curvedFaces<<",\"edges\":"<<cad.edges
            <<",\"vertices\":"<<cad.vertices<<",\"nodes\":"<<mesh.nodeCount()<<",\"cells\":"<<mesh.cellCount()
            <<",\"auxiliary_faces\":"<<cad.auxiliaryFaces<<",\"auxiliary_edges\":"<<cad.auxiliaryEdges
            <<",\"auxiliary_vertices\":"<<cad.auxiliaryVertices<<",\"auxiliary_geometry_retained\":"<<(cad.auxiliaryGeometryRetained?"true":"false")
            <<",\"cad_volume\":"<<volumeReference<<",\"legacy_import_cad_volume\":"<<cad.cadVolume
            <<",\"cad_volume_reference\":\""<<(auditedReference?(options.assemblyPolicy==geo::NativeCadAssemblyPolicy::IndependentParts?"audited_part_sum":"audited_body_sum"):"legacy_import")<<"\""
            <<",\"mesh_volume\":"<<result.quality.totalVolume
            <<",\"minimum_quality\":"<<result.quality.minimumQuality<<",\"import_seconds\":"<<importSeconds
            <<",\"quality_percentile_1\":"<<result.quality.percentile1
            <<",\"quality_percentile_5\":"<<result.quality.percentile5
            <<",\"quality_median\":"<<result.quality.medianQuality
            <<",\"quality_mean\":"<<result.quality.meanQuality
            <<",\"quality_warning\":"<<(result.qualityWarning?"true":"false")
            <<",\"audit_valid\":"<<(result.audit.valid?"true":"false")
            <<",\"boundary_missing\":"<<result.audit.missingBoundary
            <<",\"boundary_unexpected\":"<<result.audit.unexpectedBoundary
            <<",\"nonmanifold_faces\":"<<result.audit.nonmanifoldFaces
            <<",\"inconsistent_face_orientation\":"<<result.audit.inconsistentOrientation
            <<",\"separate_contact_pairs\":"<<result.separateContactPairs
            <<",\"contact_triangle_pairs\":"<<result.contactTrianglePairs
            <<",\"separate_box_face_contact_pairs\":"<<result.separateBoxFaceContactPairs
            <<",\"matched_box_face_triangles\":"<<result.matchedBoxFaceTriangles
            <<",\"native_validation_threads\":"<<options.maximumNativeValidationThreads
            <<",\"allow_centroids_within_chord_tolerance\":"<<(options.allowCentroidsWithinChordTolerance?"true":"false")
            <<",\"reuse_native_queries\":"<<(options.reuseNativeQueries?"true":"false")
            <<",\"native_projection_queries\":"<<result.nativeProjectionQueries
            <<",\"native_trimmed_projection_queries\":"<<result.nativeTrimmedProjectionQueries
            <<",\"native_trimmed_projection_memo_hits\":"<<result.nativeTrimmedProjectionMemoHits
            <<",\"native_projection_initializations\":"<<result.nativeProjectionInitializations
            <<",\"native_projection_reuses\":"<<result.nativeProjectionReuses
            <<",\"defer_discarded_native_facet_checks\":"<<(options.deferDiscardedNativeFacetChecks?"true":"false")
            <<",\"separate_matched_native_faces\":"<<(options.separateMatchedNativeFaceContacts?"true":"false")
            <<",\"native_kernel_admission\":"<<(options.useNativeKernelAdmission?"true":"false")
            <<",\"recovery_mode\":\""<<(options.recoveryMode==geo::NativeCadRecoveryMode::Advanced?"advanced":"standard")<<"\""
            <<",\"longitudinal_recovery\":"<<(options.useLongitudinalRecovery?"true":"false")
            <<",\"local_partition_recovery\":"<<(options.useLocalPartitionRecovery?"true":"false")
            <<",\"adaptive_tet_recovery\":"<<(options.useAdaptiveTetRecovery?"true":"false")
            <<",\"curvature_points_per_turn\":"<<options.curvaturePointsPerTurn
            <<",\"native_spatial_sizing\":"<<((options.useSpatialSizing||!options.sizeBoxes.empty()||options.advancingInteriorFrontEnabled(std::size_t(cad.solids)))?"true":"false")
            <<",\"face_size_controls\":"<<options.faceSizes.size()<<",\"size_box_controls\":"<<options.sizeBoxes.size()
            <<",\"live_material_recovery\":"<<(options.useLiveMaterialRecovery?"true":"false")
            <<",\"trimmed_ray_witnesses\":"<<(options.useTrimmedRayWitnesses?"true":"false")
            <<",\"local_delaunay_refinement\":"<<(options.useLocalDelaunayRefinement?"true":"false")
            <<",\"improve_surface_spacing\":"<<(options.improveNativeSurfaceSpacing?"true":"false")
            <<",\"improve_tet_quality\":"<<(options.improveNativeTetQuality?"true":"false")
            <<",\"advancing_interior_front\":"<<(interiorFrontAttempted?"true":"false")
            <<",\"advancing_interior_front_selected\":"<<(options.advancingInteriorFrontEnabled(std::size_t(cad.solids))?"true":"false")
            <<",\"advancing_interior_front_policy\":\""<<(!options.useAdvancingInteriorFront.has_value()?"auto":*options.useAdvancingInteriorFront?"enabled":"disabled")<<"\""
            <<",\"regular_surface_front\":"<<(options.useRegularSurfaceFront?"true":"false")
            <<",\"improve_regular_quality\":"<<(options.regularMeshQualityEnabled()?"true":"false")
            <<",\"maximum_edge_recovery_steiner_candidates\":"<<options.maximumEdgeRecoverySteinerCandidates
            <<",\"maximum_recovery_patch_states\":"<<options.maximumRecoveryPatchStates
            <<",\"recovery_seed\":"<<options.recoverySeed<<",\"maximum_recovery_random_retries\":"<<options.maximumRecoveryRandomRetries
            <<",\"native_contact_diagnostic\":\""<<result.nativeContactDiagnostic<<"\""
            <<",\"native_contact_source_guard\":"<<(result.nativeContactSourceGuard?"true":"false")
            <<",\"native_contact_source_audit_available\":"<<(result.nativeContactSourceAuditAvailable?"true":"false")
            <<",\"native_contact_maximum_working_payload_estimate\":"<<result.nativeContactMaximumWorkingPayloadEstimate
            <<",\"native_contact_pieces\":"<<result.nativeContactPieces
            <<",\"native_contact_planes\":"<<result.nativeContactPlanes
            <<",\"matched_native_contact_triangles\":"<<result.matchedNativeContactTriangles
            <<",\"preparation_seconds\":"<<result.preparationSeconds
            <<",\"generation_seconds\":"<<result.generationSeconds<<",\"validation_seconds\":"<<result.validationSeconds
            <<",\"audit_seconds\":"<<result.auditSeconds<<",\"quality_seconds\":"<<result.qualitySeconds
            <<",\"validation_work\":"<<result.validationWork<<",\"millimetres_per_model_unit\":"<<cad.millimetresPerModelUnit
            <<",\"total_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<",\"interferences\":[";
        for(std::size_t i=0;i<result.interferences.size();++i){if(i)std::cout<<',';const auto& p=result.interferences[i];std::cout<<"{\"first_solid_tag\":"<<p.firstSolidTag<<",\"second_solid_tag\":"<<p.secondSolidTag<<",\"state\":\""<<(p.state==geo::NativeCadInterferenceState::MaterialOverlap?"material_overlap":p.state==geo::NativeCadInterferenceState::NoMaterialOverlap?"no_material_overlap":"unresolved")<<"\",\"native_overlap_volume\":"<<p.nativeOverlapVolume<<'}';}
        std::cout<<"],\"initial_front_fallback\":";
        writeInitialFrontFallback(result.initialFrontFallback);
        std::cout<<",\"bodies\":[";
        for(std::size_t i=0;i<result.bodies.size();++i) {if(i)std::cout<<',';const auto& b=result.bodies[i];
            std::cout<<"{\"solid_tag\":"<<b.solidTag<<",\"route\":\""<<b.route<<"\",\"reason\":\""<<geo::nativeCadMeshReasonName(b.reason)
                <<"\",\"start_face_tag\":"<<b.startFaceTag<<",\"end_face_tag\":"<<b.endFaceTag
                <<",\"cad_volume\":"<<b.cadVolume<<",\"mesh_volume\":"<<b.meshVolume
                <<",\"relative_volume_error\":"<<b.relativeVolumeError<<",\"boundary_node_error\":"<<b.maximumBoundaryNodeError
                <<",\"boundary_sample_error\":"<<b.maximumBoundarySampleError<<",\"primitive_reason\":\""<<algo::primitiveReasonName(b.primitive.reason)
                <<"\",\"planar_reason\":\""<<algo::geometryReasonName(b.planar.reason)<<"\"";
            std::cout<<",\"initial_front_fallback\":";
            writeInitialFrontFallback(b.initialFrontFallback);
            std::cout<<",\"native_surface_reason\":\""<<b.nativeSurfaceReason<<"\""
                <<",\"native_volume_reason\":\""<<b.nativeVolumeReason<<"\""
                <<",\"surface_deflection\":"<<b.surfaceDeflection
                <<",\"surface_chord_error_limit\":"<<b.surfaceChordErrorLimit
                <<",\"surface_triangles\":"<<b.surfaceTriangles<<",\"surface_attempts\":"<<b.surfaceAttempts
                <<",\"rebuilt_surface_charts\":"<<b.rebuiltSurfaceCharts
                <<",\"native_surface_work\":"<<b.nativeSurfaceWork<<",\"native_volume_work\":"<<b.nativeVolumeWork
                <<",\"native_surface_seconds\":"<<b.nativeSurfaceSeconds<<",\"native_volume_seconds\":"<<b.nativeVolumeSeconds
                <<",\"native_partition_seconds\":"<<b.nativePartitionSeconds<<",\"native_assembly_seconds\":"<<b.nativeAssemblySeconds
                <<",\"exact_surface_pairs\":"<<b.exactSurfacePairs<<",\"kernel_trials\":"<<b.kernelTrials;
            const auto& region=b.recoveryRegions;
            std::cout<<",\"region_recovery\":{\"attempts\":"<<region.attempts<<",\"states\":"<<region.states<<",\"backtracks\":"<<region.backtracks<<",\"candidates\":"<<region.candidates<<",\"pair_cache_hits\":"<<region.pairCacheHits<<",\"patches\":"<<region.patches<<",\"bound_hits\":"<<region.boundHits
                <<",\"vertex_bound_hits\":"<<region.vertexBoundHits<<",\"state_bound_hits\":"<<region.stateBoundHits<<",\"largest_vertices\":"<<region.largestVertices<<",\"largest_faces\":"<<region.largestFaces
                <<",\"search_seconds\":"<<region.searchSeconds<<",\"preparation_seconds\":"<<region.preparationSeconds<<",\"visibility_seconds\":"<<region.visibilitySeconds<<",\"publication_seconds\":"<<region.publicationSeconds<<",\"validation_seconds\":"<<region.validationSeconds
                <<",\"local_front_tets\":"<<region.localFrontTets<<",\"local_front_closed\":"<<region.localFrontClosed
                <<",\"local_front_seconds\":"<<region.localFrontSeconds<<'}';
            std::cout<<",\"partition_reason\":\""<<b.partitionReason<<"\",\"partition_policy\":\""<<b.partitionPolicy<<"\""
                <<",\"partition_attempts\":"<<b.partitionAttempts<<",\"partition_pieces\":"<<b.partitionPieces
                <<",\"partition_interfaces\":"<<b.partitionInterfaces<<",\"partition_work\":"<<b.partitionWork
                <<",\"partition_volume_change\":"<<b.partitionVolumeChange
                <<",\"partition_volume_admission\":"<<b.partitionVolumeAdmission
                <<",\"sampled_partition_trim_witness\":"<<(b.sampledPartitionTrimWitness?"true":"false")
                <<",\"partition_trim_samples\":"<<b.partitionTrimSamples
                <<",\"maximum_partition_trim_drift\":"<<b.maximumPartitionTrimDrift
                <<",\"maximum_partition_trim_tolerance_ratio\":"<<b.maximumPartitionTrimToleranceRatio
                <<",\"surface_diagonal_repairs\":"<<b.surfaceDiagonalRepairs
                <<",\"boundary_embedding_certified\":"<<(b.boundaryEmbeddingCertified?"true":"false")
                <<",\"native_centroid_threads\":"<<b.nativeCentroidThreads
                <<",\"native_centroids_checked\":"<<b.nativeCentroidsChecked
                <<",\"native_centroid_exact_reuses\":"<<b.nativeCentroidExactReuses
                <<",\"native_centroid_reference_queries\":"<<b.nativeCentroidReferenceQueries
                <<",\"native_centroid_geometry_pairs\":"<<b.nativeCentroidGeometryPairs
                <<",\"native_centroid_source_guard\":"<<(b.nativeCentroidSourceGuard?"true":"false")
                <<",\"centroid_outside_distance_limit\":"<<b.centroidOutsideDistanceLimit
                <<",\"centroids_within_chord_tolerance\":"<<b.centroidsWithinChordTolerance
                <<",\"native_centroid_distance_queries\":"<<b.nativeCentroidDistanceQueries
                <<",\"maximum_outside_centroid_distance\":"<<b.maximumOutsideCentroidDistance
                <<",\"native_chord_samples_checked\":"<<b.nativeChordSamplesChecked
                <<",\"chord_samples_within_tolerance\":"<<b.chordSamplesWithinTolerance
                <<",\"maximum_outside_chord_sample_distance\":"<<b.maximumOutsideChordSampleDistance
                <<",\"native_candidate_queries\":"<<b.nativeCandidateQueries
                <<",\"native_candidate_ray_admissions\":"<<b.nativeCandidateRayAdmissions
                <<",\"native_final_ray_admissions\":"<<b.nativeFinalRayAdmissions
                <<",\"recovery_removed_interior_nodes\":"<<b.recoveryRemovedInteriorNodes
                <<",\"recovery_removed_incident_cells\":"<<b.recoveryRemovedIncidentCells
                <<",\"recovery_removed_centroid_cells\":"<<b.recoveryRemovedCentroidCells
                <<",\"recovery_steiner_nodes\":"<<b.recoverySteinerNodes
                <<",\"recovery_random_retries\":"<<b.recoveryRandomRetries
                <<",\"recovery_random_candidate_proposals\":"<<b.recoveryRandomCandidateProposals
                <<",\"recovery_merged_components\":"<<b.recoveryMergedComponents
                <<",\"recovery_source_blocked_candidates\":"<<b.recoverySourceBlockedCandidates
                <<",\"recovery_other_component_blocked_candidates\":"<<b.recoveryOtherComponentBlockedCandidates
                <<",\"recovery_source_facet_nodes\":"<<b.recoverySourceFacetNodes
                <<",\"live_material_initial_rejects\":"<<b.liveMaterialInitialRejects
                <<",\"live_material_remaining_rejects\":"<<b.liveMaterialRemainingRejects
                <<",\"live_material_patches\":"<<b.liveMaterialPatches
                <<",\"local_delaunay_support_nodes\":"<<b.localDelaunaySupportNodes
                <<",\"local_delaunay_refinement_nodes\":"<<b.localDelaunayRefinementNodes
                <<",\"local_delaunay_quality_edits\":"<<b.localDelaunayQualityEdits
                <<",\"initial_front\":{\"diagnostic\":\""<<b.initialFront.diagnostic
                <<"\",\"nodes\":"<<b.initialFront.nodes<<",\"candidates\":"<<b.initialFront.candidates
                <<",\"wall_nodes\":"<<b.initialFront.wallNodes<<",\"void_nodes\":"<<b.initialFront.voidNodes
                <<",\"void_candidates\":"<<b.initialFront.voidCandidates
                <<",\"budget_reason\":\""<<b.initialFront.budgetReason<<"\""
                <<",\"candidate_limit\":"<<b.initialFront.candidateLimit<<",\"site_limit\":"<<b.initialFront.siteLimit
                <<",\"query_limit\":"<<b.initialFront.queryLimit<<",\"work_limit\":"<<b.initialFront.workLimit
                <<",\"reserved_bulk_candidates\":"<<b.initialFront.reservedBulkCandidates
                <<",\"reserved_bulk_sites\":"<<b.initialFront.reservedBulkSites
                <<",\"source_candidates\":"<<b.initialFront.sourceCandidates
                <<",\"source_deferred_faces\":"<<b.initialFront.sourceDeferredFaces
                <<",\"source_phase_deferred\":"<<(b.initialFront.sourcePhaseDeferred?"true":"false")
                <<",\"source_sweep_complete\":"<<(b.initialFront.sourceSweepComplete?"true":"false")
                <<",\"void_queue_exhausted\":"<<(b.initialFront.voidQueueExhausted?"true":"false")
                <<",\"maximum_unresolved_void_ratio\":"<<JsonFinite{b.initialFront.maximumUnresolvedVoidRatio}
                <<",\"work\":"<<b.initialFront.work<<",\"queries\":"<<b.initialFront.queries
                <<",\"cones\":"<<b.initialFront.cones<<",\"reused_apices\":"<<b.initialFront.reusedApices
                <<",\"source_activations\":"<<b.initialFront.sourceActivations
                <<",\"propagated_activations\":"<<b.initialFront.propagatedActivations
                <<",\"propagated_sites\":"<<b.initialFront.propagatedSites
                <<",\"maximum_depth\":"<<b.initialFront.maximumDepth
                <<",\"remaining_faces\":"<<b.initialFront.remainingFaces
                <<",\"delaunay_evaluations\":"<<b.initialFront.delaunayEvaluations
                <<",\"delaunay_commits\":"<<b.initialFront.delaunayCommits
                <<",\"full_certifications\":"<<b.initialFront.fullCertifications
                <<",\"cavity_unavailable\":"<<b.initialFront.cavityUnavailable
                <<",\"cavity_limit\":"<<b.initialFront.cavityLimit
                <<",\"maximum_planned_cavity\":"<<b.initialFront.maximumPlannedCavity
                <<",\"maximum_accepted_cavity\":"<<b.initialFront.maximumAcceptedCavity
                <<",\"storage_rejections\":"<<b.initialFront.storageRejections
                <<",\"compactions\":"<<b.initialFront.compactions
                <<",\"retired_cells_reclaimed\":"<<b.initialFront.retiredCellsReclaimed
                <<",\"compaction_work\":"<<b.initialFront.compactionWork
                <<",\"auxiliary_stored_cells\":"<<b.initialFront.auxiliaryStoredCells
                <<",\"auxiliary_live_cells\":"<<b.initialFront.auxiliaryLiveCells
                <<",\"auxiliary_stored_peak\":"<<b.initialFront.auxiliaryStoredPeak
                <<",\"auxiliary_storage_capacity\":"<<b.initialFront.auxiliaryStorageCapacity
                <<",\"compaction_peak_modeled_bytes\":"<<b.initialFront.compactionPeakModeledBytes
                <<",\"auxiliary_working_bytes\":"<<b.initialFront.auxiliaryWorkingBytes
                <<",\"seconds\":"<<b.initialFront.seconds<<"}"
                <<",\"interior_front\":{\"diagnostic\":\""<<b.interiorFrontDiagnostic
                <<"\",\"nodes\":"<<b.interiorFrontNodes<<",\"candidates\":"<<b.interiorFrontCandidates
                <<",\"work\":"<<b.interiorFrontWork<<",\"queries\":"<<b.interiorFrontQueries
                <<",\"remaining_oversized\":"<<b.interiorFrontRemainingOversized
                <<",\"front_activations\":"<<b.interiorFrontActivations
                <<",\"front_advances\":"<<b.interiorFrontAdvances
                <<",\"void_activations\":"<<b.interiorVoidActivations
                <<",\"seconds\":"<<b.interiorFrontSeconds<<"}"
                <<",\"local_delaunay_used_global_fallback\":"<<(b.localDelaunayUsedGlobalFallback?"true":"false")
                <<",\"contact_cylinder_reason\":\""<<b.contactCylinderReason<<"\""
                <<",\"contact_source_edge_tag\":"<<b.contactSourceEdgeTag
                <<",\"contact_target_edge_tag\":"<<b.contactTargetEdgeTag
                <<",\"contact_ring_nodes\":"<<b.contactRingNodes
                <<",\"contact_cylinder_work\":"<<b.contactCylinderWork;
            std::cout<<",\"surface_spacing\":[";
            for(std::size_t faceIndex=0;faceIndex<b.surfaceSpacing.size();++faceIndex){
                if(faceIndex)std::cout<<',';
                const auto& face=b.surfaceSpacing[faceIndex];
                std::cout<<"{\"face_tag\":"<<face.faceTag<<",\"point_limit\":"<<face.pointLimit
                    <<",\"chord_point_limit\":"<<face.chordPointLimit<<",\"chart_node_limit\":"<<face.chartNodeLimit
                    <<",\"added_points\":"<<face.addedPoints<<",\"chord_added_points\":"<<face.chordAddedPoints
                    <<",\"pending_entries\":"<<face.pendingEntries<<",\"blocked_candidates\":"<<face.blockedCandidates
                    <<",\"termination\":\""<<face.termination<<"\",\"queue_exhausted\":"<<(face.queueExhausted?"true":"false")
                    <<",\"budget_limited\":"<<(face.budgetLimited?"true":"false")<<'}';
            }
            std::cout<<']';
            const auto& constraints=b.constraintRecovery;
            std::cout<<",\"constraint_recovery\":{\"diagnostic\":\""<<constraints.diagnostic
                <<"\",\"candidate_cells\":"<<constraints.candidateCells
                <<",\"missing_edges_before\":"<<constraints.missingEdgesBefore
                <<",\"missing_edges_after\":"<<constraints.missingEdgesAfter
                <<",\"missing_faces_before\":"<<constraints.missingFacesBefore
                <<",\"missing_faces_after\":"<<constraints.missingFacesAfter
                <<",\"patches\":"<<constraints.patches<<",\"seed_cells\":"<<constraints.seedCells
                <<",\"work\":"<<constraints.work<<",\"seconds\":"<<constraints.seconds
                <<",\"edge_walks\":"<<constraints.edgeWalks<<",\"facet_walks\":"<<constraints.facetWalks
                <<",\"walk_fallbacks\":"<<constraints.walkFallbacks
                <<",\"steiner_attempts\":"<<constraints.steinerAttempts
                <<",\"steiner_work\":"<<constraints.steinerWork<<",\"steiner_nodes\":"<<constraints.steinerNodes
                <<",\"retained_steiner_nodes\":"<<constraints.retainedSteinerNodes
                <<",\"piercing_edge_patches\":"<<constraints.piercingPatches<<'}';
            std::cout<<",\"box_contact_vertices_evaluated\":"<<b.boxContactVerticesEvaluated
                <<",\"box_contact_work\":"<<b.boxContactWork
                <<",\"box_contact_maximum_endpoint_correction\":"<<b.boxContactMaximumEndpointCorrection;
            std::cout<<",\"partition_attempt_reports\":[";
            for(std::size_t j=0;j<b.partitionAttemptReports.size();++j){if(j)std::cout<<',';const auto& a=b.partitionAttemptReports[j];
                std::cout<<"{\"policy\":\""<<a.policy<<"\",\"diagnostic\":\""<<a.diagnostic
                    <<"\",\"failed_piece_tag\":"<<a.failedPieceTag<<",\"pieces\":"<<a.pieces<<",\"interfaces\":"<<a.interfaces
                    <<",\"work\":"<<a.work<<",\"seconds\":"<<a.seconds
                    <<",\"surface_deflection\":"<<a.surfaceDeflection<<'}';}
            std::cout<<']';
            if(std::string_view(b.route)!="planar_core")std::cout
                <<",\"uniform_cylinder\":"<<(b.primitive.uniformCylinder?"true":"false")
                <<",\"uniform_circular\":"<<(b.primitive.uniformCircular?"true":"false")
                <<",\"profiled_volume\":"<<(b.primitive.profiledVolume?"true":"false")
                <<",\"minimum_section_nodes\":"<<b.primitive.minimumSectionNodes
                <<",\"maximum_section_nodes\":"<<b.primitive.maximumSectionNodes
                <<",\"minimum_axial_spacing\":"<<b.primitive.minimumAxialSpacing
                <<",\"maximum_axial_spacing\":"<<b.primitive.maximumAxialSpacing
                <<",\"maximum_axial_growth\":"<<b.primitive.maximumAxialGrowth
                <<",\"chord_error_bound\":"<<b.primitive.chordErrorBound;
            if(!b.planar.analysis.bodies.empty()){const auto& p=b.planar.analysis.bodies.front();
                std::cout<<",\"boundary_reason_code\":"<<p.boundaryReason<<",\"boundary_phase_code\":"<<p.boundaryPhase
                    <<",\"surface_reason_code\":"<<p.surfaceReason<<",\"surface_failure_face_id\":"<<p.surfaceFailureFacePersistentId;}
            if(b.optimizationAttempted)std::cout<<",\"optimization_trials\":"<<b.optimization.nodeTrials
                <<",\"optimization_moves\":"<<b.optimization.acceptedMoves
                <<",\"quality_before_optimization\":"<<b.optimization.before.minimumQuality
                <<",\"quality_after_optimization\":"<<b.optimization.after.minimumQuality;
            std::cout<<'}';
        }
        std::cout<<"]}\n";
    };
    for(int i=1;i<argc;++i)report|=std::string_view(argv[i])=="--report-json"||std::string_view(argv[i])=="--independent-parts";
    try {
        if(argc<4)throw std::invalid_argument("usage: paacs-cad-mesh input.step size output.paacsmesh --native-cad [--face-size tag size] [--native-sizing] [--size-box x0 y0 z0 x1 y1 z1 size] [--min-size h] [--growth g] [--report-json] [--trace-native-stages] [--independent-parts] [--local-partition-recovery] [--adaptive-tet-recovery] [--live-material-recovery] [--curvature-points-per-turn 0..64] [--trimmed-ray-witnesses] [--local-delaunay-refinement] [--advancing-interior-front|--no-advancing-interior-front] [--improve-tet-quality] [--improve-surface-spacing] [--regular-surface-front] [--no-improve-regular-quality] [--recovery-mode standard|advanced] [--recovery-seed N] [--recovery-retries 0..4] [--recovery-steiner-candidates 0..64] [--construction-effort legacy|auto|thorough] [--max-surface-points N] [--max-surface-chord-points N] [--max-interior-candidates N] [--max-interior-sites N]");
        options.geometry.targetSize=number<double>(argv[2]);
        for(int i=4;i<argc;++i) {
            const std::string_view a(argv[i]);
            if(a=="--native-cad"||a=="--report-json")continue;
            if(a=="--trace-native-stages"){options.stageObserver=traceNativeStage;continue;}
            if(a=="--planar-recovery")options.geometry.usePlanarSurfaceRecovery=true;
            else if(a=="--general-recovery")options.geometry.edgeGeneralRecovery.enabled=true;
            else if(a=="--optimize-profile")options.optimizeProfileInterior=true;
            else if(a=="--separate-contact-bodies")options.separateLowerDimensionalContacts=true;
            else if(a=="--separate-matched-box-faces")options.separateMatchedBoxFaceContacts=true;
            else if(a=="--no-native-fallback")options.useNativeSurfaceFallback=false;
            else if(a=="--no-native-partition")options.useNativePartitionFallback=false;
            else if(a=="--allow-centroids-within-chord-tolerance")options.allowCentroidsWithinChordTolerance=true;
            else if(a=="--max-native-pieces"&&i+1<argc)options.maximumNativePartitionPieces=number<std::uint64_t>(argv[++i]);
            else if(a=="--native-validation-threads"&&i+1<argc)options.maximumNativeValidationThreads=number<unsigned>(argv[++i]);
            else if(a=="--reuse-native-queries")options.reuseNativeQueries=true;
            else if(a=="--defer-discarded-native-facet-checks")options.deferDiscardedNativeFacetChecks=true;
            else if(a=="--separate-matched-native-faces")options.separateMatchedNativeFaceContacts=true;
            else if(a=="--retain-auxiliary-geometry")io.retainAuxiliaryGeometry=true;
            else if(a=="--native-kernel-admission")options.useNativeKernelAdmission=true;
            else if(a=="--longitudinal-recovery"){options.useLongitudinalRecovery=true;options.useNativeKernelAdmission=true;}
            else if(a=="--local-partition-recovery")options.useLocalPartitionRecovery=true;
            else if(a=="--adaptive-tet-recovery")options.useAdaptiveTetRecovery=true;
            else if(a=="--live-material-recovery")options.useLiveMaterialRecovery=true;
            else if(a=="--trimmed-ray-witnesses")options.useTrimmedRayWitnesses=true;
            else if(a=="--local-delaunay-refinement")options.useLocalDelaunayRefinement=true;
            else if(a=="--curvature-points-per-turn"&&i+1<argc)
                options.curvaturePointsPerTurn=number<unsigned>(argv[++i]);
            else if(a=="--improve-tet-quality")options.improveNativeTetQuality=true;
            else if(a=="--advancing-interior-front")options.useAdvancingInteriorFront=true;
            else if(a=="--no-advancing-interior-front")options.useAdvancingInteriorFront=false;
            else if(a=="--construction-effort"&&i+1<argc){
                const std::string_view effort(argv[++i]);
                if(effort=="legacy")options.construction.effort=geo::CadConstructionEffort::Legacy;
                else if(effort=="auto")options.construction.effort=geo::CadConstructionEffort::Automatic;
                else if(effort=="thorough")options.construction.effort=geo::CadConstructionEffort::Thorough;
                else throw std::invalid_argument("construction effort must be legacy, auto or thorough");
            }
            else if((a=="--max-surface-points"||a=="--max-surface-chord-points"||
                     a=="--max-interior-candidates"||a=="--max-interior-sites")&&i+1<argc){
                const auto value=number<std::uint64_t>(argv[++i]);
                if(!value)throw std::invalid_argument("construction stage ceilings must be positive");
                if(a=="--max-surface-points")options.construction.surfacePointsPerFace=value;
                else if(a=="--max-surface-chord-points")options.construction.surfaceChordPointsPerFace=value;
                else if(a=="--max-interior-candidates")options.construction.interiorCandidates=value;
                else options.construction.interiorSites=value;
            }
            else if(a=="--improve-surface-spacing")options.improveNativeSurfaceSpacing=true;
            else if(a=="--regular-surface-front")options.useRegularSurfaceFront=true;
            else if(a=="--improve-regular-quality")options.improveRegularMeshQuality=true;
            else if(a=="--no-improve-regular-quality")options.improveRegularMeshQuality=false;
            else if(a=="--independent-parts")options.assemblyPolicy=geo::NativeCadAssemblyPolicy::IndependentParts;
            else if(a=="--recovery-mode"&&i+1<argc){
                const std::string_view mode(argv[++i]);
                if(mode=="standard")options.recoveryMode=geo::NativeCadRecoveryMode::Standard;
                else if(mode=="advanced")options.recoveryMode=geo::NativeCadRecoveryMode::Advanced;
                else throw std::invalid_argument("recovery mode must be standard or advanced");
            }
            else if(a=="--recovery-seed"&&i+1<argc)options.recoverySeed=number<std::uint64_t>(argv[++i]);
            else if(a=="--recovery-retries"&&i+1<argc)options.maximumRecoveryRandomRetries=number<unsigned>(argv[++i]);
            else if(a=="--recovery-steiner-candidates"&&i+1<argc)options.maximumEdgeRecoverySteinerCandidates=number<unsigned>(argv[++i]);
            else if(a=="--recovery-patch-states"&&i+1<argc)options.maximumRecoveryPatchStates=number<unsigned>(argv[++i]);
            else if(a=="--brep")io.format=geo::CadFileFormat::Brep;
            else if(a=="--threads"&&i+1<argc)options.geometry.threads=number<int>(argv[++i]);
            else if(a=="--size-factor"&&i+1<argc){sizeFactor=number<double>(argv[++i]);if(!(sizeFactor>0))throw std::invalid_argument("size factor must be positive");}
            else if(a=="--max-file-bytes"&&i+1<argc)io.maximumFileBytes=number<std::uint64_t>(argv[++i]);
            else if(a=="--max-nodes"&&i+1<argc)options.geometry.maximumNodes=number<std::uint64_t>(argv[++i]);
            else if(a=="--max-cells"&&i+1<argc)options.geometry.maximumCells=number<std::uint64_t>(argv[++i]);
            else if(a=="--max-bytes"&&i+1<argc)options.geometry.maximumOutputBytes=number<std::uint64_t>(argv[++i]);
            else if(a=="--max-working-bytes"&&i+1<argc)options.maximumPrimitiveWorkingBytes=number<std::uint64_t>(argv[++i]);
            else if(a=="--max-audit-bytes"&&i+1<argc)options.geometry.maximumAuditTemporaryBytes=number<std::uint64_t>(argv[++i]);
            else if(a=="--max-work"&&i+1<argc) { options.maximumValidationWork=number<std::uint64_t>(argv[++i]);options.geometry.maximumMeshingWork=options.maximumValidationWork; }
            else if(a=="--max-chord-error"&&i+1<argc)options.maximumChordError=number<double>(argv[++i]);
            else if(a=="--max-volume-error"&&i+1<argc)options.maximumRelativeVolumeError=number<double>(argv[++i]);
            else if(a=="--min-quality"&&i+1<argc)options.geometry.minimumElementQuality=number<double>(argv[++i]);
            else if(a=="--warning-quality"&&i+1<argc)options.geometry.warningElementQuality=number<double>(argv[++i]);
            else if(a=="--min-size"&&i+1<argc)options.minimumSize=number<double>(argv[++i]);
            else if(a=="--face-size"&&i+2<argc) { const int tag=number<int>(argv[++i]); const double size=number<double>(argv[++i]); options.faceSizes.push_back({tag,size,growth}); }
            else if(a=="--native-sizing")options.useSpatialSizing=true;
            else if(a=="--size-box"&&i+7<argc) {
                const Vec3 lower{number<double>(argv[++i]),number<double>(argv[++i]),number<double>(argv[++i])};
                const Vec3 upper{number<double>(argv[++i]),number<double>(argv[++i]),number<double>(argv[++i])};
                options.sizeBoxes.push_back({lower,upper,number<double>(argv[++i]),growth});
            }
            else if(a=="--growth"&&i+1<argc) {
                growth=number<double>(argv[++i]);
                if(!std::isfinite(growth)||growth<1)throw std::invalid_argument("growth must be at least one");
                for(auto& c:options.faceSizes)c.growth=growth;
                for(auto& c:options.sizeBoxes)c.growth=growth;
            }
            else if(a=="--elements"&&i+1<argc) {const std::string_view e(argv[++i]);if(e=="tet4")options.elements=algo::PrimitiveElements::Tet4;
                else if(e=="native"||e=="hex")options.elements=algo::PrimitiveElements::HexDominant;else throw std::invalid_argument("native CAD supports tet4 or native linear cells");}
            else throw std::invalid_argument("unknown or incomplete native CAD option");
        }
        if(!std::isfinite(sizeFactor)||sizeFactor<0)throw std::invalid_argument("invalid size factor");
        std::error_code aliasError;
        const bool existingAlias=std::filesystem::equivalent(argv[1],argv[3],aliasError);
        if(existingAlias||std::filesystem::weakly_canonical(argv[1])==std::filesystem::weakly_canonical(argv[3]))
            throw std::invalid_argument("CAD input and mesh output must be different files");
        if(options.stageObserver)options.stageObserver({"import","",-1,true});
        const auto importing=std::chrono::steady_clock::now();cad=geo::importOcctCad(argv[1],io);
        importSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-importing).count();
        if(options.stageObserver)options.stageObserver({"import",geo::cadImportReasonName(cad.reason),-1,false});
        if(!cad.imported()){receipt("refused",geo::cadImportReasonName(cad.reason));std::cerr<<"CAD import: "<<geo::cadImportReasonName(cad.reason)<<'\n';return 3;}
        if(sizeFactor>0){Vec3 lo{1e300,1e300,1e300},hi{-1e300,-1e300,-1e300};for(auto s:cad.model->solids()){
            Vec3 a,b;s->bbox(a,b);lo={std::min(lo.x,a.x),std::min(lo.y,a.y),std::min(lo.z,a.z)};hi={std::max(hi.x,b.x),std::max(hi.y,b.y),std::max(hi.z,b.z)};}
            options.geometry.targetSize=(hi-lo).norm()*sizeFactor;}
        result=geo::meshNativeCad(*cad.model,options,mesh);
        if(!result.published){receipt("refused",geo::nativeCadMeshReasonName(result.reason));std::cerr<<"Native CAD mesh: "<<geo::nativeCadMeshReasonName(result.reason)<<'\n';return 4;}
        if(options.stageObserver)options.stageObserver({"archive","",-1,true});
        const auto saved=sdk::writeMeshArchive(mesh,argv[3]);archivePublished=saved.completed;
        if(options.stageObserver)options.stageObserver({"archive",saved.completed?"none":"archive_io",-1,false});
        if(!saved.completed){receipt("failed","archive_io");return 5;}
        receipt("succeeded","none");return 0;
    }catch(const std::bad_alloc&){receipt("failed","allocation_failed");return 6;}
     catch(const std::exception& e){receipt("failed","invalid_request");std::cerr<<e.what()<<'\n';return 2;}
     catch(...){receipt("failed","internal_error");return 6;}
}
