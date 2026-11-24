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
void writeRouteResolution(const paacs::geo::NativeCadRouteResolutionReport& r){
    std::cout<<"{\"reason\":\""<<r.reason<<"\",\"candidate\":\""<<r.candidate<<"\""
        <<",\"preflight_attempted\":"<<(r.preflightAttempted?"true":"false")
        <<",\"whole_body_verified\":"<<(r.wholeBodyVerified?"true":"false")
        <<",\"preflight_work\":"<<r.preflightWork<<",\"preflight_seconds\":"<<JsonFinite{r.preflightSeconds}
        <<",\"fallback_chord_limit\":";
    if(r.fallbackChordLimit)std::cout<<JsonFinite{*r.fallbackChordLimit};else std::cout<<"null";
    std::cout<<",\"attempts\":[";
    for(std::size_t i=0;i<std::min<std::size_t>(r.count,r.attempts.size());++i){
        if(i)std::cout<<',';const auto& a=r.attempts[i];
        std::cout<<"{\"method\":\""<<a.method<<"\",\"reason\":\""<<paacs::geo::nativeCadMeshReasonName(a.reason)
            <<"\",\"published\":"<<(a.published?"true":"false")<<",\"work\":"<<a.work<<",\"seconds\":"<<JsonFinite{a.seconds}<<'}';
    }
    std::cout<<"],\"scope\":\"method choice and charged attempts; not an element-size or solver certificate\"}";
}
void writeInitialTopologyQuality(const paacs::geo::NativeCadInitialTopologyQualityReport& r){
    std::cout<<"{\"attempted\":"<<(r.attempted?"true":"false")
        <<",\"complete\":"<<(r.completed?"true":"false")<<",\"changed\":"<<(r.changed?"true":"false")
        <<",\"cancelled\":"<<(r.cancelled?"true":"false")<<",\"budget_limited\":"<<(r.budgetLimited?"true":"false")
        <<",\"summary_available\":"<<(r.summaryAvailable?"true":"false")
        <<",\"candidate_list_truncated\":"<<(r.candidateListTruncated?"true":"false")
        <<",\"diagnostic\":\""<<r.diagnostic<<"\",\"work\":"<<r.work<<",\"work_limit\":"<<r.workLimit
        <<",\"seconds\":"<<JsonFinite{r.seconds}<<",\"peak_modeled_bytes\":"<<r.maximumModeledBytes
        <<",\"search_below_mean_ratio\":"<<JsonFinite{r.searchBelowMeanRatio}
        <<",\"candidates\":"<<r.candidates<<",\"proposals\":"<<r.proposals<<",\"certifications\":"<<r.certifications
        <<",\"staged_edits\":"<<r.stagedEdits<<",\"edits\":"<<r.edits<<",\"flips23\":"<<r.flips23
        <<",\"flips32\":"<<r.flips32<<",\"flips44\":"<<r.flips44
        <<",\"live_before\":"<<r.liveBefore<<",\"live_after\":"<<r.liveAfter
        <<",\"severe_before\":"<<r.severeBefore<<",\"severe_after\":"<<r.severeAfter
        <<",\"protected_source_faces\":"<<r.protectedSourceFaces<<",\"protected_source_edges\":"<<r.protectedSourceEdges;
    const auto measurement=[&](const char* name,double value){
        std::cout<<",\""<<name<<"\":";if(r.summaryAvailable)std::cout<<JsonFinite{value};else std::cout<<"null";
    };
    measurement("minimum_mean_ratio_before",r.minimumMeanRatioBefore);
    measurement("minimum_mean_ratio_after",r.minimumMeanRatioAfter);
    measurement("minimum_scaled_jacobian_before",r.minimumScaledJacobianBefore);
    measurement("minimum_scaled_jacobian_after",r.minimumScaledJacobianAfter);
    std::cout<<",\"delaunay_certified\":false,\"scope\":\"nonguard numerical cells including exterior; final Native audits remain required\"}";
}
void writeMaterialClassification(const paacs::geo::NativeCadMaterialClassificationReport& r){
    std::cout<<"{\"diagnostic\":\""<<r.diagnostic<<"\",\"complete\":"<<(r.completed?"true":"false")
        <<",\"source_model_verified_unchanged\":"<<(r.sourceModelUnmodified?"true":"false")
        <<",\"attempts\":"<<r.attempts<<",\"requested_points\":"<<r.requestedPoints
        <<",\"attempted_points\":"<<r.attemptedPoints<<",\"native_queries\":"<<r.nativeQueries
        <<",\"strict_reuses\":"<<r.strictReuses<<",\"ray_admissions\":"<<r.rayAdmissions
        <<",\"work\":"<<r.work<<",\"seconds\":"<<JsonFinite{r.seconds}
        <<",\"profiled_batches\":"<<r.profiledBatches<<",\"ray_witness_attempts\":"<<r.rayWitnessAttempts
        <<",\"classifier_setup_seconds\":"<<JsonFinite{r.classifierSetupSeconds}
        <<",\"native_query_seconds\":"<<JsonFinite{r.nativeQuerySeconds}
        <<",\"ray_witness_seconds\":"<<JsonFinite{r.rayWitnessSeconds}
        <<",\"scope\":\"cumulative material classification attempts; last batch status; not final publication\"}";
}
void writeSourceShell(const paacs::geo::NativeCadSourceShellReport& report){
    const auto count=std::min<std::size_t>(report.count,report.attempts.size());
    const auto positive=[](double value){if(value>0)std::cout<<JsonFinite{value};else std::cout<<"null";};
    std::cout<<"{\"available\":"<<(count?"true":"false")<<",\"accepted\":";
    if(count)std::cout<<(report.accepted?"true":"false");else std::cout<<"null";
    std::cout<<",\"native_fallback_requested\":"<<(report.nativeFallbackRequested?"true":"false")
        <<",\"native_fallback_attempted\":"<<(report.nativeFallbackAttempted?"true":"false")
        <<",\"native_fallback_diagnostic\":\""<<report.nativeFallbackDiagnostic<<"\"";
    std::cout<<",\"diagnostic\":\""<<report.diagnostic<<"\",\"requested_global_size\":";
    positive(report.requestedGlobalSize);std::cout<<",\"requested_chord_error\":";positive(report.requestedChordError);
    std::cout<<",\"maximum_relative_volume_error\":";positive(report.maximumRelativeVolumeError);
    std::cout<<",\"attempts\":[";
    for(std::size_t i=0;i<count;++i){
        if(i)std::cout<<',';const auto& a=report.attempts[i];
        std::cout<<"{\"attempt\":"<<i+1<<",\"chord_error_limit\":"<<JsonFinite{a.chordErrorLimit}
            <<",\"work\":"<<a.work<<",\"seconds\":"<<JsonFinite{a.seconds}
            <<",\"nodes\":"<<a.nodes<<",\"triangles\":"<<a.triangles
            <<",\"surface_published\":"<<(a.published?"true":"false")
            <<",\"surface_diagnostic\":\""<<a.surfaceDiagnostic<<"\""
            <<",\"source_model_verified_unchanged\":"<<(a.sourceModelUnmodified?"true":"false")
            <<",\"volume_criterion_met\":";
        if(a.published)std::cout<<(a.volumeCriterionMet?"true":"false");else std::cout<<"null";
        std::cout<<",\"signed_volume_mm3\":";
        if(a.published)std::cout<<JsonFinite{a.signedVolume};else std::cout<<"null";
        std::cout<<",\"cad_volume_mm3\":";
        if(a.published)std::cout<<JsonFinite{a.cadVolume};else std::cout<<"null";
        std::cout<<",\"relative_volume_error\":";
        if(a.published)std::cout<<JsonFinite{a.relativeVolumeError};else std::cout<<"null";
        std::cout<<",\"construction_budget_limited\":"<<(a.constructionBudgetLimited?"true":"false")
            <<",\"construction_stalled\":"<<(a.constructionStalled?"true":"false")
            <<",\"provider\":\""<<a.provider<<"\",\"global_sizing_applied\":"<<(a.globalSizingApplied?"true":"false")
            <<",\"curvature_points_per_turn\":"<<a.curvaturePointsPerTurn<<",\"angular_deflection\":";
        positive(a.angularDeflection);
        std::cout<<",\"witness\":{\"stage\":\""<<a.witness.stage<<"\",\"detail\":\""<<a.witness.detail<<"\",\"edge_tag\":";
        if(a.witness.edgeTag>=0)std::cout<<a.witness.edgeTag;else std::cout<<"null";
        std::cout<<",\"face_tag\":";if(a.witness.faceTag>=0)std::cout<<a.witness.faceTag;else std::cout<<"null";
        std::cout<<",\"use_id\":";if(a.witness.useId)std::cout<<*a.witness.useId;else std::cout<<"null";
        std::cout<<",\"geometric_refusal\":"<<(a.witness.geometricRefusal?"true":"false")<<"}}";
    }
    std::cout<<"],\"scope\":\"latest_source_strategy; geometry acceptance is not a sizing certificate\"}";
}
void writeInitialFrontFallback(const paacs::geo::NativeCadInitialFrontFallbackReport& report){
    std::cout<<"{\"attempted\":"<<(report.attempted?"true":"false")
        <<",\"reason\":\""<<paacs::geo::nativeCadMeshReasonName(report.reason)<<"\""
        <<",\"work\":"<<report.work<<",\"seconds\":"<<report.seconds
        <<",\"surface_diagnostic\":\""<<report.surfaceDiagnostic<<"\""
        <<",\"volume_diagnostic\":\""<<report.volumeDiagnostic<<"\",\"source_shell\":";
    writeSourceShell(report.sourceShell);std::cout<<'}';
}
void writePlanarBody(const paacs::algo::GeometryBodyReport& b){
    std::cout<<"{\"solid_tag\":"<<b.solidTag
        <<",\"strategy\":\""<<paacs::algo::geometryStrategyName(b.strategy)<<"\""
        <<",\"reason\":\""<<paacs::algo::geometryReasonName(b.reason)<<"\""
        <<",\"resolved_options\":{\"elements\":\""<<(b.elements==paacs::algo::GeometryElements::Tet4?"tet4":"native")
        <<"\",\"target_size\":"<<JsonFinite{b.targetSize}
        <<",\"source_face_tag\":"<<b.sourceFaceTag<<",\"target_face_tag\":"<<b.targetFaceTag
        <<",\"layers\":"<<b.layers<<",\"end_to_start_layer_width_ratio\":"<<JsonFinite{b.endToStartLayerWidthRatio}
        <<",\"used_spatial_field\":"<<(b.usedSpatialField?"true":"false")
        <<",\"used_size_profile\":"<<(b.usedSizeProfile?"true":"false")
        <<",\"spatial_minimum_bound\":"<<JsonFinite{b.spatialMinimumBound}<<'}'
        <<",\"estimates\":{\"nodes\":"<<b.estimatedNodes<<",\"total_cells\":"<<b.estimatedCells
        <<",\"output_bytes\":"<<b.estimatedOutputBytes<<",\"exact\":"<<(b.exactCountEstimate?"true":"false")<<'}'
        <<",\"engine\":{\"status_code\":"<<b.engineStatus<<",\"reason_code\":"<<b.engineReason
        <<",\"phase_code\":"<<b.enginePhase<<'}'
        <<",\"output\":{\"nodes\":"<<b.outputNodes<<",\"total_cells\":"<<b.outputCells<<",\"bytes\":"<<b.outputBytes<<'}'
        <<",\"generation_seconds\":"<<JsonFinite{b.generationSeconds};
    const auto& c=b.conversion;
    std::cout<<",\"conversion\":{\"published\":"<<(c.published?"true":"false")
        <<",\"reason\":\""<<paacs::algo::convexToTetReasonName(c.reason)<<"\""
        <<",\"tetrahedra\":"<<c.tetrahedra<<",\"boundary_triangles\":"<<c.boundaryTriangles
        <<",\"output_bytes\":"<<c.outputBytes<<",\"working_bytes\":"<<c.workingBytes
        <<",\"minimum_quality\":"<<JsonFinite{c.minimumQuality}<<",\"low_quality_cells\":"<<c.lowQualityCells<<'}'
        <<",\"boundary\":{\"available\":"<<(b.hasBoundaryDiagnostic?"true":"false")
        <<",\"status_code\":"<<b.boundaryStatus<<",\"reason_code\":"<<b.boundaryReason
        <<",\"phase_code\":"<<b.boundaryPhase<<'}'
        <<",\"surface\":{\"available\":"<<(b.hasSurfaceDiagnostic?"true":"false")
        <<",\"status_code\":"<<b.surfaceStatus<<",\"reason_code\":"<<b.surfaceReason
        <<",\"failure_face_id\":"<<b.surfaceFailureFacePersistentId
        <<",\"processed_faces\":"<<b.surfaceProcessedFaces<<",\"unattempted_faces\":"<<b.surfaceUnattemptedFaces<<"}}";
}
void writePlanar(const paacs::geo::NativeCadBodyMeshReport& body){
    std::cout<<"{\"attempted\":"<<(body.planarAttempted?"true":"false")
        <<",\"scope\":\"retained geometry-core receipt; default substage codes do not prove stage entry; separate from Native CAD certification\"";
    if(!body.planarAttempted){std::cout<<",\"published\":false,\"reason\":null}";return;}
    const auto& p=body.planar;
    std::cout<<",\"published\":"<<(p.published?"true":"false")
        <<",\"reason\":\""<<paacs::algo::geometryReasonName(p.reason)<<"\""
        <<",\"analysis\":{\"complete\":"<<(p.analysis.complete?"true":"false")
        <<",\"reason\":\""<<paacs::algo::geometryReasonName(p.analysis.reason)<<"\",\"work\":"<<p.analysis.work<<'}'
        <<",\"seconds\":{\"analysis\":"<<JsonFinite{p.analysisSeconds}
        <<",\"generation\":"<<JsonFinite{p.generationSeconds}<<",\"audit\":"<<JsonFinite{p.auditSeconds}
        <<",\"quality\":"<<JsonFinite{p.qualitySeconds}<<'}';
    constexpr std::size_t maximumBodyRows=16;
    const auto count=std::min(p.analysis.bodies.size(),maximumBodyRows);
    std::cout<<",\"bodies_total\":"<<p.analysis.bodies.size()<<",\"bodies_truncated\":"<<(p.analysis.bodies.size()>count?"true":"false")
        <<",\"bodies\":[";
    for(std::size_t i=0;i<count;++i){if(i)std::cout<<',';writePlanarBody(p.analysis.bodies[i]);}
    std::cout<<"]}";
}
}
int runNativeCadCli(int argc,char** argv) {
    using namespace paacs;
    geo::NativeCadImportOptions io;geo::NativeCadImportOutcome cad;
    auto options=geo::normalNativeCadMeshOptions();geo::NativeCadMeshOutcome result;mesh::Mesh mesh;
    bool report=false,archivePublished=false;double sizeFactor=0,importSeconds=0,growth=1.3;
    std::optional<std::uint64_t> requestedWork;
    std::optional<double> requestedChordError,requestedRelativeVolumeError;
    std::optional<unsigned> requestedCurvaturePoints;
    const auto boundedValidationWork=geo::NativeCadMeshOptions{}.maximumValidationWork;
    const auto boundedMeshingWork=geo::NativeCadMeshOptions{}.geometry.maximumMeshingWork;
    const auto resolveWork=[&] {
        const bool normal=options.construction.resourcePolicy==geo::CadResourcePolicy::Normal;
        options.maximumValidationWork=requestedWork.value_or(normal?geo::cadNormalWorkAllowance:boundedValidationWork);
        options.geometry.maximumMeshingWork=requestedWork.value_or(normal?geo::cadNormalWorkAllowance:boundedMeshingWork);
    };
    resolveWork();
    const auto started=std::chrono::steady_clock::now();
    const auto receipt=[&](const char* status,const char* reason) {
        resolveWork();
        if(!report)return;
        const auto& construction=options.construction;
        const auto cap=[](std::uint64_t value){return value?std::to_string(value):std::string("null");};
        bool frontSelected=options.advancingInteriorFrontEnabled(std::size_t(cad.solids));
        const char* frontSelectionBasis="request_eligibility";
        const char* resolvedInitial="not_reported";
        if(result.routeResolution.count){
            frontSelectionBasis="route_resolution_attempts";
            frontSelected=false;
            for(std::size_t i=0;i<std::min<std::size_t>(result.routeResolution.count,result.routeResolution.attempts.size());++i)
                frontSelected|=std::string_view(result.routeResolution.attempts[i].method)=="native_front";
            resolvedInitial=result.routeResolution.attempts[
                std::min<std::size_t>(result.routeResolution.count,result.routeResolution.attempts.size())-1].method;
        }
        if(std::string_view(result.routeResolution.reason)=="explicit_front_enabled")resolvedInitial="native_front";
        if(std::string_view(result.routeResolution.reason)=="explicit_legacy_dispatch")resolvedInitial="legacy_dispatch";
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
            <<",\"resources\":{\"policy\":\""<<(construction.resourcePolicy==geo::CadResourcePolicy::Normal?"normal":"bounded")
            <<"\",\"requested_max_work\":"<<(requestedWork?std::to_string(*requestedWork):"null")
            <<",\"effective_max_work\":"<<std::min(options.maximumValidationWork,options.geometry.maximumMeshingWork)
            <<",\"logical_work_cap_disabled\":"<<(construction.resourcePolicy==geo::CadResourcePolicy::Normal&&!requestedWork?"true":"false")
            <<",\"max_working_bytes\":"<<options.maximumPrimitiveWorkingBytes<<"}"
            <<",\"fidelity\":{\"requested_max_chord_error\":";
        if(requestedChordError)std::cout<<JsonFinite{*requestedChordError};else std::cout<<"null";
        std::cout<<",\"requested_max_relative_volume_error\":";
        if(requestedRelativeVolumeError)std::cout<<JsonFinite{*requestedRelativeVolumeError};else std::cout<<"null";
        std::cout<<",\"effective_max_relative_volume_error\":"<<JsonFinite{options.maximumRelativeVolumeError}
            <<",\"requested_curvature_points_per_turn\":";
        if(requestedCurvaturePoints)std::cout<<*requestedCurvaturePoints;else std::cout<<"null";
        std::cout<<",\"effective_curvature_points_per_turn\":"<<options.curvaturePointsPerTurn
            <<",\"chord_default\":\""<<(requestedChordError?"explicit":"route_specific")<<"\""
            <<",\"scope\":\"native geometry acceptance limits; measured errors and achieved sizes are separate\"}"
            <<",\"route_resolution\":";
        writeRouteResolution(result.routeResolution);
        std::cout<<",\"faces\":"<<cad.faces<<",\"curved_faces\":"<<cad.curvedFaces<<",\"edges\":"<<cad.edges
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
            <<",\"native_spatial_sizing\":"<<((options.useSpatialSizing||!options.sizeBoxes.empty()||frontSelected)?"true":"false")
            <<",\"face_size_controls\":"<<options.faceSizes.size()<<",\"size_box_controls\":"<<options.sizeBoxes.size()
            <<",\"live_material_recovery\":"<<(options.useLiveMaterialRecovery?"true":"false")
            <<",\"trimmed_ray_witnesses\":"<<(options.useTrimmedRayWitnesses?"true":"false")
            <<",\"local_delaunay_refinement\":"<<(options.useLocalDelaunayRefinement?"true":"false")
            <<",\"improve_surface_spacing\":"<<(options.improveNativeSurfaceSpacing?"true":"false")
            <<",\"improve_tet_quality\":"<<(options.improveNativeTetQuality?"true":"false")
            <<",\"advancing_interior_front\":"<<(interiorFrontAttempted?"true":"false")
            <<",\"advancing_interior_front_selected\":"<<(frontSelected?"true":"false")
            <<",\"advancing_interior_front_selection_basis\":\""<<frontSelectionBasis<<"\""
            <<",\"resolved_initial_construction\":\""<<resolvedInitial<<"\""
            <<",\"initial_construction_scope\":\"chosen constructor from actual resolution receipts; stage entry is reported separately\""
            <<",\"advancing_interior_front_policy\":\""<<(!options.useAdvancingInteriorFront.has_value()?"auto":*options.useAdvancingInteriorFront?"enabled":"disabled")<<"\""
            <<",\"requested_initial_topology_quality\":";
        if(options.initialTopologyQuality)std::cout<<(*options.initialTopologyQuality?"true":"false");else std::cout<<"null";
        std::cout
            <<",\"requested_native_source_fallback\":"<<(options.allowNativeSourceFallback?"true":"false")
            <<",\"requested_boundary_first\":"<<(options.recoverBoundaryBeforeInteriorSites?"true":"false")
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
            std::cout<<",\"planar\":";writePlanar(b);
            std::cout<<",\"initial_front_fallback\":";
            writeInitialFrontFallback(b.initialFrontFallback);
            std::cout<<",\"source_shell\":";writeSourceShell(b.sourceShell);
            std::cout<<",\"initial_topology_quality\":";writeInitialTopologyQuality(b.initialTopologyQuality);
            std::cout<<",\"material_classification\":";writeMaterialClassification(b.materialClassification);
            std::cout<<",\"route_resolution\":";writeRouteResolution(b.routeResolution);
            std::cout<<",\"fidelity\":{\"effective_max_chord_error\":";
            if(b.effectiveMaximumChordError)std::cout<<JsonFinite{*b.effectiveMaximumChordError};else std::cout<<"null";
            std::cout<<",\"chord_policy\":\""<<b.chordErrorPolicy<<"\""
                <<",\"scope\":\"selected construction limit; not a complete CAD distance certificate\"}";
            std::cout<<",\"volume_sizing\":{\"applied\":";
            if(b.volumeSizingApplied)std::cout<<(*b.volumeSizingApplied?"true":"false");else std::cout<<"null";
            std::cout<<",\"diagnostic\":\""<<b.volumeSizingDiagnostic<<"\",\"requested_global_size\":"<<JsonFinite{options.geometry.targetSize}
                <<",\"sizing_completion_certified\":false,\"basis\":\"reported_native_construction\"}";
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
                <<",\"auxiliary_imported_cells\":"<<b.initialFront.auxiliaryImportedCells
                <<",\"auxiliary_import_work\":"<<b.initialFront.auxiliaryImportWork
                <<",\"auxiliary_import_modeled_bytes\":"<<b.initialFront.auxiliaryImportModeledBytes
                <<",\"auxiliary_initial_capacity\":"<<b.initialFront.auxiliaryInitialCapacity
                <<",\"auxiliary_hard_capacity\":"<<b.initialFront.auxiliaryHardCapacity
                <<",\"auxiliary_growth_attempts\":"<<b.initialFront.auxiliaryGrowthAttempts
                <<",\"auxiliary_growths\":"<<b.initialFront.auxiliaryGrowths
                <<",\"auxiliary_growth_work\":"<<b.initialFront.auxiliaryGrowthWork
                <<",\"maximum_combined_modeled_bytes\":"<<b.initialFront.maximumCombinedModeledBytes
                <<",\"maximum_heap_capacity\":"<<b.initialFront.maximumHeapCapacity
                <<",\"resource_retries\":"<<b.initialFront.resourceRetries
                <<",\"geometry_refusal_cache_hits\":"<<b.initialFront.geometryRefusalCacheHits
                <<",\"cached_geometry_refusals\":"<<b.initialFront.cachedGeometryRefusals
                <<",\"auxiliary_resource_diagnostic\":\""<<b.initialFront.auxiliaryResourceDiagnostic<<"\""
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
                <<",\"piercing_edge_patches\":"<<constraints.piercingPatches
                <<",\"initial_topology_used\":"<<(constraints.initialTopologyUsed?"true":"false")
                <<",\"initial_topology_diagnostic\":\""<<constraints.initialTopologyDiagnostic<<"\""
                <<",\"initial_topology_cells\":"<<constraints.initialTopologyCells
                <<",\"initial_topology_import_work\":"<<constraints.initialTopologyImportWork
                <<",\"initial_topology_peak_modeled_bytes\":"<<constraints.initialTopologyPeakModeledBytes<<'}';
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
        if(argc<4)throw std::invalid_argument("usage: paacs-cad-mesh input.step size output.paacsmesh --native-cad [--face-size tag size] [--native-sizing] [--size-box x0 y0 z0 x1 y1 z1 size] [--min-size h] [--growth g] [--report-json] [--trace-native-stages] [--independent-parts] [--local-partition-recovery] [--adaptive-tet-recovery] [--live-material-recovery] [--curvature-points-per-turn 0..64] [--trimmed-ray-witnesses] [--local-delaunay-refinement] [--advancing-interior-front|--no-advancing-interior-front] [--initial-topology-quality|--no-initial-topology-quality] [--native-source-fallback|--no-native-source-fallback] [--boundary-first|--no-boundary-first] [--improve-tet-quality] [--improve-surface-spacing] [--regular-surface-front] [--no-improve-regular-quality] [--recovery-mode standard|advanced] [--recovery-seed N] [--recovery-retries 0..4] [--recovery-steiner-candidates 0..64] [--resource-policy normal|bounded] [--max-work N] [--max-working-bytes N] [--construction-effort legacy|auto|thorough] [--max-surface-points N] [--max-surface-chord-points N] [--max-interior-candidates N] [--max-interior-sites N]");
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
            else if(a=="--curvature-points-per-turn"&&i+1<argc){
                requestedCurvaturePoints=number<unsigned>(argv[++i]);options.curvaturePointsPerTurn=*requestedCurvaturePoints;
            }
            else if(a=="--improve-tet-quality")options.improveNativeTetQuality=true;
            else if(a=="--advancing-interior-front")options.useAdvancingInteriorFront=true;
            else if(a=="--no-advancing-interior-front")options.useAdvancingInteriorFront=false;
            else if(a=="--initial-topology-quality")options.initialTopologyQuality=true;
            else if(a=="--no-initial-topology-quality")options.initialTopologyQuality=false;
            else if(a=="--native-source-fallback")options.allowNativeSourceFallback=true;
            else if(a=="--no-native-source-fallback")options.allowNativeSourceFallback=false;
            else if(a=="--boundary-first")options.recoverBoundaryBeforeInteriorSites=true;
            else if(a=="--no-boundary-first")options.recoverBoundaryBeforeInteriorSites=false;
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
            else if(a=="--resource-policy"&&i+1<argc){
                const std::string_view policy(argv[++i]);
                if(policy=="normal")options.construction.resourcePolicy=geo::CadResourcePolicy::Normal;
                else if(policy=="bounded")options.construction.resourcePolicy=geo::CadResourcePolicy::Bounded;
                else throw std::invalid_argument("resource policy must be normal or bounded");
                resolveWork();
            }
            else if(a=="--max-work"&&i+1<argc) {
                const auto work=number<std::uint64_t>(argv[++i]);
                if(!work)throw std::invalid_argument("maximum work must be positive");
                requestedWork=work;
                resolveWork();
            }
            else if(a=="--max-chord-error"&&i+1<argc){
                requestedChordError=number<double>(argv[++i]);options.maximumChordError=*requestedChordError;
            }
            else if((a=="--max-volume-error"||a=="--max-relative-volume-error")&&i+1<argc){
                requestedRelativeVolumeError=number<double>(argv[++i]);options.maximumRelativeVolumeError=*requestedRelativeVolumeError;
            }
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
        resolveWork();
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
