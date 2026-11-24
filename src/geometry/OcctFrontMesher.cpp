// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctStageTrace.hpp"
#include "OcctFrontMesher.hpp"
#include "OcctBoundaryInternal.hpp"
#include "OcctModelInternal.hpp"
#include "OcctSolidClassifier.hpp"
#include "OcctBoundaryInternal.hpp"
#include "OcctStarMesher.hpp"
#include "OcctWireRefinement.hpp"
#include "OcctFacetRefinement.hpp"
#include "OcctNativeWireAuthority.hpp"
#include "OcctBoundaryDiagonalClosure.hpp"
#include "OcctPlcMaterialFilter.hpp"
#include "OcctSeedSampling.hpp"
#include "OcctTetQuality.hpp"
#include "OcctLiveCavity.hpp"
#include "OcctLocalDelaunay.hpp"
#include "OcctDelaunayMeshQuality.hpp"
#include "OcctTrimmedRayWitness.hpp"
#include "OcctInwardSites.hpp"
#include "OcctShellProximity.hpp"
#include "OcctSpatialSizing.hpp"
#include "OcctInteriorSites.hpp"
#include "OcctInteriorFront.hpp"
#include "OcctInitialPointFront.hpp"
#include "OcctInitialTopologyQuality.hpp"
#include "paacs/algo/BoundaryConformingMesher.hpp"
#include "paacs/bnd/TetPlcCavityRecovery.hpp"
#include "../bnd/TetPlcRetainedRecovery.hpp"
#include "../bnd/TetPlcGiftWrapping.hpp"
#include "../bnd/TetConstraintRecovery.hpp"
#include <Standard_OutOfMemory.hxx>
#include <Standard_Failure.hxx>
#include <Precision.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <bit>
#include <functional>
#include <chrono>
#include <stop_token>

namespace paacs::geo::detail {
namespace {
using Clock=std::chrono::steady_clock;
struct Stop{NativeCadMeshReason reason;const char* diagnostic;};
void need(bool ok,NativeCadMeshReason reason,const char* why){if(!ok)throw Stop{reason,why};}
}
OcctFrontOutcome buildOcctFrontVolume(const OcctModel& model,int tag,
    const NativeCadMeshOptions& options,std::uint64_t maximumWork,mesh::Mesh& destination,
    std::uint64_t earlyConstraintWork,OcctFrontSurface surfaceMode,bool liveTrial) noexcept {
    // Public dispatch resolves AUTO once; direct private callers opt in explicitly.
    const bool advancingInteriorFront=options.useAdvancingInteriorFront.value_or(false);
    const bool spaceSurface=surfaceMode!=OcctFrontSurface::Standard;
    const bool spatialSizing=surfaceMode==OcctFrontSurface::SpatialSizing;
    const bool improveQuality=surfaceMode==OcctFrontSurface::RegularQuality||
        (spatialSizing&&options.regularMeshQualityEnabled());
    const bool curvatureSurface=surfaceMode==OcctFrontSurface::Curvature;
    const bool regularSurface=surfaceMode==OcctFrontSurface::Regular||improveQuality||curvatureSurface||spatialSizing;
    const bool completeDiscreteSeed=liveTrial||curvatureSurface;
    const bool progressiveResources=options.construction.resourcePolicy==CadResourcePolicy::Normal;
    const bool customConstruction=progressiveResources||options.construction.effort!=CadConstructionEffort::Legacy||
        options.construction.surfacePointsPerFace||options.construction.surfaceChordPointsPerFace||
        options.construction.interiorCandidates||options.construction.interiorSites;
    const char* stage="surface";OcctStageTrace trace(stage,tag);OcctFrontOutcome out;const auto began=Clock::now();const auto& o=options.geometry;std::stop_source stop;
    const auto account=[&](std::uint64_t n){
        need(out.workUsed<=maximumWork&&n<=maximumWork-out.workUsed,NativeCadMeshReason::BudgetExceeded,"native-plc-work-budget");out.workUsed+=n;};
    const auto tick=[&](std::uint64_t n=1){if(o.cancellation&&o.cancellation->load()){stop.request_stop();throw Stop{NativeCadMeshReason::Cancelled,"cancelled"};}
        account(n);};
    // Independent stages return the work they have already performed, even
    // when they observed cancellation. Retain that receipt before stopping;
    // incremental callbacks still use tick before starting another operation.
    const auto accountCompleted=[&](std::uint64_t n){account(n);tick(0);};
    const auto remaining=[&](){tick(0);need(out.workUsed<maximumWork,NativeCadMeshReason::BudgetExceeded,"native-plc-work-budget");return maximumWork-out.workUsed;};
    try {
        tick();need(model.solid(tag)&&maximumWork>0&&options.maximumRecoveryRandomRetries<=4&&options.maximumRecoveryPatchStates<=65536,NativeCadMeshReason::InvalidRequest,"native-plc-invalid-request");Vec3 lo,hi;model.solid(tag)->bbox(lo,hi);
        std::optional<OcctSpatialSizing> sizeField;
        if(spatialSizing){
            sizing::RegionSizeSpecification specification;
            specification.baseSize=o.targetSize;
            specification.minimumSize=options.minimumSize;
            need(options.sizeBoxes.size()<=specification.maximumControls&&
                 options.faceSizes.size()<=specification.maximumControls,
                 NativeCadMeshReason::BudgetExceeded,"too-many-native-sizing-controls");
            for(const auto& box:options.sizeBoxes) {
                tick();
                need(box.lower.x<box.upper.x&&box.lower.y<box.upper.y&&box.lower.z<box.upper.z,
                     NativeCadMeshReason::InvalidRequest,"invalid-native-sizing-box");
                specification.controls.push_back({sizing::SizeRegionKind::Box,box.lower,box.upper,
                    1.,box.size,box.growth});
            }
            auto compiled=sizing::compileRegionSizeField(specification);
            need(compiled.reason!=sizing::RegionSizeReason::AllocationFailed,
                 NativeCadMeshReason::AllocationFailed,"allocation-failed");
            need(compiled.compiled,NativeCadMeshReason::InvalidRequest,"invalid-native-sizing-field");
            tick(options.faceSizes.size());
            sizeField.emplace(model,options,std::move(compiled.field));
        }
        out.chord=options.maximumChordError>0?options.maximumChordError:std::min(o.targetSize*.05,(hi-lo).norm()*.0025);
        const bool largerRecoveryBudget=options.maximumPrimitiveWorkingBytes>=1024ULL*1024*1024&&options.maximumRecoveryPatchStates>=8192;
        const unsigned regionVertices=256,expandedRegionVertices=largerRecoveryBudget?512:256;
        OcctBoundaryOptions bo;bo.linearDeflection=out.chord*.4;bo.maximumChordError=out.chord;bo.maximumNodes=std::min<std::uint64_t>(largerRecoveryBudget?50000:10000,o.maximumNodes);bo.maximumTriangles=std::min<std::uint64_t>(largerRecoveryBudget?100000:20000,o.maximumCells/2);bo.maximumWork=remaining();bo.maximumOutputBytes=o.maximumOutputBytes;bo.cancellation=o.cancellation;
        bo.construction=options.construction;
        if(customConstruction){bo.maximumNodes=o.maximumNodes;bo.maximumTriangles=o.maximumCells/2;}
        mesh::Mesh boundary;
        OcctBoundaryOutcome surface;
        if(spaceSurface&&!customConstruction){
            bo.maximumNodes=std::min<std::uint64_t>(bo.maximumNodes,curvatureSurface||spatialSizing?6'000:1'200);
            bo.maximumTriangles=std::min<std::uint64_t>(bo.maximumTriangles,curvatureSurface||spatialSizing?10'000:1'800);
        }
        out.sourceShell.requestedGlobalSize=o.targetSize;
        out.sourceShell.requestedChordError=options.maximumChordError;
        out.sourceShell.maximumRelativeVolumeError=options.maximumRelativeVolumeError;
        // Local chord witnesses alone do not establish the global source-volume
        // bound. Normal can tighten geometry resolution before volume placement;
        // the requested size field and historical Bounded dispatch stay intact.
        const unsigned sourceAttempts=progressiveResources?unsigned(out.sourceShell.attempts.size()):1;
        try {
        for(unsigned attempt=0;attempt<sourceAttempts;++attempt){
            bo.maximumWork=remaining();
            bo.maximumChordError=out.chord;bo.linearDeflection=out.chord*.4;
            auto& receipt=out.sourceShell.attempts[out.sourceShell.count++];
            receipt.chordErrorLimit=out.chord;
            const auto sourceStarted=Clock::now();
            if(spaceSurface){
                OcctBoundaryAliasReceipt spacing;
                surface=buildSpacedOcctBoundary(model,tag,bo,o.targetSize,options.maximumPrimitiveWorkingBytes/4,
                    boundary,spacing,regularSurface,improveQuality,
                    curvatureSurface||spatialSizing?options.curvaturePointsPerTurn:0,sizeField?&*sizeField:nullptr);
                out.rebuiltSurfaceCharts+=spacing.rebuiltCharts;
                out.surfaceSpacing=std::move(spacing.surfaceSpacing);
            }else if(options.recoveryMode==NativeCadRecoveryMode::Advanced){
                OcctBoundaryOutcome preparation;
                auto session=OcctBoundarySession::prepare(model,bo,preparation);
                if(session){
                    OcctBoundaryAliasOptions repair;
                    repair.maximumFlips=0;repair.maximumRebuiltCharts=32;
                    repair.maximumWorkingBytes=options.maximumPrimitiveWorkingBytes/4;
                    OcctBoundaryAliasReceipt repairs;
                    surface=session->boundary(tag,boundary,repair,repairs);
                    const auto retryRoom=bo.maximumWork-session->work();
                    if(options.useAdaptiveTetRecovery&&model.solids().size()==1&&
                        retryRoom>0&&retryRoom>=surface.work&&
                        (surface.reason==OcctBoundaryReason::NonconformingEdge||
                         surface.reason==OcctBoundaryReason::OpenBoundary)){
                        const auto firstWork=surface.work;
                        repair.curveMetricSizeCap=o.targetSize;
                        surface=session->boundary(tag,boundary,repair,repairs,retryRoom);
                        surface.work+=firstWork;
                    }
                    surface.work+=preparation.work;
                    out.rebuiltSurfaceCharts+=repairs.rebuiltCharts;
                    out.surfaceSpacing=std::move(repairs.surfaceSpacing);
                }else surface=std::move(preparation);
            }else surface=buildOcctBoundary(model,tag,bo,boundary);
            receipt.seconds=std::chrono::duration<double>(Clock::now()-sourceStarted).count();
            receipt.work=surface.work;receipt.nodes=surface.nodes;receipt.triangles=surface.triangles;
            receipt.published=surface.published;receipt.sourceModelUnmodified=surface.sourceModelUnmodified;
            receipt.surfaceDiagnostic=occtBoundaryReasonName(surface.reason);
            receipt.signedVolume=surface.signedSurfaceVolume;receipt.cadVolume=surface.cadVolume;
            receipt.relativeVolumeError=surface.relativeVolumeError;
            receipt.volumeCriterionMet=surface.published&&surface.relativeVolumeError<=options.maximumRelativeVolumeError;
            for(const auto& face:out.surfaceSpacing){
                receipt.constructionBudgetLimited|=face.budgetLimited;
                receipt.constructionStalled|=!face.queueExhausted&&!face.budgetLimited;
            }
            out.surfaceWork+=surface.work;out.surfaceTriangles=surface.triangles;
            out.maximumNodeError=surface.maximumNodeProjectionError;
            out.maximumFacetError=surface.maximumFacetSampleError;
            out.surfaceDiagnostic=receipt.surfaceDiagnostic;
            out.sourceShell.diagnostic=surface.published?"native-plc-source-volume-error":receipt.surfaceDiagnostic;
            out.surfaceSeconds=std::chrono::duration<double>(Clock::now()-began).count();
            accountCompleted(surface.work);
            need(surface.reason!=OcctBoundaryReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");
            need(surface.reason!=OcctBoundaryReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");
            need(!spaceSurface||surface.reason!=OcctBoundaryReason::SourceChanged,NativeCadMeshReason::InternalError,"spaced-surface-source-changed");
            need(!spaceSurface||surface.published||surface.reason==OcctBoundaryReason::BudgetExceeded,
                 NativeCadMeshReason::GeometryError,out.surfaceDiagnostic);
            need(surface.reason!=OcctBoundaryReason::BudgetExceeded,NativeCadMeshReason::BudgetExceeded,"native-plc-surface-budget");
            need(surface.published,NativeCadMeshReason::GeometryError,"native-plc-surface-refused");
            if(receipt.volumeCriterionMet){out.sourceShell.diagnostic="volume_criterion_met";break;}
            need(!progressiveResources||!receipt.constructionBudgetLimited,NativeCadMeshReason::BudgetExceeded,
                 "native-plc-source-construction-budget");
            need(attempt+1<sourceAttempts,NativeCadMeshReason::GeometryError,"native-plc-source-volume-error");
            const auto tighter=out.chord*.5;
            need(tighter>0&&std::isfinite(tighter)&&(spaceSurface||tighter*.4>=Precision::Confusion()),
                 NativeCadMeshReason::GeometryError,"native-plc-source-resolution-limit");
            // Failed private shells are not retained alongside the next build.
            mesh::Mesh empty;boundary.swap(empty);surface=OcctBoundaryOutcome{};
            std::vector<CadSurfaceSpacingReport>().swap(out.surfaceSpacing);
            out.chord=tighter;
        }
        }catch(const Stop& failure){out.sourceShell.diagnostic=failure.diagnostic;throw;}
        OcctStarOptions certificateOptions;certificateOptions.maximumWork=remaining();certificateOptions.maximumWorkingBytes=options.maximumPrimitiveWorkingBytes;certificateOptions.cancellation=o.cancellation;certificateOptions.maximumBoundaryNodes=bo.maximumNodes;certificateOptions.maximumBoundaryTriangles=bo.maximumTriangles;certificateOptions.maximumNodes=o.maximumNodes;certificateOptions.maximumCells=o.maximumCells;certificateOptions.maximumOutputBytes=o.maximumOutputBytes;
        OcctStageTrace shellCheck("exact_shell_certificate",tag);
        auto certificate=certifyOcctBoundarySurface(boundary,certificateOptions);
        shellCheck.finish(certificate.diagnostic);
        accountCompleted(certificate.workUsed);
        out.sourceShell.diagnostic=certificate.certified?"accepted":certificate.diagnostic;
        need(certificate.reason!=OcctStarReason::WorkBudget&&certificate.reason!=OcctStarReason::OutputBudget,NativeCadMeshReason::BudgetExceeded,"native-plc-source-certificate-budget");need(certificate.reason!=OcctStarReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");need(certificate.reason!=OcctStarReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");need(certificate.certified,NativeCadMeshReason::GeometryError,"native-plc-source-not-embedded");
        out.sourceShell.accepted=true;
        stage="seeds";trace.next(stage);out.surfaceSeconds=std::chrono::duration<double>(Clock::now()-began).count();
        const auto ns=boundary.nodeCount(),nf=boundary.cellCount();
        auto frontCells=std::min<std::uint64_t>(o.maximumCells-nf,
            progressiveResources?std::uint64_t(UINT_MAX):std::uint64_t(nf)*16);
        const auto apexStorage=completeDiscreteSeed?0:
            std::min<std::uint64_t>(256ULL*1024*1024,options.maximumPrimitiveWorkingBytes/4);
        std::uint64_t retained=progressiveResources?
            std::uint64_t(ns)*1024+std::uint64_t(nf)*4096+96ULL*1024*1024:
            std::uint64_t(ns)*1024+std::uint64_t(nf)*16384+frontCells*512+apexStorage+96ULL*1024*1024;
        need(retained<options.maximumPrimitiveWorkingBytes,NativeCadMeshReason::BudgetExceeded,"native-plc-front-storage-admission");
        // Initial placement releases its heap and proposal front before volume
        // recovery. Normal retains its checked numerical connectivity through
        // a separately admitted old-plus-new remap, then releases that owner.
        // Retain the boundary/keys/size field and existing kernel slack here;
        // the later volume stages keep their full, separate admission below.
        const std::uint64_t placementRetained=progressiveResources?
            std::uint64_t(ns)*1024+std::uint64_t(nf)*1024+96ULL*1024*1024:retained;
        if(progressiveResources){
            using occt_initial_point_front_detail::saturatedAdd;
            using occt_initial_point_front_detail::saturatedMultiply;
            const auto capacity=boundary.capacitySnapshot();
            need(capacity.complete,NativeCadMeshReason::BudgetExceeded,"native-plc-boundary-capacity-overflow");
            std::uint64_t bytes=0;
            for(const auto [count,unit]:std::array<std::pair<std::uint64_t,std::uint64_t>,9>{{
                {capacity.nodeCoordinateCapacity,sizeof(Vec3)},
                {capacity.nodeClassificationCapacity,sizeof(mesh::Classification)},
                {capacity.cellCapacity,sizeof(mesh::Cell)},
                {capacity.connectivityCapacity,sizeof(Index)},
                {capacity.dedupBucketCount,16},{capacity.dedupEntryCount,256},
                {capacity.dedupIndexCapacity,sizeof(Index)},
                {surface.nodeKeys.capacity(),sizeof(OcctBoundaryNodeKey)},
                {out.surfaceSpacing.capacity(),sizeof(CadSurfaceSpacingReport)}}})
                bytes=saturatedAdd(bytes,saturatedMultiply(count,unit));
            need(bytes<=placementRetained-96ULL*1024*1024,
                 NativeCadMeshReason::BudgetExceeded,"native-plc-placement-retained-capacity");
        }
        mesh::Mesh pendingRecovery;
        bnd::detail::TetNumericalTopology initialTopology;
        const auto buildConstraintSeed=[&](const mesh::Mesh& sites,mesh::Mesh& seed,
                                          std::uint64_t maximumSeedTets=UINT64_MAX,
                                          bool reserveRetry=false){
            OcctStageTrace recovering("constraint-recovery",tag);
            const auto priorWork=out.constraintRecovery.work;
            const auto priorSeconds=out.constraintRecovery.seconds;
            auto constraintWork=earlyConstraintWork
                ?std::min(remaining(),earlyConstraintWork)
                :progressiveResources?remaining()/3:std::min<std::uint64_t>(remaining()/3,450'000'000);
            if(reserveRetry)constraintWork=std::min(constraintWork,remaining()/3);
            const auto constraintRetained=placementRetained+std::uint64_t(sites.nodeCount())*128+
                std::uint64_t(sites.cellCount())*128;
            need(!progressiveResources||constraintRetained<options.maximumPrimitiveWorkingBytes,
                 NativeCadMeshReason::BudgetExceeded,"native-plc-constraint-input-memory");
            const auto constraintBytes=progressiveResources?
                options.maximumPrimitiveWorkingBytes-constraintRetained:options.maximumPrimitiveWorkingBytes/2;
            const auto constraints=bnd::detail::recoverTetConstraints(sites,{3,tag},
                constraintWork,constraintBytes,
                std::min(frontCells,maximumSeedTets),o.cancellation,seed,
                std::min<std::uint64_t>(64,o.maximumNodes-std::uint64_t(sites.nodeCount())),
                false,true,advancingInteriorFront,progressiveResources,
                initialTopology.available?&initialTopology:nullptr);
            const auto priorTopologyUsed=out.constraintRecovery.initialTopologyUsed;
            const auto priorTopologyDiagnostic=out.constraintRecovery.initialTopologyDiagnostic;
            const auto priorTopologyCells=out.constraintRecovery.initialTopologyCells;
            const auto priorTopologyWork=out.constraintRecovery.initialTopologyImportWork;
            const auto priorTopologyPeak=out.constraintRecovery.initialTopologyPeakModeledBytes;
            out.constraintRecovery={constraints.diagnostic,constraints.candidateCells,
                constraints.missingEdgesBefore,constraints.missingFacesBefore,
                constraints.missingEdgesAfter,constraints.missingFacesAfter,constraints.patches,
                constraints.seedCells,priorWork+constraints.work,priorSeconds+constraints.seconds,
                constraints.edgeWalks,constraints.facetWalks,constraints.walkFallbacks,
                constraints.steinerAttempts,constraints.steinerWork,constraints.steinerNodes,
                constraints.retainedSteinerNodes,constraints.piercingPatches};
            out.constraintRecovery.initialTopologyUsed=priorTopologyUsed||constraints.initialTopologyUsed;
            out.constraintRecovery.initialTopologyDiagnostic=
                std::string_view(constraints.initialTopologyDiagnostic)=="not_available"?
                    priorTopologyDiagnostic:constraints.initialTopologyDiagnostic;
            out.constraintRecovery.initialTopologyCells=std::max(priorTopologyCells,constraints.initialTopologyCells);
            out.constraintRecovery.initialTopologyImportWork=priorTopologyWork+constraints.initialTopologyImportWork;
            out.constraintRecovery.initialTopologyPeakModeledBytes=std::max(priorTopologyPeak,constraints.initialTopologyPeakModeledBytes);
            // Any refused handoff remains visible. A later boundary-only
            // fallback may not retain its old private allocation or callbacks.
            initialTopology=bnd::detail::TetNumericalTopology{};
            recovering.finish(constraints.diagnostic);
            // The private recovery has already performed this work, including
            // the last cancelled attempt. Record it before observing the stop.
            accountCompleted(constraints.work);
            need(std::string_view(constraints.diagnostic)!="cancelled",NativeCadMeshReason::Cancelled,"cancelled");
            need(std::string_view(constraints.diagnostic)!="allocation_failed",NativeCadMeshReason::AllocationFailed,"allocation-failed");
            const std::string_view topologyDiagnostic=constraints.initialTopologyDiagnostic;
            need(!topologyDiagnostic.starts_with("topology_")||topologyDiagnostic.ends_with("_budget"),
                 NativeCadMeshReason::InternalError,"native-plc-topology-handoff-invalid");
        };
        OcctSolidClassifier classifier(OcctModelAccess::solid(model,tag));
        const auto tol=model.tolerance(3,tag);
        if(spatialSizing&&(options.improveNativeTetQuality||advancingInteriorFront)){
            // Construct the first interior point set before constraint recovery.
            // The prescribed surface remains the same input constraint complex.
            OcctStageTrace placement("interior_sites",tag);
            const auto placementStarted=Clock::now();
            struct SiteLimit {};
            mesh::Mesh sites;
            const auto& construction=options.construction;
            const bool fairConstruction=construction.effort!=CadConstructionEffort::Legacy;
            const bool expandedConstruction=progressiveResources||fairConstruction||construction.interiorCandidates||construction.interiorSites;
            const auto siteBudget=advancingInteriorFront?
                (expandedConstruction?(construction.effort==CadConstructionEffort::Thorough?
                    remaining()-remaining()/3:remaining()/2):std::min<std::uint64_t>(42'000'000,remaining()/2)):
                std::min<std::uint64_t>(1'000'000,remaining()/8);
            std::uint64_t siteWork=0,siteQueries=0;
            const auto siteCharge=[&](std::uint64_t amount=1){
                if(siteWork>siteBudget||amount>siteBudget-siteWork)throw SiteLimit{};
                tick(amount);siteWork+=amount;
            };
            // The initial-front helper checks its stop token before ordinary
            // operations. Its bootstrap also sends completed constructor work
            // through this callback before propagating its observed stop.
            const auto siteAccount=[&](std::uint64_t amount){
                if(siteWork>siteBudget||amount>siteBudget-siteWork)throw SiteLimit{};
                account(amount);siteWork+=amount;
            };
            constexpr auto nodeBytes=sizeof(Vec3)+sizeof(mesh::Classification);
            constexpr auto tetBytes=sizeof(mesh::Cell)+4*sizeof(Index);
            const auto boundaryBytes=std::uint64_t(ns)*nodeBytes+
                std::uint64_t(nf)*(sizeof(mesh::Cell)+3*sizeof(Index));
            const auto outputRoom=boundaryBytes<o.maximumOutputBytes?
                o.maximumOutputBytes-boundaryBytes:0;
            using occt_initial_point_front_detail::saturatedAdd;
            using occt_initial_point_front_detail::saturatedMultiply;
            const auto maximumSites=std::min(o.maximumNodes-std::uint64_t(ns),
                outputRoom/(nodeBytes+16*tetBytes));
            std::uint64_t requestedSites=512;
            if(progressiveResources)requestedSites=construction.interiorSites?construction.interiorSites:maximumSites;
            else if(expandedConstruction){
                if(construction.interiorSites)requestedSites=construction.interiorSites;
                else if(!fairConstruction){
                    const auto factor=construction.interiorCandidates/16'384+
                        (construction.interiorCandidates%16'384!=0);
                    requestedSites=saturatedMultiply(512,std::max<std::uint64_t>(1,factor));
                }else{
                    // Budget estimates never place points. Integrate bounded
                    // actual field samples rather than imposing one local fine
                    // control's minimum over the entire body's bounding box.
                    const auto effort=construction.effort==CadConstructionEffort::Thorough?4.:2.;
                    const auto densityAt=[&](double h){
                        if(!maximumSites)return static_cast<long double>(0.);
                        long double logarithm=std::log(effort);
                        for(const auto span:std::array{hi.x-lo.x,hi.y-lo.y,hi.z-lo.z}){
                            if(!(span>0.))return static_cast<long double>(0.);
                            logarithm+=std::log(static_cast<long double>(span))-std::log(static_cast<long double>(h));
                        }
                        return logarithm>=std::log(static_cast<long double>(maximumSites))?
                            static_cast<long double>(maximumSites):std::exp(logarithm);
                    };
                    auto demand=densityAt(o.targetSize);
                    const auto queryWork=saturatedAdd(sizeField->queryWork(),1);
                    const auto sampleWork=saturatedMultiply(queryWork,27);
                    if((!options.faceSizes.empty()||!options.sizeBoxes.empty())&&sampleWork<=siteBudget/64){
                        long double integrated=0.;
                        for(unsigned x=0;x<3;++x)for(unsigned y=0;y<3;++y)for(unsigned z=0;z<3;++z){
                            const auto coordinate=[](double low,double high,unsigned i){
                                const auto f=(double(i)+.5)/3.;return low*(1-f)+high*f;
                            };
                            const Vec3 p{coordinate(lo.x,hi.x,x),coordinate(lo.y,hi.y,y),coordinate(lo.z,hi.z,z)};
                            siteCharge(queryWork);const auto h=sizeField->eval(p);
                            need(h>0.&&std::isfinite(h),NativeCadMeshReason::GeometryError,"invalid-initial-density-field");
                            integrated+=densityAt(h)/27.;
                        }
                        demand=integrated;
                    }
                    const auto density=std::isfinite(demand)&&demand<static_cast<long double>(maximumSites)?
                        std::uint64_t(std::ceil(demand)):maximumSites;
                    requestedSites=std::max({std::uint64_t(512),density,saturatedMultiply(ns,2)});
                    const auto memory=options.maximumPrimitiveWorkingBytes-placementRetained;
                    requestedSites=std::min({requestedSites,maximumSites,(siteBudget-siteWork)/8192});
                    const auto modeledBytes=[&](std::uint64_t sites){
                        const auto cones=std::min(o.maximumCells,std::max<std::uint64_t>(2048,saturatedMultiply(sites,4)));
                        const auto states=saturatedAdd(nf,saturatedMultiply(cones,3));
                        const auto front=occt_initial_point_front_detail::initialFrontStorageBytes(ns,nf,sites,cones,states);
                        const auto cells=std::min(o.maximumCells,
                            saturatedMultiply(saturatedAdd(saturatedAdd(ns,sites),8),8));
                        return saturatedAdd(front,saturatedAdd(32ULL*1024*1024,saturatedMultiply(cells,8320)));
                    };
                    std::uint64_t low=0,high=requestedSites;
                    while(low<high){const auto middle=low+(high-low)/2+1;
                        if(modeledBytes(middle)<=memory)low=middle;else high=middle-1;}
                    requestedSites=low;
                }
            }
            const auto siteCap=std::min(requestedSites,maximumSites);
            // Bound the copied shell, point/hash arrays and cached edge sizes.
            const auto preparationBytes=saturatedMultiply(saturatedAdd(saturatedAdd(ns,nf),progressiveResources?std::uint64_t(0):siteCap),1024);
            const auto candidateLimit=progressiveResources?(construction.interiorCandidates?construction.interiorCandidates:siteBudget):expandedConstruction?(construction.interiorCandidates?
                construction.interiorCandidates:
                std::max<std::uint64_t>(16'384,saturatedAdd(saturatedMultiply(siteCap,32),saturatedMultiply(nf,12)))):16'384;
            const auto queryLimit=expandedConstruction?std::min(siteBudget,
                std::max<std::uint64_t>(12'000,saturatedMultiply(candidateLimit,16))):12'000;
            if(advancingInteriorFront){
                auto& report=out.initialFront;
                report.diagnostic="budget_limited";report.siteLimit=siteCap;
                report.candidateLimit=candidateLimit;report.queryLimit=queryLimit;report.workLimit=siteBudget;
                report.budgetReason=siteBudget<=siteWork?"work_limit":!siteCap?
                    (o.maximumNodes<=std::uint64_t(ns)?"node_limit":!maximumSites?"output_limit":"modeled_memory_or_work"):
                    preparationBytes>options.maximumPrimitiveWorkingBytes-placementRetained?"modeled_memory":"initial_admission";
            }
            OcctInitialPointFrontReport built;
            const auto recordInitialPlacement=[&](){
                    auto& report=out.initialFront;
                    report.budgetReason=built.budgetReason;
                    report.candidateLimit=built.candidateLimit;report.siteLimit=built.siteLimit;
                    report.queryLimit=built.queryLimit;report.workLimit=siteBudget;
                    report.reservedBulkCandidates=built.reservedBulkCandidates;
                    report.reservedBulkSites=built.reservedBulkSites;
                    report.sourceCandidates=built.sourceCandidates;
                    report.sourceDeferredFaces=built.sourceDeferredFaces;
                    report.sourcePhaseDeferred=built.sourcePhaseDeferred;
                    report.diagnostic=built.budgetLimited?"budget_limited":!built.initialized?
                        "unavailable":built.frontClosed?"front_closed":"stalled";
                    report.nodes=built.sites;report.candidates=built.candidates;
                    report.wallNodes=built.sourcePhaseSites;report.voidNodes=built.bulkSites;
                    report.voidCandidates=built.bulkCandidates;
                    report.sourceSweepComplete=built.sourcePhaseCompleted;
                    report.voidQueueExhausted=built.bulkQueueExhausted;
                    report.maximumUnresolvedVoidRatio=built.bulkMaximumUnresolvedRatio;
                    report.cones=built.seedCones;report.reusedApices=built.reusedApices;
                    report.sourceActivations=built.originalFaceActivations;
                    report.propagatedActivations=built.propagatedFaceActivations;
                    report.propagatedSites=built.deepSites;report.maximumDepth=built.maximumDepth;
                    report.remainingFaces=built.remainingFrontFaces;
                    report.delaunayEvaluations=built.auxiliaryScorerEvaluations;
                    report.delaunayCommits=built.auxiliaryCommits;
                    report.fullCertifications=built.auxiliaryFullCertifications;
                    report.cavityUnavailable=built.auxiliaryScoreCavityUnavailable;
                    report.cavityLimit=built.auxiliaryCavityLimit;
                    report.maximumPlannedCavity=built.auxiliaryMaximumPlannedCavity;
                    report.maximumAcceptedCavity=built.auxiliaryMaximumAcceptedCavity;
                    report.storageRejections=built.auxiliaryScoreStorageRejects+built.auxiliaryCommitStorageRejects;
                    report.compactions=built.auxiliaryCompactions;
                    report.retiredCellsReclaimed=built.auxiliaryDeadCellsReclaimed;
                    report.compactionWork=built.auxiliaryCompactionWork;
                    report.auxiliaryStoredCells=built.auxiliaryStoredCells;
                    report.auxiliaryLiveCells=built.auxiliaryLiveCells;
                    report.auxiliaryStoredPeak=built.auxiliaryStoredPeak;
                    report.auxiliaryStorageCapacity=built.auxiliaryStorageCapacity;
                    report.compactionPeakModeledBytes=built.auxiliaryMaximumCompactionModeledBytes;
                    report.auxiliaryWorkingBytes=built.auxiliaryWorkingBytes;
                    report.auxiliaryImportedCells=built.auxiliaryImportedCells;
                    report.auxiliaryImportWork=built.auxiliaryImportWork;
                    report.auxiliaryImportModeledBytes=built.auxiliaryImportModeledBytes;
                    report.auxiliaryInitialCapacity=built.auxiliaryInitialCapacity;
                    report.auxiliaryHardCapacity=built.auxiliaryHardCapacity;
                    report.auxiliaryGrowthAttempts=built.auxiliaryGrowthAttempts;
                    report.auxiliaryGrowths=built.auxiliaryGrowths;
                    report.auxiliaryGrowthWork=built.auxiliaryGrowthWork;
                    report.maximumCombinedModeledBytes=built.maximumCombinedModeledBytes?
                        saturatedAdd(placementRetained,built.maximumCombinedModeledBytes):0;
                    report.maximumHeapCapacity=built.maximumHeapCapacity;
                    report.resourceRetries=built.resourceRetries;
                    report.auxiliaryResourceDiagnostic=built.auxiliaryResourceDiagnostic;
                    report.geometryRefusalCacheHits=built.geometryRefusalCacheHits;
                    report.cachedGeometryRefusals=built.cachedGeometryRefusals;
            };
            try {if(siteCap&&siteBudget&&preparationBytes<=options.maximumPrimitiveWorkingBytes-placementRetained){
                if(advancingInteriorFront){
                    OcctInitialPointFrontLimits limits;
                    limits.work=siteBudget-siteWork;limits.queries=queryLimit;limits.sites=std::size_t(siteCap);
                    limits.balanceSourceAndBulk=fairConstruction;
                    limits.progressiveResources=progressiveResources;
                    limits.cacheGeometryRefusals=progressiveResources;
                    limits.cancellation=o.cancellation;
                    limits.workingBytes=std::min<std::uint64_t>(128ULL*1024*1024,
                        options.maximumPrimitiveWorkingBytes-placementRetained);
                    if(expandedConstruction){
                        limits.candidates=std::size_t(std::min<std::uint64_t>(candidateLimit,SIZE_MAX));
                        limits.workingBytes=options.maximumPrimitiveWorkingBytes-placementRetained;
                        limits.cones=std::size_t(std::min(o.maximumCells,
                            std::max<std::uint64_t>(2048,saturatedMultiply(siteCap,4))));
                        limits.frontStates=std::size_t(std::min<std::uint64_t>(SIZE_MAX,
                            saturatedAdd(nf,saturatedMultiply(limits.cones,3))));
                        limits.maximumAuxiliaryCells=std::size_t(std::min<std::uint64_t>(
                            std::min(o.maximumCells,limits.work/24),SIZE_MAX));
                        if(progressiveResources){
                            limits.maximumAuxiliaryCells=std::size_t(std::min<std::uint64_t>(o.maximumCells,SIZE_MAX));
                            limits.cones=std::size_t(std::min<std::uint64_t>(o.maximumCells,SIZE_MAX));
                            limits.frontStates=std::size_t(std::min<std::uint64_t>(SIZE_MAX,
                                saturatedAdd(nf,saturatedMultiply(limits.cones,3))));
                        }
                    }
                    sites=appendOcctInitialPointFrontSites(boundary,*sizeField,o.targetSize,
                        mesh::Classification{3,tag},limits,[&](Vec3 p){
                            // The builder charges and bounds both point and
                            // cached material-witness queries before this call.
                            ++siteQueries;++out.nativeQueries;
                            classifier.Perform(gp_Pnt(p.x,p.y,p.z),tol);
                            return classifier.State()==TopAbs_IN;
                        },siteAccount,&built,
                        (progressiveResources||options.initialTopologyQuality.value_or(false))?&initialTopology:nullptr);
                    recordInitialPlacement();
                    need(!built.cancelled,NativeCadMeshReason::Cancelled,"cancelled");
                }else{
                    sites=appendOcctInteriorSites(boundary,*sizeField,mesh::Classification{3,tag},
                        std::size_t(siteCap),[&](Vec3 p){
                            siteCharge();
                            if(siteQueries>=2048)throw SiteLimit{};
                            ++siteQueries;++out.nativeQueries;
                            classifier.Perform(gp_Pnt(p.x,p.y,p.z),tol);
                            return classifier.State()==TopAbs_IN;
                        },siteCharge);
                }
            }
            }catch(const OcctInitialConstructionCancelled&){
                recordInitialPlacement();
                out.initialFront.diagnostic="cancelled";
                out.initialFront.budgetReason="cancelled";
                out.initialFront.work=siteWork;out.initialFront.queries=siteQueries;
                out.initialFront.seconds=std::chrono::duration<double>(Clock::now()-placementStarted).count();
                throw Stop{NativeCadMeshReason::Cancelled,"cancelled"};
            }catch(const SiteLimit&){
                // This optional preparation leaves the original shell intact.
                initialTopology=bnd::detail::TetNumericalTopology{};
                if(advancingInteriorFront)out.initialFront.budgetReason="work_limit";
            }catch(const std::bad_alloc&){
                // Preparation is transactional; retain the ordinary seed route.
                // Returning/assigning Mesh may allocate after the no-throw
                // topology move. Never pair that cloud with a different source.
                initialTopology=bnd::detail::TetNumericalTopology{};
                if(advancingInteriorFront)out.initialFront.budgetReason="allocation_failed";
            }
            if(advancingInteriorFront){
                out.initialFront.work=siteWork;out.initialFront.queries=siteQueries;
                out.initialFront.seconds=std::chrono::duration<double>(Clock::now()-placementStarted).count();
            }
            placement.finish(sites.nodeCount()>ns?"seeded":"unchanged");
            if(initialTopology.available&&options.initialTopologyQuality.value_or(progressiveResources)){
                OcctStageTrace finalizing("initial_topology_quality",tag);
                const auto qualityStarted=Clock::now();
                OcctInitialTopologyQualityLimits limits;
                limits.maximumWork=progressiveResources?remaining()/8:
                    std::min<std::uint64_t>(1'000'000,remaining()/8);
                const auto held=placementRetained+std::uint64_t(sites.nodeCount())*128+
                    std::uint64_t(sites.cellCount())*128;
                limits.maximumWorkingBytes=held<options.maximumPrimitiveWorkingBytes?
                    options.maximumPrimitiveWorkingBytes-held:0;
                limits.cancellation=o.cancellation;
                const auto finished=finalizeOcctInitialTopologyQuality(initialTopology,limits);
                out.initialTopologyQuality={finished.attempted,finished.completed,finished.changed,
                    finished.cancelled,finished.budgetLimited,finished.summaryAvailable,finished.candidateListTruncated,
                    finished.diagnostic,finished.work,finished.maximumModeledBytes,finished.candidates,
                    finished.proposals,finished.certifications,finished.stagedEdits,finished.edits,
                    finished.flips23,finished.flips32,finished.flips44,finished.liveBefore,finished.liveAfter,
                    finished.severeBefore,finished.severeAfter,finished.protectedSourceFaces,finished.protectedSourceEdges,
                    finished.minimumMeanRatioBefore,finished.minimumMeanRatioAfter,
                    finished.minimumScaledJacobianBefore,finished.minimumScaledJacobianAfter,
                    std::chrono::duration<double>(Clock::now()-qualityStarted).count()};
                out.initialTopologyQuality.searchBelowMeanRatio=limits.searchBelowMeanRatio;
                out.initialTopologyQuality.workLimit=limits.maximumWork;
                finalizing.finish(finished.diagnostic);
                accountCompleted(finished.work);
                need(!finished.cancelled,NativeCadMeshReason::Cancelled,"cancelled");
                need(std::string_view(finished.diagnostic)!="allocation_failed",NativeCadMeshReason::AllocationFailed,"allocation-failed");
                need(std::string_view(finished.diagnostic)!="invalid_topology",NativeCadMeshReason::InternalError,"native-initial-topology-invalid");
                bnd::detail::observeTetConstruction({"initial_topology_finalized",tag,&initialTopology.complex,
                    nullptr,{},std::size_t(ns),8,0,initialTopology.complex.cells.size(),{},true,
                    initialTopology.originalCoordinateRetryUsed});
            }
            if(sites.nodeCount()>ns){
                const auto siteBytes=std::uint64_t(sites.nodeCount())*nodeBytes+
                    std::uint64_t(nf)*(sizeof(mesh::Cell)+3*sizeof(Index));
                const auto steinerBytes=std::min<std::uint64_t>(64,
                    o.maximumNodes-std::uint64_t(sites.nodeCount()))*nodeBytes;
                const auto seedTets=siteBytes<o.maximumOutputBytes&&
                    steinerBytes<o.maximumOutputBytes-siteBytes?
                    (o.maximumOutputBytes-siteBytes-steinerBytes)/tetBytes:0;
                if(seedTets)buildConstraintSeed(sites,pendingRecovery,seedTets,true);
                if(!pendingRecovery.cellCount()){
                    OcctStageTrace retry("interior_seed_fallback",tag);
                    buildConstraintSeed(boundary,pendingRecovery);
                    retry.finish("boundary_only");
                }
            }else buildConstraintSeed(boundary,pendingRecovery);
        }else if(spaceSurface||options.recoveryMode==NativeCadRecoveryMode::Advanced)
            buildConstraintSeed(boundary,pendingRecovery);
        const bool usesConstraintSeed=pendingRecovery.cellCount()!=0;
        if(progressiveResources){
            if(usesConstraintSeed){
                // This route never initializes the outer gift-wrapping front.
                // Keep the actual seed and its copies/filters while recovery
                // admits its own staging workspace; no future-cell quote.
                retained+=std::uint64_t(std::max<Index>(0,pendingRecovery.nodeCount()-ns))*1024+
                    std::uint64_t(pendingRecovery.cellCount())*1024;
            }else{
                // The legacy frontier can coexist with a recovery attempt.
                // Derive its ceiling from the byte allowance, rather than a
                // fixed surface-to-volume element-count ratio.
                const std::uint64_t fixed=std::uint64_t(ns)*1024+std::uint64_t(nf)*16384+apexStorage+96ULL*1024*1024;
                need(fixed<options.maximumPrimitiveWorkingBytes,NativeCadMeshReason::BudgetExceeded,
                     "native-plc-front-storage-admission");
                frontCells=std::min(frontCells,(options.maximumPrimitiveWorkingBytes-fixed)/1024);
                retained=fixed+frontCells*512;
            }
            need(retained<options.maximumPrimitiveWorkingBytes,NativeCadMeshReason::BudgetExceeded,
                 "native-plc-seed-storage-admission");
        }
        if(advancingInteriorFront&&!usesConstraintSeed){
            const std::string_view diagnostic=out.constraintRecovery.diagnostic;
            need(diagnostic!="construction_budget"&&diagnostic!="work_budget"&&
                 diagnostic!="storage_budget"&&diagnostic!="admission_budget"&&
                 diagnostic!="memory_budget"&&diagnostic!="cell_budget"&&
                 diagnostic!="node_index_budget",
                 NativeCadMeshReason::BudgetExceeded,"native-plc-constraint-construction-budget");
        }
        need(!spaceSurface||usesConstraintSeed,NativeCadMeshReason::MeshingFailed,"spaced-surface-constraint-seed-refused");
        const bool reuseStrictAdmissions=(options.useAdaptiveTetRecovery||spaceSurface)&&usesConstraintSeed&&
            options.maximumNativeValidationThreads==1&&!options.allowCentroidsWithinChordTolerance;
        OcctNativeStrictAdmissions strictAdmissions(std::min<std::uint64_t>(
            16ULL*1024*1024,(options.maximumPrimitiveWorkingBytes-retained)/16));
        stage="seeds";trace.next(stage);
        bnd::detail::PlcGiftWrapping front;
        front.admitApexStorage=[&](std::uint64_t bytes){need(bytes<=apexStorage,NativeCadMeshReason::BudgetExceeded,"native-plc-apex-storage-budget");};
        front.cancellation=stop.get_token();front.chargeWork=tick;
        // Preserve the established global ordering when its candidate arrays
        // fit. Restrict the first scan only for surfaces needing more storage.
        front.useLocalCandidates=(std::uint64_t(ns)+64)*std::uint64_t(nf)>
            apexStorage/sizeof(std::pair<double,bnd::TetFacetVertexId>);
        const auto pointKey=[](Vec3 p){return std::array{std::bit_cast<std::uint64_t>(p.x),
            std::bit_cast<std::uint64_t>(p.y),std::bit_cast<std::uint64_t>(p.z)};};
        const std::uint64_t maximumSpacedCandidateQueries=customConstruction?maximumWork:
            curvatureSurface||spatialSizing?(options.useLocalDelaunayRefinement&&liveTrial?40'000:24'000):3'000;
        const auto admitQueries=[&](std::uint64_t n) {
            need(!spaceSurface||(out.nativeQueries<=maximumSpacedCandidateQueries&&
                 n<=maximumSpacedCandidateQueries-out.nativeQueries),
                 NativeCadMeshReason::BudgetExceeded,"spaced-surface-candidate-query-budget");
        };
        std::map<std::array<std::uint64_t,3>,bool> nativeCache;
        const auto native=[&](Vec3 p){tick();const auto key=pointKey(p);auto known=nativeCache.find(key);if(known!=nativeCache.end())return known->second;
            need(progressiveResources||nativeCache.size()<1000000,NativeCadMeshReason::BudgetExceeded,"native-plc-query-cache-budget");admitQueries(1);++out.nativeQueries;classifier.Perform(gp_Pnt(p.x,p.y,p.z),tol);const bool yes=classifier.State()==TopAbs_IN||classifier.State()==TopAbs_ON;
            if(nativeCache.size()<(progressiveResources?131072U:1000000U))nativeCache.emplace(key,yes);
            return yes;};
        std::map<std::array<std::uint64_t,3>,bool> centroidCache;
        const auto nativeCentroid=[&](Vec3 p){
            const auto key=pointKey(p);
            if(usesConstraintSeed){
                tick();
                if(auto known=centroidCache.find(key);known!=centroidCache.end())return known->second;
            }
            if(native(p))return true;
            if(!options.allowCentroidsWithinChordTolerance)return false;
            if(auto known=centroidCache.find(key);known!=centroidCache.end())return known->second;
            need(progressiveResources||centroidCache.size()<1000000,NativeCadMeshReason::BudgetExceeded,"native-plc-query-cache-budget");tick();
            admitQueries(1);classifier.Perform(gp_Pnt(p.x,p.y,p.z),tol);++out.nativeQueries;
            if(classifier.State()!=TopAbs_OUT){
                if(centroidCache.size()<(progressiveResources?131072U:1000000U))centroidCache.emplace(key,false);
                return false;
            }
            tick();BRepExtrema_DistShapeShape distance(BRepBuilderAPI_MakeVertex(gp_Pnt(p.x,p.y,p.z)).Vertex(),OcctModelAccess::solid(model,tag));
            need(distance.IsDone()&&std::isfinite(distance.Value())&&distance.Value()>=0,NativeCadMeshReason::GeometryError,"native-plc-centroid-distance-query");
            const bool admitted=distance.Value()<=out.chord;
            if(centroidCache.size()<(progressiveResources?131072U:1000000U))centroidCache.emplace(key,admitted);
            return admitted;
        };
        if(!usesConstraintSeed){
        for(Index i=0;i<ns;++i){tick();front.vertices.push_back({boundary.node(i),std::uint64_t(i)+1});}
        for(Index i=0;i<nf;++i){tick();const auto* f=boundary.cellNodes(i);front.boundary.push_back({unsigned(f[0]),unsigned(f[1]),unsigned(f[2])});}
        // Seeds are proposed by PAACS's own Delaunay builder, never a reference mesher.
        // Seed generation does not recover the boundary. Bound its Delaunay
        // input on dense tessellations; the complete source stays in front.
        const auto seedInputs=occtSeedBoundarySample(boundary.nodes(),lo,hi,tick);
        std::vector<Vec3> points;points.reserve(seedInputs.size()+8);for(auto i:seedInputs)points.push_back(boundary.node(i));
        const auto seedBoundaryCount=points.size();const auto span=std::max({hi.x-lo.x,hi.y-lo.y,hi.z-lo.z});
        for(unsigned i=0;i<8;++i)points.push_back({(i&1)?hi.x+span:lo.x-span,(i&2)?hi.y+span:lo.y-span,(i&4)?hi.z+span:lo.z-span});
        auto seedCandidate=algo::buildTetCandidate(points);accountCompleted(seedCandidate.metrics.construction.conflictInsphereTests+seedCandidate.metrics.construction.conflictOrientationTests+seedCandidate.metrics.validationFaceUses+seedCandidate.metrics.validationHullTests+seedCandidate.metrics.candidateCells);
        need(seedCandidate.reason!=algo::TetCandidateReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");need(seedCandidate.built(),NativeCadMeshReason::MeshingFailed,"native-plc-seed-candidate-refused");
        struct Proposal{double volume;Vec3 point;std::size_t ordinal;};std::vector<Proposal> proposals;
        for(std::size_t i=0;i<seedCandidate.candidate->cells.size();++i){tick();const auto c=seedCandidate.candidate->cells[i];bool guard=false;std::array<Vec3,4> p;for(unsigned j=0;j<4;++j){guard|=std::size_t(seedCandidate.candidate->pointToInput[c.vertices[j]])>=seedBoundaryCount;p[j]=seedCandidate.candidate->points[c.vertices[j]];}if(!guard)proposals.push_back({(p[1]-p[0]).cross(p[2]-p[0]).dot(p[3]-p[0]),(p[0]+p[1]+p[2]+p[3])*.25,i});}
        std::sort(proposals.begin(),proposals.end(),[](const auto& a,const auto& b){return a.volume!=b.volume?a.volume>b.volume:a.ordinal<b.ordinal;});unsigned added=0;
        for(const auto& p:proposals){if(added==64)break;if(!native(p.point))continue;bool duplicate=false;for(const auto& v:front.vertices){tick();duplicate|=(v.position-p.point).norm()<=tol;}if(!duplicate){front.vertices.push_back({p.point,front.vertices.size()+1});++added;}}
        need(front.vertices.size()<=o.maximumNodes,NativeCadMeshReason::BudgetExceeded,"native-plc-seed-node-budget");need(front.initialize(),NativeCadMeshReason::GeometryError,"native-plc-front-initialize");
        }
        bnd::TetPlcCavityOptions recovery;recovery.region={3,tag};recovery.maximumNodes=o.maximumNodes;recovery.maximumCells=frontCells;recovery.maximumSourceTriangles=bo.maximumTriangles;recovery.maximumWorkingBytes=options.maximumPrimitiveWorkingBytes-retained;recovery.maximumPatches=512;recovery.maximumCandidatePoints=128;recovery.maximumExpansionCells=512;recovery.maximumPatchCells=2048;recovery.maximumComponentFaces=2048;recovery.warningQuality=o.warningElementQuality;recovery.cancellation=stop.get_token();
        recovery.maximumLocalFrontTets=usesConstraintSeed?512:0;
        const auto partialFrom=[&](const bnd::detail::PlcGiftWrapping& value){mesh::Mesh partial;for(Index i=0;i<Index(value.vertices.size());++i){tick();partial.addNode(value.vertices[i].position,i<boundary.nodeCount()?boundary.nodeClass(i):mesh::Classification{3,tag});}for(auto c:value.cells){tick();std::array<Index,4> ids;for(unsigned j=0;j<4;++j)ids[j]=c.vertices[j];partial.addCell(mesh::CellType::Tet4,ids,{3,tag});}return partial;};
        mesh::Mesh* evolvingPartial=nullptr;
        const auto recover=[&](const mesh::Mesh& partial,mesh::Mesh& result,bool diversified=false,
                               bool combineRefill=false){
            mesh::Mesh continuation;
            const auto collect=[&](const bnd::TetPlcCavityOutcome& recovered){
                auto& profile=out.regionRecovery;
                profile.attempts+=recovered.regionSearchAttempts;profile.states+=recovered.regionSearchStates;profile.backtracks+=recovered.regionSearchBacktracks;profile.candidates+=recovered.regionSearchCandidates;profile.pairCacheHits+=recovered.regionSearchPairCacheHits;profile.patches+=recovered.regionSearchPatches;profile.boundHits+=recovered.regionSearchBoundHits;
                profile.vertexBoundHits+=recovered.regionVertexBoundHits;profile.stateBoundHits+=recovered.regionStateBoundHits;profile.largestVertices=std::max(profile.largestVertices,recovered.largestRegionVertices);profile.largestFaces=std::max(profile.largestFaces,recovered.largestRegionFaces);
                profile.searchSeconds+=recovered.regionSearchSeconds;profile.preparationSeconds+=recovered.preparationSeconds;profile.visibilitySeconds+=recovered.visibilitySeconds;profile.publicationSeconds+=recovered.publicationSeconds;profile.validationSeconds+=recovered.validationSeconds;
                profile.localFrontTets+=recovered.localFrontTets;profile.localFrontClosed+=recovered.localFrontClosed;
                profile.localFrontSeconds+=recovered.localFrontSeconds;
                out.randomCandidateProposals+=recovered.randomCandidatePointsProposed;out.mergedComponents+=recovered.mergedComponents;out.sourceBlockedCandidates+=recovered.sourceBlockedCandidates;out.otherComponentBlockedCandidates+=recovered.otherComponentBlockedCandidates;
                accountCompleted(recovered.workUsed);
            };
            recovery.randomSeed=bnd::detail::PlcGiftWrapping::mix(options.recoverySeed+out.randomRetries);
            if(diversified)++out.randomRetries;
            const auto attempt=[&](const mesh::Mesh& input,bool regions,unsigned vertices=256,
                                   bool includeVisibility=false){
                recovery.maximumWork=includeVisibility?std::max<std::uint64_t>(1,remaining()/2):remaining();
                // Keep the established inexpensive fill before combinatorial
                // repair, and reserve work for the latter when it is enabled.
                if(!regions&&options.maximumRecoveryPatchStates)
                    recovery.maximumWork=std::min<std::uint64_t>(recovery.maximumWork,20000000);
                recovery.patchSearch.maximumStates=regions?
                    (std::string_view(stage)=="region-search"?options.maximumRecoveryPatchStates:std::min(options.maximumRecoveryPatchStates,512U)):0;
                recovery.patchSearch.maximumVertices=vertices;recovery.patchSearch.maximumCandidateTets=32768;recovery.patchSearch.maximumShells=8;
                recovery.useVisibilityExpansion=!regions||includeVisibility;
                recovery.maximumRandomCandidatePoints=diversified?64:0;
                recovery.maximumMergedComponents=(diversified||regions)?4:0;
                OcctStageTrace refill("cavity_refill",tag);
                mesh::Mesh advanced;auto recovered=bnd::recoverTetPlcCavitiesRetaining(boundary,input,recovery,result,&advanced);
                refill.finish(recovered.diagnostic);
                collect(recovered);
                if(!recovered.published&&advanced.cellCount()&&!includeVisibility)continuation.swap(advanced);
                const auto why=std::string_view(recovered.diagnostic);
                if(recovered.reason==bnd::TetPlcCavityReason::BudgetExceeded&&out.workUsed<maximumWork&&
                    ((regions&&recovered.regionSearchBoundHits&&why=="no-admissible-region-refill")||
                     (!regions&&options.maximumRecoveryPatchStates&&why=="work-or-patch-budget")))
                    recovered.reason=bnd::TetPlcCavityReason::Unavailable;
                return recovered;
            };
            const bool finalRegions=std::string_view(stage)=="region-search"&&options.maximumRecoveryPatchStates;
            // A complete constraint seed with a few CAD removals already
            // localizes the cavities. Use one certified front for both
            // visibility and region refill before retrying the usual order.
            auto recovered=attempt(partial,finalRegions||combineRefill,regionVertices,combineRefill);
            if(combineRefill&&!recovered.published&&out.workUsed<maximumWork&&
                (recovered.reason==bnd::TetPlcCavityReason::Unavailable||
                 recovered.reason==bnd::TetPlcCavityReason::BudgetExceeded||
                 recovered.reason==bnd::TetPlcCavityReason::InternalError)) {
                recovered=attempt(partial,finalRegions);
            }
            // Preserve the established bounded search first. Enlarge only
            // when that final search actually reaches its vertex ceiling;
            // admitting a larger search everywhere can spend the fallback's
            // budget exploring regions that the smaller route closes quickly.
            if(finalRegions&&!recovered.published&&recovered.regionVertexBoundHits&&expandedRegionVertices>regionVertices&&
                out.workUsed<maximumWork&&(recovered.reason==bnd::TetPlcCavityReason::Unavailable||recovered.reason==bnd::TetPlcCavityReason::BudgetExceeded)){
                const auto closed=recovered.cavitiesClosed;
                recovered=attempt(continuation.cellCount()?continuation:partial,true,expandedRegionVertices);
                recovered.cavitiesClosed+=closed;
            }
            if(!recovered.published&&!finalRegions&&options.maximumRecoveryPatchStates&&
                recovered.reason==bnd::TetPlcCavityReason::Unavailable&&out.workUsed<maximumWork){
                const auto closed=recovered.cavitiesClosed;
                recovered=attempt(continuation.cellCount()?continuation:partial,true);
                recovered.cavitiesClosed+=closed;
            }
            if(!recovered.published&&continuation.cellCount()){
                if(evolvingPartial&&evolvingPartial==&partial)evolvingPartial->swap(continuation);
                else pendingRecovery.swap(continuation);
            }
            recovery.maximumRandomCandidatePoints=0;return recovered;
        };
        mesh::Mesh completed;
        const auto fillDiscreteSeed=[&](const mesh::Mesh& seed,mesh::Mesh& filled){
            if(!seed.cellCount()||seed.cellCount()>12000)return false;
            OcctStageTrace trial("discrete_fill",tag);
            auto limits=recovery;
            limits.maximumWork=std::min<std::uint64_t>(25'000'000,remaining()/2);
            limits.maximumLocalFrontTets=256;
            limits.patchSearch={8192,32768,256,4};
            const auto closed=bnd::recoverTetPlcCavities(boundary,seed,limits,filled);
            trial.finish(closed.diagnostic);
            accountCompleted(closed.workUsed);
            need(closed.reason!=bnd::TetPlcCavityReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");
            need(closed.reason!=bnd::TetPlcCavityReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");
            return closed.published;
        };
        if(completeDiscreteSeed){
            need(fillDiscreteSeed(pendingRecovery,completed),
                 NativeCadMeshReason::MeshingFailed,"native-plc-discrete-fill-refused");
            pendingRecovery=mesh::Mesh{};
        }
        std::size_t lastVertices=SIZE_MAX,lastCells=SIZE_MAX;
        stage="gift-wrapping";trace.next(stage);const unsigned frontSlices=24;for(unsigned slice=0;slice<frontSlices&&!completeDiscreteSeed&&!pendingRecovery.cellCount();++slice){tick();const auto testsBefore=front.metrics.apexTests;front.run(static_cast<unsigned>(std::min<std::uint64_t>(frontCells,UINT_MAX)),5,stop.get_token(),1000000);
            if(front.front.size()>512)continue;
            // Construction only appends nodes/cells. Revalidating the same
            // partial PLC cannot produce a different deterministic recovery.
            if(lastVertices==front.vertices.size()&&lastCells==front.cells.size()){if(testsBefore==front.metrics.apexTests)break;continue;}lastVertices=front.vertices.size();lastCells=front.cells.size();
            auto partial=partialFrom(front);auto result=recover(partial,completed);
            if(result.published){out.cavitiesClosed+=result.cavitiesClosed;break;}
            if(pendingRecovery.cellCount())break;
            need(result.reason!=bnd::TetPlcCavityReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");need(result.reason!=bnd::TetPlcCavityReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");
            need(result.reason!=bnd::TetPlcCavityReason::BudgetExceeded,NativeCadMeshReason::BudgetExceeded,result.diagnostic);
        }
        if(completed.cellCount()<=nf&&!pendingRecovery.cellCount()&&out.randomRetries<options.maximumRecoveryRandomRetries){auto partial=partialFrom(front);const auto result=recover(partial,completed,true);need(result.reason!=bnd::TetPlcCavityReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");need(result.reason!=bnd::TetPlcCavityReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");}
        // Propose local clearance points only after the established front is
        // exhausted. Exact PLC/native admission still guards every new cell.
        if(completed.cellCount()<=nf&&!pendingRecovery.cellCount()){front.interiorAdmission=native;front.centroidAdmission=nativeCentroid;const unsigned seedPasses=4;for(unsigned pass=0;pass<seedPasses&&completed.cellCount()<=nf;++pass){tick();const auto room=o.maximumNodes>front.vertices.size()?o.maximumNodes-front.vertices.size():0;if(!room)break;const auto seeded=front.seedUnfilledFront(static_cast<unsigned>(std::min<std::uint64_t>(8,room)));out.localFrontSeeds+=seeded;if(!seeded)break;front.run(static_cast<unsigned>(std::min<std::uint64_t>(frontCells,UINT_MAX)),5,stop.get_token(),front.metrics.apexTests+100000);auto partial=partialFrom(front);const auto result=recover(partial,completed);need(result.reason!=bnd::TetPlcCavityReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");need(result.reason!=bnd::TetPlcCavityReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");}}
        stage="material-filter";trace.next(stage);const bool initiallyComplete=completed.cellCount()>nf;
        const auto availableBytes=options.maximumPrimitiveWorkingBytes-retained;
        const bool useTrimmedRayWitnesses=options.useTrimmedRayWitnesses&&liveTrial;
        const bool useLocalDelaunay=options.useLocalDelaunayRefinement&&
            (liveTrial||regularSurface)&&usesConstraintSeed&&!options.allowCentroidsWithinChordTolerance&&
            o.minimumElementQuality==0.;
        const auto seedCentroids=[&](const mesh::Mesh& seed){
            const auto count=std::uint64_t(seed.countCells(mesh::CellType::Tet4));
            const auto maximumCenters=progressiveResources?
                std::min({o.maximumCells,std::uint64_t(UINT32_MAX),
                    std::uint64_t(SIZE_MAX/sizeof(Vec3)),availableBytes/128}):
                std::min<std::uint64_t>(1'000'000,availableBytes/128);
            need(count<=maximumCenters,
                 NativeCadMeshReason::BudgetExceeded,"native-plc-constraint-classifier-storage");
            std::vector<Vec3> centers;
            centers.reserve(count);
            for(Index i=0;i<seed.cellCount();++i){
                if(seed.cell(i).type!=mesh::CellType::Tet4)continue;
                tick();
                Vec3 center{};
                const auto* p=seed.cellNodes(i);
                for(unsigned j=0;j<4;++j)center=center+seed.node(p[j])*.25;
                centers.push_back(center);
            }
            return centers;
        };
        const auto classifyCenters=[&](std::span<const Vec3> centers,
                                       std::vector<unsigned char>& admitted,
                                       const OcctNativeStrictAdmissions* prior=nullptr){
            const auto count=std::uint64_t(centers.size());
            admitQueries(count);
            OcctNativeValidationOptions admission;
            admission.maximumThreads=options.maximumNativeValidationThreads;
            admission.maximumWork=remaining();admission.maximumPoints=count;
            need(!progressiveResources||(count<=availableBytes/128&&availableBytes-count*128>=4),
                 NativeCadMeshReason::BudgetExceeded,"native-plc-constraint-classifier-memory");
            admission.maximumWorkingBytes=availableBytes-count*128;
            admission.maximumSnapshotBytes=std::min<std::uint64_t>(admission.maximumSnapshotBytes,admission.maximumWorkingBytes/4);
            admission.cancellation=o.cancellation;
            admission.maximumOutsideDistance=options.allowCentroidsWithinChordTolerance?out.chord:0;
            admission.useTrimmedRayWitnesses=useTrimmedRayWitnesses;
            OcctStageTrace batchTrace("original_cad_centroid_batch",tag);
            const auto batch=classifyOcctNativeInteriorPoints(model,tag,centers,admission,admitted,
                reuseStrictAdmissions?&strictAdmissions:nullptr,prior);
            batchTrace.finish(occtNativeValidationReasonName(batch.reason));
            out.nativeQueries+=batch.nativeQueries+batch.trimmedRayAdmissions;
            out.nativeRayAdmissions+=batch.trimmedRayAdmissions;
            accountCompleted(batch.workUsed);
            need(batch.reason!=OcctNativeValidationReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");
            need(batch.reason!=OcctNativeValidationReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");
            need(batch.reason!=OcctNativeValidationReason::BudgetExceeded,NativeCadMeshReason::BudgetExceeded,"native-plc-constraint-classifier-budget");
            need(batch.completed&&batch.sourceModelUnmodified&&admitted.size()==centers.size()&&
                 batch.nativeQueries+batch.strictAdmissionReuses+batch.trimmedRayAdmissions==centers.size(),
                 NativeCadMeshReason::GeometryError,"native-plc-constraint-classifier-refused");
            for(std::size_t i=0;i<centers.size();++i){
                tick();
                if(!progressiveResources||centroidCache.size()<131072U)
                    centroidCache.emplace(pointKey(centers[i]),admitted[i]!=0);
            }
        };
        std::optional<OcctTrimmedRayWitness> localWitness;
        const auto localInsert=[&](mesh::Mesh& mesh,std::span<const Vec3> sites,
                                   std::span<const unsigned char> material={},unsigned* qualityEdits=nullptr,
                                   bool deferMaterial=false){
            const auto used=strictAdmissions.retainedPayloadBytes();
            constexpr std::uint64_t witnessBytes=128*1024;
            if(used+witnessBytes>=availableBytes) return std::size_t(0);
            const auto room=spaceSurface&&out.nativeQueries>=maximumSpacedCandidateQueries?0:
                spaceSurface?maximumSpacedCandidateQueries-out.nativeQueries:6'000;
            const auto queries=std::min<std::uint64_t>(sites.empty()?12'000:1'400,room/2);
            if(!queries&&!deferMaterial) return std::size_t(0);
            if(!deferMaterial&&useTrimmedRayWitnesses&&options.maximumNativeValidationThreads==1&&!localWitness){
                localWitness.emplace();
                localWitness->load(OcctModelAccess::solid(model,tag),tol);
            }
            const auto admitLocalCentroid=[&](Vec3 p){
                if(!localWitness) return nativeCentroid(p);
                tick();
                const auto key=pointKey(p);
                if(auto known=centroidCache.find(key);known!=centroidCache.end()) return known->second;
                if(!localWitness->admits(gp_Pnt(p.x,p.y,p.z))) return nativeCentroid(p);
                admitQueries(1);
                ++out.nativeQueries;
                ++out.nativeRayAdmissions;
                if(!progressiveResources||centroidCache.size()<131072U)centroidCache.emplace(key,true);
                return true;
            };
            const auto sizeOfMesh=std::uint64_t(mesh.nodeCount())*(sizeof(Vec3)+sizeof(mesh::Classification))+
                std::uint64_t(mesh.cellCount())*sizeof(mesh::Cell)+std::uint64_t(mesh.connectivityCount())*sizeof(Index);
            if(sizeOfMesh>=o.maximumOutputBytes) return std::size_t(0);
            // One insertion adds one vertex and at most 2k+2 cells for a
            // connected k-cell cavity. Reserve an output share for vertices.
            const auto outputGrowth=(o.maximumOutputBytes-sizeOfMesh)/(
                sizeof(Vec3)+sizeof(mesh::Classification)+sizeof(mesh::Cell)+4*sizeof(Index));
            return insertOcctDelaunaySites(mesh,sites,
                {std::min<std::uint64_t>(sites.empty()?42'000'000:5'000'000,remaining()/2),
                 deferMaterial?65'000:queries,
                 std::min(o.maximumNodes,std::uint64_t(mesh.nodeCount())+outputGrowth),
                 std::min(frontCells,std::uint64_t(mesh.countCells(mesh::CellType::Tet4))+outputGrowth),
                 availableBytes-used-witnessBytes,o.targetSize,
                 spatialSizing&&sites.empty()&& !deferMaterial?&*sizeField:nullptr,
                 options.improveNativeTetQuality},
                 [&](Vec3 p){return deferMaterial||native(p);},
                 [&](Vec3 p){return deferMaterial||admitLocalCentroid(p);},tick,material,qualityEdits);
        };
        bool refinedBeforeMaterial=false;
        if(usesConstraintSeed){
            std::vector<unsigned char> admitted;
            bool inwardRetry=false;
            bool localSupportAttempted=false;
            const auto repairLocal=[&](mesh::Mesh& mesh,std::span<const unsigned char> mask,bool afterSupport){
                LiveCavityReceipt patched;
                patched.initialRejected=patched.remainingRejected=std::ranges::count(mask,0);
                if(!liveTrial||!patched.initialRejected||
                   availableBytes<std::uint64_t(mesh.cellCount()+4096)*4096) return patched;
                OcctStageTrace local("live_material_repair",tag);
                patched=repairLiveCavities(mesh,mask,nativeCentroid,tick,
                    std::min<std::uint64_t>(afterSupport?5'000'000:2'500'000,remaining()/4),2000,
                    std::min(o.maximumNodes,std::uint64_t(mesh.nodeCount())+96),frontCells,curvatureSurface);
                local.finish(!patched.remainingRejected?"complete":patched.patched?"reduced":"unchanged");
                return patched;
            };
            const auto recordLocalRepair=[&](const LiveCavityReceipt& patched){
                out.liveMaterialInitialRejects=patched.initialRejected;
                out.liveMaterialRemainingRejects=patched.remainingRejected;
                out.liveMaterialPatches=patched.patched;
            };
            bool materialComplete=false;
            // Try the complete discrete placement and optimization before
            // native material repair. This private mesh cannot be published
            // until every retained node and centroid has a CAD witness.
            if(useLocalDelaunay&&curvatureSurface&&initiallyComplete&&liveTrial){
                OcctStageTrace placement("local_delaunay_before_material",tag);
                mesh::Mesh candidate=completed;
                unsigned edits=0;
                const auto inserted=localInsert(candidate,{}, {}, &edits,true);
                const auto centers=seedCentroids(candidate);
                const auto candidateQueries=std::uint64_t(centers.size()+candidate.nodeCount()-ns);
                if((inserted||edits)&&candidateQueries<=12'000&&
                   out.nativeQueries+candidateQueries+4000<=maximumSpacedCandidateQueries){
                    std::vector<unsigned char> mask;
                    classifyCenters(centers,mask);
                    bool nodesAdmitted=true;
                    for(Index i=ns;i<candidate.nodeCount()&&nodesAdmitted;++i)
                        nodesAdmitted=native(candidate.node(i));
                    if(nodesAdmitted){
                        const auto patched=repairLocal(candidate,mask,true);
                        if(!patched.remainingRejected){
                            completed.swap(candidate);
                            out.localDelaunayRefinementNodes=inserted;
                            out.localDelaunayQualityEdits=edits;
                            recordLocalRepair(patched);
                            refinedBeforeMaterial=materialComplete=true;
                        }
                    }
                }
                placement.finish(refinedBeforeMaterial?"accepted":"fallback");
            }
            // A shell-only fill with many rejected cells may have no nearby
            // interior apex to move. Look near the boundary for repair
            // proposals before committing to an audit of the entire fill.
            const bool canRetryInward=liveTrial&&curvatureSurface&&completed.nodeCount()==ns&&
                o.maximumNodes>std::uint64_t(ns);
            const auto tryInward=[&](std::span<const unsigned char> rejectedMask){
                if(std::ranges::count(rejectedMask,0)<=16)return false;
                const auto admitInterior=[&](Vec3 p){
                    tick();admitQueries(1);++out.nativeQueries;
                    classifier.Perform(gp_Pnt(p.x,p.y,p.z),tol);
                    return classifier.State()==TopAbs_IN;
                };
                const auto points=occtInwardSites(boundary,completed,rejectedMask,out.chord,
                    std::min<std::uint64_t>(96,o.maximumNodes-std::uint64_t(ns)),
                    admitInterior,nativeCentroid,tick);
                if(points.empty())return false;
                if(useLocalDelaunay&&!localSupportAttempted){
                    localSupportAttempted=true;
                    mesh::Mesh candidate=completed;
                    OcctStageTrace insertion("local_delaunay_support",tag);
                    const auto inserted=localInsert(candidate,points,rejectedMask);
                    insertion.finish(inserted?"inserted":"refused");
                    if(inserted){
                        std::vector<unsigned char> candidateMask;
                        classifyCenters(seedCentroids(candidate),candidateMask,&strictAdmissions);
                        // Keep the complete cavity result only if every cell
                        // has a native witness. Otherwise retain the original
                        // fill and try the established global retry.
                        const auto patched=repairLocal(candidate,candidateMask,true);
                        if(!patched.remainingRejected){
                            completed.swap(candidate);
                            out.localDelaunaySupportNodes=inserted;
                            recordLocalRepair(patched);
                            materialComplete=true;
                            return true;
                        }
                    }
                    out.localDelaunayUsedGlobalFallback=true;
                }
                mesh::Mesh sites=boundary;
                for(const auto point:points)sites.addNode(point,{3,tag});
                mesh::Mesh seed,filled;
                buildConstraintSeed(sites,seed);
                if(!fillDiscreteSeed(seed,filled))return false;
                completed.swap(filled);
                return true;
            };
            if(canRetryInward&&!materialComplete){
                OcctStageTrace proximity("curvature_shell_proximity",tag);
                const auto centers=seedCentroids(completed);
                const auto near=occtCentroidsNearShell(boundary,centers,out.chord,
                    std::min<std::uint64_t>(2'000'000,remaining()/8),tick);
                // A poorly localized shell would make this nearly a full
                // scan. Use the existing complete classification in that case.
                if(!near.empty()&&near.size()<=centers.size()/2){
                    std::vector<Vec3> candidates;
                    candidates.reserve(near.size());
                    for(const auto i:near)candidates.push_back(centers[i]);
                    std::vector<unsigned char> sampled;
                    classifyCenters(candidates,sampled);
                    // This mask is for selecting sites only. Unsampled rows
                    // stay out of the CAD cache and the admission receipt.
                    std::vector<unsigned char> proposals(centers.size(),1);
                    for(std::size_t i=0;i<near.size();++i)proposals[near[i]]=sampled[i];
                    inwardRetry=tryInward(proposals);
                }
                proximity.finish(inwardRetry?"filled":"fallback");
            }
            if(!materialComplete){
                centroidCache.clear();
                classifyCenters(seedCentroids(completeDiscreteSeed?completed:pendingRecovery),
                                admitted,&strictAdmissions);
                if(canRetryInward&&!inwardRetry){
                    OcctStageTrace retry("curvature_inward_sites",tag);
                    inwardRetry=tryInward(admitted);
                    if(inwardRetry&&!materialComplete){
                        centroidCache.clear();
                        classifyCenters(seedCentroids(completed),admitted,&strictAdmissions);
                    }
                    retry.finish(inwardRetry?"filled":"refused");
                }
            }
            recovery.maximumWorkingBytes-=strictAdmissions.retainedPayloadBytes();
            if(!materialComplete&&liveTrial)
                recordLocalRepair(repairLocal(completed,admitted,inwardRetry));
        }
        // The explicit approximation policy applies only to cell centroids.
        // Interior vertices and all recovery proposals retain strict admission.
        mesh::Mesh incomplete;if(!initiallyComplete){if(pendingRecovery.cellCount())incomplete.swap(pendingRecovery);else incomplete=partialFrom(front);}mesh::Mesh inside;const auto material=retainOcctPlcMaterial(initiallyComplete?completed:incomplete,ns,native,nativeCentroid,tick,inside);
        out.removedInteriorNodes=material.removedInteriorNodes;out.removedIncidentCells=material.removedIncidentCells;out.removedCentroidCells=material.removedCentroidCells;
        if(!initiallyComplete||material.removedInteriorNodes||material.removedIncidentCells||material.removedCentroidCells){
            need(!curvatureSurface,NativeCadMeshReason::MeshingFailed,"native-curvature-material-refused");
            for(unsigned pass=0;pass<4;++pass){
                mesh::Mesh candidate;const auto changed=closeOcctBoundaryDiagonals(model,boundary,inside,128,candidate,[&](std::uint64_t n){tick(n);});if(!changed.flips)break;
                bo.maximumWork=remaining();const auto facets=auditOcctBoundaryFacets(model,tag,candidate,bo);accountCompleted(facets.work);certificateOptions.maximumWork=remaining();const auto cert=certifyOcctBoundarySurface(candidate,certificateOptions);accountCompleted(cert.workUsed);
                
                if(!(facets.completed&&facets.sourceModelUnmodified&&cert.certified))break;boundary.swap(candidate);
            }
            recovery.interiorAdmission=native;recovery.centroidAdmission=nativeCentroid;evolvingPartial=&inside;mesh::Mesh repaired;
            const bool localCadHoles=(liveTrial||(regularSurface&&
                std::string_view(out.constraintRecovery.diagnostic)=="complete_seed"))&&options.maximumRecoveryPatchStates&&
                !material.removedInteriorNodes&&material.removedCentroidCells>0&&material.removedCentroidCells<=16;
            auto result=recover(inside,repaired,false,localCadHoles);
            if(!spaceSurface&&!result.published&&result.reason==bnd::TetPlcCavityReason::Unavailable&&out.randomRetries<options.maximumRecoveryRandomRetries)result=recover(inside,repaired,true);
            bool localRefinementBound=false;
            for(unsigned refinement=0;!spaceSurface&&!result.published&&refinement<4&&(!result.unresolvedSourceTriangles.empty()||!result.blockedSourceTriangles.empty());++refinement){
                need(result.reason!=bnd::TetPlcCavityReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");need(result.reason!=bnd::TetPlcCavityReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");
                {
                    mesh::Mesh candidate;const auto changed=closeOcctBoundaryDiagonals(model,boundary,inside,16,candidate,[&](std::uint64_t n){tick(n);});
                    if(changed.flips){bo.maximumWork=remaining();const auto facets=auditOcctBoundaryFacets(model,tag,candidate,bo);accountCompleted(facets.work);certificateOptions.maximumWork=remaining();const auto cert=certifyOcctBoundarySurface(candidate,certificateOptions);accountCompleted(cert.workUsed);
                        
                        if(facets.completed&&facets.sourceModelUnmodified&&cert.certified){boundary.swap(candidate);result=recover(inside,repaired);continue;}
                    }
                }
                using Edge=std::array<Index,2>;std::map<Edge,std::set<int>> owners;const auto edgeKey=[](Index a,Index b){if(b<a)std::swap(a,b);return Edge{a,b};};
                for(Index i=0;i<boundary.cellCount();++i){const auto* p=boundary.cellNodes(i);for(unsigned j=0;j<3;++j){tick();owners[edgeKey(p[j],p[(j+1)%3])].insert(boundary.cell(i).onEntity.tag);}}
                std::set<int> edges;for(auto i:result.unresolvedSourceTriangles){need(i>=0&&i<boundary.cellCount(),NativeCadMeshReason::InternalError,"native-plc-diagnostic-face-index");const auto* p=boundary.cellNodes(i);for(unsigned j=0;j<3;++j){const auto a=p[j],b=p[(j+1)%3];const auto edge=actualOcctWire(model,boundary.nodeClass(a),boundary.nodeClass(b),owners.at(edgeKey(a,b)),[&]{tick();});if(edge>0&&model.edge(edge)->kind()!=CurveKind::Line)edges.insert(edge);}}
                if(edges.size()>16){auto end=edges.begin();std::advance(end,16);edges.erase(end,edges.end());}const auto previousNodes=boundary.nodeCount();bool changed=false;
                // A refused first pass can terminate this loop. Admit checked
                // closed cycles here too, so their recovery does not depend
                // on an unrelated facet subdivision enabling a later pass.
                for(auto edge:edges){OcctWireRefinementOptions wire;wire.edgeTag=edge;wire.allowNonlinearWireCandidates=true;wire.allowClosedWireCandidates=true;wire.limits=certificateOptions;wire.limits.maximumWork=remaining();wire.limits.maximumBoundaryNodes=bo.maximumNodes;wire.limits.maximumBoundaryTriangles=bo.maximumTriangles;wire.limits.maximumNodes=o.maximumNodes;wire.limits.maximumCells=o.maximumCells;wire.limits.maximumOutputBytes=o.maximumOutputBytes;wire.limits.maximumWorkingBytes=options.maximumPrimitiveWorkingBytes;
                    const std::array input{OcctWireBoundary{tag,&boundary}};const auto refined=refineOcctWireRefinement(model,input,wire);accountCompleted(refined.workUsed);need(refined.reason!=OcctWireRefinementReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");need(refined.reason!=OcctWireRefinementReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");
                    if(refined.reason==OcctWireRefinementReason::BudgetExceeded){
                        // A refused local subdivision leaves the boundary
                        // unchanged. Other bounded strategies can still fit.
                        localRefinementBound=true;need(out.workUsed<maximumWork,NativeCadMeshReason::BudgetExceeded,"native-plc-adaptive-wire-budget");continue;
                    }
                    changed|=refined.published;
                }
                if(!changed){std::vector<OcctFacetSelection> selected;const auto& candidates=result.blockedSourceTriangles.empty()?result.unresolvedSourceTriangles:result.blockedSourceTriangles;for(auto i:candidates){need(i>=0&&i<boundary.cellCount(),NativeCadMeshReason::InternalError,"native-plc-diagnostic-face-index");selected.push_back({tag,i});if(selected.size()==8)break;}
                    if(!selected.empty()){OcctFacetRefinementOptions facet;facet.selections=selected;facet.maximumProjectionDistance=out.chord;facet.maximumAddedNodes=8;facet.limits=certificateOptions;facet.limits.maximumWork=remaining();const std::array input{OcctFacetBoundary{tag,&boundary}};const auto refined=refineOcctFacets(model,input,facet);out.sourceFacetNodes+=refined.addedNodes;accountCompleted(refined.workUsed);need(refined.reason!=OcctFacetRefinementReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");need(refined.reason!=OcctFacetRefinementReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");if(refined.reason==OcctFacetRefinementReason::BudgetExceeded){localRefinementBound=true;need(out.workUsed<maximumWork,NativeCadMeshReason::BudgetExceeded,"native-plc-adaptive-facet-budget");}changed=refined.published;}}
                if(!changed)break;
                // Reclassify the retained partial mesh against the refined PLC.
                // Any crossed cell is discarded before another recovery attempt.
                bnd::detail::PlcGiftWrapping checked;checked.chargeWork=tick;checked.cancellation=stop.get_token();for(Index i=0;i<boundary.nodeCount();++i)checked.vertices.push_back({boundary.node(i),checked.vertices.size()+1});
                for(Index i=previousNodes;i<inside.nodeCount();++i)checked.vertices.push_back({inside.node(i),checked.vertices.size()+1});for(Index i=0;i<boundary.cellCount();++i){const auto* p=boundary.cellNodes(i);checked.boundary.push_back({unsigned(p[0]),unsigned(p[1]),unsigned(p[2])});}need(checked.initialize(),NativeCadMeshReason::GeometryError,"native-plc-refined-front-invalid");
                const auto offset=boundary.nodeCount()-previousNodes;for(Index i=0;i<inside.cellCount();++i){bnd::TetFacetCell c;const auto* p=inside.cellNodes(i);for(unsigned j=0;j<4;++j)c.vertices[j]=unsigned(p[j]>=previousNodes?p[j]+offset:p[j]);bnd::TetFacetCell valid;const auto face=bnd::tet_facet_detail::outwardFace(c,0);if(checked.valid(bnd::detail::PlcGiftWrapping::Front{face},c.vertices[0],valid))checked.insert(c);}
                checked.interiorAdmission=native;checked.centroidAdmission=nativeCentroid;checked.run(static_cast<unsigned>(std::min<std::uint64_t>(frontCells,UINT_MAX)),10,stop.get_token(),100000);
                mesh::Mesh next;for(Index i=0;i<Index(checked.vertices.size());++i)next.addNode(checked.vertices[i].position,i<boundary.nodeCount()?boundary.nodeClass(i):mesh::Classification{3,tag});for(auto c:checked.cells){std::array<Index,4> ids;for(unsigned j=0;j<4;++j)ids[j]=c.vertices[j];next.addCell(mesh::CellType::Tet4,ids,{3,tag});}inside.swap(next);
                bo.maximumWork=remaining();const auto audited=auditOcctBoundaryFacets(model,tag,boundary,bo);out.maximumFacetError=std::max(out.maximumFacetError,audited.maximumFacetSampleError);accountCompleted(audited.work);need(audited.completed&&audited.sourceModelUnmodified,NativeCadMeshReason::GeometryError,"native-plc-adaptive-facet-error");
                result=recover(inside,repaired);if(!result.published&&result.reason==bnd::TetPlcCavityReason::Unavailable&&out.randomRetries<options.maximumRecoveryRandomRetries)result=recover(inside,repaired,true);
            }
            // Only a final refusal reaches this path: previously accepted
            // meshes pay no region-search cost. Keep the refined CAD surface
            // fixed and retry the complete affected cell region, not its points.
            if(!spaceSurface&&!result.published&&result.reason==bnd::TetPlcCavityReason::Unavailable&&options.maximumRecoveryPatchStates){
                stage="region-search";trace.next(stage);recovery.patchSearch.maximumStates=options.maximumRecoveryPatchStates;recovery.patchSearch.maximumVertices=regionVertices;recovery.patchSearch.maximumCandidateTets=32768;recovery.patchSearch.maximumShells=8;recovery.useVisibilityExpansion=false;recovery.maximumMergedComponents=4;
                result=recover(inside,repaired);recovery.patchSearch.maximumStates=0;recovery.useVisibilityExpansion=true;
            }
            need(result.reason!=bnd::TetPlcCavityReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");need(result.reason!=bnd::TetPlcCavityReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");need(result.reason!=bnd::TetPlcCavityReason::BudgetExceeded,NativeCadMeshReason::BudgetExceeded,result.diagnostic);
            need(result.published,localRefinementBound?NativeCadMeshReason::BudgetExceeded:NativeCadMeshReason::GeometryError,localRefinementBound?"native-plc-refinement-budget":"native-plc-cavity-material-refusal");out.cavitiesClosed+=result.cavitiesClosed;completed.swap(repaired);out.surfaceTriangles=boundary.cellCount();
        }
        bnd::detail::observeTetMesh("volume_before_sizing",tag,completed);
        // Uniform requests use only the initial point construction above.
        // Explicit local controls can additionally require adaptive sizing.
        if(advancingInteriorFront&&(!options.faceSizes.empty()||!options.sizeBoxes.empty())){
            OcctStageTrace construction("interior_front",tag);
            const auto started=Clock::now();
            const auto used=strictAdmissions.retainedPayloadBytes();
            const auto bytes=std::uint64_t(completed.nodeCount())*(sizeof(Vec3)+sizeof(mesh::Classification))+
                std::uint64_t(completed.cellCount())*sizeof(mesh::Cell)+
                std::uint64_t(completed.connectivityCount())*sizeof(Index);
            const auto reserve=1024*(std::uint64_t(completed.nodeCount())+
                std::uint64_t(completed.countCells(mesh::CellType::Tet4)));
            const auto room=remaining();
            const auto budget=room>reserve?(progressiveResources?(room-reserve)/2:
                std::min<std::uint64_t>(42'000'000,(room-reserve)/2)):0;
            const auto queries=out.nativeQueries<maximumSpacedCandidateQueries?
                (progressiveResources?(maximumSpacedCandidateQueries-out.nativeQueries)/2:
                 std::min<std::uint64_t>(12'000,(maximumSpacedCandidateQueries-out.nativeQueries)/2)):0;
            out.interiorFrontDiagnostic="budget_limited";
            if(budget&&queries&&used<availableBytes&&bytes<o.maximumOutputBytes){
                const auto growth=(o.maximumOutputBytes-bytes)/(
                    sizeof(Vec3)+sizeof(mesh::Classification)+sizeof(mesh::Cell)+4*sizeof(Index));
                OcctInteriorFrontLimits limits;
                limits.work=budget;limits.queries=queries;
                limits.nodes=std::min(o.maximumNodes,std::uint64_t(completed.nodeCount())+growth);
                limits.tets=std::min(frontCells,std::uint64_t(completed.countCells(mesh::CellType::Tet4))+growth);
                limits.workingBytes=availableBytes-used;
                const auto built=advanceOcctInteriorFront(completed,*sizeField,limits,
                    [&](Vec3 p){
                        tick();admitQueries(1);++out.nativeQueries;
                        classifier.Perform(gp_Pnt(p.x,p.y,p.z),tol);
                        return classifier.State()==TopAbs_IN;
                    },nativeCentroid,tick);
                out.interiorFrontDiagnostic=built.budgetLimited?"budget_limited":
                    !built.initialized?"scaffold_refused":built.targetsSatisfied?"targets_satisfied":"targets_unmet";
                out.interiorFrontNodes=built.sites;
                out.interiorFrontCandidates=built.candidates;
                out.interiorFrontWork=built.workUsed;
                out.interiorFrontQueries=built.nativeQueries;
                out.interiorFrontRemainingOversized=built.remainingOversized;
                out.interiorFrontActivations=built.frontActivations;
                out.interiorVoidActivations=built.voidActivations;
                out.interiorFrontAdvances=built.frontAdvances;
            }
            out.interiorFrontSeconds=std::chrono::duration<double>(Clock::now()-started).count();
            construction.finish(out.interiorFrontDiagnostic);
        }
        bnd::detail::observeTetMesh("volume_after_sizing",tag,completed);
        if(improveQuality&&out.nativeQueries<maximumSpacedCandidateQueries){
            // Reserve work and queries for the strict validation of the
            // candidate. Swaps retain the existing boundary and its nodes.
            const auto budget=std::min<std::uint64_t>(200'000,remaining()/8);
            auto queries=std::min<std::uint64_t>(96,(maximumSpacedCandidateQueries-out.nativeQueries)/2);
            if(budget&&queries){
                trace.next("tet_quality");
                improveOcctTetPiece(completed,{budget,o.maximumNodes,
                    std::min(frontCells,std::uint64_t(completed.countCells(mesh::CellType::Tet4))),0.,queries},
                    native,tick);
            }
        }
        if(!advancingInteriorFront&&(spatialSizing||(useLocalDelaunay&&!refinedBeforeMaterial))){
            OcctStageTrace refinement("local_delaunay_refinement",tag);
            out.localDelaunayRefinementNodes=localInsert(completed,{}, {},
                !spatialSizing||improveQuality?&out.localDelaunayQualityEdits:nullptr);
            refinement.finish(out.localDelaunayRefinementNodes||out.localDelaunayQualityEdits?"refined":"unchanged");
        }
        if(spatialSizing&&options.improveNativeTetQuality){
            OcctStageTrace optimization("volume_quality",tag);
            const auto queries=out.nativeQueries<maximumSpacedCandidateQueries?
                std::min<std::uint64_t>(3000,(maximumSpacedCandidateQueries-out.nativeQueries)/4):0;
            // Leave a conservative work share for the mandatory final audits.
            // Include native callback charges in this optional stage's share.
            const auto reserve=1024*(std::uint64_t(completed.nodeCount())+
                std::uint64_t(completed.countCells(mesh::CellType::Tet4)));
            const auto room=remaining();
            const auto budget=room>reserve?std::min<std::uint64_t>(2'000'000,(room-reserve)/8):0;
            const auto beganQuality=out.workUsed;
            const auto qualityTick=[&](std::uint64_t n){
                const auto used=out.workUsed-beganQuality;
                if(used>budget||n>budget-used)throw OcctDelaunayQualityLimit{};
                tick(n);
            };
            const auto bytes=std::uint64_t(completed.nodeCount())*(sizeof(Vec3)+sizeof(mesh::Classification))+
                std::uint64_t(completed.cellCount())*sizeof(mesh::Cell)+
                std::uint64_t(completed.connectivityCount())*sizeof(Index);
            const auto outputGrowth=bytes<o.maximumOutputBytes?
                (o.maximumOutputBytes-bytes)/(sizeof(mesh::Cell)+4*sizeof(Index)):0;
            const auto retainedQuality=strictAdmissions.retainedPayloadBytes();
            out.localDelaunayQualityEdits+=improveOcctDelaunayMeshQuality(completed,
                budget,queries,std::min(frontCells,
                    std::uint64_t(completed.countCells(mesh::CellType::Tet4))+outputGrowth),
                retainedQuality<availableBytes?availableBytes-retainedQuality:0,
                native,nativeCentroid,qualityTick);
            optimization.finish(out.localDelaunayQualityEdits?"improved":"unchanged");
        }
        bnd::detail::observeTetMesh("volume_before_cad_validation",tag,completed);
        trace.next("cad_validation");
        for(Index i=0;i<completed.nodeCount();++i)if(completed.nodeClass(i).dim==3)need(native(completed.node(i)),NativeCadMeshReason::GeometryError,"native-plc-final-interior-node-outside");
        for(Index i=0;i<completed.cellCount();++i)if(completed.cell(i).type==mesh::CellType::Tet4){const auto* p=completed.cellNodes(i);Vec3 center{};for(unsigned j=0;j<4;++j)center=center+completed.node(p[j])*.25;need(nativeCentroid(center),NativeCadMeshReason::GeometryError,"native-plc-final-centroid-outside");}
        const auto quality=algo::analyzeMixedMeshQuality(completed);tick(completed.cellCount());need(!quality.allocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");need(quality.valid,NativeCadMeshReason::MeshingFailed,"native-plc-final-quality-invalid");
        need(!spaceSurface||quality.minimumQuality>=o.minimumElementQuality,NativeCadMeshReason::QualityRequirementNotMet,"spaced-surface-quality-requirement");
        need(std::abs(quality.totalVolume-model.volume(tag))<=model.volume(tag)*options.maximumRelativeVolumeError,NativeCadMeshReason::GeometryError,"native-plc-final-volume-error");
        std::vector<Vec3> centroids;for(Index i=0;i<completed.cellCount();++i)if(completed.cell(i).type==mesh::CellType::Tet4){tick();const auto* p=completed.cellNodes(i);Vec3 center{};for(unsigned j=0;j<4;++j)center=center+completed.node(p[j])*.25;centroids.push_back(center);}
        need(!progressiveResources||availableBytes>=4,NativeCadMeshReason::BudgetExceeded,
             "native-plc-final-native-memory");
        OcctNativeValidationOptions validation;validation.maximumThreads=options.maximumNativeValidationThreads;validation.maximumWork=remaining();validation.maximumPoints=o.maximumCells;validation.maximumWorkingBytes=progressiveResources?availableBytes:options.maximumPrimitiveWorkingBytes;validation.cancellation=o.cancellation;validation.maximumOutsideDistance=options.allowCentroidsWithinChordTolerance?out.chord:0;
        if(progressiveResources)validation.maximumSnapshotBytes=std::min(validation.maximumSnapshotBytes,validation.maximumWorkingBytes/4);
        validation.useTrimmedRayWitnesses=useTrimmedRayWitnesses;
        out.nativeValidation=validateOcctNativeInteriorPoints(model,tag,centroids,validation,
            strictAdmissions.retainedPayloadBytes()?&strictAdmissions:nullptr);accountCompleted(out.nativeValidation.workUsed);
        need(out.nativeValidation.reason!=OcctNativeValidationReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");need(out.nativeValidation.reason!=OcctNativeValidationReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");need(out.nativeValidation.reason!=OcctNativeValidationReason::BudgetExceeded,NativeCadMeshReason::BudgetExceeded,"native-plc-final-native-budget");need(out.nativeValidation.completed&&out.nativeValidation.sourceModelUnmodified,NativeCadMeshReason::GeometryError,"native-plc-final-native-validation");
        need(out.nativeValidation.nativeQueries+out.nativeValidation.strictAdmissionReuses+
             out.nativeValidation.trimmedRayAdmissions==centroids.size(),
             NativeCadMeshReason::InternalError,"native-plc-incomplete-centroid-coverage");
        if(options.allowCentroidsWithinChordTolerance){
            std::vector<Vec3> samples;const auto maximumSamples=(options.maximumPrimitiveWorkingBytes-retained)/sizeof(Vec3);
            for(Index i=0;i<completed.cellCount();++i)if(completed.cell(i).type==mesh::CellType::Tet4){const auto* p=completed.cellNodes(i);Vec3 center{};for(unsigned j=0;j<4;++j)center=center+completed.node(p[j])*.25;
                need(samples.size()<=maximumSamples&&4<=maximumSamples-samples.size(),NativeCadMeshReason::BudgetExceeded,"native-plc-chord-sample-storage");
                for(unsigned j=0;j<4;++j){tick();samples.push_back(center*.5+completed.node(p[j])*.5);}}
            validation.maximumPoints=std::uint64_t(completed.cellCount())*4;validation.maximumWork=remaining();
            out.chordSampleValidation=validateOcctNativeInteriorPoints(model,tag,samples,validation);accountCompleted(out.chordSampleValidation.workUsed);
            need(out.chordSampleValidation.reason!=OcctNativeValidationReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");
            need(out.chordSampleValidation.reason!=OcctNativeValidationReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");
            need(out.chordSampleValidation.reason!=OcctNativeValidationReason::BudgetExceeded,NativeCadMeshReason::BudgetExceeded,"native-plc-chord-sample-budget");
            need(out.chordSampleValidation.completed&&out.chordSampleValidation.sourceModelUnmodified,NativeCadMeshReason::GeometryError,"native-plc-chord-sample-validation");
        }
        // The curvature builder checked each original trimmed-face centroid
        // and three edge midpoints. Reuse its centroid check at assembly
        // only if all exterior coordinates, tags and directed facets match.
        if(curvatureSurface&&surface.sourceModelUnmodified&&completed.nodeCount()>=ns){
            bool identical=true;
            for(Index i=0;i<ns&&identical;++i){
                tick();
                identical=pointKey(completed.node(i))==pointKey(boundary.node(i))&&
                    completed.nodeClass(i)==boundary.nodeClass(i);
            }
            Index facet=0;
            for(Index i=0;i<completed.cellCount()&&identical;++i){
                const auto& cell=completed.cell(i);
                if(cell.onEntity.dim!=2)continue;
                tick();
                const auto* ids=completed.cellNodes(i);
                identical=cell.type==mesh::CellType::Tri3&&facet<nf&&cell.onEntity==boundary.cell(facet).onEntity&&
                    std::equal(ids,ids+3,boundary.cellNodes(facet));
                ++facet;
            }
            out.originalFacetCentersCertified=identical&&facet==nf;
        }
        need(std::uint64_t(completed.nodeCount())*(sizeof(Vec3)+sizeof(mesh::Classification))+std::uint64_t(completed.cellCount())*sizeof(mesh::Cell)+std::uint64_t(completed.connectivityCount())*sizeof(Index)<=o.maximumOutputBytes,NativeCadMeshReason::BudgetExceeded,"native-plc-output-bytes");
        need(std::uint64_t(completed.nodeCount())<=o.maximumNodes&&std::uint64_t(completed.cellCount())<=o.maximumCells,NativeCadMeshReason::BudgetExceeded,"native-plc-final-count-budget");tick(0);out.meshVolume=quality.totalVolume;destination.swap(completed);out.published=true;out.reason=NativeCadMeshReason::None;out.diagnostic="none";
    }catch(const Stop& s){out.reason=s.reason;out.diagnostic=s.diagnostic;}catch(const Standard_OutOfMemory&){out.reason=NativeCadMeshReason::AllocationFailed;out.diagnostic="allocation-failed";}catch(const std::bad_alloc&){out.reason=NativeCadMeshReason::AllocationFailed;out.diagnostic="allocation-failed";}catch(const Standard_Failure&){out.reason=NativeCadMeshReason::GeometryError;out.diagnostic="native-plc-query-failed";}catch(...){out.reason=NativeCadMeshReason::InternalError;out.diagnostic="native-plc-internal-error";}
    trace.finish(out.diagnostic);
    out.volumeSeconds=std::max(0.,std::chrono::duration<double>(Clock::now()-began).count()-out.surfaceSeconds);return out;
}
}
