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
#include "paacs/algo/BoundaryConformingMesher.hpp"
#include "paacs/bnd/TetPlcCavityRecovery.hpp"
#include "../bnd/TetPlcRetainedRecovery.hpp"
#include "../bnd/TetPlcGiftWrapping.hpp"
#include "../bnd/TetConstraintRecovery.hpp"
#include <Standard_OutOfMemory.hxx>
#include <Standard_Failure.hxx>
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
    const bool spaceSurface=surfaceMode!=OcctFrontSurface::Standard;
    const bool spatialSizing=surfaceMode==OcctFrontSurface::SpatialSizing;
    const bool improveQuality=surfaceMode==OcctFrontSurface::RegularQuality||
        (spatialSizing&&options.regularMeshQualityEnabled());
    const bool curvatureSurface=surfaceMode==OcctFrontSurface::Curvature;
    const bool regularSurface=surfaceMode==OcctFrontSurface::Regular||improveQuality||curvatureSurface||spatialSizing;
    const bool completeDiscreteSeed=liveTrial||curvatureSurface;
    const char* stage="surface";OcctStageTrace trace(stage,tag);OcctFrontOutcome out;const auto began=Clock::now();const auto& o=options.geometry;std::stop_source stop;
    const auto tick=[&](std::uint64_t n=1){if(o.cancellation&&o.cancellation->load()){stop.request_stop();throw Stop{NativeCadMeshReason::Cancelled,"cancelled"};}
        need(out.workUsed<=maximumWork&&n<=maximumWork-out.workUsed,NativeCadMeshReason::BudgetExceeded,"native-plc-work-budget");out.workUsed+=n;};
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
        mesh::Mesh boundary;
        OcctBoundaryOutcome surface;
        if(spaceSurface){
            bo.maximumNodes=std::min<std::uint64_t>(bo.maximumNodes,curvatureSurface||spatialSizing?6'000:1'200);
            bo.maximumTriangles=std::min<std::uint64_t>(bo.maximumTriangles,curvatureSurface||spatialSizing?10'000:1'800);
            OcctBoundaryAliasReceipt spacing;
            surface=buildSpacedOcctBoundary(model,tag,bo,o.targetSize,options.maximumPrimitiveWorkingBytes/4,
                boundary,spacing,regularSurface,improveQuality,
                curvatureSurface||spatialSizing?options.curvaturePointsPerTurn:0,sizeField?&*sizeField:nullptr);
            out.rebuiltSurfaceCharts=spacing.rebuiltCharts;
        }else if(options.recoveryMode==NativeCadRecoveryMode::Advanced){
            OcctBoundaryOutcome preparation;
            auto session=OcctBoundarySession::prepare(model,bo,preparation);
            if(session){
                OcctBoundaryAliasOptions repair;
                repair.maximumFlips=0;
                repair.maximumRebuiltCharts=32;
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
                out.rebuiltSurfaceCharts=repairs.rebuiltCharts;
            }else surface=std::move(preparation);
        }else surface=buildOcctBoundary(model,tag,bo,boundary);
        tick(surface.work);
        out.surfaceWork=surface.work;
        out.surfaceTriangles=surface.triangles;
        out.maximumNodeError=surface.maximumNodeProjectionError;
        out.maximumFacetError=surface.maximumFacetSampleError;
        out.surfaceDiagnostic=occtBoundaryReasonName(surface.reason);
        out.surfaceSeconds=std::chrono::duration<double>(Clock::now()-began).count();
        need(surface.reason!=OcctBoundaryReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");need(surface.reason!=OcctBoundaryReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");
        need(!spaceSurface||surface.reason!=OcctBoundaryReason::SourceChanged,NativeCadMeshReason::InternalError,"spaced-surface-source-changed");
        need(!spaceSurface||surface.published||surface.reason==OcctBoundaryReason::BudgetExceeded,
             NativeCadMeshReason::GeometryError,out.surfaceDiagnostic);
        need(surface.reason!=OcctBoundaryReason::BudgetExceeded,NativeCadMeshReason::BudgetExceeded,"native-plc-surface-budget");need(surface.published,NativeCadMeshReason::GeometryError,"native-plc-surface-refused");
        need(surface.relativeVolumeError<=options.maximumRelativeVolumeError,NativeCadMeshReason::GeometryError,"native-plc-source-volume-error");
        OcctStarOptions certificateOptions;certificateOptions.maximumWork=remaining();certificateOptions.maximumWorkingBytes=options.maximumPrimitiveWorkingBytes;certificateOptions.cancellation=o.cancellation;certificateOptions.maximumBoundaryNodes=bo.maximumNodes;certificateOptions.maximumBoundaryTriangles=bo.maximumTriangles;certificateOptions.maximumNodes=o.maximumNodes;certificateOptions.maximumCells=o.maximumCells;certificateOptions.maximumOutputBytes=o.maximumOutputBytes;
        OcctStageTrace shellCheck("exact_shell_certificate",tag);
        auto certificate=certifyOcctBoundarySurface(boundary,certificateOptions);
        shellCheck.finish(certificate.diagnostic);
        tick(certificate.workUsed);
        need(certificate.reason!=OcctStarReason::WorkBudget&&certificate.reason!=OcctStarReason::OutputBudget,NativeCadMeshReason::BudgetExceeded,"native-plc-source-certificate-budget");need(certificate.reason!=OcctStarReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");need(certificate.reason!=OcctStarReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");need(certificate.certified,NativeCadMeshReason::GeometryError,"native-plc-source-not-embedded");
        stage="seeds";trace.next(stage);out.surfaceSeconds=std::chrono::duration<double>(Clock::now()-began).count();
        const auto ns=boundary.nodeCount(),nf=boundary.cellCount();const auto frontCells=std::min<std::uint64_t>(o.maximumCells-nf,std::uint64_t(nf)*16);
        const auto apexStorage=completeDiscreteSeed?0:
            std::min<std::uint64_t>(256ULL*1024*1024,options.maximumPrimitiveWorkingBytes/4);
        const auto retained=std::uint64_t(ns)*1024+std::uint64_t(nf)*16384+frontCells*512+apexStorage+96ULL*1024*1024;
        need(retained<options.maximumPrimitiveWorkingBytes,NativeCadMeshReason::BudgetExceeded,"native-plc-front-storage-admission");
        mesh::Mesh pendingRecovery;
        const auto buildConstraintSeed=[&](const mesh::Mesh& sites,mesh::Mesh& seed,
                                          std::uint64_t maximumSeedTets=UINT64_MAX,
                                          bool reserveRetry=false){
            OcctStageTrace recovering("constraint-recovery",tag);
            const auto priorWork=out.constraintRecovery.work;
            const auto priorSeconds=out.constraintRecovery.seconds;
            auto constraintWork=earlyConstraintWork
                ?std::min(remaining(),earlyConstraintWork)
                :std::min<std::uint64_t>(remaining()/3,450'000'000);
            if(reserveRetry)constraintWork=std::min(constraintWork,remaining()/3);
            const auto constraints=bnd::detail::recoverTetConstraints(sites,{3,tag},
                constraintWork,options.maximumPrimitiveWorkingBytes/2,
                std::min(frontCells,maximumSeedTets),o.cancellation,seed,
                std::min<std::uint64_t>(64,o.maximumNodes-std::uint64_t(sites.nodeCount())));
            out.constraintRecovery={constraints.diagnostic,constraints.candidateCells,
                constraints.missingEdgesBefore,constraints.missingFacesBefore,
                constraints.missingEdgesAfter,constraints.missingFacesAfter,constraints.patches,
                constraints.seedCells,priorWork+constraints.work,priorSeconds+constraints.seconds,
                constraints.edgeWalks,constraints.facetWalks,constraints.walkFallbacks,
                constraints.steinerAttempts,constraints.steinerWork,constraints.steinerNodes,
                constraints.retainedSteinerNodes,constraints.piercingPatches};
            recovering.finish(constraints.diagnostic);
            tick(constraints.work);
            need(std::string_view(constraints.diagnostic)!="cancelled",NativeCadMeshReason::Cancelled,"cancelled");
            need(std::string_view(constraints.diagnostic)!="allocation_failed",NativeCadMeshReason::AllocationFailed,"allocation-failed");
        };
        OcctSolidClassifier classifier(OcctModelAccess::solid(model,tag));
        const auto tol=model.tolerance(3,tag);
        if(spatialSizing&&(options.improveNativeTetQuality||options.useAdvancingInteriorFront)){
            // Put size-aware interior sites into the first Delaunay point set.
            // The prescribed surface remains the same input constraint complex.
            OcctStageTrace placement("interior_sites",tag);
            struct SiteLimit {};
            mesh::Mesh sites;
            const auto siteBudget=std::min<std::uint64_t>(1'000'000,remaining()/8);
            constexpr auto nodeBytes=sizeof(Vec3)+sizeof(mesh::Classification);
            constexpr auto tetBytes=sizeof(mesh::Cell)+4*sizeof(Index);
            const auto boundaryBytes=std::uint64_t(ns)*nodeBytes+
                std::uint64_t(nf)*(sizeof(mesh::Cell)+3*sizeof(Index));
            const auto outputRoom=boundaryBytes<o.maximumOutputBytes?
                o.maximumOutputBytes-boundaryBytes:0;
            const auto siteCap=std::min({std::uint64_t(512),
                o.maximumNodes-std::uint64_t(ns),outputRoom/(nodeBytes+16*tetBytes)});
            // Bound the copied shell, point/hash arrays and cached edge sizes.
            const auto preparationBytes=(std::uint64_t(ns)+std::uint64_t(nf)+siteCap)*1024;
            std::uint64_t siteWork=0,siteQueries=0;
            const auto siteCharge=[&](std::uint64_t amount=1){
                if(amount>siteBudget-siteWork)throw SiteLimit{};
                tick(amount);siteWork+=amount;
            };
            try {if(siteCap&&siteBudget&&preparationBytes<=options.maximumPrimitiveWorkingBytes-retained){
                sites=appendOcctInteriorSites(boundary,*sizeField,mesh::Classification{3,tag},
                    std::size_t(siteCap),
                    [&](Vec3 p){
                        siteCharge();
                        if(siteQueries>=2048)throw SiteLimit{};
                        ++siteQueries;++out.nativeQueries;
                        classifier.Perform(gp_Pnt(p.x,p.y,p.z),tol);
                        return classifier.State()==TopAbs_IN;
                    },siteCharge,nullptr,!options.useAdvancingInteriorFront);
            }
            }catch(const SiteLimit&){
                // This optional preparation leaves the original shell intact.
            }catch(const std::bad_alloc&){
                // Preparation is transactional; retain the ordinary seed route.
            }
            placement.finish(sites.nodeCount()>ns?"seeded":"unchanged");
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
        const std::uint64_t maximumSpacedCandidateQueries=curvatureSurface||spatialSizing?
            (options.useLocalDelaunayRefinement&&liveTrial?40'000:24'000):3'000;
        const auto admitQueries=[&](std::uint64_t n) {
            need(!spaceSurface||(out.nativeQueries<=maximumSpacedCandidateQueries&&
                 n<=maximumSpacedCandidateQueries-out.nativeQueries),
                 NativeCadMeshReason::BudgetExceeded,"spaced-surface-candidate-query-budget");
        };
        std::map<std::array<std::uint64_t,3>,bool> nativeCache;
        const auto native=[&](Vec3 p){tick();const auto key=pointKey(p);auto known=nativeCache.find(key);if(known!=nativeCache.end())return known->second;
            need(nativeCache.size()<1000000,NativeCadMeshReason::BudgetExceeded,"native-plc-query-cache-budget");admitQueries(1);++out.nativeQueries;classifier.Perform(gp_Pnt(p.x,p.y,p.z),tol);const bool yes=classifier.State()==TopAbs_IN||classifier.State()==TopAbs_ON;nativeCache.emplace(key,yes);return yes;};
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
            need(centroidCache.size()<1000000,NativeCadMeshReason::BudgetExceeded,"native-plc-query-cache-budget");tick();
            admitQueries(1);classifier.Perform(gp_Pnt(p.x,p.y,p.z),tol);++out.nativeQueries;
            if(classifier.State()!=TopAbs_OUT){centroidCache.emplace(key,false);return false;}
            tick();BRepExtrema_DistShapeShape distance(BRepBuilderAPI_MakeVertex(gp_Pnt(p.x,p.y,p.z)).Vertex(),OcctModelAccess::solid(model,tag));
            need(distance.IsDone()&&std::isfinite(distance.Value())&&distance.Value()>=0,NativeCadMeshReason::GeometryError,"native-plc-centroid-distance-query");
            const bool admitted=distance.Value()<=out.chord;centroidCache.emplace(key,admitted);return admitted;
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
        auto seedCandidate=algo::buildTetCandidate(points);tick(seedCandidate.metrics.construction.conflictInsphereTests+seedCandidate.metrics.construction.conflictOrientationTests+seedCandidate.metrics.validationFaceUses+seedCandidate.metrics.validationHullTests+seedCandidate.metrics.candidateCells);
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
                tick(recovered.workUsed);out.randomCandidateProposals+=recovered.randomCandidatePointsProposed;out.mergedComponents+=recovered.mergedComponents;out.sourceBlockedCandidates+=recovered.sourceBlockedCandidates;out.otherComponentBlockedCandidates+=recovered.otherComponentBlockedCandidates;
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
            tick(closed.workUsed);
            trial.finish(closed.diagnostic);
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
            need(count<=std::min<std::uint64_t>(1'000'000,availableBytes/128),
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
            admission.maximumWorkingBytes=availableBytes-count*128;
            admission.maximumSnapshotBytes=std::min<std::uint64_t>(admission.maximumSnapshotBytes,admission.maximumWorkingBytes/4);
            admission.cancellation=o.cancellation;
            admission.maximumOutsideDistance=options.allowCentroidsWithinChordTolerance?out.chord:0;
            admission.useTrimmedRayWitnesses=useTrimmedRayWitnesses;
            OcctStageTrace batchTrace("original_cad_centroid_batch",tag);
            const auto batch=classifyOcctNativeInteriorPoints(model,tag,centers,admission,admitted,
                reuseStrictAdmissions?&strictAdmissions:nullptr,prior);
            batchTrace.finish(occtNativeValidationReasonName(batch.reason));
            tick(batch.workUsed);out.nativeQueries+=batch.nativeQueries+batch.trimmedRayAdmissions;
            out.nativeRayAdmissions+=batch.trimmedRayAdmissions;
            need(batch.reason!=OcctNativeValidationReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");
            need(batch.reason!=OcctNativeValidationReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");
            need(batch.reason!=OcctNativeValidationReason::BudgetExceeded,NativeCadMeshReason::BudgetExceeded,"native-plc-constraint-classifier-budget");
            need(batch.completed&&batch.sourceModelUnmodified&&admitted.size()==centers.size()&&
                 batch.nativeQueries+batch.strictAdmissionReuses+batch.trimmedRayAdmissions==centers.size(),
                 NativeCadMeshReason::GeometryError,"native-plc-constraint-classifier-refused");
            for(std::size_t i=0;i<centers.size();++i){
                tick();
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
                centroidCache.emplace(key,true);
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
                bo.maximumWork=remaining();const auto facets=auditOcctBoundaryFacets(model,tag,candidate,bo);tick(facets.work);certificateOptions.maximumWork=remaining();const auto cert=certifyOcctBoundarySurface(candidate,certificateOptions);tick(cert.workUsed);
                
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
                    if(changed.flips){bo.maximumWork=remaining();const auto facets=auditOcctBoundaryFacets(model,tag,candidate,bo);tick(facets.work);certificateOptions.maximumWork=remaining();const auto cert=certifyOcctBoundarySurface(candidate,certificateOptions);tick(cert.workUsed);
                        
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
                    const std::array input{OcctWireBoundary{tag,&boundary}};const auto refined=refineOcctWireRefinement(model,input,wire);tick(refined.workUsed);need(refined.reason!=OcctWireRefinementReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");need(refined.reason!=OcctWireRefinementReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");
                    if(refined.reason==OcctWireRefinementReason::BudgetExceeded){
                        // A refused local subdivision leaves the boundary
                        // unchanged. Other bounded strategies can still fit.
                        localRefinementBound=true;need(out.workUsed<maximumWork,NativeCadMeshReason::BudgetExceeded,"native-plc-adaptive-wire-budget");continue;
                    }
                    changed|=refined.published;
                }
                if(!changed){std::vector<OcctFacetSelection> selected;const auto& candidates=result.blockedSourceTriangles.empty()?result.unresolvedSourceTriangles:result.blockedSourceTriangles;for(auto i:candidates){need(i>=0&&i<boundary.cellCount(),NativeCadMeshReason::InternalError,"native-plc-diagnostic-face-index");selected.push_back({tag,i});if(selected.size()==8)break;}
                    if(!selected.empty()){OcctFacetRefinementOptions facet;facet.selections=selected;facet.maximumProjectionDistance=out.chord;facet.maximumAddedNodes=8;facet.limits=certificateOptions;facet.limits.maximumWork=remaining();const std::array input{OcctFacetBoundary{tag,&boundary}};const auto refined=refineOcctFacets(model,input,facet);tick(refined.workUsed);need(refined.reason!=OcctFacetRefinementReason::Cancelled,NativeCadMeshReason::Cancelled,"cancelled");need(refined.reason!=OcctFacetRefinementReason::AllocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");if(refined.reason==OcctFacetRefinementReason::BudgetExceeded){localRefinementBound=true;need(out.workUsed<maximumWork,NativeCadMeshReason::BudgetExceeded,"native-plc-adaptive-facet-budget");}changed=refined.published;out.sourceFacetNodes+=refined.addedNodes;}}
                if(!changed)break;
                // Reclassify the retained partial mesh against the refined PLC.
                // Any crossed cell is discarded before another recovery attempt.
                bnd::detail::PlcGiftWrapping checked;checked.chargeWork=tick;checked.cancellation=stop.get_token();for(Index i=0;i<boundary.nodeCount();++i)checked.vertices.push_back({boundary.node(i),checked.vertices.size()+1});
                for(Index i=previousNodes;i<inside.nodeCount();++i)checked.vertices.push_back({inside.node(i),checked.vertices.size()+1});for(Index i=0;i<boundary.cellCount();++i){const auto* p=boundary.cellNodes(i);checked.boundary.push_back({unsigned(p[0]),unsigned(p[1]),unsigned(p[2])});}need(checked.initialize(),NativeCadMeshReason::GeometryError,"native-plc-refined-front-invalid");
                const auto offset=boundary.nodeCount()-previousNodes;for(Index i=0;i<inside.cellCount();++i){bnd::TetFacetCell c;const auto* p=inside.cellNodes(i);for(unsigned j=0;j<4;++j)c.vertices[j]=unsigned(p[j]>=previousNodes?p[j]+offset:p[j]);bnd::TetFacetCell valid;const auto face=bnd::tet_facet_detail::outwardFace(c,0);if(checked.valid(bnd::detail::PlcGiftWrapping::Front{face},c.vertices[0],valid))checked.insert(c);}
                checked.interiorAdmission=native;checked.centroidAdmission=nativeCentroid;checked.run(static_cast<unsigned>(std::min<std::uint64_t>(frontCells,UINT_MAX)),10,stop.get_token(),100000);
                mesh::Mesh next;for(Index i=0;i<Index(checked.vertices.size());++i)next.addNode(checked.vertices[i].position,i<boundary.nodeCount()?boundary.nodeClass(i):mesh::Classification{3,tag});for(auto c:checked.cells){std::array<Index,4> ids;for(unsigned j=0;j<4;++j)ids[j]=c.vertices[j];next.addCell(mesh::CellType::Tet4,ids,{3,tag});}inside.swap(next);
                bo.maximumWork=remaining();const auto audited=auditOcctBoundaryFacets(model,tag,boundary,bo);tick(audited.work);need(audited.completed&&audited.sourceModelUnmodified,NativeCadMeshReason::GeometryError,"native-plc-adaptive-facet-error");out.maximumFacetError=std::max(out.maximumFacetError,audited.maximumFacetSampleError);
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
        if(options.useAdvancingInteriorFront){
            OcctStageTrace construction("interior_front",tag);
            const auto started=Clock::now();
            const auto used=strictAdmissions.retainedPayloadBytes();
            const auto bytes=std::uint64_t(completed.nodeCount())*(sizeof(Vec3)+sizeof(mesh::Classification))+
                std::uint64_t(completed.cellCount())*sizeof(mesh::Cell)+
                std::uint64_t(completed.connectivityCount())*sizeof(Index);
            const auto reserve=1024*(std::uint64_t(completed.nodeCount())+
                std::uint64_t(completed.countCells(mesh::CellType::Tet4)));
            const auto room=remaining();
            const auto budget=room>reserve?std::min<std::uint64_t>(42'000'000,(room-reserve)/2):0;
            const auto queries=out.nativeQueries<maximumSpacedCandidateQueries?
                std::min<std::uint64_t>(12'000,(maximumSpacedCandidateQueries-out.nativeQueries)/2):0;
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
                out.interiorFrontDiagnostic=!built.initialized?"scaffold_refused":
                    built.targetsSatisfied?"targets_satisfied":built.budgetLimited?"budget_limited":"targets_unmet";
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
        if(!options.useAdvancingInteriorFront&&(spatialSizing||(useLocalDelaunay&&!refinedBeforeMaterial))){
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
        trace.next("cad_validation");
        for(Index i=0;i<completed.nodeCount();++i)if(completed.nodeClass(i).dim==3)need(native(completed.node(i)),NativeCadMeshReason::GeometryError,"native-plc-final-interior-node-outside");
        for(Index i=0;i<completed.cellCount();++i)if(completed.cell(i).type==mesh::CellType::Tet4){const auto* p=completed.cellNodes(i);Vec3 center{};for(unsigned j=0;j<4;++j)center=center+completed.node(p[j])*.25;need(nativeCentroid(center),NativeCadMeshReason::GeometryError,"native-plc-final-centroid-outside");}
        const auto quality=algo::analyzeMixedMeshQuality(completed);tick(completed.cellCount());need(!quality.allocationFailed,NativeCadMeshReason::AllocationFailed,"allocation-failed");need(quality.valid,NativeCadMeshReason::MeshingFailed,"native-plc-final-quality-invalid");
        need(!spaceSurface||quality.minimumQuality>=o.minimumElementQuality,NativeCadMeshReason::QualityRequirementNotMet,"spaced-surface-quality-requirement");
        need(std::abs(quality.totalVolume-model.volume(tag))<=model.volume(tag)*options.maximumRelativeVolumeError,NativeCadMeshReason::GeometryError,"native-plc-final-volume-error");
        std::vector<Vec3> centroids;for(Index i=0;i<completed.cellCount();++i)if(completed.cell(i).type==mesh::CellType::Tet4){tick();const auto* p=completed.cellNodes(i);Vec3 center{};for(unsigned j=0;j<4;++j)center=center+completed.node(p[j])*.25;centroids.push_back(center);}
        OcctNativeValidationOptions validation;validation.maximumThreads=options.maximumNativeValidationThreads;validation.maximumWork=remaining();validation.maximumPoints=o.maximumCells;validation.maximumWorkingBytes=options.maximumPrimitiveWorkingBytes;validation.cancellation=o.cancellation;validation.maximumOutsideDistance=options.allowCentroidsWithinChordTolerance?out.chord:0;
        validation.useTrimmedRayWitnesses=useTrimmedRayWitnesses;
        out.nativeValidation=validateOcctNativeInteriorPoints(model,tag,centroids,validation,
            strictAdmissions.retainedPayloadBytes()?&strictAdmissions:nullptr);tick(out.nativeValidation.workUsed);
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
            out.chordSampleValidation=validateOcctNativeInteriorPoints(model,tag,samples,validation);tick(out.chordSampleValidation.workUsed);
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
