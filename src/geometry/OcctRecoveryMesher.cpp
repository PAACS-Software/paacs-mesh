// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/BoundaryConformingMesher.hpp"
#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/topo/Workspace.hpp"
#include <numeric>
#include <queue>
#include <set>
#include "OcctRecoveryMesher.hpp"
#include "OcctPlcFrontRecovery.hpp"
#include "../bnd/TetEdgeGeneralCavityCadPolicy.hpp"
#include "OcctModelInternal.hpp"
#include "OcctNativeWireAuthority.hpp"
#include "OcctSolidClassifier.hpp"
#include <Standard_OutOfMemory.hxx>
#include <Standard_Failure.hxx>
#include <algorithm>
#include <cmath>
#include <map>
#include <new>

namespace paacs::geo::detail {
namespace {
struct CoreStop {const char* reason;OcctPlcRecoveryReason code=OcctPlcRecoveryReason::InvalidMesh;};
void coreNeed(bool ok,const char* why,OcctPlcRecoveryReason code=OcctPlcRecoveryReason::InvalidMesh){if(!ok)throw CoreStop{why,code};}
struct NativeBoundaryOwner final:algo::VolumeBoundaryOwnershipProvider {
 int body=1;std::set<algo::FacePersistentId> allowedFaces;
 algo::VolumeBoundaryProviderStamp stamp()const noexcept override{return {1,1};}
 algo::VolumeFaceOwnership ownership(algo::FacePersistentId f)const override{
  coreNeed(allowedFaces.contains(f),"unknown-native-artifact-face");return {stamp(),f,algo::VolumeProviderAvailability::Available,{{body,1}}};}
 algo::ProviderRecoveryEnvelope supplied;algo::ProviderRecoveryEnvelope recoveryEnvelope(std::span<const Vec3>)const override{return supplied;}
};
struct FixedCoreOutcome {OcctPlcRecoveryReason reason=OcctPlcRecoveryReason::InternalError;bool published=false;const char* diagnostic="internal-error";std::uint64_t workUsed=0,generalReservedWork=0,nativeSeedQueries=0,steinerNodes=0;int edgeReason=-1,facetReason=-1,candidateReason=-1;std::array<Index,2> failedEdge{-1,-1};bool failedEdgeAvailable=false;};
FixedCoreOutcome buildFixedPlcCore(const mesh::Mesh& boundary,int bodyTag,const OcctStarOptions& options,mesh::Mesh& destination,const OcctModel* nativePieceModel,const OcctPartitionRecoveryOptions& recovery) noexcept {
 bnd::detail::ScopedCadGeneralPruning pruning(true);
 FixedCoreOutcome receipt;const char* why="internal-error";bool complete=false;OcctPlcRecoveryReason code=OcctPlcRecoveryReason::InternalError;
 const std::size_t nodes=boundary.nodeCount(),triangles=boundary.cellCount();std::size_t cells=0;double plc=0;
 const auto charge=[&](std::initializer_list<std::uint64_t> values){for(const auto n:values){coreNeed(n<=UINT64_MAX-receipt.workUsed,"work-overflow",OcctPlcRecoveryReason::BudgetExceeded);receipt.workUsed+=n;}coreNeed(receipt.workUsed<=options.maximumWork,"core-phase-work-budget",OcctPlcRecoveryReason::BudgetExceeded);if(options.cancellation&&options.cancellation->load())throw CoreStop{"cancelled",OcctPlcRecoveryReason::Cancelled};};
 const auto tick=[&](std::uint64_t n=1){charge({n});};
 try {coreNeed(bodyTag>0&&nodes>=4&&nodes<=10000&&triangles>=4&&triangles<=20000,"fixed-plc-input-budget",OcctPlcRecoveryReason::BudgetExceeded);
  coreNeed(std::uint64_t(nodes+triangles)*16384<=options.maximumWorkingBytes,"working-payload-budget",OcctPlcRecoveryReason::BudgetExceeded);NativeBoundaryOwner owner;owner.body=bodyTag;

  coreNeed(nativePieceModel&&nativePieceModel->solid(bodyTag),"native-artifact-authority-required",OcctPlcRecoveryReason::InvalidRequest);
  mesh::Mesh original=boundary,surface=boundary;algo::SurfaceSetResult set;set.setStamp={1,1};
  const auto* solid=nativePieceModel->solid(bodyTag);
  std::set<int> nativeFaces;for(int tag:solid->boundaryFaceTags()){tick();nativeFaces.insert(tag);owner.allowedFaces.insert(tag);}
  std::map<int,std::unique_ptr<algo::SurfaceFaceResult>> faces;
  std::map<int,std::map<Index,std::size_t>> localByGlobal;
  std::map<std::array<Index,2>,std::set<int>> edgeOwners;
  for(Index ci=0;ci<boundary.cellCount();++ci){tick();const auto& cell=boundary.cell(ci);
    coreNeed(cell.type==mesh::CellType::Tri3&&cell.onEntity.dim==2&&nativeFaces.contains(cell.onEntity.tag),"native-boundary-face-ownership");
    const int tag=cell.onEntity.tag;auto& face=faces[tag];
    if(!face){face=std::make_unique<algo::SurfaceFaceResult>();face->faceStamp={std::uint64_t(tag),1};face->faceClass={2,tag};}
    const auto* ids=boundary.cellNodes(ci);std::array<std::size_t,3> local;
    for(unsigned j=0;j<3;++j){tick();coreNeed(ids[j]>=0&&std::uint64_t(ids[j])<nodes,"facet-index");
      const auto p=boundary.node(ids[j]);const auto cls=boundary.nodeClass(ids[j]);
      coreNeed(cls.dim>=0&&cls.dim<=2&&cls.tag>0,"native-boundary-node-class");
      auto [where,inserted]=localByGlobal[tag].emplace(ids[j],face->localNodes.size());local[j]=where->second;
      if(inserted){face->localNodes.push_back({p,cls});face->publishedGlobalByLocalNode.push_back(ids[j]);}
      auto edge=std::array<Index,2>{ids[j],ids[(j+1)%3]};std::sort(edge.begin(),edge.end());edgeOwners[edge].insert(tag);
    }
    face->localCells.push_back({local,{2,tag}});face->publishedGlobalCellByLocalCell.push_back(ci);
    plc+=boundary.node(ids[0]).dot(boundary.node(ids[1]).cross(boundary.node(ids[2])))/6.;
  }
  coreNeed(plc>0,"plc-volume");
  for(const auto& [endpoints,ownedBy]:edgeOwners){tick();
    const auto ca=boundary.nodeClass(endpoints[0]),cb=boundary.nodeClass(endpoints[1]);
    const int edgeTag=actualOcctWire(*nativePieceModel,ca,cb,ownedBy,[&]{tick();});
    if(ownedBy.size()>1)coreNeed(edgeTag>0,"ambiguous-native-shared-edge");
    if(edgeTag<=0)continue;
    for(const auto id:endpoints){tick();const algo::SharedCurveNodeKey key{std::uint64_t(edgeTag),1,std::uint64_t(id)};
      set.sharedLedger.emplace(key,algo::SharedCurveNodeBinding{key,id,boundary.nodeClass(id)});}
    for(const int tag:ownedBy){tick();algo::PersistentConstraintEdge edge;
      edge.faceStamp=faces.at(tag)->faceStamp;edge.classification={1,edgeTag};edge.chainPersistentId=edgeTag;
      edge.publishedGlobalEndpoints=endpoints;faces.at(tag)->persistentEdges.push_back(edge);}
  }
  for(auto& [tag,face]:faces){tick();set.faces.push_back(std::move(face));}
  
  auto built=algo::buildClosedTri3Artifact(set,surface,set.setStamp,owner.stamp(),owner);
  charge({built.metrics.edgeEntries,built.metrics.facetEntries,built.metrics.treeNodes,built.metrics.axisSelectionVisits,built.metrics.sortComparisons,built.metrics.nodePairVisits,built.metrics.aabbTests,built.metrics.edgeEdgePairs,built.metrics.edgeFacetPairs,built.metrics.facetFacetPairs,built.metrics.exactEdgeEdgeTests,built.metrics.exactEdgeFacetTests,built.metrics.exactFacetFacetTests,built.metrics.legalSharedContacts,built.metrics.rejectedContacts});
  coreNeed(built.failure.reason!=algo::BoundaryConformingV3Reason::AllocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);
  coreNeed(bool(built.artifact),"artifact-rejected");auto& artifact=*built.artifact;

  mesh::Mesh volumePoints=surface;
  if(recovery.nativeKernelSeed){
    coreNeed(std::uint64_t(nodes+triangles+1)*16384<=options.maximumWorkingBytes,"seed-working-payload-budget",OcctPlcRecoveryReason::BudgetExceeded);
    auto seedLimits=options;seedLimits.maximumWork=options.maximumWork-receipt.workUsed;
    seedLimits.nativeAdmissionModel=nativePieceModel;seedLimits.nativeAdmissionBodyTag=bodyTag;
    mesh::Mesh seedVolume;const auto star=buildOcctStarVolume(surface,bodyTag,seedLimits,seedVolume);
    receipt.nativeSeedQueries+=star.nativeCandidateQueries;tick(star.workUsed);
    coreNeed(star.reason!=OcctStarReason::Cancelled,"cancelled",OcctPlcRecoveryReason::Cancelled);
    coreNeed(star.reason!=OcctStarReason::AllocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);
    coreNeed(star.reason!=OcctStarReason::NativeQueryFailed,"native-query-failed",OcctPlcRecoveryReason::NativeQueryFailed);
    coreNeed(star.reason!=OcctStarReason::WorkBudget&&star.reason!=OcctStarReason::OutputBudget,"seed-budget",OcctPlcRecoveryReason::BudgetExceeded);
    coreNeed(star.reason!=OcctStarReason::InvalidBoundary,"invalid-seed-boundary",OcctPlcRecoveryReason::InvalidMesh);
    coreNeed(star.reason!=OcctStarReason::InternalError,"seed-internal-error",OcctPlcRecoveryReason::InternalError);
    if(star.published){
      paacs::geo::detail::OcctSolidClassifier strict(OcctModelAccess::solid(*nativePieceModel,bodyTag));tick();++receipt.nativeSeedQueries;
      strict.Perform(gp_Pnt(star.kernel.x,star.kernel.y,star.kernel.z),nativePieceModel->tolerance(3,bodyTag));
      if(strict.State()==TopAbs_IN){
        coreNeed(std::uint64_t(volumePoints.nodeCount())<options.maximumNodes,"seed-node-budget",OcctPlcRecoveryReason::BudgetExceeded);
        volumePoints.addNode(star.kernel,{3,bodyTag});
      }
    }
  }

  const auto volumeNodes=volumePoints.nodeCount();
  auto env=algo::buildTetRecoveryEnvelope(algo::RecoveryEnvelopeMode::DerivedFiniteBox,surface.nodes(),owner.stamp(),owner);
  coreNeed(env.failure.reason!=algo::BoundaryConformingV3Reason::AllocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);coreNeed(bool(env.envelope),"envelope");
  {
    constexpr double fraction=1.;coreNeed(std::isfinite(fraction)&&fraction>0&&fraction<=16,"guard-span-profile");
    Vec3 lo=surface.node(0),hi=lo;for(const auto q:surface.nodes()){tick();lo={std::min(lo.x,q.x),std::min(lo.y,q.y),std::min(lo.z,q.z)};hi={std::max(hi.x,q.x),std::max(hi.y,q.y),std::max(hi.z,q.z)};}
    const auto margin=std::max({hi.x-lo.x,hi.y-lo.y,hi.z-lo.z})*fraction;coreNeed(std::isfinite(margin)&&margin>0,"guard-margin");
    const Vec3 low{lo.x-margin,lo.y-margin,lo.z-margin},high{hi.x+margin,hi.y+margin,hi.z+margin};
    owner.supplied.providerStamp=owner.stamp();owner.supplied.outwardHullFaces=env.envelope->outwardHullFaces;
    for(std::size_t i=0;i<8;++i){tick();owner.supplied.guardCoordinates.push_back({(i&1U)?high.x:low.x,(i&2U)?high.y:low.y,(i&4U)?high.z:low.z});}
    env=algo::buildTetRecoveryEnvelope(algo::RecoveryEnvelopeMode::ProviderSuppliedConvexHull,volumePoints.nodes(),owner.stamp(),owner);
    coreNeed(env.failure.reason!=algo::BoundaryConformingV3Reason::AllocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);coreNeed(bool(env.envelope),"wider-envelope-rejected");
  }
  mesh::Mesh construction=volumePoints;
  for(auto p:env.envelope->guardCoordinates)construction.addNode(p,{3,0});
  auto candidate=algo::buildTetCandidate(construction.nodes());receipt.candidateReason=int(candidate.reason);
  charge({candidate.metrics.inputPoints,candidate.metrics.uniquePoints,candidate.metrics.candidateCells,candidate.metrics.interiorFaces,candidate.metrics.boundaryFaces,candidate.metrics.validationFaceUses,candidate.metrics.validationHullTests,candidate.metrics.discardedCandidateWork,candidate.metrics.construction.originalOrientationAuditTests,candidate.metrics.construction.discardedConstructionWork});charge({candidate.metrics.construction.inputPoints,candidate.metrics.construction.pointsInserted,candidate.metrics.construction.pointsDuplicate,candidate.metrics.construction.pointsSkipped,candidate.metrics.construction.walkFallbackScans,candidate.metrics.construction.cavityTrimResurrections,candidate.metrics.construction.bootstrapDistinctTests,candidate.metrics.construction.bootstrapOrientationTests,candidate.metrics.construction.walkCellVisits,candidate.metrics.construction.walkOrientationTests,candidate.metrics.construction.conflictQueries,candidate.metrics.construction.conflictOrientationTests,candidate.metrics.construction.conflictInsphereTests,candidate.metrics.construction.fallbackCellTests,candidate.metrics.construction.cavityCellVisits,candidate.metrics.construction.cavityNeighbourVisits,candidate.metrics.construction.trimCellVisits,candidate.metrics.construction.trimOrientationTests,candidate.metrics.construction.shellFacesEmitted,candidate.metrics.construction.adjacencyAuditFaceVisits,candidate.metrics.construction.outputCells});coreNeed(candidate.reason!=algo::TetCandidateReason::AllocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);if(!candidate.built()){
    Vec3 lo=surface.node(0),hi=lo;for(auto p:surface.nodes()){tick();lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}
    const auto span=std::max({hi.x-lo.x,hi.y-lo.y,hi.z-lo.z});const auto oldLo=env.envelope->guardCoordinates.front();
    for(const double fraction:{.125,.5,1.}){tick(8);owner.supplied={};owner.supplied.providerStamp=owner.stamp();owner.supplied.outwardHullFaces=env.envelope->outwardHullFaces;
      for(auto p:env.envelope->guardCoordinates)owner.supplied.guardCoordinates.push_back({p.x==oldLo.x?lo.x-span*fraction:hi.x+span*fraction,p.y==oldLo.y?lo.y-span*fraction:hi.y+span*fraction,p.z==oldLo.z?lo.z-span*fraction:hi.z+span*fraction});
      auto wider=algo::buildTetRecoveryEnvelope(algo::RecoveryEnvelopeMode::ProviderSuppliedConvexHull,surface.nodes(),owner.stamp(),owner);
      coreNeed(wider.failure.reason!=algo::BoundaryConformingV3Reason::AllocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);if(!wider.envelope)continue;construction=volumePoints;for(auto p:wider.envelope->guardCoordinates)construction.addNode(p,{3,0});
      candidate=algo::buildTetCandidate(construction.nodes());receipt.candidateReason=int(candidate.reason);
      charge({candidate.metrics.inputPoints,candidate.metrics.uniquePoints,candidate.metrics.candidateCells,candidate.metrics.interiorFaces,candidate.metrics.boundaryFaces,candidate.metrics.validationFaceUses,candidate.metrics.validationHullTests,candidate.metrics.discardedCandidateWork,candidate.metrics.construction.originalOrientationAuditTests,candidate.metrics.construction.discardedConstructionWork});charge({candidate.metrics.construction.inputPoints,candidate.metrics.construction.pointsInserted,candidate.metrics.construction.pointsDuplicate,candidate.metrics.construction.pointsSkipped,candidate.metrics.construction.walkFallbackScans,candidate.metrics.construction.cavityTrimResurrections,candidate.metrics.construction.bootstrapDistinctTests,candidate.metrics.construction.bootstrapOrientationTests,candidate.metrics.construction.walkCellVisits,candidate.metrics.construction.walkOrientationTests,candidate.metrics.construction.conflictQueries,candidate.metrics.construction.conflictOrientationTests,candidate.metrics.construction.conflictInsphereTests,candidate.metrics.construction.fallbackCellTests,candidate.metrics.construction.cavityCellVisits,candidate.metrics.construction.cavityNeighbourVisits,candidate.metrics.construction.trimCellVisits,candidate.metrics.construction.trimOrientationTests,candidate.metrics.construction.shellFacesEmitted,candidate.metrics.construction.adjacencyAuditFaceVisits,candidate.metrics.construction.outputCells});
      coreNeed(candidate.reason!=algo::TetCandidateReason::AllocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);if(candidate.built())break;
    }
  }coreNeed(candidate.reason!=algo::TetCandidateReason::AllocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);
  coreNeed(candidate.built(),"candidate",OcctPlcRecoveryReason::UnsupportedGeometry);
  std::vector<Index> selected(candidate.candidate->cells.size());std::iota(selected.begin(),selected.end(),Index(0));
  std::vector<Index> globals;std::vector<mesh::Classification> classes;
  for(auto input:candidate.candidate->pointToInput){globals.push_back(input);classes.push_back(construction.nodeClass(input));}
  std::vector<algo::TetImportProtectedEdge> edges;std::vector<algo::TetImportProtectedFace> facets;
  for(const auto& e:artifact.edges)edges.push_back({e.key.canonicalGlobalNodes,e.classification});
  for(const auto& f:artifact.facets)facets.push_back({f.directedGlobalNodes,f.classification});
  // Only candidate entities proven exact-present are protected at import.
  // The complete immutable artifact still schedules every prescribed entity.
  std::set<std::array<Index,2>> presentEdges;std::map<std::array<Index,3>,unsigned> presentFaces;
  for(const auto& tet:candidate.candidate->cells){tick(10);std::array<Index,4> ids;for(unsigned j=0;j<4;++j)ids[j]=globals[tet.vertices[j]];
    for(unsigned a=0;a<4;++a)for(unsigned b=a+1;b<4;++b){std::array<Index,2> key{ids[a],ids[b]};std::sort(key.begin(),key.end());presentEdges.insert(key);}
    for(unsigned skip=0;skip<4;++skip){std::array<Index,3> key;unsigned at=0;for(unsigned j=0;j<4;++j)if(j!=skip)key[at++]=ids[j];std::sort(key.begin(),key.end());++presentFaces[key];}}
  std::erase_if(edges,[&](const auto& row){tick();return !presentEdges.contains(row.existingGlobalNodes);});
  std::erase_if(facets,[&](const auto& row){tick();auto key=row.existingGlobalNodes;std::sort(key.begin(),key.end());const auto it=presentFaces.find(key);return it==presentFaces.end()||it->second!=2;});
  auto imported=algo::importTetCandidateWithProtectedBoundary(*candidate.candidate,selected,globals,classes,construction,{3,0},edges,facets);
  coreNeed(imported.reason!=algo::TetCandidateReason::AllocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);
  coreNeed(imported.imported(),"import");
  charge({imported.result->metrics.selectedCells,imported.result->metrics.importedNodes,imported.result->metrics.importedCells,imported.result->metrics.importEdgeValidationCellVisits,imported.result->metrics.importEdgeBuildCellVisits,imported.result->metrics.importEdgeUses,imported.result->metrics.importUniqueEdges,imported.result->metrics.importFaceValidationCellVisits,imported.result->metrics.importFaceBuildCellVisits,imported.result->metrics.importFaceUses,imported.result->metrics.importUniqueFaces,imported.result->metrics.importEdgeAudits,imported.result->metrics.importFaceAudits,imported.result->metrics.importPerCellGlobalAudits,imported.result->metrics.importCandidateValidationCalls,imported.result->metrics.importCandidateValidationPointRows,imported.result->metrics.importCandidateValidationCellRows,imported.result->metrics.importCandidateValidationFaceUses,imported.result->metrics.importCandidateValidationFaceRows,imported.result->metrics.importCandidateValidationBoundaryRows,imported.result->metrics.importCandidateValidationDualNeighbourVisits,imported.result->metrics.importCandidateValidationBoundaryEdgeUses,imported.result->metrics.importCandidateValidationBoundaryNeighbourVisits,imported.result->metrics.importCandidateValidationHullTests,imported.result->metrics.importCandidateValidationSlotsInitialized,imported.result->metrics.importSelectionRowsCopied,imported.result->metrics.importSelectionSortComparisons,imported.result->metrics.importSelectionDuplicateComparisons,imported.result->metrics.importSelectionCellVisits,imported.result->metrics.importSelectionVertexVisits,imported.result->metrics.importSelectionPointSlotsInitialized,imported.result->metrics.importGlobalMapRows,imported.result->metrics.importGlobalMapComparisons,imported.result->metrics.importGlobalMapInsertions,imported.result->metrics.importLocalPointSlotsInitialized,imported.result->metrics.importLocalGlobalSlotsInitialized,imported.result->metrics.importLocalCellSlotsInitialized,imported.result->metrics.importPersistentNodeSlotsInitialized,imported.result->metrics.partialImportRowsCopied,imported.result->metrics.partialImportGlobalSlotsInitialized,imported.result->metrics.partialImportCandidateMappingVisits,imported.result->metrics.partialImportEdgeValidationRows,imported.result->metrics.partialImportFaceValidationRows,imported.result->metrics.partialImportProtectedEndpointVisits,imported.result->metrics.partialImportAuthorityPointValidations,imported.result->metrics.partialImportDuplicateKeyComparisons,imported.result->metrics.partialImportResolutionNodeVisits,imported.result->metrics.partialImportEdgeResolutionLookups,imported.result->metrics.partialImportFaceResolutionLookups,imported.result->metrics.partialImportEdgeDataReads,imported.result->metrics.partialImportFaceDataCopies,imported.result->metrics.partialImportFaceIncidenceRowsCopied,imported.result->metrics.partialImportSortComparisons});candidate.candidate.reset();auto& workspace=*imported.result->workspace;
  std::map<Index,topo::NodeHandle> local;std::map<topo::NodeHandle,Index> global;
  for(const auto& b:imported.result->nodeBindings){local.emplace(b.existingGlobalIndex,b.local);global.emplace(b.local,b.existingGlobalIndex);}
  auto contextBuilt=bnd::buildTetRecoveryContext(workspace,{3,0});
  coreNeed(contextBuilt.reason!=bnd::TetRecoveryContextReason::AllocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);coreNeed(contextBuilt.built(),"context");
  charge({contextBuilt.context->metrics.nodes,contextBuilt.context->metrics.cells,contextBuilt.context->metrics.edges,contextBuilt.context->metrics.faces,contextBuilt.context->metrics.auditVisits,contextBuilt.context->metrics.indexNodes});auto& context=*contextBuilt.context;
  bnd::TetFacetSetState state;state.workspaceStamp=topo::WorkspaceStamp::of(workspace);state.artifactStamp=artifact.surfaceSetStamp;
  bnd::TetEdgeGeneralRecoveryPolicy policy;policy.enabled=true;const auto generalCap=recovery.maximumGeneralWorkPerAttempt;policy.maximumTotalWorkUnits=std::min<std::uint64_t>(generalCap,options.maximumWork-receipt.workUsed);policy.limits.maxWorkUnits=policy.maximumTotalWorkUnits;policy.limits.maxCavityCells=128;policy.limits.maxCavityVertices=32;policy.limits.maxEnumeratedTets=100000;policy.limits.maxSearchStates=100000;policy.limits.maximumSteinerCandidates=recovery.maximumSteinerCandidates;
  coreNeed(policy.maximumTotalWorkUnits>0,"core-phase-work-budget",OcctPlcRecoveryReason::BudgetExceeded);
  // Reserve a documented strict upper bound for every general-cavity solve
  // in this invocation; the current public scheduler does not export its used count.
  receipt.generalReservedWork=policy.maximumTotalWorkUnits;charge({receipt.generalReservedWork});
  bnd::TetFacetSetMixedPresentBatchSink mixed;
  const auto recovered=bnd::recoverTetFacetSetIndexed(artifact,local,{3,0},state,workspace,context,policy,nullptr,nullptr,&mixed);
  charge({recovered.schedulerMetrics.routeAssessments,recovered.schedulerMetrics.dependencyEdges,recovered.schedulerMetrics.dependencyInvalidations,recovered.schedulerMetrics.selectedRebuilds,recovered.schedulerMetrics.batchPreflightAssessments,recovered.schedulerMetrics.alreadyPresentBatchEntities,recovered.schedulerMetrics.protectedConflicts,recovered.schedulerMetrics.conflictSccVisits,recovered.schedulerMetrics.fullStageRescans});charge({mixed.metrics.presentBindingBatches,mixed.metrics.presentBindingEntities,mixed.metrics.stateCopies,mixed.metrics.copiedRows,mixed.metrics.schedulerBindingAuditRows,mixed.metrics.pendingRowsVisited,mixed.metrics.dependencyPairsTested,mixed.metrics.exactPresentSelectedRebuilds,mixed.metrics.exactPresentDependencyTests,mixed.metrics.nontrivialRecoveries,mixed.metrics.modeledBatchPreparationAllocationRequests});charge({mixed.metrics.assessmentCounters.cellsScanned,mixed.metrics.assessmentCounters.crossedEdgeStars,mixed.metrics.assessmentCounters.closureInsertions,mixed.metrics.assessmentCounters.fixedPointPasses,mixed.metrics.assessmentCounters.apexAttempts,mixed.metrics.assessmentCounters.tetTetIntersectionTests});charge({mixed.metrics.recoveryCounters.cellsScanned,mixed.metrics.recoveryCounters.crossedEdgeStars,mixed.metrics.recoveryCounters.closureInsertions,mixed.metrics.recoveryCounters.fixedPointPasses,mixed.metrics.recoveryCounters.apexAttempts,mixed.metrics.recoveryCounters.tetTetIntersectionTests});charge({mixed.metrics.assessmentBindingWork.edgeBindingAuditCalls,mixed.metrics.assessmentBindingWork.edgeBindingRowsVisited,mixed.metrics.assessmentBindingWork.edgeAliasRowsVisited,mixed.metrics.assessmentBindingWork.facetBindingAuditCalls,mixed.metrics.assessmentBindingWork.facetBindingRowsVisited,mixed.metrics.assessmentBindingWork.facetBindingRequestRowsVisited,mixed.metrics.assessmentBindingWork.facetBindingDuplicateRowsCompared,mixed.metrics.assessmentBindingWork.facetBindingReceiptBuilds,mixed.metrics.assessmentBindingWork.facetBindingReceiptReuses,mixed.metrics.assessmentBindingWork.facetBindingReceiptRowsCopied,mixed.metrics.assessmentBindingWork.facetBindingReceiptLookups,mixed.metrics.assessmentBindingWork.facetBindingReceiptReserveCalls,mixed.metrics.assessmentBindingWork.edgeDataCopyRows,mixed.metrics.assessmentBindingWork.faceDataCopyRows,mixed.metrics.assessmentBindingWork.faceIncidenceRowsCopied,mixed.metrics.assessmentBindingWork.stateCopies,mixed.metrics.assessmentBindingWork.stateCopyRows});charge({mixed.metrics.recoveryBindingWork.edgeBindingAuditCalls,mixed.metrics.recoveryBindingWork.edgeBindingRowsVisited,mixed.metrics.recoveryBindingWork.edgeAliasRowsVisited,mixed.metrics.recoveryBindingWork.facetBindingAuditCalls,mixed.metrics.recoveryBindingWork.facetBindingRowsVisited,mixed.metrics.recoveryBindingWork.facetBindingRequestRowsVisited,mixed.metrics.recoveryBindingWork.facetBindingDuplicateRowsCompared,mixed.metrics.recoveryBindingWork.facetBindingReceiptBuilds,mixed.metrics.recoveryBindingWork.facetBindingReceiptReuses,mixed.metrics.recoveryBindingWork.facetBindingReceiptRowsCopied,mixed.metrics.recoveryBindingWork.facetBindingReceiptLookups,mixed.metrics.recoveryBindingWork.facetBindingReceiptReserveCalls,mixed.metrics.recoveryBindingWork.edgeDataCopyRows,mixed.metrics.recoveryBindingWork.faceDataCopyRows,mixed.metrics.recoveryBindingWork.faceIncidenceRowsCopied,mixed.metrics.recoveryBindingWork.stateCopies,mixed.metrics.recoveryBindingWork.stateCopyRows});
  if(recovery.maximumSteinerCandidates)receipt.steinerNodes=workspace.nodes().activeBySequence().size()-construction.nodeCount();
  coreNeed(recovered.reason!=bnd::TetFacetSetRecoveryReason::AllocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);
  if(!recovered.processed.empty()){const auto& row=recovered.processed.back();receipt.edgeReason=int(row.edgeOutcome.reason);receipt.facetReason=int(row.facetOutcome.reason);
    coreNeed(row.edgeOutcome.reason!=bnd::TetEdgeRecoveryReason::AllocationFailed&&row.facetOutcome.reason!=bnd::TetFacetRecoveryReason::AllocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);
    if(recovered.status!=bnd::TetFacetSetRecoveryStatus::Complete&&row.kind==bnd::TetFacetSetWorkKind::Edge){receipt.failedEdge=row.edgeKey.canonicalGlobalNodes;receipt.failedEdgeAvailable=true;}
  }
    
coreNeed(recovered.status==bnd::TetFacetSetRecoveryStatus::Complete,"recovery-refused",receipt.edgeReason==int(bnd::TetEdgeRecoveryReason::ResourceBudgetReached)?OcctPlcRecoveryReason::BudgetExceeded:OcctPlcRecoveryReason::UnsupportedGeometry);
  // Exact side seeds and immutable recovered-facet barriers, as in accepted V3 extraction.
  constexpr int unknown=-2;std::map<topo::CellHandle,int> labels;std::set<topo::FaceHandle> barriers;
  const auto seed=[&](topo::CellHandle cell,int value){auto& old=labels.at(cell);coreNeed(old==unknown||old==value,"side-label-conflict");old=value;};
  for(const auto& [cell,record]:context.cells){labels.emplace(cell,unknown);for(auto id:record.nodes)
   if(const auto known=global.find(context.denseToNode[id]);known!=global.end()&&std::size_t(known->second)>=std::size_t(volumeNodes)){seed(cell,-1);break;}}
  for(const auto& facet:artifact.facets){const auto binding=state.bindings.find(facet.key);coreNeed(binding!=state.bindings.end(),"missing-binding");
   const auto h=binding->second.recoveredFace;coreNeed(barriers.insert(h).second,"duplicate-barrier");const auto f=workspace.faces().get(h);coreNeed(f.incident.size()==2,"barrier-incidence");
   std::set<topo::NodeHandle> corners;for(auto id:facet.directedGlobalNodes)corners.insert(local.at(id));
   for(const auto& inc:f.incident){const topo::CellHandle cell{inc.cellSlot,inc.cellGen,workspace.cells().poolId()};const auto data=workspace.cells().get(cell);
    bool found=false;for(std::size_t j=0;j<4;++j)if(!corners.contains(data.nodes[j])){coreNeed(!found,"opposite-vertex-count");found=true;
     const auto sign=predicates::orient3d(surface.node(facet.directedGlobalNodes[0]),surface.node(facet.directedGlobalNodes[1]),surface.node(facet.directedGlobalNodes[2]),workspace.nodes().get(data.nodes[j]).coord);
     coreNeed(sign!=predicates::Sign::Zero,"zero-barrier-side");seed(cell,sign==predicates::Sign::Negative?owner.body:-1);}coreNeed(found,"missing-opposite-vertex");}}
  tick(context.cells.size()+artifact.facets.size()*16);std::queue<topo::CellHandle> queue;for(const auto& [cell,label]:labels)if(label!=unknown)queue.push(cell);
  constexpr std::array<std::array<int,3>,4> tf{{{{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};
  while(!queue.empty()){tick(5);const auto cell=queue.front();queue.pop();const auto& record=context.cells.at(cell);
   for(const auto& v:tf){const auto key=bnd::DenseTetFaceKey::of({record.nodes[v[0]],record.nodes[v[1]],record.nodes[v[2]]});const auto& row=context.faces.at(key);
    if(barriers.contains(row.handle)||row.owners.size()!=2)continue;const auto next=row.owners[0].cell==cell?row.owners[1].cell:row.owners[0].cell;
    auto& label=labels.at(next);if(label==unknown){label=labels.at(cell);queue.push(next);}else coreNeed(label==labels.at(cell),"flood-label-conflict");}}
  mesh::Mesh result=original;
  for(Index i=static_cast<Index>(nodes);i<volumeNodes;++i){tick();result.addNode(volumePoints.node(i),volumePoints.nodeClass(i));}
  std::set<Index> used;std::map<topo::NodeHandle,Index> generated;
  for(const auto& [cell,label]:labels){tick(5);coreNeed(label!=unknown,"unseeded-component");if(label<0)continue;const auto data=workspace.cells().get(cell);std::array<Index,4> ids;
   for(unsigned j=0;j<4;++j){const auto known=global.find(data.nodes[j]);
    if(known==global.end()){auto published=generated.find(data.nodes[j]);if(published==generated.end()){coreNeed(std::uint64_t(result.nodeCount())<options.maximumNodes,"recovery-node-budget",OcctPlcRecoveryReason::BudgetExceeded);const auto id=result.addNode(workspace.nodes().get(data.nodes[j]).coord,{3,bodyTag});published=generated.emplace(data.nodes[j],id).first;}ids[j]=published->second;}
    else {ids[j]=known->second;coreNeed(ids[j]>=0&&ids[j]<volumeNodes,"guard-retained");}used.insert(ids[j]);}
   coreNeed(predicates::orient3d(result.node(ids[0]),result.node(ids[1]),result.node(ids[2]),result.node(ids[3]))==predicates::Sign::Positive,"inverted-tet");
   result.addCell(mesh::CellType::Tet4,ids,{3,owner.body});++cells;}
  
  coreNeed(cells>0,"empty-recovered-volume");for(Index i=0;i<volumeNodes;++i)coreNeed(used.contains(i),"missing-original-or-seed-point");
  tick(result.connectivityCount()+result.nodeCount());auto auditLimits=options.auditLimits;auditLimits.cancellation=options.cancellation;const auto audit=mesh::auditLinearVolumeBoundary(result,auditLimits);coreNeed(!audit.allocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);coreNeed(!audit.cancelled,"cancelled",OcctPlcRecoveryReason::Cancelled);coreNeed(!audit.budgetExceeded,"core-phase-work-budget",OcctPlcRecoveryReason::BudgetExceeded);coreNeed(audit.valid,"exact-boundary-incidence");
  tick(result.cellCount()+result.nodeCount());const auto q=algo::analyzeMixedMeshQuality(result);coreNeed(!q.allocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);coreNeed(q.valid,"quality");const double volume=q.totalVolume;
  coreNeed(std::abs(volume-plc)<=plc*1e-10,"plc-volume-difference");

  // Prefix is byte/value-identical by construction; verify to make certificate explicit.
  for(Index i=0;i<Index(nodes);++i){const auto a=result.node(i),b=original.node(i);coreNeed(a.x==b.x&&a.y==b.y&&a.z==b.z&&result.nodeClass(i)==original.nodeClass(i),"original-node-changed");}
  for(Index i=0;i<Index(triangles);++i){coreNeed(result.cell(i).onEntity==original.cell(i).onEntity,"native-face-class-changed");for(unsigned j=0;j<3;++j)coreNeed(result.cellNodes(i)[j]==original.cellNodes(i)[j],"original-directed-facet-changed");}
  coreNeed(result.nodeCount()<=options.maximumNodes&&result.cellCount()<=options.maximumCells,"output-budget",OcctPlcRecoveryReason::BudgetExceeded);
  coreNeed(std::uint64_t(result.nodeCount())*(sizeof(Vec3)+sizeof(mesh::Classification))+std::uint64_t(result.cellCount())*sizeof(mesh::Cell)+std::uint64_t(result.connectivityCount())*sizeof(Index)<=options.maximumOutputBytes,"output-budget",OcctPlcRecoveryReason::BudgetExceeded);
  tick();destination.swap(result);complete=true;why="none";code=OcctPlcRecoveryReason::None;

 }catch(const CoreStop& s){why=s.reason;code=s.code;}catch(const Standard_OutOfMemory&){why="allocation-failed";code=OcctPlcRecoveryReason::AllocationFailed;}catch(const std::bad_alloc&){why="allocation-failed";code=OcctPlcRecoveryReason::AllocationFailed;}catch(const Standard_Failure&){why="native-query-failed";code=OcctPlcRecoveryReason::NativeQueryFailed;}catch(...){why="internal-error";}
 receipt.published=complete;receipt.diagnostic=why;receipt.reason=code;return receipt;
}
} // namespace
} // namespace paacs::geo::detail

namespace paacs::geo::detail {
namespace {
struct Stop{const char* why;OcctPlcRecoveryReason reason=OcctPlcRecoveryReason::UnsupportedGeometry;};
using Tri=std::array<Index,3>;using Edge=std::array<Index,2>;
Edge edge(Index a,Index b){if(b<a)std::swap(a,b);return {a,b};}
bool same(Vec3 a,Vec3 b){return a.x==b.x&&a.y==b.y&&a.z==b.z;}
}
OcctPlcRecoveryOutcome recoverOcctPlcVolume(mesh::Mesh& source,int pieceTag,
 const OcctModel& pieceModel,const OcctModel& original,int originalBody,
 const std::vector<mesh::Classification>& parentByFace,const OcctStarOptions& options,
 const OcctPartitionRecoveryOptions& recovery,mesh::Mesh& destination) noexcept{
 if(recovery.enabled&&recovery.useLiveConstraints)
  return recoverOcctPlcLiveVolume(source,pieceTag,options,recovery,destination);
 if(recovery.enabled&&recovery.useAdvancingFront)
  return recoverOcctPlcFrontVolume(source,pieceTag,pieceModel,original,originalBody,options,recovery,destination);
 OcctPlcRecoveryOutcome out;
 const auto collect=[&](std::uint64_t count){if(count>UINT64_MAX-out.workUsed)throw Stop{"work-overflow",OcctPlcRecoveryReason::BudgetExceeded};out.workUsed+=count;if(out.workUsed>options.maximumWork)throw Stop{"core-phase-work-budget",OcctPlcRecoveryReason::BudgetExceeded};if(options.cancellation&&options.cancellation->load())throw Stop{"cancelled",OcctPlcRecoveryReason::Cancelled};};
 const auto tick=[&](std::uint64_t count=1){if(options.cancellation&&options.cancellation->load())throw Stop{"cancelled",OcctPlcRecoveryReason::Cancelled};collect(count);};
 const auto need=[](bool yes,const char* why,OcctPlcRecoveryReason reason=OcctPlcRecoveryReason::UnsupportedGeometry){if(!yes)throw Stop{why,reason};};
 try{
  need(recovery.enabled&&recovery.maximumCoreAttempts>0&&recovery.maximumCoreAttempts<=16&&recovery.maximumQuadFlipsPerAttempt>0&&recovery.maximumQuadFlipsPerAttempt<=64&&recovery.maximumGeneralWorkPerAttempt>0&&recovery.maximumSteinerCandidates<=64,"invalid-recovery-options",OcctPlcRecoveryReason::InvalidRequest);
  out.attempted=true;out.piecesAttempted=1;out.failedPieceTag=pieceTag;
  tick();need(pieceTag>0&&originalBody>0&&pieceModel.solid(pieceTag)&&original.solid(originalBody),"native-body-missing",OcctPlcRecoveryReason::InvalidRequest);
  need(source.nodeCount()>=4&&source.cellCount()>=4,"empty-or-incomplete-boundary",OcctPlcRecoveryReason::InvalidRequest);
  need(source.nodeCount()<=10000&&source.cellCount()<=20000,"core-count-bound",OcctPlcRecoveryReason::BudgetExceeded);
  // Phase payload admission is conservative and excludes OCCT heap storage.
  need(std::uint64_t(source.nodeCount()+source.cellCount()+1)*32768<=options.maximumWorkingBytes,"core-working-admission",OcctPlcRecoveryReason::BudgetExceeded);
  mesh::Mesh boundary=source;std::set<std::vector<Tri>> visited;
  paacs::geo::detail::OcctSolidClassifier pieceClassifier(OcctModelAccess::solid(pieceModel,pieceTag));
  paacs::geo::detail::OcctSolidClassifier originalClassifier(OcctModelAccess::solid(original,originalBody));
  const auto pieceTol=pieceModel.tolerance(3,pieceTag),originalTol=original.tolerance(3,originalBody);
  const auto inside=[&](Vec3 p){tick(2);out.nativeQueries+=2;pieceClassifier.Perform(gp_Pnt(p.x,p.y,p.z),pieceTol);originalClassifier.Perform(gp_Pnt(p.x,p.y,p.z),originalTol);return (pieceClassifier.State()==TopAbs_IN||pieceClassifier.State()==TopAbs_ON)&&(originalClassifier.State()==TopAbs_IN||originalClassifier.State()==TopAbs_ON);};
  for(unsigned attempt=0;attempt<recovery.maximumCoreAttempts;++attempt){
   std::vector<Tri> signature;signature.reserve(boundary.cellCount());
   for(Index i=0;i<boundary.cellCount();++i){tick();need(boundary.cell(i).type==mesh::CellType::Tri3&&boundary.cell(i).onEntity.dim==2,"core-boundary-kind");const auto* p=boundary.cellNodes(i);signature.push_back({p[0],p[1],p[2]});}
   need(visited.insert(signature).second,"core-quad-cycle");
   auto limited=options;limited.maximumWork=options.maximumWork-out.workUsed;
   const auto certificate=certifyOcctBoundarySurface(boundary,limited);collect(certificate.workUsed);
   if(!certificate.certified){
    auto reason=OcctPlcRecoveryReason::UnsupportedGeometry;
    if(certificate.reason==OcctStarReason::WorkBudget||certificate.reason==OcctStarReason::OutputBudget)reason=OcctPlcRecoveryReason::BudgetExceeded;
    else if(certificate.reason==OcctStarReason::Cancelled)reason=OcctPlcRecoveryReason::Cancelled;
    else if(certificate.reason==OcctStarReason::AllocationFailed)reason=OcctPlcRecoveryReason::AllocationFailed;
    else if(certificate.reason==OcctStarReason::InvalidBoundary)reason=OcctPlcRecoveryReason::InvalidMesh;
    throw Stop{certificate.diagnostic,reason};
   }
   limited.maximumWork=options.maximumWork-out.workUsed;mesh::Mesh result;
   const auto core=buildFixedPlcCore(boundary,pieceTag,limited,result,&pieceModel,recovery);++out.coreAttempts;out.generalWorkReserved+=core.generalReservedWork;out.nativeQueries+=core.nativeSeedQueries;out.steinerNodes+=core.steinerNodes;collect(core.workUsed);
   
   std::vector<std::array<Index,2>> exteriorPairs;
   if(core.published){need(result.nodeCount()>=boundary.nodeCount(),"core-node-prefix-size");
   for(Index i=0;i<boundary.nodeCount();++i){tick();need(same(result.node(i),source.node(i))&&result.nodeClass(i)==source.nodeClass(i),"core-original-node-prefix-changed");}
   for(Index i=0;i<boundary.cellCount();++i){tick();const auto* a=boundary.cellNodes(i);const auto* b=result.cellNodes(i);need(result.cell(i).type==mesh::CellType::Tri3&&result.cell(i).onEntity==boundary.cell(i).onEntity&&a[0]==b[0]&&a[1]==b[1]&&a[2]==b[2],"core-directed-facet-prefix-changed");}
   // Immutable native boundary nodes retain their source classification. The
   // outer consumer validates each against its owning Original entity; a
   // boundary solid-classifier result is not an entity projection certificate.
   // Newly introduced interior nodes require BOTH native solid inclusions.
   for(Index i=boundary.nodeCount();i<result.nodeCount();++i)need(inside(result.node(i)),"core-native-interior-node-outside");
   std::vector<Index> rejected;
   for(Index i=0;i<result.cellCount();++i)if(result.cell(i).type==mesh::CellType::Tet4){const auto* p=result.cellNodes(i);if(!inside((result.node(p[0])+result.node(p[1])+result.node(p[2])+result.node(p[3]))*.25))rejected.push_back(i);}
   if(rejected.empty()){
    // Full whole-original projected facet/chord/volume/native checks remain in
    // the outer partition consumer, after exact interface assembly.
    auto auditLimits=options.auditLimits;auditLimits.cancellation=options.cancellation;
    const auto audit=mesh::auditLinearVolumeBoundary(result,auditLimits);need(!audit.allocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);need(!audit.cancelled,"cancelled",OcctPlcRecoveryReason::Cancelled);need(!audit.budgetExceeded,"core-phase-work-budget",OcctPlcRecoveryReason::BudgetExceeded);need(audit.valid,"core-final-boundary-audit");
    tick(result.cellCount());const auto quality=algo::analyzeMixedMeshQuality(result);need(!quality.allocationFailed,"allocation-failed",OcctPlcRecoveryReason::AllocationFailed);need(quality.valid,"core-final-positive-quality");
    need(std::abs(quality.totalVolume-certificate.signedVolume)<=quality.totalVolume*1e-10,"core-final-plc-volume");
    tick(0);source.swap(boundary);destination.swap(result);out.published=true;out.piecesRecovered=1;out.failedPieceTag=-1;out.reason=OcctPlcRecoveryReason::None;out.diagnostic="none";return out;
   }
   for(const auto cell:rejected){const auto* ids=result.cellNodes(cell);std::vector<Index> exterior;
    for(Index f=0;f<boundary.cellCount();++f){tick();const auto* t=boundary.cellNodes(f);if(std::find(ids,ids+4,t[0])!=ids+4&&std::find(ids,ids+4,t[1])!=ids+4&&std::find(ids,ids+4,t[2])!=ids+4)exterior.push_back(f);}
    if(exterior.size()==2)exteriorPairs.push_back({exterior[0],exterior[1]});
   }
   }else{
    need(core.failedEdgeAvailable&&(core.edgeReason==int(bnd::TetEdgeRecoveryReason::NonCycleLink)||core.edgeReason==int(bnd::TetEdgeRecoveryReason::RecoveryConflictCycle)),core.diagnostic,core.reason);
    std::vector<Index> exterior;
    for(Index f=0;f<boundary.cellCount();++f){tick();const auto* t=boundary.cellNodes(f);if(std::find(t,t+3,core.failedEdge[0])!=t+3&&std::find(t,t+3,core.failedEdge[1])!=t+3)exterior.push_back(f);}
    need(exterior.size()==2,"core-refused-edge-is-not-surface-diagonal");exteriorPairs.push_back({exterior[0],exterior[1]});
    
   }
   std::vector<Tri> proposed=signature;std::set<Index> changed;std::set<Edge> proposedDiagonals;std::map<Edge,unsigned> allEdges;
   std::map<Edge,std::set<int>> nativeOwners;
   for(Index i=0;i<boundary.cellCount();++i)for(unsigned j=0;j<3;++j){tick();const auto key=edge(signature[i][j],signature[i][(j+1)%3]);
    ++allEdges[key];nativeOwners[key].insert(boundary.cell(i).onEntity.tag);}
   const auto nativeSegment=[&](Index a,Index b,int face){const auto found=nativeOwners.find(edge(a,b));
    if(found!=nativeOwners.end())return actualOcctWire(pieceModel,boundary.nodeClass(a),boundary.nodeClass(b),found->second,[&]{tick();})!=0;
    const std::array faces{face};return actualOcctWire(pieceModel,boundary.nodeClass(a),boundary.nodeClass(b),faces,[&]{tick();})!=0;};
   unsigned flipped=0;
   for(const auto& exterior:exteriorPairs){if(flipped==recovery.maximumQuadFlipsPerAttempt)break;
    if(changed.contains(exterior[0])||changed.contains(exterior[1]))continue;
    const auto cls=boundary.cell(exterior[0]).onEntity;if(cls!=boundary.cell(exterior[1]).onEntity||cls.tag<=0||static_cast<std::size_t>(cls.tag)>=parentByFace.size()||parentByFace[cls.tag].dim!=2)continue;
    const auto a=signature[exterior[0]],b=signature[exterior[1]];Index u=-1,v=-1;
    for(unsigned j=0;j<3;++j)if(std::find(b.begin(),b.end(),a[j])!=b.end()&&std::find(b.begin(),b.end(),a[(j+1)%3])!=b.end()){u=a[j];v=a[(j+1)%3];break;}
    if(u<0||allEdges[edge(u,v)]!=2)continue;
    if(nativeSegment(u,v,cls.tag)){continue;}
    Index c=-1,d=-1;for(auto x:a)if(x!=u&&x!=v)c=x;for(auto x:b)if(x!=u&&x!=v)d=x;
    if(c<0||d<0||c==d||allEdges.contains(edge(c,d))||proposedDiagonals.contains(edge(c,d))){continue;}
    bool opposite=false;for(unsigned j=0;j<3;++j)opposite|=b[j]==v&&b[(j+1)%3]==u;need(opposite,"core-quad-winding");
    proposed[exterior[0]]={c,u,d};proposed[exterior[1]]={c,d,v};changed.insert(exterior[0]);changed.insert(exterior[1]);proposedDiagonals.insert(edge(c,d));++flipped;
   }
   
   need(flipped>0,"core-no-eligible-exterior-quad");
   mesh::Mesh next;next.reserveAppend(boundary.nodeCount(),boundary.cellCount(),boundary.connectivityCount());
   for(Index i=0;i<boundary.nodeCount();++i){tick();next.addNode(boundary.node(i),boundary.nodeClass(i));}
   std::set<Edge> beforeNative,afterNative;
   for(Index i=0;i<boundary.cellCount();++i){tick();const auto cls=boundary.cell(i).onEntity;next.addCell(mesh::CellType::Tri3,proposed[i],cls);
    if(changed.contains(i))need(parentByFace[cls.tag].dim==2,"core-interface-retriangulated");
    for(unsigned j=0;j<3;++j){const auto a=signature[i][j],b=signature[i][(j+1)%3];if(nativeSegment(a,b,cls.tag))beforeNative.insert(edge(a,b));const auto c=proposed[i][j],d=proposed[i][(j+1)%3];if(nativeSegment(c,d,cls.tag))afterNative.insert(edge(c,d));}}
   
   need(beforeNative==afterNative,"core-native-curve-segments-changed");boundary.swap(next);out.quadFlips+=flipped;
  }
  throw Stop{"core-quad-attempt-bound"};
 }catch(const Stop& stop){out.diagnostic=stop.why;out.reason=stop.reason;}
 catch(const Standard_OutOfMemory&){out.diagnostic="allocation-failed";out.reason=OcctPlcRecoveryReason::AllocationFailed;}
 catch(const std::bad_alloc&){out.diagnostic="allocation-failed";out.reason=OcctPlcRecoveryReason::AllocationFailed;}
 catch(const Standard_Failure&){out.diagnostic="native-query-failed";out.reason=OcctPlcRecoveryReason::NativeQueryFailed;}
 catch(...){out.diagnostic="core-internal-error";out.reason=OcctPlcRecoveryReason::InternalError;}
 return out;
}
}
