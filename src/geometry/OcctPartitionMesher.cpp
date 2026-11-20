// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctPartitionMesher.hpp"
#include "OcctStageTrace.hpp"
#include "OcctBoundaryInternal.hpp"
#include "OcctBoundaryDiagonalClosure.hpp"
#include "OcctModelInternal.hpp"
#include "OcctPlcMaterialFilter.hpp"
#include "OcctTetQuality.hpp"
#include "paacs/bnd/TetPlcCavityRecovery.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "OcctSolidClassifier.hpp"
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <Standard_OutOfMemory.hxx>
#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <map>
#include <set>
#include <string_view>

namespace paacs::geo::detail {
namespace {
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point start){return std::chrono::duration<double>(Clock::now()-start).count();}
struct Timer {double& seconds;Clock::time_point start=Clock::now();~Timer(){seconds+=elapsed(start);}};
struct Stop {OcctPartitionMeshReason reason;const char* diagnostic;};
void need(bool ok,OcctPartitionMeshReason reason,const char* diagnostic){if(!ok)throw Stop{reason,diagnostic};}
bool finite(Vec3 p){return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);}
bool identical(Vec3 a,Vec3 b){return std::bit_cast<std::uint64_t>(a.x)==std::bit_cast<std::uint64_t>(b.x)&&std::bit_cast<std::uint64_t>(a.y)==std::bit_cast<std::uint64_t>(b.y)&&std::bit_cast<std::uint64_t>(a.z)==std::bit_cast<std::uint64_t>(b.z);}
using Tri=std::array<Index,3>;
Tri sorted(Tri t){std::sort(t.begin(),t.end());return t;}
bool sameWinding(Tri a,Tri b){for(unsigned i=0;i<3;++i)if(a[0]==b[i]&&a[1]==b[(i+1)%3]&&a[2]==b[(i+2)%3])return true;return false;}
std::uint64_t add(std::uint64_t a,std::uint64_t b){need(b<=UINT64_MAX-a,OcctPartitionMeshReason::BudgetExceeded,"payload-overflow");return a+b;}
std::uint64_t multiply(std::uint64_t a,std::uint64_t b){need(b==0||a<=UINT64_MAX/b,OcctPartitionMeshReason::BudgetExceeded,"payload-overflow");return a*b;}
std::uint64_t payload(std::uint64_t n,std::uint64_t c,std::uint64_t connectivity){return add(add(multiply(n,sizeof(Vec3)+sizeof(mesh::Classification)),multiply(c,sizeof(mesh::Cell))),multiply(connectivity,sizeof(Index)));}
struct Context {
    const OcctPartitionMeshOptions& options;OcctPartitionMeshOutcome& out;
    std::uint64_t remaining()const{return options.maximumWork-out.workUsed;}
    void charge(std::uint64_t n,std::uint64_t& phase){
        need(!options.cancellation||!options.cancellation->load(std::memory_order_relaxed),OcctPartitionMeshReason::Cancelled,"cancelled");
        need(n<=remaining(),OcctPartitionMeshReason::BudgetExceeded,"cumulative-work-budget");out.workUsed+=n;phase+=n;
    }
    void collect(std::uint64_t n,std::uint64_t& phase){
        // A completed child attempt already consumed these visits, including
        // work reported by a cancelled attempt. Preserve its receipt first.
        need(n<=remaining(),OcctPartitionMeshReason::BudgetExceeded,"cumulative-work-budget");out.workUsed+=n;phase+=n;
        charge(0,phase);
    }
    void output(std::uint64_t n,std::uint64_t c,std::uint64_t cn){need(n<=options.maximumNodes&&c<=options.maximumCells&&payload(n,c,cn)<=options.maximumOutputBytes,OcctPartitionMeshReason::BudgetExceeded,"output-budget");}
    void working(std::uint64_t bytes){out.maximumWorkingPayloadEstimate=std::max(out.maximumWorkingPayloadEstimate,bytes);need(bytes<=options.maximumWorkingBytes,OcctPartitionMeshReason::BudgetExceeded,"working-payload-admission");}
};
void boundaryFailure(OcctBoundaryReason reason){
    if(reason==OcctBoundaryReason::BudgetExceeded)throw Stop{OcctPartitionMeshReason::BudgetExceeded,"surface-budget"};
    if(reason==OcctBoundaryReason::Cancelled)throw Stop{OcctPartitionMeshReason::Cancelled,"surface-cancelled"};
    if(reason==OcctBoundaryReason::AllocationFailed)throw Stop{OcctPartitionMeshReason::AllocationFailed,"surface-allocation-failed"};
    throw Stop{OcctPartitionMeshReason::SurfaceFailed,occtBoundaryReasonName(reason)};
}
void starFailure(OcctStarReason reason,const char* diagnostic){
    if(reason==OcctStarReason::WorkBudget||reason==OcctStarReason::OutputBudget)throw Stop{OcctPartitionMeshReason::BudgetExceeded,diagnostic};
    if(reason==OcctStarReason::Cancelled)throw Stop{OcctPartitionMeshReason::Cancelled,diagnostic};
    if(reason==OcctStarReason::AllocationFailed)throw Stop{OcctPartitionMeshReason::AllocationFailed,diagnostic};
    if(reason==OcctStarReason::NativeQueryFailed)throw Stop{OcctPartitionMeshReason::GeometryError,diagnostic};
    throw Stop{OcctPartitionMeshReason::StarFailed,diagnostic};
}
void nativeValidationFailure(OcctNativeValidationReason reason){
    using R=OcctNativeValidationReason;
    const auto* diagnostic=occtNativeValidationReasonName(reason);
    if(reason==R::BudgetExceeded)throw Stop{OcctPartitionMeshReason::BudgetExceeded,diagnostic};
    if(reason==R::Cancelled)throw Stop{OcctPartitionMeshReason::Cancelled,diagnostic};
    if(reason==R::AllocationFailed)throw Stop{OcctPartitionMeshReason::AllocationFailed,diagnostic};
    if(reason==R::OutsideMaterial)throw Stop{OcctPartitionMeshReason::GeometryError,"native-tet-centroid-outside"};
    if(reason==R::CopyGeometryMismatch||reason==R::ReferenceStateMismatch||reason==R::SourceChanged)
        throw Stop{OcctPartitionMeshReason::ProvenanceLost,diagnostic};
    if(reason==R::NativeQueryFailed)throw Stop{OcctPartitionMeshReason::GeometryError,diagnostic};
    if(reason==R::InvalidRequest)throw Stop{OcctPartitionMeshReason::InvalidRequest,diagnostic};
    throw Stop{OcctPartitionMeshReason::InternalError,diagnostic};
}
mesh::Classification parent(const OcctPartitionOutcome& p,mesh::Classification c){
    const std::vector<mesh::Classification>* map=c.dim==0?&p.parentByVertex:c.dim==1?&p.parentByEdge:c.dim==2?&p.parentByFace:nullptr;
    need(map&&c.tag>0&&static_cast<std::size_t>(c.tag)<map->size(),OcctPartitionMeshReason::ProvenanceLost,"derived-entity-has-no-parent");
    const auto cls=(*map)[static_cast<std::size_t>(c.tag)];
    need(cls.dim>=0&&cls.dim<=3&&cls.tag>=0,OcctPartitionMeshReason::ProvenanceLost,"invalid-original-parent");return cls;
}
struct FaceUse {Tri ids;int piece=-1,derivedFace=-1;mesh::Classification original;};
struct ExteriorDiagonalWorkLimit {};

template<class Admit>
bool repairPartitionMaterial(const OcctPartitionOutcome& partition,int tag,mesh::Mesh& boundary,
    const OcctStarOptions& options,Context& context,Admit&& originalAdmits,mesh::Mesh& volume,
    bool allowAdaptiveRefill) {
    auto& out=context.out;
    OcctStageTrace materialTrace("material-filter",tag);
    mesh::Mesh retained;
    const auto charge=[&](std::uint64_t count){context.charge(count,out.validationWork);};
    const auto material=retainOcctPlcMaterial(volume,boundary.nodeCount(),originalAdmits,charge,retained);
    if(!material.removedInteriorNodes&&!material.removedIncidentCells&&!material.removedCentroidCells)return false;

    const auto volumePayload=payload(volume.nodeCount(),volume.cellCount(),volume.connectivityCount());
    const auto cavityLimits=[&](const mesh::Mesh& shell,std::uint64_t work,std::uint64_t workingBytes){
        bnd::TetPlcCavityOptions limits;
        limits.region={3,tag};
        limits.maximumNodes=options.maximumNodes;
        limits.maximumCells=options.maximumCells-std::uint64_t(shell.cellCount());
        limits.maximumSourceTriangles=4096;
        limits.maximumWorkingBytes=workingBytes;
        limits.maximumWork=work;
        limits.maximumComponentFaces=4096;
        limits.maximumPatchCells=2048;
        limits.maximumExpansionCells=256;
        limits.maximumCandidatePoints=64;
        limits.maximumRandomCandidatePoints=24;
        limits.maximumMergedComponents=4;
        limits.maximumPatches=256;
        limits.maximumLocalFrontTets=256;
        limits.patchSearch={8192,32768,256,4};
        limits.interiorAdmission=originalAdmits;
        limits.centroidAdmission=originalAdmits;
        return limits;
    };
    const auto tryExteriorDiagonal=[&](){
        const auto workCap=std::min<std::uint64_t>(5'000'000,context.remaining()/4);
        const auto start=context.out.workUsed;
        if(!material.removedCentroidCells||workCap<50'000)return false;
        const auto working=add(add(volumePayload,payload(retained.nodeCount(),retained.cellCount(),
            retained.connectivityCount())),multiply(add(boundary.nodeCount(),boundary.cellCount()),2048));
        if(working>=options.maximumWorkingBytes)return false;
        const auto workLeft=[&](){return workCap-std::min(workCap,context.out.workUsed-start);};
        const auto trialCharge=[&](std::uint64_t n){
            if(n>workLeft())throw ExteriorDiagonalWorkLimit{};
            context.charge(n,out.starWork);
        };
        const auto originalExterior=[&](mesh::Classification face){
            return parent(partition,face).dim==2&&partition.model->faceOwners(face.tag).size()==1;
        };
        OcctStageTrace stage("exterior_material_diagonal",tag);
        try {
            std::map<Tri,Index> exterior;
            for(Index i=0;i<boundary.cellCount();++i){
                trialCharge(1);
                if(!originalExterior(boundary.cell(i).onEntity))continue;
                const auto* p=boundary.cellNodes(i);
                exterior.emplace(sorted({p[0],p[1],p[2]}),i);
            }
            using Pair=std::array<Index,2>;
            std::set<Pair> rejectedPairs;
            for(Index i=0;i<volume.cellCount();++i){
                trialCharge(1);
                if(volume.cell(i).type!=mesh::CellType::Tet4)continue;
                const auto* p=volume.cellNodes(i);
                std::array<Index,4> owners{};
                unsigned count=0;
                for(unsigned side=0;side<4;++side){
                    const auto at=exterior.find(sorted({p[(side+1)%4],p[(side+2)%4],p[(side+3)%4]}));
                    if(at!=exterior.end())owners[count++]=at->second;
                }
                if(count!=2||boundary.cell(owners[0]).onEntity!=boundary.cell(owners[1]).onEntity)continue;
                Vec3 centroid{};
                for(unsigned k=0;k<4;++k)centroid=centroid+volume.node(p[k])*.25;
                if(!originalAdmits(centroid))
                    rejectedPairs.insert({std::min(owners[0],owners[1]),std::max(owners[0],owners[1])});
            }
            if(rejectedPairs.empty())return false;
            // A rejected exterior sliver can expose the opposite pair of
            // triangles on the same CAD face. Prefer that existing tetrahedral
            // boundary when the diagonal closure certifies both exact owners.
            mesh::Mesh staged;
            const auto diagonal=closeOcctBoundaryDiagonals(*partition.model,boundary,retained,64,staged,
                trialCharge,[&](Index a,Index b,mesh::Classification){
                    return rejectedPairs.contains({std::min(a,b),std::max(a,b)});
                });
            if(!diagonal.flips)return false;
            for(Index i=0;i<staged.cellCount();++i){
                trialCharge(1);
                const auto* p=staged.cellNodes(i);
                if(std::equal(p,p+3,boundary.cellNodes(i)))continue;
                const auto face=staged.cell(i).onEntity;
                if(!originalExterior(face))return false;
                const auto a=staged.node(p[0]),b=staged.node(p[1]),c=staged.node(p[2]);
                const auto chord=context.options.boundary.maximumChordError+partition.model->tolerance(2,face.tag);
                for(const auto sample:std::array{(a+b+c)/3.,(a+b)*.5,(b+c)*.5,(c+a)*.5}){
                    trialCharge(1);
                    Vec3 closest;
                    if(!closestPointOnOcctFace(*partition.model,face.tag,sample,closest)||
                       !((sample-closest).norm()<=chord))return false;
                }
            }
            // All old nodes and cut interfaces are identical. Re-certify the
            // revised shell and admit a volume transaction before using it.
            auto certificateLimits=options;
            certificateLimits.maximumWork=workLeft();
            certificateLimits.maximumWorkingBytes-=working;
            const auto certificate=certifyOcctBoundarySurface(staged,certificateLimits);
            trialCharge(certificate.workUsed);
            out.exactSurfacePairs=add(out.exactSurfacePairs,certificate.exactBoundaryPairs);
            if(certificate.reason==OcctStarReason::Cancelled||certificate.reason==OcctStarReason::AllocationFailed)
                starFailure(certificate.reason,certificate.diagnostic);
            if(!certificate.certified||workLeft()<1000)return false;
            const auto solverWork=workLeft()/2;
            auto admissionWork=workLeft()-solverWork;
            bool admissionWorkExhausted=false;
            Stop admissionFailure{OcctPartitionMeshReason::None,nullptr};
            const auto admitWithinTrial=[&](Vec3 point){
                if(!admissionWork){admissionWorkExhausted=true;throw ExteriorDiagonalWorkLimit{};}
                --admissionWork;
                try{return originalAdmits(point);}
                catch(const Stop& failure){admissionFailure=failure;throw;}
                catch(const Standard_OutOfMemory&){throw std::bad_alloc{};}
            };
            auto limits=cavityLimits(staged,solverWork,options.maximumWorkingBytes-working);
            limits.interiorAdmission=admitWithinTrial;
            limits.centroidAdmission=admitWithinTrial;
            mesh::Mesh closed;
            bnd::TetPlcCavityOutcome refill;
            {Timer timer{out.starSeconds};refill=bnd::recoverTetPlcCavities(staged,retained,limits,closed);}
            context.collect(refill.workUsed,out.starWork);
            if(admissionFailure.reason!=OcctPartitionMeshReason::None)throw admissionFailure;
            if(refill.reason==bnd::TetPlcCavityReason::Cancelled)
                throw Stop{OcctPartitionMeshReason::Cancelled,refill.diagnostic};
            if(refill.reason==bnd::TetPlcCavityReason::AllocationFailed)
                throw Stop{OcctPartitionMeshReason::AllocationFailed,refill.diagnostic};
            if(admissionWorkExhausted)throw ExteriorDiagonalWorkLimit{};
            if(!refill.published){stage.finish(refill.diagnostic);return false;}
            if(context.out.workUsed-start>workCap)throw ExteriorDiagonalWorkLimit{};
            out.surfaceDiagonalRepairs=add(out.surfaceDiagonalRepairs,diagonal.flips);
            boundary.swap(staged);
            volume.swap(closed);
            stage.finish("complete");
            return true;
        }catch(const ExteriorDiagonalWorkLimit&){stage.finish("diagonal-work-cap");return false;}
    };
    if(allowAdaptiveRefill&&tryExteriorDiagonal())return true;

    // Prefer native-admitted interior points near the exposed CAD triangles.
    // Original curve nodes and shared-cut triangles are never changed here.
    OcctSolidClassifier pieceCheck(OcctModelAccess::solid(*partition.model,tag));
    const auto pieceTolerance=partition.model->tolerance(3,tag);
    std::set<Index> affectedFaces;
    for(Index ci=0;ci<volume.cellCount()&&affectedFaces.size()<16;++ci){
        charge(1);
        if(volume.cell(ci).type!=mesh::CellType::Tet4)continue;
        const auto* ids=volume.cellNodes(ci);
        const auto center=(volume.node(ids[0])+volume.node(ids[1])+volume.node(ids[2])+volume.node(ids[3]))*.25;
        if(originalAdmits(center))continue;
        double closest=(std::numeric_limits<double>::max)();
        Index face=-1;
        for(Index fi=0;fi<boundary.cellCount();++fi){
            charge(1);
            if(parent(partition,boundary.cell(fi).onEntity).dim!=2)continue;
            const auto* vertices=boundary.cellNodes(fi);
            const auto offset=(boundary.node(vertices[0])+boundary.node(vertices[1])+boundary.node(vertices[2]))/3.-center;
            const auto distance=offset.dot(offset);
            if(distance<closest){closest=distance;face=fi;}
        }
        if(face>=0)affectedFaces.insert(face);
    }
    for(const auto face:affectedFaces){
        charge(1);
        const auto* vertices=boundary.cellNodes(face);
        const auto a=boundary.node(vertices[0]),b=boundary.node(vertices[1]),c=boundary.node(vertices[2]);
        const auto normal=(b-a).cross(c-a);
        const auto length=normal.norm();
        if(!(length>0))continue;
        const auto edge=std::min({(b-a).norm(),(c-b).norm(),(a-c).norm()});
        for(const double scale:{.15,.45}){
            if(std::uint64_t(retained.nodeCount())>=options.maximumNodes)break;
            const auto point=(a+b+c)/3.-normal*(edge*scale/length);
            if(!originalAdmits(point))continue;
            charge(1);
            pieceCheck.Perform(gp_Pnt(point.x,point.y,point.z),pieceTolerance);
            if(pieceCheck.State()==TopAbs_IN)retained.addNode(point,{3,tag});
        }
    }
    need(volumePayload<options.maximumWorkingBytes,OcctPartitionMeshReason::BudgetExceeded,"local-recovery-working-budget");
    const auto limits=cavityLimits(boundary,std::min<std::uint64_t>(12'000'000,context.remaining()/2),
                                  options.maximumWorkingBytes-volumePayload);
    mesh::Mesh fixed;
    bnd::TetPlcCavityOutcome repair;
    materialTrace.finish();
    {OcctStageTrace stage("cavity_refill",tag);Timer timer{out.starSeconds};repair=bnd::recoverTetPlcCavities(boundary,retained,limits,fixed);}
    context.collect(repair.workUsed,out.starWork);
    if(repair.reason==bnd::TetPlcCavityReason::AllocationFailed)
        throw Stop{OcctPartitionMeshReason::AllocationFailed,repair.diagnostic};
    if(repair.reason==bnd::TetPlcCavityReason::BudgetExceeded)
        throw Stop{OcctPartitionMeshReason::BudgetExceeded,repair.diagnostic};
    need(repair.published,OcctPartitionMeshReason::StarFailed,repair.diagnostic);
    volume.swap(fixed);
    return true;
}
} // namespace

static OcctPartitionMeshOutcome buildOcctPartitionVolumeAttempt(const OcctModel& original,int solidTag,
    const OcctPartitionMeshOptions& options,mesh::Mesh& destination,
    const OcctPartitionRecoveryOptions& recovery,OcctPartitionRecoveryReceipt* receipt) noexcept {
    const auto started=Clock::now();OcctPartitionMeshOutcome out;Context context{options,out};
    try {
        const auto* solid=original.solid(solidTag);
        need(solid&&std::isfinite(options.maximumRelativeVolumeError)&&options.maximumRelativeVolumeError>0&&options.maximumRelativeVolumeError<1&&
             std::isfinite(options.boundary.maximumChordError)&&options.boundary.maximumChordError>0,
             OcctPartitionMeshReason::InvalidRequest,"invalid-request");
        need(options.maximumNativeValidationThreads>=1&&options.maximumNativeValidationThreads<=8,
             OcctPartitionMeshReason::InvalidRequest,"invalid-native-validation-threads");
        need(options.maximumWork&&options.maximumNodes&&options.maximumCells&&options.maximumOutputBytes&&options.maximumWorkingBytes,
             OcctPartitionMeshReason::BudgetExceeded,"zero-admission-budget");
        context.charge(0,out.assemblyWork);out.policy=options.partition.policy;out.surfaceChordLimit=options.boundary.maximumChordError;
        auto po=options.partition;po.sourceBodyTags={solidTag};po.cancellation=options.cancellation;
        po.maximumHistoryWork=std::min(po.maximumHistoryWork,context.remaining());
        OcctPartitionOutcome partition;
        {OcctStageTrace stage("partition");Timer timer{out.partitionSeconds};partition=partitionOcctModel(original,po);}
        out.partitionReason=partition.reason;out.planes=partition.planes;out.pieces=partition.pieces;out.interfaces=partition.interfaces.size();
        out.cadVolume=partition.sourceVolume;out.partitionVolume=partition.partitionVolume;out.partitionVolumeChange=partition.relativeVolumeChange;
        out.trimTolerance=partition.trimTolerance;
        out.sourceIntegrationError=partition.sourceIntegrationError;out.maximumPieceIntegrationError=partition.maximumPieceIntegrationError;
        out.sourceImmutabilityVerified=partition.sourceImmutabilityVerified;
        context.collect(partition.historyWork,out.partitionWork);
        if(!partition.partitioned()){
            if(partition.reason==OcctPartitionReason::BudgetExceeded)throw Stop{OcctPartitionMeshReason::BudgetExceeded,"partition-budget"};
            if(partition.reason==OcctPartitionReason::Cancelled)throw Stop{OcctPartitionMeshReason::Cancelled,"partition-cancelled"};
            if(partition.reason==OcctPartitionReason::AllocationFailed)throw Stop{OcctPartitionMeshReason::AllocationFailed,"partition-allocation-failed"};
            throw Stop{OcctPartitionMeshReason::PartitionFailed,occtPartitionReasonName(partition.reason)};
        }
        need(out.sourceImmutabilityVerified&&out.cadVolume>0&&std::isfinite(out.cadVolume),OcctPartitionMeshReason::ProvenanceLost,"uncertified-partition-source");
        std::set<int> ownedFaces,ownedEdges,ownedVertices;
        for(int face:solid->boundaryFaceTags()){context.charge(1,out.assemblyWork);ownedFaces.insert(face);
            for(const auto& wire:original.boundaryUses(face))for(const auto& use:wire.edges){context.charge(1,out.assemblyWork);ownedEdges.insert(use.edgeTag);if(use.startVertexTag>=0)ownedVertices.insert(use.startVertexTag);if(use.endVertexTag>=0)ownedVertices.insert(use.endVertexTag);}}
        std::map<int,OcctPartitionInterface> interfaces;
        for(const auto& cut:partition.interfaces){context.charge(1,out.assemblyWork);
            need(cut.faceTag>0&&cut.originalBodyTag==solidTag&&cut.pieceTags[0]>0&&cut.pieceTags[1]>0&&cut.pieceTags[0]!=cut.pieceTags[1]&&
                 cut.senses[0]*cut.senses[1]==-1&&interfaces.emplace(cut.faceTag,cut).second,
                 OcctPartitionMeshReason::InterfaceMismatch,"invalid-native-interface");
            need(parent(partition,{2,cut.faceTag})==mesh::Classification{3,solidTag},OcctPartitionMeshReason::ProvenanceLost,"cut-face-parent-mismatch");}
        auto bo=options.boundary;bo.cancellation=options.cancellation;bo.maximumWork=std::min(bo.maximumWork,context.remaining());
        need(bo.maximumWork>0,OcctPartitionMeshReason::BudgetExceeded,"surface-work-budget");
        bo.maximumNodes=std::min(bo.maximumNodes,options.maximumNodes);bo.maximumTriangles=std::min(bo.maximumTriangles,options.maximumCells);
        bo.maximumOutputBytes=std::min(bo.maximumOutputBytes,options.maximumOutputBytes);
        OcctBoundaryOutcome preparation;std::unique_ptr<OcctBoundarySession> session;
        {OcctStageTrace stage("surface");Timer timer{out.surfaceSeconds};session=OcctBoundarySession::prepare(*partition.model,bo,preparation,OcctBoundaryAuditPolicy::Deferred);}
        out.boundaryReason=preparation.reason;context.collect(preparation.work,out.surfaceWork);
        if(!session)boundaryFailure(preparation.reason);
        mesh::Mesh assembled;
        std::map<OcctBoundaryNodeKey,Index> nodes;
        std::map<Tri,std::vector<FaceUse>> faceUses;
        constexpr std::uint64_t materialCacheLimit=32768;
        const bool improveQuality=options.tetQualityTargetSize>0&&out.pieces<=24;
        const bool useMaterialCache=recovery.useLiveConstraints||improveQuality;
        std::uint64_t qualityNativeQueries=1200;
        const auto materialCacheBytes=useMaterialCache?materialCacheLimit*128:0;
        const auto baseWorking=add(materialCacheBytes,add(
            multiply(partition.parentByVertex.size()+partition.parentByEdge.size()+partition.parentByFace.size(),sizeof(mesh::Classification)),
            add(multiply(partition.parentBodyByPiece.size(),sizeof(int)),multiply(partition.interfaces.size(),sizeof(OcctPartitionInterface)+128))));
        const auto workingNow=[&](){return add(baseWorking,add(payload(assembled.nodeCount(),assembled.cellCount(),assembled.connectivityCount()),add(multiply(nodes.size(),128),multiply(faceUses.size(),512))));};
        const auto pieces=partition.model->solids();need(pieces.size()==out.pieces,OcctPartitionMeshReason::ProvenanceLost,"piece-count-mismatch");
        std::uint64_t exteriorCount=0;
        OcctSolidClassifier originalPieceClassifier;
        std::map<std::array<std::uint64_t,3>,bool> originalMaterialCache;
        const auto originalTolerance=original.tolerance(3,solidTag);
        if(useMaterialCache){
            context.working(workingNow());
            originalPieceClassifier.Load(OcctModelAccess::solid(original,solidTag));
        }
        const auto pointKey=[](Vec3 point){return std::array{std::bit_cast<std::uint64_t>(point.x),
            std::bit_cast<std::uint64_t>(point.y),std::bit_cast<std::uint64_t>(point.z)};};
        const auto originalAdmits=[&](Vec3 point){
            context.charge(1,out.validationWork);
            const auto key=pointKey(point);
            if(const auto known=originalMaterialCache.find(key);known!=originalMaterialCache.end())return known->second;
            ++out.nativeCandidateQueries;
            originalPieceClassifier.Perform(gp_Pnt(point.x,point.y,point.z),originalTolerance);
            const auto state=originalPieceClassifier.State();
            const bool inside=state==TopAbs_IN||state==TopAbs_ON;
            if(originalMaterialCache.size()<materialCacheLimit)originalMaterialCache.emplace(key,inside);
            return inside;
        };
        for(const auto* piece:pieces){
            const int tag=piece->tag();out.failedPieceTag=tag;context.charge(1,out.assemblyWork);
            if(options.maximumAdaptiveCuts){
                context.charge(1,out.assemblyWork);
                piece->bbox(out.failedPieceBBoxLow,out.failedPieceBBoxHigh);
                out.failedPieceBoundsAvailable=true;
            }
            need(tag>0&&static_cast<std::size_t>(tag)<partition.parentBodyByPiece.size()&&partition.parentBodyByPiece[tag]==solidTag,
                 OcctPartitionMeshReason::ProvenanceLost,"piece-body-parent-mismatch");
            mesh::Mesh boundary;OcctBoundaryOutcome surface;
            {OcctStageTrace stage("surface");Timer timer{out.surfaceSeconds};
                if(options.repairSurfaceAliases){
                    std::vector<int> allowed;
                    for(int face:piece->boundaryFaceTags()){context.charge(1,out.surfaceWork);if(parent(partition,{2,face}).dim==2)allowed.push_back(face);}
                    context.working(workingNow());need(workingNow()<options.maximumWorkingBytes,OcctPartitionMeshReason::BudgetExceeded,"alias-working-budget");
                    OcctBoundaryAliasOptions aliases;aliases.faceTags=allowed;aliases.maximumWorkingBytes=options.maximumWorkingBytes-workingNow();OcctBoundaryAliasReceipt repaired;
                    surface=session->boundary(tag,boundary,aliases,repaired,context.remaining());out.surfaceDiagonalRepairs+=repaired.flips;
                    context.working(add(workingNow(),repaired.maximumWorkingPayloadEstimate));
                }else surface=session->boundary(tag,boundary,context.remaining());
            }
            out.boundaryReason=surface.reason;context.collect(surface.work,out.surfaceWork);
            if(!surface.published)boundaryFailure(surface.reason);
            // This retained private session stages pieces without repeating the
            // whole-source byte audit. The final audit below precedes the only
            // publication; public/default sessions continue auditing each call.
            need(surface.nodeKeys.size()==static_cast<std::size_t>(boundary.nodeCount()),OcctPartitionMeshReason::ProvenanceLost,"surface-identity-receipt-missing");
            const auto before=add(workingNow(),payload(boundary.nodeCount(),boundary.cellCount(),boundary.connectivityCount()));
            context.working(add(before,multiply(add(boundary.nodeCount(),boundary.cellCount()),1024)));
            auto so=options.star;so.maximumWork=std::min({so.maximumWork,context.remaining(),options.maximumPieceStarWork});so.cancellation=options.cancellation;
            so.maximumNodes=std::min(so.maximumNodes,options.maximumNodes);so.maximumCells=std::min(so.maximumCells,options.maximumCells);
            so.maximumOutputBytes=std::min(so.maximumOutputBytes,options.maximumOutputBytes);
            so.maximumWorkingBytes=std::min(so.maximumWorkingBytes,options.maximumWorkingBytes-before);
            if(options.useNativeKernelAdmission&&!recovery.useLiveConstraints){
                so.nativeAdmissionModel=&original;so.nativeAdmissionBodyTag=solidTag;
            }
            mesh::Mesh fan;OcctStarOutcome star;
            {OcctStageTrace stage("volume_fill");Timer timer{out.starSeconds};star=buildOcctStarVolume(boundary,tag,so,fan);}
            out.starReason=star.reason;out.nativeCandidateQueries+=star.nativeCandidateQueries;
            context.collect(star.workUsed,out.starWork);out.kernelTrials+=star.kernelTrials;out.exactSurfacePairs+=star.exactBoundaryPairs;
            const bool nonzeroGenus=star.reason==OcctStarReason::UnsupportedTopology&&
                std::string_view(star.diagnostic)=="boundary-genus-not-zero";
            if(options.maximumAdaptiveCuts&&nonzeroGenus){
                // The star topology guard precedes its embedding audit. A
                // genus refusal can trigger a cut only after the independent
                // closed-surface certificate establishes exact embedding.
                so.maximumWork=std::min({context.remaining(),so.maximumWork-star.workUsed,
                                        options.maximumPieceStarWork-star.workUsed});
                OcctBoundaryCertificate certificate;
                {OcctStageTrace stage("volume_fill");Timer timer{out.starSeconds};certificate=certifyOcctBoundarySurface(boundary,so);}
                context.collect(certificate.workUsed,out.starWork);
                out.exactSurfacePairs+=certificate.exactBoundaryPairs;
                if(!certificate.certified){out.starReason=certificate.reason;starFailure(certificate.reason,certificate.diagnostic);}
                out.failedPieceBoundaryCertified=true;
            }
            bool recoveredByCore=false;
            if(recovery.enabled&&!star.published&&(star.reason==OcctStarReason::KernelNotFound||
               (recovery.useLiveConstraints&&nonzeroGenus))){
                context.working(add(before,multiply(add(add(boundary.nodeCount(),boundary.cellCount()),1),32768)));
                so.maximumWork=std::min(context.remaining(),recovery.maximumPieceWork);
                OcctPlcRecoveryOutcome recovered;
                {OcctStageTrace stage("piece_recovery",tag);Timer timer{out.starSeconds};recovered=recoverOcctPlcVolume(boundary,tag,*partition.model,original,solidTag,partition.parentByFace,so,recovery,fan);}
                if(receipt){
                    receipt->attempted|=recovered.attempted;receipt->reason=recovered.reason;receipt->diagnostic=recovered.diagnostic;receipt->failedPieceTag=recovered.failedPieceTag;
                    receipt->piecesAttempted=add(receipt->piecesAttempted,recovered.piecesAttempted);receipt->piecesRecovered=add(receipt->piecesRecovered,recovered.piecesRecovered);
                    receipt->steinerNodes=add(receipt->steinerNodes,recovered.steinerNodes);receipt->coreAttempts=add(receipt->coreAttempts,recovered.coreAttempts);receipt->quadFlips=add(receipt->quadFlips,recovered.quadFlips);
                    receipt->workUsed=add(receipt->workUsed,recovered.workUsed);receipt->generalWorkReserved=add(receipt->generalWorkReserved,recovered.generalWorkReserved);receipt->nativeQueries=add(receipt->nativeQueries,recovered.nativeQueries);
                }
                // Retained core stages can report an admission overrun after
                // the phase. Keep completed work even when refusing it.
                out.workUsed=add(out.workUsed,recovered.workUsed);out.starWork=add(out.starWork,recovered.workUsed);
                need(out.workUsed<=options.maximumWork,OcctPartitionMeshReason::BudgetExceeded,"recovery-cumulative-work-budget");context.charge(0,out.starWork);
                if(!recovered.published){
                    if(recovered.reason==OcctPlcRecoveryReason::BudgetExceeded)throw Stop{OcctPartitionMeshReason::BudgetExceeded,recovered.diagnostic};
                    if(recovered.reason==OcctPlcRecoveryReason::Cancelled)throw Stop{OcctPartitionMeshReason::Cancelled,recovered.diagnostic};
                    if(recovered.reason==OcctPlcRecoveryReason::AllocationFailed)throw Stop{OcctPartitionMeshReason::AllocationFailed,recovered.diagnostic};
                    if(recovered.reason==OcctPlcRecoveryReason::NativeQueryFailed)throw Stop{OcctPartitionMeshReason::GeometryError,recovered.diagnostic};
                    if(recovered.reason==OcctPlcRecoveryReason::InvalidRequest)throw Stop{OcctPartitionMeshReason::InvalidRequest,recovered.diagnostic};
                    if(recovered.reason==OcctPlcRecoveryReason::InternalError)throw Stop{OcctPartitionMeshReason::InternalError,recovered.diagnostic};
                    throw Stop{OcctPartitionMeshReason::StarFailed,recovered.diagnostic};
                }
                recoveredByCore=true;
            }
            if(!star.published&&!recoveredByCore)starFailure(star.reason,star.diagnostic);
            if(recovery.useLiveConstraints){
                if(options.maximumNativeValidationThreads>1&&fan.cellCount()>600&&
                   originalMaterialCache.size()<materialCacheLimit){
                    OcctStageTrace stage("cad_validation");Timer timer{out.validationSeconds};
                    const auto reserved=add(before,add(payload(fan.nodeCount(),fan.cellCount(),fan.connectivityCount()),
                        multiply(fan.cellCount(),sizeof(Vec3))));
                    context.working(reserved);
                    std::vector<Vec3> centroids;
                    centroids.reserve(fan.cellCount());
                    for(Index i=0;i<fan.cellCount();++i)if(fan.cell(i).type==mesh::CellType::Tet4){
                        context.charge(1,out.validationWork);
                        const auto* ids=fan.cellNodes(i);
                        centroids.push_back((fan.node(ids[0])+fan.node(ids[1])+fan.node(ids[2])+fan.node(ids[3]))*.25);
                    }
                    OcctNativeValidationOptions validation;
                    validation.maximumThreads=options.maximumNativeValidationThreads;
                    validation.maximumWork=context.remaining();
                    validation.maximumWorkingBytes=std::min(validation.maximumWorkingBytes,options.maximumWorkingBytes-reserved);
                    validation.maximumSnapshotBytes=std::min(validation.maximumSnapshotBytes,validation.maximumWorkingBytes/4);
                    validation.cancellation=options.cancellation;
                    std::vector<unsigned char> admitted;
                    const auto batch=classifyOcctNativeInteriorPoints(original,solidTag,centroids,validation,admitted);
                    context.collect(batch.workUsed,out.validationWork);
                    if(!batch.completed)nativeValidationFailure(batch.reason);
                    need(batch.sourceModelUnmodified&&admitted.size()==centroids.size(),
                         OcctPartitionMeshReason::ProvenanceLost,"incomplete-local-centroid-validation");
                    out.nativeCandidateQueries+=batch.nativeQueries+batch.originalReferenceQueries;
                    for(std::size_t i=0;i<centroids.size()&&originalMaterialCache.size()<materialCacheLimit;++i)
                        originalMaterialCache.try_emplace(pointKey(centroids[i]),admitted[i]!=0);
                }
                const auto repaired=repairPartitionMaterial(partition,tag,boundary,so,context,originalAdmits,fan,
                                                            recovery.useAdaptiveRefill);
                recoveredByCore|=repaired;
            }
            if(improveQuality&&qualityNativeQueries&&fan.cellCount()<=8192&&fan.nodeCount()<=4096&&
               add(before,16ULL*1024*1024)<options.maximumWorkingBytes&&context.remaining()>=10'000){
                OcctStageTrace stage("tet_quality_cavities",tag);
                Timer timer{out.starSeconds};
                const auto shapeBytes=payload(fan.nodeCount(),fan.cellCount(),fan.connectivityCount());
                const auto stepBytes=sizeof(mesh::Cell)+4*sizeof(Index)+sizeof(Vec3)+sizeof(mesh::Classification);
                const auto room=options.maximumOutputBytes>shapeBytes?
                    (options.maximumOutputBytes-shapeBytes)/stepBytes:0;
                const auto tets=std::uint64_t(fan.countCells(mesh::CellType::Tet4));
                const bool improved=improveOcctTetPiece(fan,
                    {std::min<std::uint64_t>(500'000,context.remaining()/4),
                     std::min(options.maximumNodes,std::uint64_t(fan.nodeCount())+room),
                     std::min(options.maximumCells-boundary.cellCount(),tets+room),
                     options.tetQualityTargetSize,qualityNativeQueries},originalAdmits,
                    [&](std::uint64_t n){context.charge(n,out.starWork);});
                recoveredByCore|=improved;
                stage.finish(improved?"improved":"unchanged");
            }
            {OcctStageTrace stage("partition_assembly");Timer timer{out.assemblySeconds};
                std::vector<Index> remap(static_cast<std::size_t>(fan.nodeCount()),-1);
                const auto activeWorking=add(payload(boundary.nodeCount(),boundary.cellCount(),boundary.connectivityCount()),
                    add(payload(fan.nodeCount(),fan.cellCount(),fan.connectivityCount()),multiply(remap.size(),sizeof(Index))));
                for(Index i=0;i<boundary.nodeCount();++i){context.charge(1,out.assemblyWork);const auto cls=boundary.nodeClass(i);const auto& key=surface.nodeKeys[i];
                    need(cls.dim==key.dimension&&cls.tag==key.tag,OcctPartitionMeshReason::ProvenanceLost,"key-classification-mismatch");
                    const auto mapped=parent(partition,cls);const auto p=boundary.node(i);
                    auto found=nodes.find(key);
                    if(found!=nodes.end()){const auto id=found->second;need(identical(assembled.node(id),p)&&assembled.nodeClass(id)==mapped,OcctPartitionMeshReason::InterfaceMismatch,"repeated-node-key-payload-mismatch");remap[i]=id;}
                    else{context.output(assembled.nodeCount()+1,add(assembled.cellCount(),exteriorCount),add(assembled.connectivityCount(),multiply(exteriorCount,3)));
                        context.working(add(add(workingNow(),activeWorking),1024));const auto id=assembled.addNode(p,mapped);nodes.emplace(key,id);remap[i]=id;}}
                if(!recoveredByCore){
                const auto kernelIndex=boundary.nodeCount();need(fan.nodeCount()==kernelIndex+1&&fan.nodeClass(kernelIndex)==mesh::Classification{3,tag},OcctPartitionMeshReason::InvalidMesh,"unexpected-star-interior-nodes");
                context.output(assembled.nodeCount()+1,add(assembled.cellCount(),exteriorCount),add(assembled.connectivityCount(),multiply(exteriorCount,3)));
                remap[kernelIndex]=assembled.addNode(fan.node(kernelIndex),{3,solidTag});
                }else{
                    for(Index k=boundary.nodeCount();k<fan.nodeCount();++k){context.charge(1,out.assemblyWork);
                        need(fan.nodeClass(k)==mesh::Classification{3,tag},OcctPartitionMeshReason::InvalidMesh,"unexpected-recovery-interior-class");
                        context.output(assembled.nodeCount()+1,add(assembled.cellCount(),exteriorCount),add(assembled.connectivityCount(),multiply(exteriorCount,3)));
                        remap[k]=assembled.addNode(fan.node(k),{3,solidTag});
                    }
                }
                for(Index i=0;i<fan.cellCount();++i){context.charge(1,out.assemblyWork);const auto& cell=fan.cell(i);const auto* ids=fan.cellNodes(i);
                    if(cell.type==mesh::CellType::Tet4&&cell.onEntity==mesh::Classification{3,tag}){
                        const std::array<Index,4> t{remap[ids[0]],remap[ids[1]],remap[ids[2]],remap[ids[3]]};
                        need(predicates::orient3d(assembled.node(t[0]),assembled.node(t[1]),assembled.node(t[2]),assembled.node(t[3]))==predicates::Sign::Positive,OcctPartitionMeshReason::InvalidMesh,"nonpositive-joined-tet");
                        context.output(assembled.nodeCount(),add(add(assembled.cellCount(),1),exteriorCount),add(add(assembled.connectivityCount(),4),multiply(exteriorCount,3)));
                        assembled.addCell(mesh::CellType::Tet4,t,{3,solidTag});
                    }else{
                        need(cell.type==mesh::CellType::Tri3&&cell.onEntity.dim==2,OcctPartitionMeshReason::InvalidMesh,"unexpected-star-cell");
                        const Tri t{remap[ids[0]],remap[ids[1]],remap[ids[2]]};const auto mapped=parent(partition,cell.onEntity);
                        if(mapped.dim==2){++exteriorCount;context.output(assembled.nodeCount(),add(assembled.cellCount(),exteriorCount),add(assembled.connectivityCount(),multiply(exteriorCount,3)));}
                        else need(mapped==mesh::Classification{3,solidTag}&&interfaces.contains(cell.onEntity.tag),OcctPartitionMeshReason::ProvenanceLost,"undeclared-cut-face");
                        context.working(add(add(workingNow(),activeWorking),1024));auto& uses=faceUses[sorted(t)];need(uses.size()<2,OcctPartitionMeshReason::InterfaceMismatch,"nonmanifold-joined-facet");uses.push_back({t,tag,cell.onEntity.tag,mapped});
                    }
                }
            }
        }
        out.failedPieceTag=-1;
        out.failedPieceBoundsAvailable=false;
        out.failedPieceBoundaryCertified=false;
        mesh::Mesh exterior;std::vector<Index> exteriorNodes(static_cast<std::size_t>(assembled.nodeCount()),-1);
        std::map<int,std::uint64_t> pairedByFace;
        std::set<int> retainedOriginalFaces;
        {OcctStageTrace stage("partition_assembly");Timer timer{out.assemblySeconds};
            for(const auto& [key,uses]:faceUses){context.charge(1,out.assemblyWork);
                if(uses.size()==2){const auto& a=uses[0];const auto& b=uses[1];auto native=interfaces.find(a.derivedFace);
                    need(native!=interfaces.end()&&a.derivedFace==b.derivedFace&&a.original==mesh::Classification{3,solidTag}&&b.original==a.original&&
                         !sameWinding(a.ids,b.ids),OcctPartitionMeshReason::InterfaceMismatch,"undeclared-or-miswound-paired-facet");
                    auto owners=std::array{a.piece,b.piece},expected=native->second.pieceTags;std::sort(owners.begin(),owners.end());std::sort(expected.begin(),expected.end());
                    need(owners==expected,OcctPartitionMeshReason::InterfaceMismatch,"paired-facet-owner-mismatch");
                    ++out.pairedInterfaceTriangles;++pairedByFace[a.derivedFace];
                }else{
                    need(uses.size()==1&&uses[0].original.dim==2&&ownedFaces.contains(uses[0].original.tag),OcctPartitionMeshReason::InterfaceMismatch,"unpaired-cut-or-foreign-exterior-facet");
                    const auto& use=uses[0];Tri compact;
                    for(unsigned j=0;j<3;++j){const auto id=use.ids[j];auto& compactId=exteriorNodes[id];if(compactId<0){const auto cls=assembled.nodeClass(id);need(cls.dim>=0&&cls.dim<=2,OcctPartitionMeshReason::ProvenanceLost,"interior-node-on-original-exterior");compactId=exterior.addNode(assembled.node(id),cls);}compact[j]=compactId;}
                    assembled.addCell(mesh::CellType::Tri3,use.ids,use.original);exterior.addCell(mesh::CellType::Tri3,compact,use.original);
                    retainedOriginalFaces.insert(use.original.tag);
                }
            }
            for(const auto& [tag,cut]:interfaces){context.charge(1,out.assemblyWork);need(pairedByFace.contains(tag)&&pairedByFace[tag]>0,OcctPartitionMeshReason::InterfaceMismatch,"native-interface-not-paired");}
            need(retainedOriginalFaces==ownedFaces&&static_cast<std::uint64_t>(exterior.cellCount())==exteriorCount,
                 OcctPartitionMeshReason::ProvenanceLost,"original-face-coverage-mismatch");
        }
        out.nodes=assembled.nodeCount();out.cells=assembled.cellCount();out.surfaceTriangles=exterior.cellCount();
        context.output(out.nodes,out.cells,assembled.connectivityCount());
        const auto beforeCertificate=add(workingNow(),add(payload(exterior.nodeCount(),exterior.cellCount(),exterior.connectivityCount()),multiply(exteriorNodes.size(),sizeof(Index))));
        context.working(add(beforeCertificate,multiply(add(exterior.nodeCount(),exterior.cellCount()),1024)));
        {OcctStageTrace stage("cad_validation");Timer timer{out.validationSeconds};
            auto so=options.star;so.maximumWork=std::min(so.maximumWork,context.remaining());so.cancellation=options.cancellation;
            so.maximumWorkingBytes=std::min(so.maximumWorkingBytes,options.maximumWorkingBytes-beforeCertificate);
            out.exteriorCertificate=certifyOcctBoundarySurface(exterior,so);context.collect(out.exteriorCertificate.workUsed,out.validationWork);
            if(!out.exteriorCertificate.certified)starFailure(out.exteriorCertificate.reason,out.exteriorCertificate.diagnostic);
            out.exactSurfacePairs+=out.exteriorCertificate.exactBoundaryPairs;
            // Loading the native solid builds its classifier search structure.
            // Keep that structure local to this call and perform every original
            // volume-node and Tet4-centroid query with NativeSolid's tolerance.
            // No shared mutable model cache or inclusion sample is removed.
            paacs::geo::detail::OcctSolidClassifier classifier;
            try{if(!useMaterialCache)classifier.Load(OcctModelAccess::solid(original,solidTag));}
            catch(const Standard_OutOfMemory&){throw;}
            catch(const Standard_Failure&){throw Stop{OcctPartitionMeshReason::GeometryError,"original-solid-classifier-failed"};}
            const auto solidTolerance=original.tolerance(3,solidTag);
            const auto nativeInside=[&](const Vec3& p){
                if(!finite(p))return false;
                if(useMaterialCache)return originalAdmits(p);
                try{classifier.Perform(gp_Pnt(p.x,p.y,p.z),solidTolerance);
                    return classifier.State()==TopAbs_IN||classifier.State()==TopAbs_ON;}
                catch(const Standard_OutOfMemory&){throw;}
                catch(const Standard_Failure&){return false;}
            };
            OcctStageTrace witnesses("cad_nodes",solidTag);
            for(Index i=0;i<assembled.nodeCount();++i){context.charge(1,out.validationWork);const auto p=assembled.node(i);const auto cls=assembled.nodeClass(i);need(finite(p),OcctPartitionMeshReason::GeometryError,"nonfinite-original-node");
                if(cls.dim==3){need(cls.tag==solidTag&&nativeInside(p),OcctPartitionMeshReason::GeometryError,"interior-node-outside-original-solid");continue;}
                Vec3 q;
                if(cls.dim==0){need(ownedVertices.contains(cls.tag)&&original.vertex(cls.tag),OcctPartitionMeshReason::ProvenanceLost,"foreign-original-vertex");q=original.vertex(cls.tag)->position();}
                else if(cls.dim==1){need(ownedEdges.contains(cls.tag)&&original.edge(cls.tag),OcctPartitionMeshReason::ProvenanceLost,"foreign-original-edge");q=original.edge(cls.tag)->closestPoint(p);}
                else{need(cls.dim==2&&ownedFaces.contains(cls.tag)&&closestPointOnOcctFace(original,cls.tag,p,q),OcctPartitionMeshReason::ProvenanceLost,"foreign-or-unprojectable-original-face");}
                const auto error=(p-q).norm();out.maximumBoundaryNodeError=std::max(out.maximumBoundaryNodeError,error);
                const auto tolerance=original.tolerance(cls.dim,cls.tag);
                need(finite(q)&&std::isfinite(error)&&std::isfinite(tolerance)&&tolerance>=0&&error<=std::max(8*tolerance,1e-9),OcctPartitionMeshReason::GeometryError,"original-boundary-node-error");
            }
            witnesses.next("cad_centroids");
            if(options.maximumNativeValidationThreads==1||recovery.useLiveConstraints){
                // Recheck every final centroid. The local route can reuse
                // prior classifications at bit-identical coordinates.
                for(Index i=0;i<assembled.cellCount();++i)if(assembled.cell(i).type==mesh::CellType::Tet4){context.charge(1,out.validationWork);const auto* ids=assembled.cellNodes(i);
                    const auto centroid=(assembled.node(ids[0])+assembled.node(ids[1])+assembled.node(ids[2])+assembled.node(ids[3]))*.25;
                    need(nativeInside(centroid),OcctPartitionMeshReason::GeometryError,"native-tet-centroid-outside");}
            }else{
                // Gather every Tet4 centroid in stable cell order. Native copy
                // and classifier heap storage remain outside this logical
                // payload admission; the helper owns all mutable query state.
                const auto centroidPayload=multiply(assembled.cellCount(),sizeof(Vec3));
                const auto retained=add(beforeCertificate,centroidPayload);
                context.working(retained);
                std::vector<Vec3> centroids;centroids.reserve(static_cast<std::size_t>(assembled.cellCount()));
                for(Index i=0;i<assembled.cellCount();++i)if(assembled.cell(i).type==mesh::CellType::Tet4){
                    context.charge(1,out.validationWork);const auto* ids=assembled.cellNodes(i);
                    centroids.push_back((assembled.node(ids[0])+assembled.node(ids[1])+assembled.node(ids[2])+assembled.node(ids[3]))*.25);
                }
                need(!centroids.empty(),OcctPartitionMeshReason::InvalidMesh,"assembled-volume-has-no-tets");
                OcctNativeValidationOptions validation;
                validation.maximumThreads=options.maximumNativeValidationThreads;
                validation.maximumPoints=options.maximumCells;validation.maximumWork=context.remaining();
                validation.maximumWorkingBytes=options.maximumWorkingBytes-retained;
                validation.maximumSnapshotBytes=std::min(options.partition.maximumSnapshotAuditBytes,
                                                          validation.maximumWorkingBytes/8);
                validation.cancellation=options.cancellation;
                need(validation.maximumWork&&validation.maximumWorkingBytes&&validation.maximumSnapshotBytes,
                     OcctPartitionMeshReason::BudgetExceeded,"native-validation-admission");
                out.nativeCentroidValidation=validateOcctNativeInteriorPoints(original,solidTag,centroids,validation);
                context.collect(out.nativeCentroidValidation.workUsed,out.validationWork);
                context.working(add(retained,out.nativeCentroidValidation.maximumWorkingPayloadEstimate));
                if(!out.nativeCentroidValidation.completed)nativeValidationFailure(out.nativeCentroidValidation.reason);
                need(out.nativeCentroidValidation.sourceModelUnmodified&&
                     out.nativeCentroidValidation.points==centroids.size()&&
                     out.nativeCentroidValidation.nativeQueries==centroids.size(),
                     OcctPartitionMeshReason::ProvenanceLost,"incomplete-native-centroid-validation");
            }
            witnesses.next("cad_facets");
            const auto validateFacets=[&]{
                constexpr std::uint64_t midpointEntryBytes=128;
                const auto midpointLimit=recovery.useAdaptiveRefill?std::min({
                    std::uint64_t{16384},multiply(exterior.cellCount(),3),
                    (options.maximumWorkingBytes-beforeCertificate)/midpointEntryBytes}):0;
                context.working(add(beforeCertificate,multiply(midpointLimit,midpointEntryBytes)));
                std::set<std::array<Index,3>> acceptedMidpoints;
                const auto chord=options.boundary.maximumChordError;
                for(Index i=0;i<exterior.cellCount();++i){
                    const auto tag=exterior.cell(i).onEntity.tag;
                    const auto* ids=exterior.cellNodes(i);
                    const auto a=exterior.node(ids[0]),b=exterior.node(ids[1]),c=exterior.node(ids[2]);
                    const auto* face=recovery.useAdaptiveRefill?original.face(tag):nullptr;
                    const auto limit=chord+original.tolerance(2,tag);
                    const std::array samples{(a+b+c)/3.,(a+b)*.5,(b+c)*.5,(c+a)*.5};
                    for(unsigned k=0;k<samples.size();++k){
                        context.charge(1,out.validationWork);++out.sampleWork;
                        const auto sample=samples[k];
                        std::array<Index,3> midpointKey{};
                        if(k&&midpointLimit){
                            const auto u=ids[k-1],v=ids[k%3];
                            midpointKey={Index(tag),std::min(u,v),std::max(u,v)};
                            // Each accepted witness already contributed to the maximum.
                            if(acceptedMidpoints.contains(midpointKey))continue;
                        }
                        // An on-trim point gives an upper bound on the distance.
                        // If it cannot prove the chord bound, use the full solve.
                        Vec2 uv;
                        Vec3 closest;
                        double error=std::numeric_limits<double>::infinity();
                        if(face&&face->projectPoint(sample,uv,closest))error=(sample-closest).norm();
                        if(!(error<=chord)){
                            need(closestPointOnOcctFace(original,tag,sample,closest),
                                 OcctPartitionMeshReason::GeometryError,"original-facet-projection-failed");
                            error=(sample-closest).norm();
                        }
                        need(finite(closest)&&std::isfinite(error)&&error<=limit,
                             OcctPartitionMeshReason::GeometryError,"original-facet-chord-error");
                        out.maximumBoundarySampleError=std::max(out.maximumBoundarySampleError,error);
                        if(k&&acceptedMidpoints.size()<midpointLimit)acceptedMidpoints.insert(midpointKey);
                    }
                }
            };
            validateFacets();
            witnesses.finish();
            auto auditLimits=options.star.auditLimits;auditLimits.cancellation=options.cancellation;
            out.audit=mesh::auditLinearVolumeBoundary(assembled,auditLimits);
            if(out.audit.cancelled)throw Stop{OcctPartitionMeshReason::Cancelled,"final-audit-cancelled"};
            if(out.audit.budgetExceeded)throw Stop{OcctPartitionMeshReason::BudgetExceeded,"final-audit-budget"};
            if(out.audit.allocationFailed)throw Stop{OcctPartitionMeshReason::AllocationFailed,"final-audit-allocation-failed"};
            need(out.audit.valid,OcctPartitionMeshReason::InvalidMesh,"final-boundary-incidence-invalid");
            context.charge(assembled.cellCount(),out.validationWork);out.quality=algo::analyzeMixedMeshQuality(assembled);
            need(!out.quality.allocationFailed,OcctPartitionMeshReason::AllocationFailed,"final-quality-allocation-failed");
            need(out.quality.valid,OcctPartitionMeshReason::InvalidMesh,"final-positive-quality-invalid");out.meshVolume=out.quality.totalVolume;
            need(std::abs(out.meshVolume-out.exteriorCertificate.signedVolume)<=out.meshVolume*1e-10,OcctPartitionMeshReason::InvalidMesh,"joined-exterior-volume-mismatch");
            out.relativeVolumeError=std::abs(out.meshVolume-out.cadVolume)/out.cadVolume;
            need(std::isfinite(out.relativeVolumeError)&&out.relativeVolumeError<=options.maximumRelativeVolumeError,OcctPartitionMeshReason::VolumeMismatch,"original-cad-volume-error");
        }
        OcctBoundaryOutcome sourceAudit;
        {OcctStageTrace stage("surface");Timer timer{out.surfaceSeconds};sourceAudit=session->auditSource(context.remaining());}
        out.boundaryReason=sourceAudit.reason;context.collect(sourceAudit.work,out.surfaceWork);
        if(sourceAudit.reason!=OcctBoundaryReason::None)boundaryFailure(sourceAudit.reason);
        need(sourceAudit.sourceModelUnmodified,OcctPartitionMeshReason::ProvenanceLost,"surface-source-audit-missing");
        context.charge(0,out.assemblyWork);destination.swap(assembled);out.published=true;out.destinationUnchanged=false;out.reason=OcctPartitionMeshReason::None;out.diagnostic="none";
    }catch(const Stop& stop){out.reason=stop.reason;out.diagnostic=stop.diagnostic;}
     catch(const std::bad_alloc&){out.reason=OcctPartitionMeshReason::AllocationFailed;out.diagnostic="allocation-failed";}
     catch(const Standard_OutOfMemory&){out.reason=OcctPartitionMeshReason::AllocationFailed;out.diagnostic="native-allocation-failed";}
     catch(...){out.reason=OcctPartitionMeshReason::InternalError;out.diagnostic="internal-error";}
    out.totalSeconds=elapsed(started);return out;
}

OcctPartitionMeshOutcome buildOcctPartitionVolume(const OcctModel& original,int solidTag,
    const OcctPartitionMeshOptions& options,mesh::Mesh& destination) noexcept {
    return buildOcctPartitionVolume(original,solidTag,options,OcctPartitionRecoveryOptions{},destination,nullptr);
}
OcctPartitionMeshOutcome buildOcctPartitionVolume(const OcctModel& original,int solidTag,
    const OcctPartitionMeshOptions& options,const OcctPartitionRecoveryOptions& recovery,
    mesh::Mesh& destination,OcctPartitionRecoveryReceipt* receipt) noexcept {
    if(receipt)*receipt={};
    if(recovery.enabled&&(!options.useNativeKernelAdmission||!recovery.maximumPieceWork||!recovery.maximumCoreAttempts||recovery.maximumCoreAttempts>16||!recovery.maximumQuadFlipsPerAttempt||recovery.maximumQuadFlipsPerAttempt>64||!recovery.maximumGeneralWorkPerAttempt)){
        OcctPartitionMeshOutcome invalid;invalid.reason=OcctPartitionMeshReason::InvalidRequest;invalid.diagnostic="invalid-recovery-options";return invalid;
    }
    if(!options.maximumAdaptiveCuts)
        return buildOcctPartitionVolumeAttempt(original,solidTag,options,destination,recovery,receipt);
    const auto started=Clock::now();OcctPartitionMeshOutcome out;
    OcctPartitionMeshOutcome totals;
    std::vector<OcctPartitionAdaptiveCut> cuts;
    unsigned attempts=0;int principalFace=-1;
    try {
        need(options.maximumAdaptiveCuts<=16&&
             options.partition.policy==OcctPartitionPolicy::NativeConcaveFeatureAxes,
             OcctPartitionMeshReason::InvalidRequest,"invalid-adaptive-policy");
        const auto* solid=original.solid(solidTag);
        need(solid,OcctPartitionMeshReason::InvalidRequest,"invalid-adaptive-solid");
        const auto retained=multiply(options.partition.extraPlanes.size()+options.maximumAdaptiveCuts,
                                     sizeof(OcctPartitionPlane));
        const auto heldPayload=add(retained,multiply(options.maximumAdaptiveCuts,sizeof(OcctPartitionAdaptiveCut)));
        need(heldPayload<options.maximumWorkingBytes,OcctPartitionMeshReason::BudgetExceeded,"adaptive-working-payload-admission");
        cuts.reserve(options.maximumAdaptiveCuts);
        auto current=options;
        current.partition.extraPlanes.reserve(options.partition.extraPlanes.size()+options.maximumAdaptiveCuts);
        Vec3 principalAxis;double largest=-1;
        const auto checkpoint=[&](std::uint64_t count=0){
            need(!options.cancellation||!options.cancellation->load(std::memory_order_relaxed),
                 OcctPartitionMeshReason::Cancelled,"cancelled");
            need(count<=options.maximumWork-totals.workUsed,
                 OcctPartitionMeshReason::BudgetExceeded,"adaptive-cumulative-work-budget");
            totals.workUsed+=count;totals.adaptiveWork+=count;
        };
        for(const auto& use:original.faceUses(solidTag)){
            checkpoint(1);const auto* descriptor=original.surfaceDescriptor(use.faceTag);
            need(descriptor,OcctPartitionMeshReason::ProvenanceLost,"adaptive-face-descriptor-missing");
            if(descriptor->kind!=SurfaceKind::Cylinder)continue;
            checkpoint(1);GProp_GProps area;
            BRepGProp::SurfaceProperties(OcctModelAccess::face(original,use.faceTag),area);
            need(std::isfinite(area.Mass())&&area.Mass()>=0,
                 OcctPartitionMeshReason::PartitionFailed,"adaptive-native-area-invalid");
            if(area.Mass()>largest){largest=area.Mass();principalFace=use.faceTag;principalAxis=descriptor->axis;}
        }
        need(principalFace>0&&finite(principalAxis),OcctPartitionMeshReason::PartitionFailed,"adaptive-principal-axis-missing");
        mesh::Mesh staged;
        for(;;){
            checkpoint();current.maximumWork=options.maximumWork-totals.workUsed;
            current.maximumWorkingBytes=options.maximumWorkingBytes-heldPayload;
            auto admittedRecovery=recovery;
            if(recovery.onlyFinalAdaptiveAttempt&&cuts.size()!=options.maximumAdaptiveCuts)admittedRecovery.enabled=false;
            auto child=buildOcctPartitionVolumeAttempt(original,solidTag,current,staged,admittedRecovery,receipt);
            const bool completed=child.published;
            ++attempts;
            totals.workUsed=add(totals.workUsed,child.workUsed);
            totals.partitionWork=add(totals.partitionWork,child.partitionWork);
            totals.surfaceWork=add(totals.surfaceWork,child.surfaceWork);
            totals.starWork=add(totals.starWork,child.starWork);
            totals.assemblyWork=add(totals.assemblyWork,child.assemblyWork);
            totals.validationWork=add(totals.validationWork,child.validationWork);
            totals.sampleWork=add(totals.sampleWork,child.sampleWork);
            totals.kernelTrials=add(totals.kernelTrials,child.kernelTrials);
            totals.nativeCandidateQueries=add(totals.nativeCandidateQueries,child.nativeCandidateQueries);
            totals.surfaceDiagonalRepairs=add(totals.surfaceDiagonalRepairs,child.surfaceDiagonalRepairs);
            totals.exactSurfacePairs=add(totals.exactSurfacePairs,child.exactSurfacePairs);
            totals.maximumWorkingPayloadEstimate=std::max(totals.maximumWorkingPayloadEstimate,
                                                         add(child.maximumWorkingPayloadEstimate,heldPayload));
            totals.partitionSeconds+=child.partitionSeconds;totals.surfaceSeconds+=child.surfaceSeconds;
            totals.starSeconds+=child.starSeconds;totals.assemblySeconds+=child.assemblySeconds;
            totals.validationSeconds+=child.validationSeconds;
            out=std::move(child);out.published=false;out.destinationUnchanged=true;
            need(totals.workUsed<=options.maximumWork,OcctPartitionMeshReason::BudgetExceeded,"adaptive-cumulative-work-budget");
            if(completed){
                checkpoint();destination.swap(staged);out.published=true;out.destinationUnchanged=false;break;
            }
            const bool certifiedGenus=out.starReason==OcctStarReason::UnsupportedTopology&&
                std::string_view(out.diagnostic)=="boundary-genus-not-zero"&&out.failedPieceBoundaryCertified;
            if(cuts.size()==options.maximumAdaptiveCuts||!out.failedPieceBoundsAvailable||out.failedPieceTag<1||
               (!certifiedGenus&&out.starReason!=OcctStarReason::KernelNotFound&&out.starReason!=OcctStarReason::WorkBudget)||
               (out.reason!=OcctPartitionMeshReason::StarFailed&&out.reason!=OcctPartitionMeshReason::BudgetExceeded))break;
            const auto low=out.failedPieceBBoxLow,high=out.failedPieceBBoxHigh;
            need(finite(low)&&finite(high)&&high.x>low.x&&high.y>low.y&&high.z>low.z,
                 OcctPartitionMeshReason::InvalidMesh,"adaptive-failed-piece-bounds-invalid");
            const auto extent=high-low,midpoint=low*.5+high*.5;
            const std::array<Vec3,3> axes{Vec3{1,0,0},Vec3{0,1,0},Vec3{0,0,1}};
            const std::array<double,3> spans{extent.x,extent.y,extent.z};
            unsigned chosen=0;double longest=-1;
            for(unsigned i=0;i<3;++i){
                checkpoint(1);const auto score=spans[i]*(1-std::abs(axes[i].dot(principalAxis)));
                if(score>longest){longest=score;chosen=i;}
            }
            need(longest>0&&finite(midpoint),OcctPartitionMeshReason::InvalidMesh,"adaptive-transverse-span-invalid");
            const OcctPartitionPlane proposed{midpoint,axes[chosen]};
            bool duplicate=false;
            for(const auto& prior:current.partition.extraPlanes){
                checkpoint(1);const double norm=prior.normal.norm();
                need(std::isfinite(norm)&&norm>0,OcctPartitionMeshReason::InvalidRequest,"adaptive-prior-plane-invalid");
                const auto n=prior.normal/norm;const auto dot=n.dot(proposed.normal);
                if(std::abs(dot)>1-1e-10&&std::abs(proposed.normal.dot(proposed.origin)-
                    (dot>0?n.dot(prior.origin):-n.dot(prior.origin)))<=1e-8){duplicate=true;break;}
            }
            if(duplicate)break;
            need(current.partition.extraPlanes.size()<current.partition.maximumPlanes,
                 OcctPartitionMeshReason::BudgetExceeded,"adaptive-plane-budget");
            checkpoint(1);
            cuts.push_back({out.failedPieceTag,low,high,proposed,out.starReason,out.workUsed,totals.workUsed});
            current.partition.extraPlanes.push_back(proposed);
        }
    }catch(const Stop& stop){out.reason=stop.reason;out.diagnostic=stop.diagnostic;}
     catch(const std::bad_alloc&){out.reason=OcctPartitionMeshReason::AllocationFailed;out.diagnostic="allocation-failed";}
     catch(const Standard_OutOfMemory&){out.reason=OcctPartitionMeshReason::AllocationFailed;out.diagnostic="native-allocation-failed";}
     catch(...){out.reason=OcctPartitionMeshReason::InternalError;out.diagnostic="internal-error";}
    out.workUsed=totals.workUsed;out.partitionWork=totals.partitionWork;out.surfaceWork=totals.surfaceWork;
    out.starWork=totals.starWork;out.assemblyWork=totals.assemblyWork;out.validationWork=totals.validationWork;
    out.sampleWork=totals.sampleWork;out.kernelTrials=totals.kernelTrials;out.exactSurfacePairs=totals.exactSurfacePairs;out.surfaceDiagonalRepairs=totals.surfaceDiagonalRepairs;
    out.nativeCandidateQueries=totals.nativeCandidateQueries;
    out.maximumWorkingPayloadEstimate=totals.maximumWorkingPayloadEstimate;
    out.partitionSeconds=totals.partitionSeconds;out.surfaceSeconds=totals.surfaceSeconds;
    out.starSeconds=totals.starSeconds;out.assemblySeconds=totals.assemblySeconds;
    out.validationSeconds=totals.validationSeconds;out.totalSeconds=elapsed(started);
    out.adaptiveAttempts=attempts;out.adaptivePrincipalFaceTag=principalFace;out.adaptiveWork=totals.adaptiveWork;
    out.adaptiveCutReports=std::move(cuts);return out;
}
} // namespace paacs::geo::detail
