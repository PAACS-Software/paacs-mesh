// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctPlcFrontRecovery.hpp"
#include "OcctStageTrace.hpp"
#include "OcctModelInternal.hpp"
#include "OcctSolidClassifier.hpp"
#include "../bnd/TetConstraintRecovery.hpp"
#include "../bnd/TetPlcGiftWrapping.hpp"
#include "../bnd/TetPlcRetainedRecovery.hpp"
#include "paacs/bnd/TetPlcCavityRecovery.hpp"
#include <Standard_OutOfMemory.hxx>
#include <Standard_Failure.hxx>
#include <cmath>
#include <optional>
#include <stop_token>
#include <string_view>

namespace paacs::geo::detail {
namespace {
struct Stop {OcctPlcRecoveryReason reason;const char* diagnostic;};
void require(bool yes,OcctPlcRecoveryReason reason,const char* why){if(!yes)throw Stop{reason,why};}
}
OcctPlcRecoveryOutcome recoverOcctPlcLiveVolume(const mesh::Mesh& source,int pieceTag,
    const OcctStarOptions& options,const OcctPartitionRecoveryOptions& recovery,
    mesh::Mesh& destination) noexcept {
    using R=OcctPlcRecoveryReason;
    OcctPlcRecoveryOutcome out;
    out.attempted=true;
    out.piecesAttempted=1;
    out.failedPieceTag=pieceTag;
    const auto cap=std::min({options.maximumWork,recovery.maximumPieceWork,std::uint64_t{12'000'000}});
    const auto charge=[&](std::uint64_t count){
        require(count<=cap-out.workUsed,R::BudgetExceeded,"live-piece-work-budget");
        out.workUsed+=count;
        require(!options.cancellation||!options.cancellation->load(std::memory_order_relaxed),R::Cancelled,"cancelled");
    };
    try {
        charge(0);
        const auto nodes=std::uint64_t(source.nodeCount()),faces=std::uint64_t(source.cellCount());
        require(nodes>=4&&faces>=4&&pieceTag>0&& &source!=&destination,R::InvalidRequest,"live-piece-invalid-source");
        require(nodes<=2048&&faces<=4096,R::UnsupportedGeometry,"live-piece-count-limit");
        require(nodes<=options.maximumNodes&&faces<options.maximumCells,R::BudgetExceeded,"live-piece-output-budget");
        const auto retained=(nodes+faces)*4096+18ULL*1024*1024;
        require(retained<options.maximumWorkingBytes,R::BudgetExceeded,"live-piece-working-budget");
        auto certificateOptions=options;
        certificateOptions.maximumWork=cap;
        const auto certificate=certifyOcctBoundarySurface(source,certificateOptions);
        charge(certificate.workUsed);
        if(!certificate.certified){
            const auto reason=certificate.reason==OcctStarReason::AllocationFailed?R::AllocationFailed:
                certificate.reason==OcctStarReason::WorkBudget||certificate.reason==OcctStarReason::OutputBudget?
                R::BudgetExceeded:R::UnsupportedGeometry;
            throw Stop{reason,certificate.diagnostic};
        }
        const auto maxTets=std::min<std::uint64_t>(10'000,options.maximumCells-faces);
        mesh::Mesh seed;
        OcctStageTrace stage("constraint-recovery",pieceTag);
        const auto constraints=bnd::detail::recoverTetConstraints(source,{3,pieceTag},
            (cap-out.workUsed)*2/3,(options.maximumWorkingBytes-retained)/2,maxTets,
            options.cancellation,seed,std::min<std::uint64_t>(32,options.maximumNodes-nodes));
        stage.finish(constraints.diagnostic);
        charge(constraints.work);
        out.steinerNodes=constraints.retainedSteinerNodes;
        require(std::string_view(constraints.diagnostic)!="allocation_failed",R::AllocationFailed,"allocation-failed");
        require(seed.cellCount()>0,R::UnsupportedGeometry,"live-piece-no-seed");
        bnd::TetPlcCavityOptions limits;
        limits.region={3,pieceTag};
        limits.maximumNodes=std::min<std::uint64_t>(4096,options.maximumNodes);
        limits.maximumCells=maxTets;
        limits.maximumSourceTriangles=4096;
        limits.maximumWorkingBytes=options.maximumWorkingBytes-retained;
        limits.maximumWork=(cap-out.workUsed)-(cap-out.workUsed)/4;
        limits.maximumComponentFaces=4096;
        limits.maximumPatchCells=2048;
        limits.maximumExpansionCells=256;
        limits.maximumCandidatePoints=64;
        limits.maximumRandomCandidatePoints=32;
        limits.maximumMergedComponents=8;
        limits.maximumPatches=256;
        limits.patchSearch={32768,32768,512,4};
        limits.maximumLocalFrontTets=256;
        if(options.cancellation)limits.interiorAdmission=[&](Vec3){charge(0);return true;};
        // A new snapshot can coexist with the previous progress while its
        // flat buffers grow. Leave room for both plus the index remap.
        const auto progressPayload=limits.maximumNodes*(sizeof(Vec3)+sizeof(mesh::Classification))+
            maxTets*(sizeof(mesh::Cell)+4*sizeof(Index));
        const auto progressReserve=6*progressPayload+limits.maximumNodes*(sizeof(Index)+sizeof(unsigned char));
        const bool adaptive=recovery.useAdaptiveRefill&&progressReserve<limits.maximumWorkingBytes;
        if(adaptive)limits.maximumWorkingBytes-=progressReserve;
        mesh::Mesh complete,progress;
        stage.next("cavity_refill");
        auto closed=adaptive
            ?bnd::recoverTetPlcCavitiesRetaining(source,seed,limits,complete,&progress)
            :bnd::recoverTetPlcCavities(source,seed,limits,complete);
        stage.finish(closed.diagnostic);
        charge(closed.workUsed);
        using C=bnd::TetPlcCavityReason;
        const auto refusal=std::string_view(closed.diagnostic);
        const bool exhaustedRegion=refusal=="no-admissible-region-refill";
        const bool interruptedPatch=progress.cellCount()>0&&
            (refusal=="work-or-patch-budget"||refusal=="patch-count-budget");
        if(adaptive&&!closed.published&&out.workUsed<cap&&
           (closed.reason==C::Unavailable||closed.reason==C::BudgetExceeded)&&
           (exhaustedRegion||interruptedPatch)){
            if(progress.cellCount()>0){
                seed.swap(progress);
                mesh::Mesh{}.swap(progress);
            }
            // Reopen a few more layers around the unresolved front. The
            // region solver treats source triangles as barriers on both passes.
            limits.maximumWork=cap-out.workUsed;
            limits.maximumExpansionCells=384;
            limits.patchSearch.maximumShells=8;
            if(closed.regionStateBoundHits)limits.patchSearch.maximumStates=65536;
            limits.useVisibilityExpansion=false;
            limits.maximumLocalFrontTets=64;
            stage.next("cavity_refill");
            closed=bnd::recoverTetPlcCavities(source,seed,limits,complete);
            stage.finish(closed.diagnostic);
            charge(closed.workUsed);
        }
        if(!closed.published){
            const auto reason=closed.reason==C::AllocationFailed?R::AllocationFailed:
                closed.reason==C::Cancelled?R::Cancelled:
                closed.reason==C::BudgetExceeded?R::BudgetExceeded:R::UnsupportedGeometry;
            throw Stop{reason,closed.diagnostic};
        }
        const auto bytes=std::uint64_t(complete.nodeCount())*(sizeof(Vec3)+sizeof(mesh::Classification))+
            std::uint64_t(complete.cellCount())*sizeof(mesh::Cell)+std::uint64_t(complete.connectivityCount())*sizeof(Index);
        require(bytes<=options.maximumOutputBytes,R::BudgetExceeded,"live-piece-output-bytes");
        destination.swap(complete);
        out.published=true;
        out.piecesRecovered=1;
        out.failedPieceTag=-1;
        out.reason=R::None;
        out.diagnostic="none";
    }catch(const Stop& failure){out.reason=failure.reason;out.diagnostic=failure.diagnostic;}
     catch(const std::bad_alloc&){out.reason=R::AllocationFailed;out.diagnostic="allocation-failed";}
     catch(...){out.reason=R::InternalError;out.diagnostic="live-piece-internal-error";}
    return out;
}
OcctPlcRecoveryOutcome recoverOcctPlcFrontVolume(const mesh::Mesh& source,int pieceTag,
    const OcctModel& pieces,const OcctModel& original,int originalBody,
    const OcctStarOptions& options,const OcctPartitionRecoveryOptions& recovery,
    mesh::Mesh& destination) noexcept {
    OcctPlcRecoveryOutcome out;out.attempted=true;out.piecesAttempted=1;out.failedPieceTag=pieceTag;
    std::stop_source stop;std::optional<Stop> predicateFailure;
    using R=OcctPlcRecoveryReason;
    try {
        const auto cap=std::min(options.maximumWork,recovery.maximumPieceWork);
        const auto charge=[&](std::uint64_t count){
            if(options.cancellation&&options.cancellation->load(std::memory_order_relaxed)){
                stop.request_stop();throw Stop{R::Cancelled,"cancelled"};
            }
            require(out.workUsed<=cap&&count<=cap-out.workUsed,R::BudgetExceeded,"front-piece-work-budget");
            out.workUsed+=count;
        };
        charge(0);
        require(&source!=&destination&&cap&&pieces.solid(pieceTag)&&original.solid(originalBody),
                R::InvalidRequest,"front-piece-invalid-request");
        require(std::isfinite(recovery.maximumFrontSeconds)&&recovery.maximumFrontSeconds>0&&
                recovery.maximumFrontSeconds<=30.,R::InvalidRequest,"front-piece-invalid-time-limit");
        require(source.nodeCount()>=4&&source.cellCount()>=4,R::InvalidRequest,"front-piece-empty-source");
        require(source.nodeCount()<=512&&source.cellCount()<=4096,R::UnsupportedGeometry,"front-piece-count-limit");
        require(std::uint64_t(source.nodeCount())<=options.maximumNodes&&
                std::uint64_t(source.cellCount())<options.maximumCells,R::BudgetExceeded,"front-piece-output-budget");
        // This route deliberately stays local. Larger native pieces need more
        // subdivision, not unbounded apex arrays or a global workspace copy.
        const auto retained=std::uint64_t(source.nodeCount()+source.cellCount())*4096+16ULL*1024*1024;
        require(retained<options.maximumWorkingBytes,R::BudgetExceeded,"front-piece-working-budget");
        auto certificateOptions=options;certificateOptions.maximumWork=cap-out.workUsed;
        const auto certificate=certifyOcctBoundarySurface(source,certificateOptions);
        charge(certificate.workUsed);
        if(!certificate.certified){
            const auto code=certificate.reason==OcctStarReason::Cancelled?R::Cancelled:
                certificate.reason==OcctStarReason::AllocationFailed?R::AllocationFailed:
                certificate.reason==OcctStarReason::WorkBudget||certificate.reason==OcctStarReason::OutputBudget?R::BudgetExceeded:R::UnsupportedGeometry;
            throw Stop{code,certificate.diagnostic};
        }
        OcctSolidClassifier pieceClassifier(OcctModelAccess::solid(pieces,pieceTag));
        OcctSolidClassifier originalClassifier(OcctModelAccess::solid(original,originalBody));
        const auto pieceTolerance=pieces.tolerance(3,pieceTag),originalTolerance=original.tolerance(3,originalBody);
        bool inCavity=false;std::uint64_t predicateAllowance=0;
        const auto inside=[&](Vec3 point){
            try {
                if(inCavity){require(predicateAllowance>=2,R::BudgetExceeded,"front-piece-native-query-budget");predicateAllowance-=2;}
                charge(2);++out.nativeQueries;
                pieceClassifier.Perform(gp_Pnt(point.x,point.y,point.z),pieceTolerance);
                require(pieceClassifier.State()!=TopAbs_UNKNOWN,R::NativeQueryFailed,"front-piece-native-query-failed");
                if(pieceClassifier.State()!=TopAbs_IN&&pieceClassifier.State()!=TopAbs_ON)return false;
                ++out.nativeQueries;originalClassifier.Perform(gp_Pnt(point.x,point.y,point.z),originalTolerance);
                require(originalClassifier.State()!=TopAbs_UNKNOWN,R::NativeQueryFailed,"front-original-native-query-failed");
                return originalClassifier.State()==TopAbs_IN||originalClassifier.State()==TopAbs_ON;
            }catch(const Stop& failure){predicateFailure=failure;throw;}
             catch(const Standard_OutOfMemory&){predicateFailure=Stop{R::AllocationFailed,"allocation-failed"};throw *predicateFailure;}
             catch(const std::bad_alloc&){predicateFailure=Stop{R::AllocationFailed,"allocation-failed"};throw *predicateFailure;}
             catch(const Standard_Failure&){predicateFailure=Stop{R::NativeQueryFailed,"native-query-failed"};throw *predicateFailure;}
        };
        const auto maximumTets=std::min<std::uint64_t>(10000,options.maximumCells-source.cellCount());
        mesh::Mesh partial;
        {
            bnd::detail::PlcGiftWrapping front;front.chargeWork=charge;front.interiorAdmission=inside;
            for(Index i=0;i<source.nodeCount();++i)front.vertices.push_back({source.node(i),std::uint64_t(i)+1});
            for(Index i=0;i<source.cellCount();++i){const auto* f=source.cellNodes(i);front.boundary.push_back({unsigned(f[0]),unsigned(f[1]),unsigned(f[2])});}
            require(front.initialize(),R::UnsupportedGeometry,"front-piece-not-closed-oriented");
            front.seedUnfilledFront(unsigned(std::min<std::uint64_t>(16,options.maximumNodes-source.nodeCount())));
            front.run(unsigned(std::min<std::uint64_t>(4096,maximumTets)),recovery.maximumFrontSeconds,stop.get_token(),50000);
            require(front.vertices.size()<=options.maximumNodes,R::BudgetExceeded,"front-piece-node-budget");
            for(Index i=0;i<source.nodeCount();++i)partial.addNode(source.node(i),source.nodeClass(i));
            for(std::size_t i=std::size_t(source.nodeCount());i<front.vertices.size();++i)
                partial.addNode(front.vertices[i].position,{3,pieceTag});
            for(auto cell:front.cells){std::array<Index,4> ids{Index(cell.vertices[0]),Index(cell.vertices[1]),Index(cell.vertices[2]),Index(cell.vertices[3])};partial.addCell(mesh::CellType::Tet4,ids,{3,pieceTag});}
        }
        require(out.workUsed<cap,R::BudgetExceeded,"front-piece-work-budget");
        bnd::TetPlcCavityOptions limits;limits.region={3,pieceTag};
        limits.maximumNodes=std::min<std::uint64_t>(4096,options.maximumNodes);limits.maximumCells=maximumTets;
        limits.maximumSourceTriangles=4096;limits.maximumWorkingBytes=options.maximumWorkingBytes-retained;
        // The PLC solver accounts for its own geometric work. Reserve native
        // callback work separately, so the combined receipt cannot overrun cap.
        predicateAllowance=(cap-out.workUsed)/4;limits.maximumWork=cap-out.workUsed-predicateAllowance;
        require(predicateAllowance>=2&&limits.maximumWork,R::BudgetExceeded,"front-piece-cavity-budget");
        limits.maximumComponentFaces=4096;limits.maximumPatchCells=2048;limits.maximumExpansionCells=256;
        limits.maximumCandidatePoints=64;limits.maximumRandomCandidatePoints=32;limits.maximumMergedComponents=8;
        limits.maximumPatches=256;limits.patchSearch={32768,32768,512,4};limits.cancellation=stop.get_token();limits.interiorAdmission=inside;
        mesh::Mesh complete;inCavity=true;
        const auto closed=bnd::recoverTetPlcCavities(source,partial,limits,complete);
        inCavity=false;charge(closed.workUsed);
        if(predicateFailure)throw *predicateFailure;
        if(!closed.published){
            const auto code=closed.reason==bnd::TetPlcCavityReason::Cancelled?R::Cancelled:
                closed.reason==bnd::TetPlcCavityReason::AllocationFailed?R::AllocationFailed:
                closed.reason==bnd::TetPlcCavityReason::BudgetExceeded?R::BudgetExceeded:R::UnsupportedGeometry;
            throw Stop{code,closed.diagnostic};
        }
        for(Index i=source.nodeCount();i<complete.nodeCount();++i)
            require(inside(complete.node(i)),R::InvalidMesh,"front-piece-native-node-outside");
        for(Index i=0;i<complete.cellCount();++i)if(complete.cell(i).type==mesh::CellType::Tet4){const auto* f=complete.cellNodes(i);
            require(inside((complete.node(f[0])+complete.node(f[1])+complete.node(f[2])+complete.node(f[3]))*.25),R::InvalidMesh,"front-piece-native-centroid-outside");}
        require(std::uint64_t(complete.nodeCount())<=options.maximumNodes&&std::uint64_t(complete.cellCount())<=options.maximumCells,
                R::BudgetExceeded,"front-piece-output-budget");
        const auto payload=std::uint64_t(complete.nodeCount())*(sizeof(Vec3)+sizeof(mesh::Classification))+
            std::uint64_t(complete.cellCount())*sizeof(mesh::Cell)+std::uint64_t(complete.connectivityCount())*sizeof(Index);
        require(payload<=options.maximumOutputBytes,R::BudgetExceeded,"front-piece-output-bytes");
        charge(0);destination.swap(complete);out.published=true;out.piecesRecovered=1;out.reason=R::None;out.diagnostic="none";
    }catch(const Stop& failure){out.reason=failure.reason;out.diagnostic=failure.diagnostic;}
     catch(const Standard_OutOfMemory&){out.reason=R::AllocationFailed;out.diagnostic="allocation-failed";}
     catch(const std::bad_alloc&){out.reason=R::AllocationFailed;out.diagnostic="allocation-failed";}
     catch(const Standard_Failure&){out.reason=R::NativeQueryFailed;out.diagnostic="native-query-failed";}
     catch(...){out.reason=R::InternalError;out.diagnostic="front-piece-internal-error";}
    return out;
}
}
