// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctFacetRefinement.hpp"
#include <Standard_OutOfMemory.hxx>
#include <Standard_Failure.hxx>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <compare>
#include <map>
#include <new>
#include <set>

namespace paacs::geo::detail {
namespace {
struct Stop { const char* diagnostic; OcctFacetRefinementReason reason; };
void need(bool b,const char* why,OcctFacetRefinementReason r=OcctFacetRefinementReason::UnsupportedGeometry) {
    if(!b)throw Stop{why,r};
}
bool finite(Vec3 p){return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);}
using Point=std::array<std::uint64_t,3>;
Point bits(Vec3 p){return {std::bit_cast<std::uint64_t>(p.x),std::bit_cast<std::uint64_t>(p.y),std::bit_cast<std::uint64_t>(p.z)};}
struct Key {
    int face=-1;std::array<Point,3> points{};
    auto operator<=>(const Key&)const=default;
};
Key key(const mesh::Mesh& s,Index ci){
    Key k;k.face=s.cell(ci).onEntity.tag;const auto* ids=s.cellNodes(ci);
    for(unsigned j=0;j<3;++j)k.points[j]=bits(s.node(ids[j]));
    std::sort(k.points.begin(),k.points.end());return k;
}
Vec3 point(Point p){return {std::bit_cast<double>(p[0]),std::bit_cast<double>(p[1]),std::bit_cast<double>(p[2])};}
struct Patch { Vec3 inserted{};std::set<int> owners,seen; };
}
OcctFacetRefinementOutcome refineOcctFacets(const OcctModel& model,
    std::span<const OcctFacetBoundary> inputs,const OcctFacetRefinementOptions& options) noexcept {
    OcctFacetRefinementOutcome out;
    const auto tick=[&](std::uint64_t n=1){
        if(options.limits.cancellation&&options.limits.cancellation->load())
            throw Stop{"cancelled",OcctFacetRefinementReason::Cancelled};
        need(n<=options.limits.maximumWork-out.workUsed,"work-budget",OcctFacetRefinementReason::BudgetExceeded);out.workUsed+=n;
    };
    const auto certificate=[&](const mesh::Mesh& s){
        auto limits=options.limits;limits.maximumWork-=out.workUsed;
        const auto c=certifyOcctBoundarySurface(s,limits);tick(c.workUsed);
        if(!c.certified)throw Stop{c.diagnostic,
            c.reason==OcctStarReason::Cancelled?OcctFacetRefinementReason::Cancelled:
            c.reason==OcctStarReason::WorkBudget||c.reason==OcctStarReason::OutputBudget?OcctFacetRefinementReason::BudgetExceeded:
            c.reason==OcctStarReason::AllocationFailed?OcctFacetRefinementReason::AllocationFailed:OcctFacetRefinementReason::InvalidBoundary};
    };
    try {
        tick(0);
        need(!inputs.empty()&&inputs.size()<=64&&!options.selections.empty()&&options.selections.size()<=1024&&
             options.maximumCandidatePatches<=1024&&
             std::isfinite(options.maximumProjectionDistance)&&options.maximumProjectionDistance>0,
             "invalid-facet-options",OcctFacetRefinementReason::InvalidRequest);
        std::map<int,const mesh::Mesh*> sources;std::set<const mesh::Mesh*> destinations;
        std::uint64_t estimated=0;
        for(const auto& input:inputs){tick();
            need(input.boundary&&model.solid(input.bodyTag)&&sources.emplace(input.bodyTag,input.boundary).second&&
                 destinations.insert(input.boundary).second,"invalid-or-duplicate-body-boundary",OcctFacetRefinementReason::InvalidRequest);
            const auto& s=*input.boundary;
            need(std::uint64_t(s.nodeCount())<=options.limits.maximumBoundaryNodes&&std::uint64_t(s.cellCount())<=options.limits.maximumBoundaryTriangles,
                 "input-count-budget",OcctFacetRefinementReason::BudgetExceeded);
            const auto rows=std::uint64_t(s.nodeCount())+std::uint64_t(s.cellCount());
            need(rows<=(UINT64_MAX-estimated)/512,"working-overflow",OcctFacetRefinementReason::BudgetExceeded);estimated+=rows*512;
        }
        need(options.maximumAddedNodes<=(UINT64_MAX-estimated)/1024,"working-overflow",OcctFacetRefinementReason::BudgetExceeded);
        estimated+=options.maximumAddedNodes*1024;
        need(estimated<=options.limits.maximumWorkingBytes,"working-payload-budget",OcctFacetRefinementReason::BudgetExceeded);
        for(const auto& input:inputs)certificate(*input.boundary);
        std::map<Key,Patch> patches;
        for(const auto& selected:options.selections){tick();
            const auto found=sources.find(selected.bodyTag);
            need(found!=sources.end()&&selected.triangle>=0&&selected.triangle<found->second->cellCount(),"invalid-facet-selection",OcctFacetRefinementReason::InvalidRequest);
            const auto& s=*found->second;const auto& c=s.cell(selected.triangle);
            need(c.type==mesh::CellType::Tri3&&c.onEntity.dim==2&&model.face(c.onEntity.tag),"invalid-native-facet",OcctFacetRefinementReason::InvalidBoundary);
            const auto k=key(s,selected.triangle);if(patches.contains(k))continue;
            out.failedFaceTag=k.face;out.failedTriangle=selected.triangle;
            Patch patch;
            for(const auto& owner:model.faceOwners(k.face)){tick();
                need(sources.contains(owner.solidTag),"missing-shared-face-coowner");patch.owners.insert(owner.solidTag);}
            need(patch.owners.contains(selected.bodyTag),"facet-not-owned-by-body",OcctFacetRefinementReason::InvalidBoundary);
            // Coordinate bits provide a common summation order across opposite
            // face uses and differing body-local node indices.
            const auto* face=model.face(k.face);
            if(options.requireGeometricChange&&face->kind()==SurfaceKind::Plane){
                ++out.unchangedCandidates;continue;
            }
            try{
                const auto a=point(k.points[0]),b=point(k.points[1]),cpt=point(k.points[2]);
                const auto center=(a+b+cpt)/3.;Vec2 uv;Vec3 projected;
                tick();++out.nativeQueries;
                need(face->projectPoint(center,uv,projected)&&face->containsParam(uv)&&!face->singularAt(uv)&&finite(projected),"native-trim-projection-failed");
                const auto displacement=(projected-center).norm();
                need(std::isfinite(displacement)&&displacement<=options.maximumProjectionDistance,"projection-displacement-limit");
                if(options.requireGeometricChange&&displacement<=model.tolerance(2,k.face)){
                    ++out.unchangedCandidates;continue;
                }
                // Re-evaluate the accepted UV on the native surface, and require
                // agreement at its original raw face tolerance.
                tick();++out.nativeQueries;
                need((face->pointAt(uv)-projected).norm()<=model.tolerance(2,k.face),"native-evaluation-disagreement");
                const auto v=b-a,w=cpt-a,d=projected-a;
                const auto vv=v.dot(v),vw=v.dot(w),ww=w.dot(w),dv=d.dot(v),dw=d.dot(w),denom=vv*ww-vw*vw;
                need(std::isfinite(denom)&&denom>0,"projected-point-outside-source-triangle");
                const auto u=(ww*dv-vw*dw)/denom,t=(vv*dw-vw*dv)/denom;
                need(std::isfinite(u)&&std::isfinite(t)&&u>0&&t>0&&u+t<1,
                     "projected-point-outside-source-triangle");
                patch.inserted=projected;out.maximumDisplacement=std::max(out.maximumDisplacement,displacement);
            }catch(const Stop& e){
                if(!options.skipUnavailableCandidates||e.reason!=OcctFacetRefinementReason::UnsupportedGeometry)throw;
                ++out.unavailableCandidates;continue;
            }
            patches.emplace(k,std::move(patch));
            if(options.maximumCandidatePatches&&patches.size()>=options.maximumCandidatePatches)break;
        }
        need(!patches.empty(),"no-usable-facet-candidates");
        std::vector<bnd::detail::FacetSourcePatch> committed;
        committed.reserve(static_cast<std::size_t>(std::min<std::uint64_t>(options.maximumAddedNodes,options.selections.size()*inputs.size())));
        std::vector<mesh::Mesh> rebuilt;rebuilt.reserve(inputs.size());
        for(const auto& input:inputs){const auto& s=*input.boundary;mesh::Mesh candidate;std::set<int> bodyFaces;
            for(const int face:model.solid(input.bodyTag)->boundaryFaceTags()){tick();bodyFaces.insert(face);}
            for(Index i=0;i<s.nodeCount();++i){tick();candidate.addNode(s.node(i),s.nodeClass(i));}
            for(Index ci=0;ci<s.cellCount();++ci){tick();const auto& cell=s.cell(ci);const auto* ids=s.cellNodes(ci);
                need(cell.type==mesh::CellType::Tri3&&cell.onEntity.dim==2&&bodyFaces.contains(cell.onEntity.tag),"invalid-native-face-ownership",OcctFacetRefinementReason::InvalidBoundary);
                auto found=patches.find(key(s,ci));
                if(found==patches.end()){candidate.addCell(mesh::CellType::Tri3,std::span<const Index>(ids,3),cell.onEntity);continue;}
                auto& patch=found->second;
                need(patch.owners.contains(input.bodyTag)&&patch.seen.insert(input.bodyTag).second,"duplicate-or-unowned-shared-facet",OcctFacetRefinementReason::InvalidBoundary);
                need(out.addedNodes<options.maximumAddedNodes,"added-node-budget",OcctFacetRefinementReason::BudgetExceeded);
                const auto newNode=candidate.addNode(patch.inserted,{2,cell.onEntity.tag});++out.addedNodes;
                committed.push_back({input.bodyTag,cell.onEntity.tag,ci,{ids[0],ids[1],ids[2]},newNode});
                for(unsigned j=0;j<3;++j){tick();const std::array triangle{ids[j],ids[(j+1)%3],newNode};candidate.addCell(mesh::CellType::Tri3,triangle,cell.onEntity);}
                out.addedTriangles+=2;
            }
            need(std::uint64_t(candidate.nodeCount())<=options.limits.maximumBoundaryNodes&&std::uint64_t(candidate.cellCount())<=options.limits.maximumBoundaryTriangles,
                 "output-count-budget",OcctFacetRefinementReason::BudgetExceeded);
            const auto payload=std::uint64_t(candidate.nodeCount())*(sizeof(Vec3)+sizeof(mesh::Classification))+
                std::uint64_t(candidate.cellCount())*sizeof(mesh::Cell)+std::uint64_t(candidate.connectivityCount())*sizeof(Index);
            need(payload<=options.limits.maximumOutputBytes,"output-payload-budget",OcctFacetRefinementReason::BudgetExceeded);
            certificate(candidate);rebuilt.push_back(std::move(candidate));
        }
        for(const auto& [key,patch]:patches){tick();need(patch.seen==patch.owners,"missing-bitwise-shared-facet");}
        tick(0);for(std::size_t i=0;i<inputs.size();++i)inputs[i].boundary->swap(rebuilt[i]);
        out.committedPatches.swap(committed);
        out.published=true;out.reason=OcctFacetRefinementReason::None;out.diagnostic="none";
        out.failedFaceTag=-1;out.failedTriangle=-1;
    }catch(const Stop& e){out.reason=e.reason;out.diagnostic=e.diagnostic;}
    catch(const Standard_OutOfMemory&){out.reason=OcctFacetRefinementReason::AllocationFailed;out.diagnostic="allocation-failed";}
    catch(const std::bad_alloc&){out.reason=OcctFacetRefinementReason::AllocationFailed;out.diagnostic="allocation-failed";}
    catch(const Standard_Failure&){out.reason=OcctFacetRefinementReason::NativeQueryFailed;out.diagnostic="native-query-failed";}
    catch(...){out.reason=OcctFacetRefinementReason::InternalError;out.diagnostic="internal-error";}
    return out;
}
} // namespace paacs::geo::detail
