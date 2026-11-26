// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Diagnostic ordering only. Native-rejected points are never admitted here.
#include "PlcVoidGeometry.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <atomic>
#include <bit>
#include <climits>
#include <map>
#include <new>
#include <span>
#include <set>

namespace paacs::geo::detail {
struct OcctSourceMismatchStamp {
    int bodyTag=-1;
    std::uint64_t generation=0,sourceMeshHash=0;
    friend bool operator==(const OcctSourceMismatchStamp&,const OcctSourceMismatchStamp&)=default;
};
struct OcctSourceMismatchLimits {
    std::uint64_t maximumWork=2'000'000,maximumWorkingBytes=64ULL*1024*1024;
    std::uint64_t maximumNodes=1'000'000,maximumSourceTriangles=1'000'000;
    std::uint64_t maximumRejectedPoints=65536,maximumCandidates=65536;
    const std::atomic<bool>* cancellation=nullptr;
};
enum class OcctSourceMismatchReason {None,InvalidInput,BudgetExceeded,Cancelled,AllocationFailed,InternalError};
struct OcctSourceMismatchFacetPriority {
    Index triangle=-1;
    std::uint64_t insideSourceWitnesses=0;
    double nearestDistanceSquared=std::numeric_limits<double>::infinity();
};
struct OcctSourceMismatchOutcome {
    bool complete=false;
    OcctSourceMismatchReason reason=OcctSourceMismatchReason::InternalError;
    const char* diagnostic="not-attempted";
    OcctSourceMismatchStamp source;
    std::uint64_t workUsed=0,modeledBytes=0,nativeRejectedPoints=0;
    std::uint64_t insideSourceOutsideCad=0,outsideSource=0,ambiguousSource=0;
    std::uint64_t prioritizedFacets=0,insideWitnessesNearMissingFacets=0;
    // Every supplied candidate appears once. Supported facets precede the
    // unchanged ascending-ID fallback; this receipt certifies no mesh or CAD.
    std::vector<OcctSourceMismatchFacetPriority> facets;
};
namespace source_mismatch_detail {
struct Stop {OcctSourceMismatchReason reason;const char* diagnostic;};
using Front=bnd::detail::PlcGiftWrapping;
using Face=Front::Face;
inline bool finite(Vec3 p){return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);}
inline double coordinate(Vec3 p,unsigned axis){return axis==0?p.x:axis==1?p.y:p.z;}
// Ranking metric only. Outward rounding retains near ties in the BVH.
inline double distanceLower(Vec3 p,Front::Box box){
    double result=0;
    for(unsigned axis=0;axis<3;++axis){
        const auto x=coordinate(p,axis),lo=coordinate(box.lo,axis),hi=coordinate(box.hi,axis);
        if(x<lo)result=std::max(result,std::nextafter(lo-x,0.));
        if(x>hi)result=std::max(result,std::nextafter(x-hi,0.));
    }
    return result;
}
}

// Preconditions: source is the independently certified closed/embedded source
// of current.bodyTag; rejected points were refused by ORIGINAL CAD for that
// body. The helper verifies the current source hash and generation binding.
// Exact PLC parity only diagnoses the mismatch. Nearest distances only order
// candidates; original Native projection, trim, source/contact and final volume
// gates remain mandatory. No cache or point classification survives this call.
inline OcctSourceMismatchOutcome rankOcctSourceMismatchFacets(const mesh::Mesh& source,
    OcctSourceMismatchStamp current,OcctSourceMismatchStamp expected,
    std::span<const Vec3> rejected,std::span<const Index> candidates,
    const OcctSourceMismatchLimits& limits={}) noexcept {
    using namespace source_mismatch_detail;
    OcctSourceMismatchOutcome out;out.source=current;
    const auto need=[](bool ok,OcctSourceMismatchReason why,const char* diagnostic){if(!ok)throw Stop{why,diagnostic};};
    const auto tick=[&](std::uint64_t count=1){
        if(limits.cancellation&&limits.cancellation->load(std::memory_order_relaxed))
            throw Stop{OcctSourceMismatchReason::Cancelled,"cancelled"};
        need(out.workUsed<=limits.maximumWork&&count<=limits.maximumWork-out.workUsed,
            OcctSourceMismatchReason::BudgetExceeded,"source-mismatch-work-budget");
        out.workUsed+=count;
    };
    try {
        tick(0);
        need(current==expected&&current.bodyTag>=0,OcctSourceMismatchReason::InvalidInput,"source-mismatch-stale-identity");
        const auto nodes=std::uint64_t(source.nodeCount()),faces=std::uint64_t(source.cellCount());
        need(nodes>=4&&faces>=4&&nodes<bnd::kInvalidTetFacetVertex&&faces<=UINT_MAX/2,
            OcctSourceMismatchReason::InvalidInput,"source-mismatch-source-count");
        need(nodes<=limits.maximumNodes&&faces<=limits.maximumSourceTriangles&&
            rejected.size()<=limits.maximumRejectedPoints&&candidates.size()<=limits.maximumCandidates,
            OcctSourceMismatchReason::BudgetExceeded,"source-mismatch-count-budget");
        // Includes source BVH, unique-face/edge checks, returned ranks and sort
        // scratch. Points are borrowed. No complete volume mesh is copied.
        need(nodes<=UINT64_MAX/128&&faces<=(UINT64_MAX-nodes*128-65536)/1024,
            OcctSourceMismatchReason::BudgetExceeded,"source-mismatch-byte-overflow");
        out.modeledBytes=nodes*128+faces*1024+65536;
        need(candidates.size()<=(UINT64_MAX-out.modeledBytes)/128,
            OcctSourceMismatchReason::BudgetExceeded,"source-mismatch-byte-overflow");
        out.modeledBytes+=std::uint64_t(candidates.size())*128;
        need(out.modeledBytes<=limits.maximumWorkingBytes,OcctSourceMismatchReason::BudgetExceeded,"source-mismatch-working-budget");
        // The hash is a bounded identity screen, not a source certificate.
        tick(nodes*8+faces*12);
        need(mesh::meshHash(source)==current.sourceMeshHash,OcctSourceMismatchReason::InvalidInput,"source-mismatch-source-changed");
        Front index;index.vertices.reserve(nodes);index.boundary.reserve(faces);index.boundaryBoxes.reserve(faces);
        index.order.resize(faces);index.tree.reserve(faces*2);
        for(Index i=0;i<source.nodeCount();++i){tick();
            const auto p=source.node(i);const auto entity=source.nodeClass(i);
            need(finite(p)&&entity.dim>=0&&entity.dim<=2&&entity.tag>=0,
                OcctSourceMismatchReason::InvalidInput,"source-mismatch-source-node");
            index.vertices.push_back({p,std::uint64_t(i)+1});
        }
        std::set<Face> unique;
        std::map<std::array<Index,2>,std::pair<unsigned,int>> edges;
        const auto lookupWork=2*std::bit_width(faces)+1;
        for(Index i=0;i<source.cellCount();++i){tick(lookupWork+4);
            const auto& cell=source.cell(i);const auto* ids=source.cellNodes(i);
            need(cell.type==mesh::CellType::Tri3&&cell.onEntity.dim==2&&cell.onEntity.tag>=0,
                OcctSourceMismatchReason::InvalidInput,"source-mismatch-source-face");
            Face f{};for(unsigned j=0;j<3;++j){
                need(ids[j]>=0&&std::uint64_t(ids[j])<nodes,OcctSourceMismatchReason::InvalidInput,"source-mismatch-source-index");
                f[j]=static_cast<bnd::TetFacetVertexId>(ids[j]);
            }
            need(f[0]!=f[1]&&f[1]!=f[2]&&f[0]!=f[2]&&unique.insert(Front::key(f)).second&&
                bnd::tet_facet_detail::projectionForTriangle(source.node(ids[0]),source.node(ids[1]),source.node(ids[2]))!=
                    bnd::tet_facet_detail::Projection::Invalid,
                OcctSourceMismatchReason::InvalidInput,"source-mismatch-degenerate-source");
            for(unsigned j=0;j<3;++j){tick(lookupWork);auto a=ids[j],b=ids[(j+1)%3];
                auto& edge=edges[{std::min(a,b),std::max(a,b)}];++edge.first;edge.second+=a<b?1:-1;
            }
            index.boundary.push_back(f);index.boundaryBoxes.push_back(index.box(f));index.order[i]=static_cast<unsigned>(i);
        }
        for(const auto& [edge,use]:edges){(void)edge;tick();
            need(use.first==2&&use.second==0,OcctSourceMismatchReason::InvalidInput,"source-mismatch-open-source");
        }
        // Do not call Front::initialize(): this query owns no advancing front,
        // and its BVH construction must expose every scan/comparison to charge.
        const auto build=[&](auto&& self,unsigned first,unsigned last)->unsigned {
            tick();const auto id=static_cast<unsigned>(index.tree.size());index.tree.push_back({});
            auto box=index.boundaryBoxes[index.order[first]];
            for(auto i=first+1;i<last;++i){tick();const auto other=index.boundaryBoxes[index.order[i]];
                box.lo={std::min(box.lo.x,other.lo.x),std::min(box.lo.y,other.lo.y),std::min(box.lo.z,other.lo.z)};
                box.hi={std::max(box.hi.x,other.hi.x),std::max(box.hi.y,other.hi.y),std::max(box.hi.z,other.hi.z)};
            }
            index.tree[id]={box,first,last,0,0};
            if(last-first>8){
                const auto span=box.hi-box.lo;const unsigned axis=span.x>=span.y&&span.x>=span.z?0:span.y>=span.z?1:2;
                const auto middle=first+(last-first)/2;
                std::nth_element(index.order.begin()+first,index.order.begin()+middle,index.order.begin()+last,[&](unsigned a,unsigned b){
                    tick();const auto x=index.boundaryBoxes[a],y=index.boundaryBoxes[b];
                    const auto ca=coordinate(x.lo,axis)*.5+coordinate(x.hi,axis)*.5;
                    const auto cb=coordinate(y.lo,axis)*.5+coordinate(y.hi,axis)*.5;
                    return std::pair{ca,a}<std::pair{cb,b};
                });
                const auto left=self(self,first,middle),right=self(self,middle,last);
                index.tree[id].left=left;index.tree[id].right=right;
            }return id;
        };
        build(build,0,static_cast<unsigned>(faces));index.initialized=true;
        std::vector<Index> candidateRow(faces,kInvalidIndex);
        out.facets.reserve(candidates.size());
        for(auto triangle:candidates){tick();
            need(triangle>=0&&std::uint64_t(triangle)<faces&&candidateRow[triangle]==kInvalidIndex,
                OcctSourceMismatchReason::InvalidInput,"source-mismatch-candidate-index");
            candidateRow[triangle]=static_cast<Index>(out.facets.size());out.facets.push_back({triangle});
        }
        PlcVoidGeometry geometry(index,tick);
        for(auto p:rejected){tick();++out.nativeRejectedPoints;
            need(finite(p),OcctSourceMismatchReason::InvalidInput,"source-mismatch-point");
            const auto inclusion=geometry.inside(p);
            if(inclusion==PlcVoidGeometry::Inclusion::Outside){++out.outsideSource;continue;}
            if(inclusion==PlcVoidGeometry::Inclusion::Ambiguous){++out.ambiguousSource;continue;}
            ++out.insideSourceOutsideCad;
            std::array<std::pair<double,unsigned>,3> nearest{};unsigned found=0;
            std::array<unsigned,128> stack{};std::size_t count=1,visits=0;stack[0]=0;
            while(count){tick();const auto branch=stack[--count];
                need(branch<index.tree.size()&&++visits<=index.tree.size(),OcctSourceMismatchReason::InternalError,"source-mismatch-tree");
                const auto& node=index.tree[branch];
                if(found==nearest.size()&&distanceLower(p,node.bounds)>std::nextafter(std::sqrt(nearest[found-1].first),
                    std::numeric_limits<double>::infinity()))continue;
                if(node.left||node.right){
                    need(node.left&&node.right&&count+2<=stack.size(),OcctSourceMismatchReason::InternalError,"source-mismatch-stack");
                    auto near=node.left,far=node.right;
                    if(distanceLower(p,index.tree[far].bounds)<distanceLower(p,index.tree[near].bounds))std::swap(near,far);
                    stack[count++]=far;stack[count++]=near;
                }else for(unsigned j=node.first;j<node.last;++j){tick(16);const auto row=index.order[j];const auto f=index.boundary[row];
                    const auto squared=pointTriangleDistanceSquared(p,index.vertices[f[0]].position,index.vertices[f[1]].position,index.vertices[f[2]].position);
                    if(!(squared>=0)||!std::isfinite(squared))continue;
                    const auto key=std::pair{squared,row};
                    if(found<nearest.size())nearest[found++]=key;
                    else if(key<nearest[found-1])nearest[found-1]=key;else continue;
                    std::sort(nearest.begin(),nearest.begin()+found);
                }
            }
            bool supported=false;
            for(unsigned j=0;j<found;++j){tick();const auto row=candidateRow[nearest[j].second];if(row==kInvalidIndex)continue;
                auto& candidate=out.facets[row];++candidate.insideSourceWitnesses;
                candidate.nearestDistanceSquared=std::min(candidate.nearestDistanceSquared,nearest[j].first);supported=true;
            }
            out.insideWitnessesNearMissingFacets+=supported;
        }
        std::sort(out.facets.begin(),out.facets.end(),[&](const auto& a,const auto& b){tick();
            if(a.insideSourceWitnesses!=b.insideSourceWitnesses)return a.insideSourceWitnesses>b.insideSourceWitnesses;
            if(a.insideSourceWitnesses&&a.nearestDistanceSquared!=b.nearestDistanceSquared)
                return a.nearestDistanceSquared<b.nearestDistanceSquared;
            return a.triangle<b.triangle;
        });
        for(const auto& candidate:out.facets){tick();out.prioritizedFacets+=candidate.insideSourceWitnesses!=0;}
        tick(0);out.complete=true;out.reason=OcctSourceMismatchReason::None;out.diagnostic="source-mismatch-priority-only";
    }catch(const Stop& stop){out.reason=stop.reason;out.diagnostic=stop.diagnostic;}
     catch(const std::bad_alloc&){out.reason=OcctSourceMismatchReason::AllocationFailed;out.diagnostic="allocation-failed";}
     catch(...){out.reason=OcctSourceMismatchReason::InternalError;out.diagnostic="source-mismatch-internal";}
    if(!out.complete){std::vector<OcctSourceMismatchFacetPriority>{}.swap(out.facets);out.prioritizedFacets=0;}
    return out;
}
} // namespace paacs::geo::detail
