// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "../bnd/PlcSpatialIndex.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/predicates/Predicates.hpp"
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace paacs::geo::detail {
struct OcctThinRegionOrderingLimits {
    std::uint64_t maximumWork=5'000'000;
    std::uint64_t maximumWorkingBytes=32ULL*1024*1024;
    unsigned maximumSamples=64;
    // Charge this known field-oracle cost before the otherwise pure callback.
    std::uint64_t sizeSampleWork=1;
    // The initial constructor offers .85, 1 and 1.15 times its computed
    // height. The largest reach diagnoses competing fronts, not proof that
    // every possible apex is infeasible. This is a private scheduling policy.
    double proposalFactor=1.15;
};
struct OcctThinRegionOrderingReport {
    bool attempted=true,complete=false,selectBoundaryFirst=false,cancelled=false,budgetLimited=false;
    const char* diagnostic="not_attempted";
    std::uint64_t workUsed=0,samples=0,distinctSamples=0,blockedSamples=0,unknownSamples=0,maximumModeledBytes=0;
    // Opposite-boundary distance divided by the sampled proposal reach.
    // Zero means no strict intersection distance was measured.
    double minimumClearanceRatio=0;
    double proposalFactor=1.15;
};
namespace thin_region_ordering_detail {
struct Stop {const char* diagnostic;bool cancelled=false,budget=false;};
inline bool finite(Vec3 p){return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);}
inline Vec3 centroid(Vec3 a,Vec3 b,Vec3 c){return a+(b-a)/3.+(c-a)/3.;}
inline bool inTriangle(Vec3 p,Vec3 a,Vec3 b,Vec3 c,Vec3 normal) {
    const auto x=std::abs(normal.x),y=std::abs(normal.y),z=std::abs(normal.z);
    const unsigned axis=x>=y&&x>=z?0:y>=z?1:2;
    const auto project=[&](Vec3 q){return axis==0?Vec2{q.y,q.z}:axis==1?Vec2{q.x,q.z}:Vec2{q.x,q.y};};
    const std::array signs{predicates::orient2d(project(a),project(b),project(p)),
        predicates::orient2d(project(b),project(c),project(p)),
        predicates::orient2d(project(c),project(a),project(p))};
    bool positive=false,negative=false;
    for(auto sign:signs){positive|=sign==predicates::Sign::Positive;negative|=sign==predicates::Sign::Negative;}
    return !(positive&&negative);
}
}

// Advisory ONLY. The caller already owns a certified, consistently outward
// oriented, closed one-body source PLC. This helper neither certifies that
// precondition nor admits points/cells. Every Native/source/volume gate remains
// necessary. It borrows the immutable source and constructs one bounded BVH.
// Uncertain rays/resources abstain; allocation and external charge exceptions
// propagate so an enclosing hard resource/cancellation policy cannot be hidden.
template<class SampleSize,class Charge>
OcctThinRegionOrderingReport evaluateOcctThinRegionOrdering(const mesh::Mesh& boundary,
    double targetH,SampleSize sampleSize,const OcctThinRegionOrderingLimits& limits,
    Charge charge,const std::atomic<bool>* cancellation=nullptr) {
    using namespace thin_region_ordering_detail;
    using IndexTree=bnd::detail::PlcGrowingAabbIndex;
    using Sign=predicates::Sign;
    OcctThinRegionOrderingReport out;out.proposalFactor=limits.proposalFactor;
    const auto workLimit=std::min<std::uint64_t>(5'000'000,limits.maximumWork);
    const auto tick=[&](std::uint64_t n=1){
        if(cancellation&&cancellation->load(std::memory_order_relaxed))throw Stop{"cancelled",true,false};
        if(out.workUsed>workLimit||n>workLimit-out.workUsed)throw Stop{"work_budget",false,true};
        charge(n);out.workUsed+=n;
    };
    const auto need=[](bool okay,const char* why){if(!okay)throw Stop{why};};
    try {
        tick(0);
        need(std::isfinite(targetH)&&targetH>0&&std::isfinite(limits.proposalFactor)&&
             limits.proposalFactor>=.85&&limits.proposalFactor<=1.15,"invalid_request");
        const auto count=std::uint64_t(boundary.cellCount());
        need(count>=4&&boundary.nodeCount()>=4,"invalid_source");
        constexpr std::uint64_t fixed=65536,perFace=192;
        const auto memoryLimit=std::min<std::uint64_t>(32ULL*1024*1024,limits.maximumWorkingBytes);
        if(memoryLimit<fixed||count>(memoryLimit-fixed)/perFace)
            throw Stop{"working_budget",false,true};
        // BVH: at most 2N boxes+four size_t fields (160N on 64-bit), N order
        // entries (8N), N cumulative long-double areas (<=16N), plus fixed
        // sample/recursion/vector scratch. No full coordinate copy is retained.
        static_assert(sizeof(long double)<=16&&sizeof(std::size_t)<=8);
        out.maximumModeledBytes=fixed+perFace*count;
        std::vector<long double> cumulative(static_cast<std::size_t>(count));
        long double totalArea=0;
        Vec3 lower=boundary.node(0),upper=lower;
        const auto points=[&](std::size_t id){
            const auto* ids=boundary.cellNodes(static_cast<Index>(id));
            return std::array{boundary.node(ids[0]),boundary.node(ids[1]),boundary.node(ids[2])};
        };
        for(Index i=0;i<boundary.cellCount();++i){
            tick(8);
            const auto& cell=boundary.cell(i);
            need(cell.type==mesh::CellType::Tri3&&cell.onEntity.dim==2,"invalid_source");
            const auto* ids=boundary.cellNodes(i);
            for(unsigned j=0;j<3;++j)need(ids[j]>=0&&ids[j]<boundary.nodeCount(),"invalid_source");
            const auto p=points(std::size_t(i));
            for(auto point:p){
                need(finite(point),"coordinate_range");
                lower={std::min(lower.x,point.x),std::min(lower.y,point.y),std::min(lower.z,point.z)};
                upper={std::max(upper.x,point.x),std::max(upper.y,point.y),std::max(upper.z,point.z)};
            }
            const auto normal=(p[1]-p[0]).cross(p[2]-p[0]);
            const auto twiceArea=std::hypot(normal.x,normal.y,normal.z);
            need(finite(normal)&&std::isfinite(twiceArea)&&twiceArea>0,"coordinate_range");
            const auto next=totalArea+static_cast<long double>(twiceArea);
            need(std::isfinite(next)&&next>totalArea,"area_range");
            cumulative[std::size_t(i)]=totalArea=next;
        }
        const auto span=std::max({upper.x-lower.x,upper.y-lower.y,upper.z-lower.z});
        // Conservative envelope for the exact-sign kernels' normal products.
        need(std::isfinite(span)&&span>=1e-50&&span<=1e50,"coordinate_range");
        const auto bounds=[&](std::size_t id){
            tick(1);const auto p=points(id);
            return IndexTree::Box{{std::min({p[0].x,p[1].x,p[2].x}),std::min({p[0].y,p[1].y,p[2].y}),std::min({p[0].z,p[1].z,p[2].z})},
                {std::max({p[0].x,p[1].x,p[2].x}),std::max({p[0].y,p[1].y,p[2].y}),std::max({p[0].z,p[1].z,p[2].z})}};
        };
        IndexTree index;index.rebuild(std::size_t(count),bounds,tick);
        const unsigned sampleCount=unsigned(std::min({count,std::uint64_t(limits.maximumSamples),std::uint64_t(64)}));
        if(sampleCount<16){out.complete=true;out.diagnostic="insufficient_samples";return out;}
        std::size_t previous=std::size_t(count);
        bool previousBlocked=false,previousUnknown=false;
        double previousRatio=0;
        for(unsigned sample=0;sample<sampleCount;++sample){
            tick();
            const auto wanted=totalArea*(static_cast<long double>(sample)+.5L)/sampleCount;
            std::size_t a=0,b=std::size_t(count);
            while(a<b){tick();const auto middle=a+(b-a)/2;if(cumulative[middle]<wanted)a=middle+1;else b=middle;}
            need(a<count,"area_range");
            ++out.samples;
            bool blocked=false,ambiguous=false;double nearestRatio=0;
            if(a==previous){blocked=previousBlocked;ambiguous=previousUnknown;nearestRatio=previousRatio;}
            else {
                ++out.distinctSamples;
                const auto p=points(a);const auto ab=p[1]-p[0],ac=p[2]-p[0];
                const auto normal=ab.cross(ac);const auto norm=std::hypot(normal.x,normal.y,normal.z);
                const auto middle=centroid(p[0],p[1],p[2]);auto center=middle;
                if(ab.dot(ac)>0&&(p[0]-p[1]).dot(p[2]-p[1])>0&&(p[0]-p[2]).dot(p[1]-p[2])>0)
                    center=p[0]+(ac.cross(normal)*ab.norm2()+normal.cross(ab)*ac.norm2())/(2*normal.norm2());
                std::array lengths{ab.norm(),(p[2]-p[1]).norm(),ac.norm()};std::sort(lengths.begin(),lengths.end());
                tick(limits.sizeSampleWork);const auto requested=sampleSize(middle);
                need(std::isfinite(requested)&&requested>0,"invalid_size_sample");
                const auto h=std::min({targetH,requested,1.15*lengths[1]});
                const auto radiusSquared=std::max({(center-p[0]).norm2(),(center-p[1]).norm2(),(center-p[2]).norm2()});
                const auto height=radiusSquared<h*h?std::min(std::sqrt(h*h-radiusSquared),std::sqrt(2./3.)*h):std::sqrt(2./3.)*h;
                const auto reach=limits.proposalFactor*height;
                const auto q=middle-normal*(reach/norm);
                if(!finite(center)||!finite(middle)||!finite(q)||!std::isfinite(reach)||!(reach>0)||
                   (q.x==middle.x&&q.y==middle.y&&q.z==middle.z))ambiguous=true;
                else {
                    const auto down=[](double v){return std::nextafter(v,-std::numeric_limits<double>::infinity());};
                    const auto up=[](double v){return std::nextafter(v,std::numeric_limits<double>::infinity());};
                    const IndexTree::Box query{{down(std::min(middle.x,q.x)),down(std::min(middle.y,q.y)),down(std::min(middle.z,q.z))},
                        {up(std::max(middle.x,q.x)),up(std::max(middle.y,q.y)),up(std::max(middle.z,q.z))}};
                    index.all(std::size_t(count),query,bounds,tick,[&](std::size_t id){
                        tick(8);if(id==a)return true;const auto triangle=points(id);
                        const auto start=predicates::orient3d(triangle[0],triangle[1],triangle[2],middle);
                        const auto end=predicates::orient3d(triangle[0],triangle[1],triangle[2],q);
                        if(start==Sign::Zero){
                            if(end==Sign::Zero||inTriangle(middle,triangle[0],triangle[1],triangle[2],
                                (triangle[1]-triangle[0]).cross(triangle[2]-triangle[0])))ambiguous=true;
                            return true;
                        }
                        if(start==end)return true;
                        const std::array sides{predicates::orient3d(middle,q,triangle[0],triangle[1]),
                            predicates::orient3d(middle,q,triangle[1],triangle[2]),
                            predicates::orient3d(middle,q,triangle[2],triangle[0])};
                        bool positive=false,negative=false,zero=false;
                        for(auto sign:sides){positive|=sign==Sign::Positive;negative|=sign==Sign::Negative;zero|=sign==Sign::Zero;}
                        if(positive&&negative)return true;
                        if(zero||end==Sign::Zero||start!=Sign::Negative||end!=Sign::Positive){ambiguous=true;return true;}
                        const auto s=predicates::orient3dDeterminantEstimate(triangle[0],triangle[1],triangle[2],middle);
                        const auto e=predicates::orient3dDeterminantEstimate(triangle[0],triangle[1],triangle[2],q);
                        const auto ratio=s/(s-e);
                        if(!std::isfinite(ratio)||!(ratio>0&&ratio<1)){ambiguous=true;return true;}
                        blocked=true;if(nearestRatio==0||ratio<nearestRatio)nearestRatio=ratio;
                        return true;
                    });
                }
                previous=a;previousBlocked=blocked;previousUnknown=ambiguous;previousRatio=nearestRatio;
            }
            if(ambiguous)++out.unknownSamples;else if(blocked)++out.blockedSamples;
            if(nearestRatio>0&&(out.minimumClearanceRatio==0||nearestRatio<out.minimumClearanceRatio))
                out.minimumClearanceRatio=nearestRatio;
        }
        tick(0);out.complete=true;
        out.selectBoundaryFirst=out.unknownSamples==0&&out.distinctSamples>=16&&out.blockedSamples*5>=out.samples*4;
        out.diagnostic=out.unknownSamples?"ambiguous_rays":out.distinctSamples<16?"insufficient_distinct_samples":
            out.selectBoundaryFirst?"opposing_boundary_in_proposal_reach":"original_order";
    }catch(const Stop& stop){out.diagnostic=stop.diagnostic;out.cancelled=stop.cancelled;out.budgetLimited=stop.budget;}
    return out;
}
} // namespace paacs::geo::detail
