// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctThinRegionOrdering.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <array>
#include <atomic>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <new>
#include <stdexcept>
#include <vector>

using namespace paacs;
using namespace paacs::geo::detail;
namespace {
unsigned checks=0;
void check(bool value,const char* message){++checks;if(!value)throw std::runtime_error(message);}

// Owned piecewise-planar hollow prisms. Both profiles are CCW; inner walls
// reverse the material winding. Annular caps leave the bore open at both ends.
mesh::Mesh hollow(const std::vector<Vec2>& outer,const std::vector<Vec2>& inner,double length) {
    mesh::Mesh result;const auto count=outer.size();
    check(count==inner.size()&&count>=4,"invalid owned hollow profile");
    constexpr unsigned layers=4;
    std::vector<std::vector<Index>> outside(layers+1),inside(layers+1);
    for(unsigned layer=0;layer<=layers;++layer)for(unsigned ring=0;ring<2;++ring){
        auto& ids=ring?inside[layer]:outside[layer];
        for(std::size_t i=0;i<count;++i){const auto p=ring?inner[i]:outer[i];
            ids.push_back(result.addNode({p.u,p.v,length*layer/layers},{1,int(ring*count+i+1)}));}
    }
    const auto quad=[&](std::array<Index,4> ids,int face,bool reverse){
        std::array<Index,3> a{ids[0],ids[1],ids[2]},b{ids[0],ids[2],ids[3]};
        if(reverse){std::swap(a[1],a[2]);std::swap(b[1],b[2]);}
        result.addCell(mesh::CellType::Tri3,a,{2,face});result.addCell(mesh::CellType::Tri3,b,{2,face});
    };
    for(unsigned layer=0;layer<layers;++layer)for(std::size_t i=0;i<count;++i){const auto next=(i+1)%count;
        quad({outside[layer][i],outside[layer][next],outside[layer+1][next],outside[layer+1][i]},int(i+1),false);
        quad({inside[layer][i],inside[layer][next],inside[layer+1][next],inside[layer+1][i]},int(count+i+1),true);
    }
    for(std::size_t i=0;i<count;++i){const auto next=(i+1)%count;
        quad({outside[0][i],outside[0][next],inside[0][next],inside[0][i]},int(2*count+1),true);
        quad({outside[layers][i],outside[layers][next],inside[layers][next],inside[layers][i]},int(2*count+2),false);
    }
    std::map<std::array<Index,2>,std::pair<unsigned,int>> edges;
    double signedVolume=0;
    for(Index i=0;i<result.cellCount();++i){const auto* ids=result.cellNodes(i);
        for(unsigned j=0;j<3;++j){const auto a=ids[j],b=ids[(j+1)%3];
            auto& use=edges[{std::min(a,b),std::max(a,b)}];++use.first;use.second+=a<b?1:-1;}
        const auto a=result.node(ids[0]),b=result.node(ids[1]),c=result.node(ids[2]);
        signedVolume+=a.dot(b.cross(c))/6.;
        if(result.cell(i).onEntity.tag>int(2*count)){
            const std::array signs{predicates::orient2d({a.x,a.y},{b.x,b.y},{0,0}),
                predicates::orient2d({b.x,b.y},{c.x,c.y},{0,0}),
                predicates::orient2d({c.x,c.y},{a.x,a.y},{0,0})};
            bool positive=false,negative=false;
            for(auto sign:signs){positive|=sign==predicates::Sign::Positive;negative|=sign==predicates::Sign::Negative;}
            check(positive&&negative,"owned cap filled the open bore");
        }
    }
    for(const auto& [edge,use]:edges){(void)edge;check(use.first==2&&use.second==0,"owned source is not closed and consistently oriented");}
    const auto area=[](const auto& loop){double value=0;for(std::size_t i=0;i<loop.size();++i){
        const auto next=(i+1)%loop.size();value+=loop[i].u*loop[next].v-loop[next].u*loop[i].v;}return value*.5;};
    const auto expected=(area(outer)-area(inner))*length;
    check(expected>0&&std::abs(signedVolume-expected)<=expected*1e-12,"owned hollow source volume/orientation differs from its analytic polygon profile");
    return result;
}
mesh::Mesh rectangle(double width,double height,double wall,double length=40.) {
    const double x=width/2,y=height/2;
    return hollow({{-x,-y},{x,-y},{x,y},{-x,y}},
        {{wall-x,wall-y},{x-wall,wall-y},{x-wall,y-wall},{wall-x,y-wall}},length);
}
mesh::Mesh circular(double outerRadius,double wall) {
    std::vector<Vec2> outer,inner;
    for(unsigned i=0;i<32;++i){const auto angle=2*std::acos(-1.)*i/32;
        outer.push_back({outerRadius*std::cos(angle),outerRadius*std::sin(angle)});
        inner.push_back({(outerRadius-wall)*std::cos(angle),(outerRadius-wall)*std::sin(angle)});}
    return hollow(outer,inner,40.);
}
mesh::Mesh transformed(mesh::Mesh source) {
    const auto axis=Vec3{1,2,3}/std::sqrt(14.);const double c=std::cos(.61),s=std::sin(.61);
    for(Index i=0;i<source.nodeCount();++i){const auto p=source.node(i);
        source.setNode(i,p*c+axis.cross(p)*s+axis*(axis.dot(p)*(1-c))+Vec3{17,-23,11});}
    return source;
}
OcctThinRegionOrderingReport run(const mesh::Mesh& source,double size=4.,double fieldSize=4.,
    OcctThinRegionOrderingLimits limits={}) {
    const auto signature=mesh::meshHash(source);std::uint64_t charged=0,fieldCalls=0;
    limits.sizeSampleWork=7;
    const auto report=evaluateOcctThinRegionOrdering(source,size,[&](Vec3){++fieldCalls;return fieldSize;},
        limits,[&](std::uint64_t work){charged+=work;});
    check(report.workUsed==charged&&charged<=limits.maximumWork,"ordering work is not charged exactly once");
    check(report.maximumModeledBytes<=limits.maximumWorkingBytes,"ordering exceeded its admitted modeled memory");
    check(report.samples<=64&&fieldCalls<=report.distinctSamples,"ordering exceeded its sampled field allowance");
    check(mesh::meshHash(source)==signature,"advisory ordering changed source coordinates, connectivity or tags");
    return report;
}
void geometry() {
    const std::array sources{rectangle(40,40,1),rectangle(60,30,1),circular(20,1)};
    for(const auto& source:sources){
        const auto report=run(source);
        std::cout<<"thin: samples="<<report.samples<<" distinct="<<report.distinctSamples
                 <<" blocked="<<report.blockedSamples<<" unknown="<<report.unknownSamples
                 <<" ratio="<<report.minimumClearanceRatio<<" diagnostic="<<report.diagnostic<<'\n';
        check(report.complete&&report.selectBoundaryFirst&&!report.cancelled&&!report.budgetLimited,
              "owned thin hollow source did not select the bounded boundary-first schedule");
        check(report.distinctSamples>=16&&report.unknownSamples==0&&report.blockedSamples*5>=report.samples*4,
              "thin schedule lacks its complete, substantial opposite-wall evidence");
        check(report.minimumClearanceRatio>0&&report.minimumClearanceRatio<1,
              "thin schedule has no strictly nearer opposite boundary");
        const auto rotated=run(transformed(source));
        check(rotated.complete&&rotated.selectBoundaryFirst&&rotated.unknownSamples==0,
              "rotation or translation changed the geometry-based schedule");
        const auto finer=run(source,4.,.25);
        check(finer.complete&&!finer.selectBoundaryFirst&&finer.blockedSamples==0,
              "selected finer field was ignored when computing the actual proposal reach");
    }
    for(const auto& source:std::array{rectangle(40,40,4),rectangle(60,30,4),circular(20,4)}){
        const auto report=run(source);
        check(report.complete&&!report.selectBoundaryFirst&&report.blockedSamples==0,
              "thick-wall control was incorrectly scheduled as a competing thin front");
    }
    auto coarse=rectangle(40,40,1);
    auto strict=OcctThinRegionOrderingLimits{};strict.proposalFactor=.85;
    const auto lower=run(coarse,4.,4.,strict),upper=run(coarse);
    check(lower.selectBoundaryFirst&&upper.selectBoundaryFirst&&
          lower.minimumClearanceRatio>upper.minimumClearanceRatio,
          "reported clearance ratio does not describe the chosen actual proposal height");
}
void refusals() {
    const auto source=rectangle(40,40,1);const auto signature=mesh::meshHash(source);
    OcctThinRegionOrderingLimits limits;
    limits.maximumWork=1;const auto work=run(source,4,4,limits);
    check(!work.complete&&!work.selectBoundaryFirst&&work.budgetLimited,"work limit did not abstain");
    limits={};limits.maximumWorkingBytes=1;const auto memory=run(source,4,4,limits);
    check(!memory.complete&&!memory.selectBoundaryFirst&&memory.budgetLimited,"memory limit did not abstain before allocation");
    limits={};limits.maximumSamples=8;const auto samples=run(source,4,4,limits);
    check(samples.complete&&!samples.selectBoundaryFirst&&samples.samples==0,"insufficient samples did not abstain");
    std::atomic<bool> cancelled{true};std::uint64_t charged=0;
    const auto stop=evaluateOcctThinRegionOrdering(source,4.,[](Vec3){return 4.;},limits,
        [&](std::uint64_t n){charged+=n;},&cancelled);
    check(stop.cancelled&&!stop.complete&&!stop.selectBoundaryFirst&&charged==0,
          "pre-cancelled ordering started work or selected a route");
    cancelled=false;limits={};
    const auto interrupted=evaluateOcctThinRegionOrdering(source,4.,[](Vec3){return 4.;},limits,
        [&](std::uint64_t n){charged+=n;if(charged>100)cancelled=true;},&cancelled);
    check(interrupted.cancelled&&!interrupted.complete&&!interrupted.selectBoundaryFirst,
          "mid-probe cancellation was ignored");
    bool allocation=false;
    try{(void)evaluateOcctThinRegionOrdering(source,4.,[](Vec3)->double{throw std::bad_alloc{};},limits,
        [](std::uint64_t){});}catch(const std::bad_alloc&){allocation=true;}
    check(allocation,"allocation failure from field evaluation was hidden as advisory abstention");
    const auto invalid=evaluateOcctThinRegionOrdering(source,4.,[](Vec3){return std::numeric_limits<double>::quiet_NaN();},
        limits,[](std::uint64_t){});
    check(!invalid.complete&&!invalid.selectBoundaryFirst,"nonfinite field sample selected a schedule");
    mesh::Mesh duplicated=source;
    for(Index i=0;i<source.cellCount();++i){const auto* ids=source.cellNodes(i);
        duplicated.addCell(mesh::CellType::Tri3,{ids[0],ids[1],ids[2]},source.cell(i).onEntity);}
    const auto ambiguous=run(duplicated);
    check(ambiguous.complete&&ambiguous.unknownSamples>0&&!ambiguous.selectBoundaryFirst,
          "ambiguous coincident source faces were presented as thin-region evidence");
    check(mesh::meshHash(source)==signature,"refused advisory changed its immutable source");
}
}
int main(){try{geometry();refusals();std::cout<<checks<<" thin-region scheduling checks passed\n";return 0;}
catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
