// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctInitialPointFront.hpp"
#include <array>
#include <bit>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>

using namespace paacs;
using namespace paacs::geo::detail;
namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
struct Field {
    double h=.45;
    double eval(Vec3)const{return h;}
    std::uint64_t queryWork()const{return 1;}
};
mesh::Mesh boundary(){
    mesh::Mesh mesh;std::map<std::array<int,3>,Index> ids;
    const auto node=[&](std::array<int,3> p){
        if(const auto it=ids.find(p);it!=ids.end())return it->second;
        const auto id=mesh.addNode({p[0]*.5,p[1]*.5,p[2]*.5},{2,int(ids.size()+1)});
        ids.emplace(p,id);return id;
    };
    for(unsigned axis=0;axis<3;++axis)for(const int sign:std::array{-1,1})
        for(int u=-2;u<2;++u)for(int v=-2;v<2;++v){
            std::array<Index,4> q{};const std::array du{0,1,1,0},dv{0,0,1,1};
            for(unsigned j=0;j<4;++j){std::array<int,3> p{};p[axis]=sign*2;
                p[(axis+1)%3]=u+du[j];p[(axis+2)%3]=v+dv[j];q[j]=node(p);}
            std::array<Index,3> a{q[0],q[1],q[2]},b{q[0],q[2],q[3]};
            if(sign<0){std::swap(a[1],a[2]);std::swap(b[1],b[2]);}
            const mesh::Classification face{2,int(2*axis+(sign>0?2:1))};
            mesh.addCell(mesh::CellType::Tri3,a,face);mesh.addCell(mesh::CellType::Tri3,b,face);
        }
    return mesh;
}
auto bits(Vec3 p){return std::array{std::bit_cast<std::uint64_t>(p.x),
    std::bit_cast<std::uint64_t>(p.y),std::bit_cast<std::uint64_t>(p.z)};}
struct Result {mesh::Mesh mesh;OcctInitialPointFrontReport report;};
Result run(const mesh::Mesh& source,const Field& field,const OcctInitialPointFrontLimits& limits){
    Result result;std::uint64_t charged=0,queries=0;
    result.mesh=appendOcctInitialPointFrontSites(source,field,field.h,{3,7},limits,
        [&](Vec3 p){++queries;return std::abs(p.x)<1&&std::abs(p.y)<1&&std::abs(p.z)<1;},
        [&](std::uint64_t n){charged+=n;},&result.report);
    require(result.report.workUsed==charged&&charged<=limits.work,"initial work accounting escaped its allowance");
    require(result.report.nativeQueries==queries&&queries<=limits.queries,"native queries escaped their allowance");
    require(result.mesh.cellCount()==source.cellCount(),"construction exported temporary occupancy cells");
    require(result.mesh.nodeCount()==source.nodeCount()+Index(result.report.sites),"accepted site count differs");
    for(Index i=0;i<source.nodeCount();++i)
        require(bits(result.mesh.node(i))==bits(source.node(i))&&result.mesh.nodeClass(i)==source.nodeClass(i),
            "construction changed prescribed boundary coordinates or identities");
    for(Index i=0;i<source.cellCount();++i){
        require(result.mesh.cell(i).onEntity==source.cell(i).onEntity,"construction changed a face identity");
        for(unsigned j=0;j<3;++j)require(result.mesh.cellNodes(i)[j]==source.cellNodes(i)[j],"construction changed a source facet");
    }
    for(Index i=source.nodeCount();i<result.mesh.nodeCount();++i){
        const auto p=result.mesh.node(i);
        require(std::abs(p.x)<1&&std::abs(p.y)<1&&std::abs(p.z)<1&&
            result.mesh.nodeClass(i)==mesh::Classification{3,7},"accepted point is not strict same-body material");
    }
    require(result.report.sites<=limits.sites&&result.report.candidates<=limits.candidates,"point/candidate limits escaped");
    require(result.report.minimumConeMeanRatio>=.01&&result.report.minimumConeScaledJacobian>=1e-6,
        "phase scheduling weakened the temporary cone shape gate");
    return result;
}
void exerciseCandidateFairness(){
    const auto source=boundary();const Field field;
    OcctInitialPointFrontLimits limits;limits.sites=64;limits.cones=512;
    limits.candidates=64;limits.work=16'000'000;limits.queries=20'000;
    const auto legacy=run(source,field,limits);
    require(legacy.report.sourceCandidates>0&&legacy.report.bulkCandidates==0,
        "fixture did not expose source-wave candidate starvation");
    limits.balanceSourceAndBulk=true;
    const auto fair=run(source,field,limits);
    require(fair.report.sourcePhaseDeferred&&fair.report.sourceDeferredFaces>0,
        "first source wave did not defer a resumable face");
    require(fair.report.reservedBulkCandidates>0&&fair.report.reservedBulkSites>0&&fair.report.bulkCandidates>0,
        "a candidate-limited source wave starved the reserved bulk attempt");
    require(fair.report.bulkSites>0,"reserved bulk attempt produced no material point");
    const auto repeated=run(source,field,limits);
    require(repeated.mesh.nodeCount()==fair.mesh.nodeCount(),"fair scheduling is not deterministic");
    for(Index i=0;i<fair.mesh.nodeCount();++i)require(bits(fair.mesh.node(i))==bits(repeated.mesh.node(i)),
        "fair scheduling changed an accepted point on repeat");
    std::cout<<"candidate fairness: legacy="<<legacy.report.sites<<" fair="<<fair.report.sites
        <<" bulk="<<fair.report.bulkSites<<" deferred="<<fair.report.sourceDeferredFaces<<'\n';
}
void exerciseWorkAndQueryFairness(){
    const auto source=boundary();const Field field;
    OcctInitialPointFrontLimits limits;limits.sites=64;limits.cones=512;
    limits.balanceSourceAndBulk=true;limits.work=16'000'000;limits.queries=20'000;
    limits.reservedBulkWork=limits.work;
    const auto work=run(source,field,limits);
    require(work.report.sourcePhaseDeferred&&work.report.sourceCandidates==0&&work.report.bulkCandidates>0,
        "source work refusal did not hand remaining work to bulk");
    limits.reservedBulkWork=0;limits.reservedBulkQueries=limits.queries;
    const auto query=run(source,field,limits);
    require(query.report.sourcePhaseDeferred&&query.report.bulkCandidates>0,
        "source native-query refusal prevented the reserved bulk attempt");
}
void exerciseExpandedSiteAllowance(){
    const auto source=boundary();const Field field{.1};
    OcctInitialPointFrontLimits limits;limits.sites=1536;limits.cones=6144;
    limits.frontStates=source.cellCount()+3*limits.cones;
    limits.candidates=131'072;limits.work=250'000'000;limits.queries=250'000;
    limits.workingBytes=512ULL*1024*1024;
    const auto result=run(source,field,limits);
    require(result.report.siteLimit==1536&&result.report.sites>1024,
        "explicit larger site allowance still stopped at the old 1024-site ceiling");
    std::cout<<"expanded initial cloud: sites="<<result.report.sites<<" bulk="<<result.report.bulkSites
        <<" candidates="<<result.report.candidates<<" work="<<result.report.workUsed<<'\n';
}
void exerciseProgressiveUnusedCeilings(){
    const auto source=boundary();const Field field;
    OcctInitialPointFrontLimits limits;
    limits.progressiveResources=true;limits.work=16'000'000;limits.queries=20'000;
    limits.workingBytes=512ULL*1024*1024;limits.candidates=20'000;
    const auto small=run(source,field,limits);
    limits.sites=1'000'000;limits.cones=5'000'000;limits.frontStates=15'000'192;
    const auto large=run(source,field,limits);
    require(large.report.initialized&&large.report.auxiliaryInitialized&&large.report.sites>0,
        "unused future front ceilings prevented progressive construction");
    require(large.report.maximumCombinedModeledBytes>0&&
        large.report.maximumCombinedModeledBytes<=limits.workingBytes,
        "front, heap and numerical scratch exceeded combined admission");
    require(large.mesh.nodeCount()==small.mesh.nodeCount(),"inactive future ceilings changed initial point count");
    for(Index i=0;i<large.mesh.nodeCount();++i)
        require(bits(large.mesh.node(i))==bits(small.mesh.node(i)),"inactive future ceilings changed initial point positions");
    std::cout<<"progressive front: sites="<<large.report.sites<<" peak_model="
        <<large.report.maximumCombinedModeledBytes<<" retries="<<large.report.resourceRetries<<'\n';
}
}
int main(){try{
    exerciseCandidateFairness();exerciseWorkAndQueryFairness();exerciseExpandedSiteAllowance();
    exerciseProgressiveUnusedCeilings();
    std::cout<<"initial construction allowances and fair phase scheduling passed\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
