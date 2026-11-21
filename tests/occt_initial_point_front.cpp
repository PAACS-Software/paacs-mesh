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
#include <numbers>
#include <stdexcept>

using namespace paacs;
using namespace paacs::geo::detail;
namespace {
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
auto bits(Vec3 p){return std::array{std::bit_cast<std::uint64_t>(p.x),
    std::bit_cast<std::uint64_t>(p.y),std::bit_cast<std::uint64_t>(p.z)};}
struct UniformField {
    double h=.45;
    double eval(Vec3)const{return h;}
    std::uint64_t queryWork()const{return 1;}
};
mesh::Mesh cubeSurface(){
    mesh::Mesh mesh;std::map<std::array<int,3>,Index> ids;
    const auto vertex=[&](std::array<int,3> p){
        if(const auto it=ids.find(p);it!=ids.end())return it->second;
        unsigned on=0;for(auto x:p)on+=std::abs(x)==2;
        const auto id=mesh.addNode({p[0]*.5,p[1]*.5,p[2]*.5},{3-int(on),int(ids.size()+1)});
        ids.emplace(p,id);return id;
    };
    for(unsigned axis=0;axis<3;++axis)for(int sign:std::array{-1,1})
        for(int u=-2;u<2;++u)for(int v=-2;v<2;++v){
            std::array<Index,4> q{};const std::array du{0,1,1,0},dv{0,0,1,1};
            for(unsigned j=0;j<4;++j){std::array<int,3> p{};p[axis]=sign*2;
                p[(axis+1)%3]=u+du[j];p[(axis+2)%3]=v+dv[j];q[j]=vertex(p);}
            std::array<Index,3> a{q[0],q[1],q[2]},b{q[0],q[2],q[3]};
            if(sign<0){std::swap(a[1],a[2]);std::swap(b[1],b[2]);}
            const mesh::Classification face{2,int(2*axis+(sign>0?2:1))};
            mesh.addCell(mesh::CellType::Tri3,a,face);mesh.addCell(mesh::CellType::Tri3,b,face);
        }
    return mesh;
}
void originalPrefix(const mesh::Mesh& a,const mesh::Mesh& b){
    check(b.nodeCount()>=a.nodeCount()&&b.cellCount()==a.cellCount(),"source counts changed");
    for(Index i=0;i<a.nodeCount();++i)
        check(bits(a.node(i))==bits(b.node(i))&&a.nodeClass(i)==b.nodeClass(i),"fixed source node changed");
    for(Index i=0;i<a.cellCount();++i){
        check(a.cell(i).type==b.cell(i).type&&a.cell(i).onEntity==b.cell(i).onEntity,"source cell identity changed");
        for(unsigned j=0;j<3;++j)check(a.cellNodes(i)[j]==b.cellNodes(i)[j],"source oriented triangle changed");
    }
}
void exerciseUniformFront(){
    const auto source=cubeSurface();const UniformField field;
    OcctInitialPointFrontLimits limits;limits.sites=96;limits.cones=384;
    OcctInitialPointFrontReport report;std::uint64_t charged=0,native=0;
    const auto inside=[&](Vec3 p){++native;return std::abs(p.x)<1&&std::abs(p.y)<1&&std::abs(p.z)<1;};
    const auto result=appendOcctInitialPointFrontSites(source,field,field.h,{3,7},limits,inside,
        [&](std::uint64_t n){charged+=n;},&report);
    originalPrefix(source,result);
    check(report.initialized&&report.auxiliaryInitialized&&report.sites>0,"uniform front created no initial sites");
    check(report.auxiliaryCommits==report.sites&&report.auxiliaryFullCertifications>=report.sites,
        "initial site bypassed full auxiliary patch admission");
    check(report.auxiliaryComparisonCommits<=report.auxiliaryRankedCommits&&
        report.auxiliaryRankedCommits<=report.sites,"auxiliary quality availability accounting failed");
    check(report.deepSites>0&&report.propagatedFaceAccepted>0&&report.maximumDepth>1,
        "initial construction did not propagate beyond source faces");
    check(report.reusedApices>0&&report.seedCones<=limits.cones&&report.sites<=limits.sites&&
        report.candidates<=limits.candidates,"front closure or construction limits failed");
    check(report.workUsed==charged&&charged<=limits.work&&report.nativeQueries==native&&native<=limits.queries,
        "initial construction escaped work or native accounting");
    check(report.minimumConeMeanRatio>=.01&&report.minimumConeScaledJacobian>=1e-6,
        "accepted temporary cone missed its angular floor");
    check(report.effectiveScales.size()==std::size_t(result.nodeCount()),"construction scales missing");
    double deepest=0.;
    for(Index i=source.nodeCount();i<result.nodeCount();++i){
        const auto p=result.node(i);
        check(std::abs(p.x)<1&&std::abs(p.y)<1&&std::abs(p.z)<1&&
            result.nodeClass(i)==mesh::Classification{3,7},"new point is not strict same-body material");
        deepest=std::max(deepest,1-std::max({std::abs(p.x),std::abs(p.y),std::abs(p.z)}));
        for(Index j=0;j<i;++j)check((p-result.node(j)).norm()+1e-13>=
            .65*std::min(report.effectiveScales[i],report.effectiveScales[j]),
            "initial point violates separation from an earlier live vertex");
    }
    check(deepest>=field.h,"topological propagation never reached one requested size into the volume");
    const auto again=appendOcctInitialPointFrontSites(source,field,field.h,{3,7},limits,
        [](Vec3 p){return std::abs(p.x)<1&&std::abs(p.y)<1&&std::abs(p.z)<1;},[](std::uint64_t){});
    check(again.nodeCount()==result.nodeCount(),"initial point count is not deterministic");
    for(Index i=0;i<result.nodeCount();++i)check(bits(again.node(i))==bits(result.node(i)),"initial point is not deterministic");
    std::cout<<"uniform initial front: sites="<<report.sites<<" deep="<<report.deepSites
        <<" physicalDepth="<<deepest<<" work="<<charged<<'\n';
}
void exerciseRefusals(){
    const auto source=cubeSurface();const UniformField field;
    OcctInitialPointFrontLimits limits;limits.sites=8;limits.cones=64;
    for(unsigned which=0;which<7;++which){
        auto cap=limits;
        switch(which){case 0:cap.sites=0;break;case 1:cap.cones=0;break;case 2:cap.work=0;break;
            case 3:cap.queries=0;break;case 4:cap.workingBytes=1;break;case 5:cap.candidates=0;break;
            default:cap.frontStates=0;break;}
        OcctInitialPointFrontReport out;
        const auto result=appendOcctInitialPointFrontSites(source,field,field.h,{3,7},cap,
            [](Vec3){return true;},[](std::uint64_t){},&out);
        originalPrefix(source,result);check(result.nodeCount()==source.nodeCount()&&out.sites==0,
            "an exhausted construction limit leaked a point");
    }
    OcctInitialPointFrontReport refused;
    const auto result=appendOcctInitialPointFrontSites(source,field,field.h,{3,7},limits,
        [](Vec3){return false;},[](std::uint64_t){},&refused);
    originalPrefix(source,result);
    check(result.nodeCount()==source.nodeCount()&&refused.nativeRejects>0,"native refusal leaked an initial point");
    bool threw=false;
    try{(void)appendOcctInitialPointFrontSites(source,field,field.h,{3,7},limits,
        [](Vec3){return true;},[](std::uint64_t){throw std::runtime_error("charge-stop");});}
    catch(const std::runtime_error&){threw=true;}
    check(threw,"external construction exception was swallowed");
}
mesh::Mesh prismSurface(){
    constexpr Index count=48;mesh::Mesh source;
    for(unsigned layer=0;layer<2;++layer)for(Index i=0;i<count;++i){
        const auto angle=2*std::numbers::pi*double(i)/count;
        source.addNode({std::cos(angle),std::sin(angle),layer?.5:-.5},{1,int(i+1+layer*count)});
    }
    for(Index i=1;i<count-1;++i){
        source.addCell(mesh::CellType::Tri3,std::array<Index,3>{0,i+1,i},{2,1});
        source.addCell(mesh::CellType::Tri3,std::array<Index,3>{count,count+i,count+i+1},{2,2});
    }
    for(Index i=0;i<count;++i){const auto j=(i+1)%count;
        source.addCell(mesh::CellType::Tri3,std::array<Index,3>{i,j,count+j},{2,3});
        source.addCell(mesh::CellType::Tri3,std::array<Index,3>{i,count+j,count+i},{2,3});
    }return source;
}
void exerciseCompleteConflictCavity(){
    using Front=bnd::detail::PlcGiftWrapping;
    const auto source=prismSurface();Front front;
    for(const auto p:source.nodes())front.vertices.push_back({p,front.vertices.size()+1});
    for(Index i=0;i<source.cellCount();++i){const auto* v=source.cellNodes(i);
        front.boundary.push_back({bnd::TetFacetVertexId(v[0]),bnd::TetFacetVertexId(v[1]),bnd::TetFacetVertexId(v[2])});}
    check(front.initialize(),"large-conflict source is not a closed oriented front");
    std::uint64_t work=0;bool refuseWork=false;
    const auto charge=[&](std::uint64_t n){if(refuseWork)throw std::runtime_error("selected-stop");
        check(n<=100'000'000-work,"large-conflict work guard exceeded");work+=n;};
    front.chargeWork=charge;
    const auto inside=[&](Vec3 p){
        for(const auto f:front.boundary)
            if(predicates::orient3dFiltered(front.vertices[f[0]].position,front.vertices[f[1]].position,
                front.vertices[f[2]].position,p)!=predicates::Sign::Negative)return false;
        return true;
    };
    OcctInitialDelaunayScorer limited,complete;
    check(limited.initialize(source,{3,7},4096,64ULL*1024*1024,charge,inside,64)&&
        complete.initialize(source,{3,7},4096,64ULL*1024*1024,charge,inside,0),
        "auxiliary Delaunay initialization refused the owned prism");
    const auto originalStored=complete.storedCells(),originalMap=complete.auxiliaryToFront().size();
    const auto stagedFront=bnd::TetFacetVertexId(front.vertices.size());
    front.vertices.push_back({{},front.vertices.size()+1});
    front.interiorAdmission=inside;front.centroidAdmission=inside;
    bnd::TetFacetCell occupancy;
    check(inside({})&&front.valid(Front::Front{front.boundary.front()},stagedFront,occupancy),
        "large-conflict point failed strict material and exact temporary occupancy admission");
    check(!limited.score({},front)&&limited.scoreCavityUnavailable==1&&limited.scoreLocateMisses==0,
        "64-cell guard did not refuse the complete conflict cavity");
    auto plan=complete.score({},front);
    check(plan&&plan->valid&&!plan->certified&&plan->patch.removed.size()>64,
        "automatic cavity admission did not retain the complete large conflict");
    check(complete.storedCells()==originalStored&&complete.acceptedSites()==0&&
        complete.auxiliaryToFront().size()==originalMap,"candidate scoring mutated auxiliary topology");
    refuseWork=true;bool threw=false;
    try{(void)complete.commit(*plan,stagedFront);}catch(const std::runtime_error&){threw=true;}
    refuseWork=false;
    check(threw&&complete.storedCells()==originalStored&&complete.acceptedSites()==0&&
        complete.auxiliaryToFront().size()==originalMap,"selected admission exception changed auxiliary state");
    check(complete.commit(*plan,stagedFront)&&complete.acceptedSites()==1&&
        complete.maximumAcceptedCavity>64&&complete.certificationRefusals==0,
        "complete selected cavity failed exact patch and survivor admission");
    front.insert(occupancy);
    check(!complete.commit(*plan,stagedFront),"stale auxiliary plan was applied twice");
    const auto& mapping=complete.auxiliaryToFront();
    for(Index i=0;i<source.nodeCount();++i)check(mapping[i]==bnd::TetFacetVertexId(i),"original auxiliary identity changed");
    for(unsigned i=0;i<8;++i)check(mapping[source.nodeCount()+i]==bnd::kInvalidTetFacetVertex,"guard became a source point");
    check(mapping.back()==stagedFront,"accepted auxiliary point mapping changed");
    std::cout<<"complete initial cavity: removed="<<plan->patch.removed.size()<<" added="<<plan->patch.added.size()
        <<" work="<<work<<'\n';
}
}
int main(){try{
    exerciseUniformFront();exerciseRefusals();exerciseCompleteConflictCavity();
    std::cout<<"initial volume points, complete conflict cavities and transactional guards passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
