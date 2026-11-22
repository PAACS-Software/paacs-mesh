// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctInitialPointFront.hpp"
#include "PlcVoidGeometry.hpp"
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <iostream>
#include <map>
#include <new>
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
    OcctInitialPointFrontLimits limits;limits.sites=192;limits.cones=384;limits.work=16'000'000;
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
    check(report.sourcePhaseCompleted&&report.sourcePhaseSites>0&&report.bulkSites>0&&
        report.sourcePhaseSites+report.bulkSites<=report.sites&&report.bulkCandidates>=report.bulkSites,
        "initial source wave did not hand its point cloud to bulk placement");
    check(std::isfinite(report.bulkMaximumUnresolvedRatio)&&report.bulkMaximumUnresolvedRatio>=0.,
        "remaining bulk-potential diagnostic is not finite");
    check(report.propagatedFaceAccepted>0,"closure did not advance a parked propagated face");
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
        for(Index j=source.nodeCount();j<i;++j)check((p-result.node(j)).norm()+1e-13>=.65*field.h,
            "two accepted interior sites violate requested-field separation");
    }
    check(deepest>=field.h,"initial points never reached one requested size into the volume");
    OcctInitialPointFrontReport repeated;
    const auto again=appendOcctInitialPointFrontSites(source,field,field.h,{3,7},limits,
        [](Vec3 p){return std::abs(p.x)<1&&std::abs(p.y)<1&&std::abs(p.z)<1;},[](std::uint64_t){},&repeated);
    check(again.nodeCount()==result.nodeCount(),"initial point count is not deterministic");
    for(Index i=0;i<result.nodeCount();++i)check(bits(again.node(i))==bits(result.node(i)),"initial point is not deterministic");
    check(repeated.sourcePhaseCompleted==report.sourcePhaseCompleted&&repeated.sourcePhaseSites==report.sourcePhaseSites&&
        repeated.bulkSites==report.bulkSites&&repeated.bulkCandidates==report.bulkCandidates&&
        repeated.bulkQueueExhausted==report.bulkQueueExhausted&&
        repeated.bulkMaximumUnresolvedRatio==report.bulkMaximumUnresolvedRatio,
        "source/bulk phase report is not deterministic");
    std::cout<<"uniform initial front: sites="<<report.sites<<" source="<<report.sourcePhaseSites<<" bulk="<<report.bulkSites
        <<" physicalDepth="<<deepest<<" work="<<charged<<'\n';
}
void exerciseFineBoundaryInteriorExclusion(){
    const auto source=cubeSurface();const UniformField field{.8};
    OcctInitialPointFrontLimits limits;limits.sites=64;limits.cones=384;
    OcctInitialPointFrontReport report;
    const auto result=appendOcctInitialPointFrontSites(source,field,field.h,{3,7},limits,
        [](Vec3 p){return std::abs(p.x)<1&&std::abs(p.y)<1&&std::abs(p.z)<1;},
        [](std::uint64_t){},&report);
    originalPrefix(source,result);
    check(report.auxiliaryInitialized&&report.sites>=2,
        "fine fixed boundary did not produce an interior-pair fixture");
    bool finerBirth=false;
    for(Index i=source.nodeCount();i<result.nodeCount();++i){
        finerBirth=finerBirth||report.effectiveScales[i]<field.h;
        for(Index j=0;j<source.nodeCount();++j)check((result.node(i)-result.node(j)).norm()+1e-13>=
            .65*std::min(report.effectiveScales[i],report.effectiveScales[j]),
            "requested packing changed the geometric boundary exclusion contract");
        for(Index j=source.nodeCount();j<i;++j)check((result.node(i)-result.node(j)).norm()+1e-13>=.65*field.h,
            "a fine boundary birth radius persisted between interior sites");
    }
    check(finerBirth,"fixture did not exercise geometry finer than the requested field");
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
using Scorer=OcctInitialDelaunayScorer;
using Complex=Scorer::Complex;
using Front=Scorer::Front;
mesh::Mesh compactionCube(){
    mesh::Mesh source;std::map<std::array<int,3>,Index> ids;
    const auto vertex=[&](std::array<int,3> p){
        if(const auto i=ids.find(p);i!=ids.end())return i->second;
        const auto id=source.addNode({double(p[0]),double(p[1]),double(p[2])},{0,int(ids.size()+1)});
        ids.emplace(p,id);return id;
    };
    for(unsigned axis=0;axis<3;++axis)for(int sign:std::array{-1,1}){
        std::array<Index,4> q{};const std::array u{-1,1,1,-1},v{-1,-1,1,1};
        for(unsigned i=0;i<4;++i){std::array<int,3> p{};p[axis]=sign;
            p[(axis+1)%3]=u[i];p[(axis+2)%3]=v[i];q[i]=vertex(p);}
        std::array<Index,3> a{q[0],q[1],q[2]},b{q[0],q[2],q[3]};
        if(sign<0){std::swap(a[1],a[2]);std::swap(b[1],b[2]);}
        source.addCell(mesh::CellType::Tri3,a,{2,int(2*axis+(sign>0?2:1))});
        source.addCell(mesh::CellType::Tri3,b,{2,int(2*axis+(sign>0?2:1))});
    }return source;
}
struct Snapshot {
    std::vector<std::array<std::uint64_t,4>> vertices;
    std::vector<std::array<std::uint64_t,13>> raw;
    std::vector<std::array<std::uint64_t,7>> live;
    std::vector<std::array<std::uint64_t,5>> faces,edges;
    std::vector<Scorer::Vertex> mapping;
    std::vector<std::size_t> anchors;
    std::uint64_t epoch=0;std::size_t accepted=0;
    bool operator==(const Snapshot&)const=default;
};
Snapshot snapshot(const Scorer& scorer){
    Snapshot out;const auto& c=scorer.complexForDiagnostics();
    const auto& relevance=scorer.relevanceForDiagnostics();
    check(relevance.size()==c.cells.size(),"relevance does not cover stored cells");
    for(const auto& vertex:c.vertices){const auto p=bits(vertex.position);
        out.vertices.push_back({p[0],p[1],p[2],vertex.sequence});}
    for(std::size_t id=0;id<c.cells.size();++id){
        const auto& cell=c.cells[id];std::array<std::uint64_t,13> raw{};unsigned at=0;
        for(auto v:cell.vertices)raw[at++]=v;
        raw[at++]=cell.sequence;raw[at++]=cell.classification.dim;raw[at++]=cell.classification.tag;
        raw[at++]=c.alive[id];raw[at++]=relevance[id];
        for(auto next:c.neighbours[id])raw[at++]=next;
        out.raw.push_back(raw);
        if(c.alive[id])out.live.push_back({cell.vertices[0],cell.vertices[1],cell.vertices[2],cell.vertices[3],
            std::uint64_t(cell.classification.dim),std::uint64_t(cell.classification.tag),relevance[id]});
    }
    for(const auto& [key,owners]:c.faces)for(auto owner:owners)
        out.faces.push_back({key[0],key[1],key[2],owner,0});
    for(const auto& [key,owners]:c.edges)for(auto owner:owners)
        out.edges.push_back({key[0],key[1],owner,0,0});
    out.mapping=scorer.auxiliaryToFront();out.anchors=c.vertexCell;
    out.epoch=scorer.planEpoch();out.accepted=scorer.acceptedSites();return out;
}
void coherent(const Scorer& scorer){
    const auto& c=scorer.complexForDiagnostics();
    check(c.cells.size()==c.neighbours.size()&&c.cells.size()==c.alive.size()&&
        c.cells.size()==c.cellBoxes.size(),"compacted parallel cell tables differ");
    std::size_t alive=0;
    for(std::size_t id=0;id<c.cells.size();++id)if(c.alive[id]){
        ++alive;const auto& cell=c.cells[id];const auto& v=cell.vertices;
        check(predicates::orient3d(c.vertices[v[0]].position,c.vertices[v[1]].position,
            c.vertices[v[2]].position,c.vertices[v[3]].position)==predicates::Sign::Positive,
            "compaction changed exact cell orientation");
        check(cell.classification==mesh::Classification{3,7},"compaction changed body classification");
        const auto bounds=c.bounds(v);
        check(bits(bounds.lo)==bits(c.cellBoxes[id].lo)&&bits(bounds.hi)==bits(c.cellBoxes[id].hi),
            "compacted cell box no longer encloses its exact vertices");
        const auto keys=Complex::faceKeys(cell);
        for(unsigned side=0;side<4;++side){const auto next=c.neighbours[id][side];
            if(next==Complex::noCell)continue;
            check(next<c.cells.size()&&c.alive[next],"compacted adjacency points to a retired cell");
            const auto other=Complex::faceKeys(c.cells[next]);
            const auto found=std::ranges::find(other,keys[side]);
            check(found!=other.end()&&c.neighbours[next][std::size_t(found-other.begin())]==id,
                "compacted adjacency is not reciprocal across the same face");
        }
    }
    check(alive==scorer.liveCells(),"compacted live count differs from live table");
    for(const auto& [key,owners]:c.faces){
        check(!owners.empty()&&owners.size()<=2,"compacted face incidence is nonmanifold");
        for(auto id:owners){check(id<c.cells.size()&&c.alive[id],"face owner retained a tombstone");
            const auto keys=Complex::faceKeys(c.cells[id]);
            check(std::ranges::find(keys,key)!=keys.end(),"face owner has different vertices");}
    }
    for(std::size_t v=0;v<c.vertexCell.size();++v){const auto id=c.vertexCell[v];
        check(id<c.cells.size()&&c.alive[id]&&std::ranges::find(c.cells[id].vertices,v)!=c.cells[id].vertices.end(),
            "compacted point anchor is not an incident live cell");}
}

mesh::Mesh tetraSurface(){
    mesh::Mesh source;
    for(const auto p:std::array{Vec3{0,0,0},Vec3{1,0,0},Vec3{0,1,0},Vec3{0,0,1}})
        source.addNode(p,{0,int(source.nodeCount()+1)});
    const Scorer::Cell cell{{0,1,2,3},0,{3,7}};
    for(unsigned side=0;side<4;++side){const auto f=bnd::tet_facet_detail::outwardFace(cell,side);
        source.addCell(mesh::CellType::Tri3,std::array<Index,3>{f[0],f[1],f[2]},{2,int(side+1)});}
    return source;
}
Front sourceFront(const mesh::Mesh& source){
    Front front;
    for(const auto p:source.nodes())front.vertices.push_back({p,front.vertices.size()+1});
    for(Index i=0;i<source.cellCount();++i){const auto* v=source.cellNodes(i);
        front.boundary.push_back({Scorer::Vertex(v[0]),Scorer::Vertex(v[1]),Scorer::Vertex(v[2])});}
    check(front.initialize(),"inheritance source front refused");return front;
}
void exercisePlcVoidGeometry(){
    const auto source=cubeSurface();auto front=sourceFront(source);
    const auto vertices=front.vertices;
    const auto boundary=front.boundary;
    const auto sourceFlags=front.sourceVertex;const auto fixed=front.fixedVertexCount;
    std::uint64_t work=0;PlcVoidGeometry plc(front,[&](std::uint64_t n){work+=n;});
    using Inclusion=PlcVoidGeometry::Inclusion;
    for(const auto p:std::array{Vec3{},Vec3{.3,-.2,.7},Vec3{-.8,.8,-.8}})
        check(plc.inside(p)==Inclusion::Inside,"exact PLC ray missed strict material");
    check(std::abs(plc.clearance({},10.)-1.)<1e-13&&plc.clearance({},.4)==.4,
        "source clearance or upper clamp is incorrect");
    check(plc.clearance({1.2,0,0},10.)>0.&&plc.inside({1.2,0,0})==Inclusion::Outside,
        "positive source clearance was mistaken for material inclusion");
    check(plc.inside({1,0,0})==Inclusion::Ambiguous,
        "an exact source-boundary point received strict material parity");
    check(std::abs(pointTriangleDistanceSquared({.25,.25,2},{0,0,0},{1,0,0},{0,1,0})-4.)<1e-13&&
        std::abs(pointTriangleDistanceSquared({2,0,0},{0,0,0},{1,0,0},{0,1,0})-1.)<1e-13,
        "source triangle distance missed a face or vertex nearest point");
    auto hollow=source;const auto inner=compactionCube();const auto offset=hollow.nodeCount();
    for(auto p:inner.nodes())hollow.addNode(p*.25,{0,20});
    for(Index i=0;i<inner.cellCount();++i){const auto* v=inner.cellNodes(i);
        hollow.addCell(mesh::CellType::Tri3,std::array<Index,3>{offset+v[0],offset+v[2],offset+v[1]},{2,20});}
    auto hollowFront=sourceFront(hollow);PlcVoidGeometry hollowPlc(hollowFront);
    check(hollowPlc.clearance({},2.)>0.&&hollowPlc.inside({})==Inclusion::Outside,
        "positive-clearance enclosed void was classified as material");
    check(hollowPlc.inside({.7,.1,.1})==Inclusion::Inside,"material around an enclosed void was rejected");
    auto tetra=sourceFront(tetraSurface());PlcVoidGeometry tetraPlc(tetra);
    check(tetraPlc.inside({.6,.6,0})==Inclusion::Outside,
        "supporting-plane extension was treated as an ambiguous source contact");
    mesh::Mesh oblique;
    const auto transform=[](Vec3 p){return Vec3{3+p.x+.25*p.y,-2+p.y+.5*p.z,5+p.z};};
    for(Index i=0;i<source.nodeCount();++i)oblique.addNode(transform(source.node(i)),source.nodeClass(i));
    for(Index i=0;i<source.cellCount();++i){const auto* v=source.cellNodes(i);
        oblique.addCell(mesh::CellType::Tri3,std::array<Index,3>{v[0],v[1],v[2]},source.cell(i).onEntity);}
    auto obliqueFront=sourceFront(oblique);PlcVoidGeometry obliquePlc(obliqueFront);
    check(obliquePlc.inside(transform({.3,-.2,.7}))==Inclusion::Inside&&
        obliquePlc.inside(transform({1.2,0,0}))==Inclusion::Outside,
        "PLC inclusion changed under an oblique affine source transform");
    unsigned calls=0;bool stopped=false;
    PlcVoidGeometry bounded(front,[&](std::uint64_t){if(++calls>=5)throw std::runtime_error("PLC-query-work-stop");});
    try{(void)bounded.inside({});}catch(const std::runtime_error&){stopped=true;}
    check(stopped,"source query swallowed an external work exception");
    check(work>0,"source query work was not charged");
    check(front.vertices.size()==vertices.size()&&front.boundary==boundary&&front.sourceVertex==sourceFlags&&
        front.fixedVertexCount==fixed&&front.cells.empty(),"source query changed source identities or occupancy");
    for(std::size_t i=0;i<vertices.size();++i)check(bits(front.vertices[i].position)==bits(vertices[i].position)&&
        front.vertices[i].sequence==vertices[i].sequence,"source query changed a stored coordinate or sequence");
}
void exerciseCanonicalOwnerScope(){
    const auto source=tetraSurface();auto front=sourceFront(source);
    check(!front.initialOwnerCanonical,"generic front enabled an owner-only canonical shortcut");
    const auto apex=Scorer::Vertex(front.vertices.size());
    front.vertices.push_back({{.125,.125,.125},front.vertices.size()+1});
    const auto f=front.boundary.front();const Scorer::Cell cone{{f[0],f[2],f[1],apex},0,{3,7}};
    const auto bounds=front.box(cone.vertices);
    check(front.boundaryAccepts(cone,bounds,0),"canonical contact fixture refused valid geometry");
    const auto sequence=front.vertices.back().sequence;
    front.vertices.back().sequence=front.vertices[apex-1].sequence;
    check(!front.boundaryAccepts(cone,bounds,0),"default contact path accepted a noncanonical vertex table");
    front.vertices.back().sequence=sequence;
    check(front.boundaryAccepts(cone,bounds,0),"restored canonical contact fixture remained invalid");
}
void exerciseMaterialInheritance(){
    const auto source=tetraSurface();auto front=sourceFront(source);
    const auto inside=[](Vec3 p){return p.x>0&&p.y>0&&p.z>0&&p.x+p.y+p.z<1;};
    front.interiorAdmission=inside;front.centroidAdmission=inside;
    Scorer scorer;check(scorer.initialize(source,{3,7},512,32ULL*1024*1024,[](std::uint64_t){},inside),
        "material inheritance initialization failed");
    const Vec3 point{.125,.125,.125};const auto apex=Scorer::Vertex(front.vertices.size());
    front.vertices.push_back({point,front.vertices.size()+1});Scorer::Cell occupancy;
    check(front.valid(Front::Front{front.boundary.front()},apex,occupancy),"material inheritance occupancy refused");
    auto plan=scorer.score(point,front);
    check(plan&&plan->materialInherited&&plan->patch.removed.size()==1,
        "interior source tetrahedron did not inherit its material classification");
    const auto before=snapshot(scorer);const auto removed=plan->patch.removed.front();
    auto& relevance=const_cast<std::vector<unsigned char>&>(scorer.relevanceForDiagnostics());
    relevance[removed]=2;const auto changed=snapshot(scorer);
    check(!scorer.commit(*plan,apex)&&snapshot(scorer)==changed,
        "changed removed material cache was trusted by inherited commit");
    relevance[removed]=1;auto stale=*plan;++stale.generation;
    check(!scorer.commit(stale,apex)&&snapshot(scorer)==before,"changed plan epoch escaped admission");
    // Mutable provisional cache metadata must not become committed truth.
    std::fill(plan->addedRelevance.begin(),plan->addedRelevance.end(),4);
    check(scorer.commit(*plan,apex),"certified inherited material patch refused");
    const auto& cached=scorer.relevanceForDiagnostics();
    for(std::size_t i=before.raw.size();i<cached.size();++i)
        check(cached[i]==1,"inherited commit trusted mutable added relevance metadata");
    coherent(scorer);front.insert(occupancy);
    check(!scorer.commit(*plan,apex),"accepted inherited plan remained current");
}
void exerciseInheritedSourceBarrier(){
    const auto source=tetraSurface();auto front=sourceFront(source);Scorer scorer;
    check(scorer.initialize(source,{3,7},512,32ULL*1024*1024,[](std::uint64_t){}),
        "source barrier scorer initialization failed");
    // Controlled fault injection through the existing const diagnostic views:
    // give both sides of one original source facet a cached material label.
    auto& live=const_cast<Complex&>(scorer.complexForDiagnostics());
    Complex replacement;replacement.vertices=live.vertices;replacement.vertices[4].position={0,0,-1};
    replacement.walkConstraints=true;replacement.maximumStoredCells=live.maximumStoredCells;
    replacement.chargeWork=live.chargeWork;
    replacement.insert({{0,1,2,3},0,{3,7}});replacement.insert({{0,2,1,4},0,{3,7}});
    live=std::move(replacement);
    auto& relevance=const_cast<std::vector<unsigned char>&>(scorer.relevanceForDiagnostics());relevance.assign(2,1);
    const Vec3 point{.25,.25,0};const auto apex=Scorer::Vertex(front.vertices.size());
    front.vertices.push_back({point,front.vertices.size()+1});auto plan=scorer.score(point,front);
    check(plan&&plan->patch.removed.size()==2&&!plan->materialInherited&&plan->guardAdded>0,
        "internal original source facet did not fall back to complete source contacts");
    auto forged=*plan;forged.materialInherited=true;
    forged.addedRelevance.assign(forged.patch.added.size(),1);forged.patch.convexBipyramid=true;
    const auto before=snapshot(scorer);const auto certifications=scorer.fullCertifications;
    check(!scorer.commit(forged,apex)&&snapshot(scorer)==before&&scorer.fullCertifications==certifications,
        "mutable inherited plan erased an internal original source facet");
    auto changedRemoved=forged;changedRemoved.patch.removed.resize(1);
    const auto validityRejects=scorer.commitValidityRejects;
    check(!scorer.commit(changedRemoved,apex)&&snapshot(scorer)==before&&
        scorer.commitValidityRejects==validityRejects+1,
        "changed removed set bypassed full selected patch certification");
    live.alive[1]=0;const auto retired=snapshot(scorer);
    check(!scorer.commit(forged,apex)&&snapshot(scorer)==retired,
        "retired removed cell was trusted by inherited commit");
    live.alive[1]=1;
}
void exerciseConicSurvivorAdmission(){
    const auto source=tetraSurface();auto front=sourceFront(source);Scorer scorer;
    check(scorer.initialize(source,{3,7},512,32ULL*1024*1024,[](std::uint64_t){}),
        "conic survivor scorer initialization failed");
    auto& live=const_cast<Complex&>(scorer.complexForDiagnostics());
    Complex replacement;replacement.vertices=live.vertices;
    // A disconnected survivor with distinct identities occupies the same tetrahedron.
    // The local positive shell check passes; exact survivor contacts must refuse it.
    for(unsigned i=0;i<4;++i)replacement.vertices[4+i].position=source.node(i);
    replacement.walkConstraints=true;replacement.maximumStoredCells=live.maximumStoredCells;
    replacement.chargeWork=live.chargeWork;
    replacement.insert({{0,1,2,3},0,{3,7}});replacement.insert({{4,5,6,7},0,{3,7}});
    live=std::move(replacement);
    auto& relevance=const_cast<std::vector<unsigned char>&>(scorer.relevanceForDiagnostics());relevance={1,3};
    const Vec3 point{.125,.125,.125};const auto apex=Scorer::Vertex(front.vertices.size());
    front.vertices.push_back({point,front.vertices.size()+1});auto plan=scorer.score(point,front);
    check(plan&&plan->materialInherited&&plan->patch.removed==std::vector<std::size_t>{0},
        "conic survivor fixture did not isolate its local replacement");
    plan->patch.convexBipyramid=true;
    const auto before=snapshot(scorer);const auto rejects=scorer.commitSurvivorRejects;
    check(!scorer.commit(*plan,apex)&&snapshot(scorer)==before&&scorer.commitSurvivorRejects==rejects+1,
        "mutable bipyramid flag bypassed conic pair/survivor certification");
}

void exerciseCompaction(){
    const auto source=compactionCube();Front front;
    for(const auto p:source.nodes())front.vertices.push_back({p,front.vertices.size()+1});
    for(Index i=0;i<source.cellCount();++i){const auto* v=source.cellNodes(i);
        front.boundary.push_back({Scorer::Vertex(v[0]),Scorer::Vertex(v[1]),Scorer::Vertex(v[2])});}
    check(front.initialize(),"owned cube source refused");
    enum class Failure{None,Work,Allocation};Failure failure=Failure::None;
    unsigned imports=0;std::uint64_t charged=0,native=0;
    const auto charge=[&](std::uint64_t n){
        if(failure!=Failure::None&&n==11&&++imports==2){
            if(failure==Failure::Allocation)throw std::bad_alloc{};
            throw std::runtime_error("partial-clone-work-stop");
        }check(n<=100'000'000-charged,"compaction test work bound exceeded");charged+=n;
    };
    const auto inside=[&](Vec3 p){++native;return std::abs(p.x)<1&&std::abs(p.y)<1&&std::abs(p.z)<1;};
    front.chargeWork=charge;front.interiorAdmission=inside;front.centroidAdmission=inside;
    Scorer scorer;check(scorer.initialize(source,{3,7},512,32ULL*1024*1024,charge,inside,0),
        "compaction scorer initialization failed");
    const auto capacity=scorer.storageCapacity();const auto admittedBytes=scorer.workingBytes();
    const auto first=Scorer::Vertex(front.vertices.size());front.vertices.push_back({{},first+1});
    bnd::TetFacetCell occupancy;
    check(front.valid(Front::Front{front.boundary[0]},first,occupancy),"first occupancy cone refused");
    auto initial=scorer.score({},front);check(initial&&scorer.commit(*initial,first),"first auxiliary edit refused");
    front.insert(occupancy);
    check(scorer.storedCells()>scorer.liveCells(),"fixture has no retired auxiliary cells");
    const Vec3 point{.35,.1,.1};const auto second=Scorer::Vertex(front.vertices.size());
    front.vertices.push_back({point,second+1});
    auto stale=scorer.score(point,front);front.vertices.pop_back();
    check(stale&&stale->generation==scorer.planEpoch(),"stale-plan fixture could not be staged");
    const auto before=snapshot(scorer);const auto queries=native;
    const auto workBefore=charged,compactionBefore=scorer.compactionWork;
    for(auto mode:std::array{Failure::Work,Failure::Allocation}){
        failure=mode;imports=0;bool threw=false;
        try{(void)scorer.prepareForBatch(true);}
        catch(const std::bad_alloc&){check(mode==Failure::Allocation,"unexpected allocation exception");threw=true;}
        catch(const std::runtime_error&){check(mode==Failure::Work,"unexpected work exception");threw=true;}
        failure=Failure::None;
        check(threw&&imports==2&&snapshot(scorer)==before,
            "partially constructed compaction clone changed the old complex/cache/epoch");
    }
    check(scorer.compactions==0&&scorer.prepareForBatch(true),"forced compaction did not complete");
    const auto after=snapshot(scorer);
    check(after.vertices==before.vertices&&after.mapping==before.mapping&&after.live==before.live,
        "compaction changed point identities, oriented live cells, body or relevance cache");
    check(after.accepted==before.accepted&&after.epoch==before.epoch+1&&
        scorer.storedCells()==scorer.liveCells()&&scorer.compactions==1&&
        scorer.deadCellsReclaimed==before.raw.size()-before.live.size(),"compaction lifetime counters are inconsistent");
    check(scorer.storageCapacity()==capacity&&scorer.workingBytes()==admittedBytes&&
        scorer.maximumCompactionModeledBytes>0&&scorer.maximumCompactionModeledBytes<=admittedBytes,
        "compaction raised capacity or exceeded its admitted modeled payload");
    check(native==queries&&charged-workBefore==scorer.compactionWork-compactionBefore,
        "compaction repeated native classification or escaped work accounting");
    coherent(scorer);
    check(!scorer.commit(*stale,second)&&snapshot(scorer)==after,
        "pre-compaction plan survived cell-ID renumbering");
    front.vertices.push_back({point,second+1});
    auto fresh=scorer.score(point,front);
    check(fresh&&fresh->generation==scorer.planEpoch()&&
        front.valid(Front::Front{front.boundary[2]},second,occupancy)&&scorer.commit(*fresh,second),
        "fresh post-compaction query/occupancy/selected admission failed");
    front.insert(occupancy);coherent(scorer);
    check(scorer.acceptedSites()==2,"compaction altered accepted-point generation");
    const auto disabled=snapshot(scorer);scorer.compactionEnabled=false;
    check(!scorer.prepareForBatch(true)&&snapshot(scorer)==disabled,"disabled compaction changed state");
    std::cout<<"compaction PASS reclaimed="<<scorer.deadCellsReclaimed<<" work="<<scorer.compactionWork
        <<" attempts="<<scorer.compactionAttempts<<" native="<<native<<'\n';
}

void exerciseBoundedBootstrapAndCancellation(){
    using Scorer=OcctInitialDelaunayScorer;
    const auto source=compactionCube(),unchanged=source;
    constexpr std::size_t cells=512;
    const auto preparation=std::uint64_t(source.nodeCount()+source.cellCount());
    const auto importReserve=std::uint64_t(source.nodeCount()+8)+12*cells;
    for(unsigned which=0;which<3;++which){
        Scorer scorer;std::uint64_t charged=0;
        const auto work=which==0?preparation+importReserve+1:1'000'000;
        const auto capacity=which==1?std::size_t(1):cells;
        const auto bytes=which==2?10ULL*1024*1024:64ULL*1024*1024;
        const bool initialized=scorer.initialize(source,{3,7},capacity,bytes,
            [&](std::uint64_t n){check(n<=work-charged,"bounded bootstrap overcharged work");charged+=n;},
            {},0,work);
        check(!initialized&&scorer.budgetLimited&&!scorer.cancelled&&
            scorer.complexForDiagnostics().vertices.empty()&&scorer.acceptedSites()==0,
            "bounded bootstrap refusal imported a partial candidate");
        check(charged==preparation+scorer.bootstrapWork&&scorer.bootstrapWork<=scorer.bootstrapWorkLimit&&
            scorer.bootstrapPeakTemporaryBytes<=scorer.bootstrapTemporaryLimitBytes,
            "bounded bootstrap work was omitted or counted twice");
        if(which==0)check(scorer.bootstrapWork<=1&&scorer.bootstrapPeakTemporaryBytes>=4ULL*1024*1024,
            "tiny bootstrap work cap did not stop the online bounded constructor");
        if(which==2)check(scorer.bootstrapTemporaryLimitBytes<4ULL*1024*1024,
            "small-memory fixture did not exercise the bounded workspace refusal");
        originalPrefix(unchanged,source);
    }
    // The final preparation charge arms cancellation. The bounded constructor
    // observes it, then the accounting callback resets the atomic before the
    // initializer can react. The observed Cancelled outcome must remain sticky.
    std::atomic_bool cancellation{false};std::uint64_t preparationCharges=0,charged=0;
    bool reset=false,caught=false;Scorer cancelled;
    try{
        (void)cancelled.initialize(source,{3,7},cells,64ULL*1024*1024,
            [&](std::uint64_t n){
                charged+=n;
                if(n==1&&++preparationCharges==preparation)cancellation.store(true);
                if(n==0&&cancellation.load()){cancellation.store(false);reset=true;}
            },{},0,1'000'000,&cancellation);
    }catch(const OcctInitialConstructionCancelled&){caught=true;}
    check(caught&&reset&&!cancellation.load()&&cancelled.cancelled&&
        !cancelled.initialized()&&cancelled.complexForDiagnostics().vertices.empty()&&charged==preparation,
        "resetting the atomic erased a bounded bootstrap cancellation");
    OcctInitialPointFrontLimits limits;limits.cancellation=&cancellation;
    cancellation.store(true);OcctInitialPointFrontReport report;std::uint64_t calls=0;
    const UniformField field;
    bool cancelledBeforeReturn=false;
    try{(void)appendOcctInitialPointFrontSites(source,field,field.h,{3,7},limits,
        [&](Vec3){++calls;return true;},[](std::uint64_t){},&report);
    }catch(const OcctInitialConstructionCancelled&){cancelledBeforeReturn=true;}
    cancellation.store(false);originalPrefix(unchanged,source);
    check(cancelledBeforeReturn&&report.cancelled&&!report.auxiliaryInitialized&&
        report.sites==0&&calls==0,
        "cancelled initial construction was mistaken for a successful fallback");
    std::cout<<"bounded initial bootstrap refusals and sticky cancellation passed\n";
}

template<class Run>void withCertifiedSource(const mesh::Mesh& source,Run run){
    auto front=sourceFront(source);std::atomic_bool cancellation{false};
    std::function<void(std::uint64_t)> hook;
    Scorer scorer;
    check(scorer.initialize(source,{3,7},2048,64ULL*1024*1024,[&](std::uint64_t n){if(hook)hook(n);},
        {},0,4'000'000,&cancellation),"certificate fixture bootstrap refused");
    check(!scorer.useCertifiedConeCavities,"generic scorer enabled a trusted-owner certificate");
    scorer.useCertifiedConeCavities=true;front.initialOwnerCanonical=true;
    run(front,scorer,hook,cancellation);
}
Scorer::Plan stageCertificatePoint(Front& front,Scorer& scorer,Vec3 point){
    front.vertices.push_back({point,front.vertices.size()+1});
    auto plan=scorer.score(point,front);check(plan&&plan->valid,"certificate point proposal unavailable");
    return std::move(*plan);
}
void acceptCertificatePoint(Front& front,Scorer& scorer,Vec3 point){
    auto plan=stageCertificatePoint(front,scorer,point);
    const auto certificates=scorer.coneCertificates,fallbacks=scorer.coneCertificateFallbacks;
    check(scorer.commit(plan,Scorer::Vertex(front.vertices.size()-1)),"certified geometric insertion refused");
    check(scorer.coneCertificates==certificates+1&&scorer.coneCertificateFallbacks==fallbacks,
        "complete cone fan did not receive its exact local certificate");
    coherent(scorer);
}
void exerciseCertificateDegeneracies(){
    withCertifiedSource(tetraSurface(),[](Front& front,Scorer& scorer,auto&,auto&){
        // These sites lie exactly on an old face and edge. Visibility growth
        // must consume every zero shell cone before the local proof succeeds.
        acceptCertificatePoint(front,scorer,{.25,.25,0});
        acceptCertificatePoint(front,scorer,{.5,0,0});
        acceptCertificatePoint(front,scorer,{.125,.125,std::ldexp(1.,-40)});
        const auto before=snapshot(scorer);const auto validityRejects=scorer.commitValidityRejects;
        auto duplicate=stageCertificatePoint(front,scorer,{0,0,0});
        const auto cached=snapshot(scorer);
        check(!scorer.commit(duplicate,Scorer::Vertex(front.vertices.size()-1))&&snapshot(scorer)==cached&&
            scorer.commitValidityRejects==validityRejects+1,
            "duplicate site replaced an original vertex identity");
        check(before.accepted==scorer.acceptedSites(),"duplicate site advanced the accepted generation");
    });
    withCertifiedSource(compactionCube(),[](Front& front,Scorer& scorer,auto&,auto&){
        // The original cube and these exact dyadic sites repeatedly expose
        // co-spherical and coplanar configurations without symbolic perturbation.
        for(const auto p:std::array{Vec3{0,0,0},Vec3{.5,0,0},Vec3{0,.5,0},Vec3{0,0,.5},
            Vec3{-.5,0,0},Vec3{0,-.5,0},Vec3{0,0,-.5}})acceptCertificatePoint(front,scorer,p);
    });
}
void exerciseActualPatchCertification(){
    for(unsigned mutation=0;mutation<5;++mutation)withCertifiedSource(tetraSurface(),[&](Front& front,Scorer& scorer,auto&,auto&){
        auto plan=stageCertificatePoint(front,scorer,{.125,.125,.125});
        switch(mutation){
        case 0:std::swap(plan.patch.added.front().vertices[1],plan.patch.added.front().vertices[2]);break;
        case 1:plan.patch.added.pop_back();plan.addedRelevance.pop_back();break;
        case 2:plan.patch.added.push_back(plan.patch.added.front());plan.addedRelevance.push_back(1);break;
        case 3:{const auto& live=scorer.complexForDiagnostics();
            for(std::size_t id=0;id<live.alive.size();++id)if(live.alive[id]&&id!=plan.patch.removed.front()){
                plan.patch.removed.push_back(id);break;}break;}
        default:plan.point={100,100,100};break;
        }
        const auto before=snapshot(scorer);
        check(!scorer.commit(plan,Scorer::Vertex(front.vertices.size()-1))&&snapshot(scorer)==before,
            "mutated actual point/base/cone/removed set escaped complete selected admission");
    });
    withCertifiedSource(tetraSurface(),[](Front& front,Scorer& scorer,auto&,auto&){
        auto plan=stageCertificatePoint(front,scorer,{.125,.125,.125});
        // A cyclically reordered base has the same geometry, but differs from
        // the exact builder quartet. It must use the full fallback certificate.
        auto& v=plan.patch.added.front().vertices;std::rotate(v.begin()+1,v.begin()+2,v.end());
        const auto certificates=scorer.coneCertificates,fallbacks=scorer.coneCertificateFallbacks;
        check(scorer.commit(plan,Scorer::Vertex(front.vertices.size()-1))&&
            scorer.coneCertificates==certificates&&scorer.coneCertificateFallbacks==fallbacks+1,
            "unsupported quartet ordering bypassed the complete fallback proof");
        coherent(scorer);
    });
}
void exerciseTouchedIncidenceFallback(){
    for(bool corruptOwner:std::array{false,true})withCertifiedSource(tetraSurface(),[&](Front& front,Scorer& scorer,auto& hook,auto&){
        auto plan=stageCertificatePoint(front,scorer,{.125,.125,.125});
        auto& live=const_cast<Complex&>(scorer.complexForDiagnostics());
        const auto id=plan.patch.removed.front();const auto key=Complex::faceKeys(live.cells[id])[0];
        const auto oldNeighbour=live.neighbours[id][0];const auto oldOwners=live.faces.at(key);
        check(oldNeighbour!=Complex::noCell&&oldOwners.size()==2,"local reciprocity fixture has no survivor neighbor");
        if(corruptOwner)live.faces.at(key).erase(id);else live.neighbours[id][0]=Complex::noCell;
        const auto before=snapshot(scorer);const auto fallbacks=scorer.coneCertificateFallbacks;
        hook=[&](std::uint64_t){if(scorer.coneCertificateFallbacks>fallbacks)throw std::runtime_error("fallback-stop");};
        bool threw=false;try{(void)scorer.commit(plan,Scorer::Vertex(front.vertices.size()-1));}
        catch(const std::runtime_error&){threw=true;}hook={};
        check(threw&&scorer.coneCertificateFallbacks==fallbacks+1&&snapshot(scorer)==before,
            "local incidence/reciprocity fault neither fell back nor preserved state on refusal");
        live.neighbours[id][0]=oldNeighbour;live.faces.at(key)=oldOwners;
        check(scorer.commit(plan,Scorer::Vertex(front.vertices.size()-1)),"restored local ownership remained invalid");
        coherent(scorer);
    });
}
void exerciseCertificateTransactions(){
    for(unsigned failure=0;failure<3;++failure)withCertifiedSource(tetraSurface(),[&](Front& front,Scorer& scorer,auto& hook,auto& cancellation){
        auto plan=stageCertificatePoint(front,scorer,{.125,.125,.125});
        const auto before=snapshot(scorer);const auto attempts=scorer.coneCertificateAttempts;
        hook=[&](std::uint64_t){if(scorer.coneCertificateAttempts==attempts)return;
            if(failure==0)throw std::runtime_error("certificate-work-stop");
            if(failure==1)throw std::bad_alloc{};
            cancellation.store(true);
        };
        bool threw=false;try{(void)scorer.commit(plan,Scorer::Vertex(front.vertices.size()-1));}
        catch(const OcctInitialConstructionCancelled&){check(failure==2,"unexpected certificate cancellation");threw=true;}
        catch(const std::bad_alloc&){check(failure==1,"unexpected certificate allocation refusal");threw=true;}
        catch(const std::runtime_error&){check(failure==0,"unexpected certificate work refusal");threw=true;}
        hook={};check(threw&&snapshot(scorer)==before&&scorer.coneCertificateAttempts==attempts+1,
            "certificate work/allocation/cancellation refusal changed auxiliary topology or cache");
        if(failure==2)check(scorer.cancelled,"certificate cancellation observation was not sticky");
    });
}
}
int main(){try{
    exerciseUniformFront();exerciseFineBoundaryInteriorExclusion();exerciseRefusals();exerciseCompleteConflictCavity();exerciseCompaction();
    exercisePlcVoidGeometry();
    exerciseCanonicalOwnerScope();exerciseMaterialInheritance();exerciseInheritedSourceBarrier();exerciseConicSurvivorAdmission();
    exerciseBoundedBootstrapAndCancellation();
    exerciseCertificateDegeneracies();exerciseActualPatchCertification();
    exerciseTouchedIncidenceFallback();exerciseCertificateTransactions();
    std::cout<<"initial volume points, complete conflict cavities and transactional guards passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
