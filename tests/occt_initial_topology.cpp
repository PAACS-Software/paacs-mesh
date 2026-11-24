// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctInitialPointFront.hpp"
#include "../src/bnd/TetConstraintRecovery.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

using namespace paacs;
using namespace paacs::geo::detail;
namespace {
using Scorer=OcctInitialDelaunayScorer;
using Carrier=bnd::detail::TetNumericalTopology;
using Complex=Scorer::Complex;
using Vertex=Scorer::Vertex;
using Quartet=std::array<Vertex,4>;
constexpr std::uint64_t MiB=1024ULL*1024;
constexpr mesh::Classification region{3,7};
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
auto bits(Vec3 p){return std::array{std::bit_cast<std::uint64_t>(p.x),
    std::bit_cast<std::uint64_t>(p.y),std::bit_cast<std::uint64_t>(p.z)};}
bool inside(Vec3 p){return std::abs(p.x)<1&&std::abs(p.y)<1&&std::abs(p.z)<1;}

mesh::Mesh cube(unsigned divisions=1){
    mesh::Mesh result;std::map<std::array<unsigned,3>,Index> ids;
    const auto node=[&](std::array<unsigned,3> p){
        if(const auto it=ids.find(p);it!=ids.end())return it->second;
        const auto id=result.addNode({-1.+2.*p[0]/divisions,-1.+2.*p[1]/divisions,
            -1.+2.*p[2]/divisions},{0,int(ids.size()+1)});
        ids.emplace(p,id);return id;
    };
    for(unsigned axis=0;axis<3;++axis)for(unsigned side=0;side<2;++side)
        for(unsigned u=0;u<divisions;++u)for(unsigned v=0;v<divisions;++v){
            std::array<Index,4> q{};const std::array du{0U,1U,1U,0U},dv{0U,0U,1U,1U};
            for(unsigned j=0;j<4;++j){std::array<unsigned,3> p{};p[axis]=side*divisions;
                p[(axis+1)%3]=u+du[j];p[(axis+2)%3]=v+dv[j];q[j]=node(p);}
            std::array<Index,3> a{q[0],q[1],q[2]},b{q[0],q[2],q[3]};
            if(!side){std::swap(a[1],a[2]);std::swap(b[1],b[2]);}
            result.addCell(mesh::CellType::Tri3,a,{2,int(2*axis+side+1)});
            result.addCell(mesh::CellType::Tri3,b,{2,int(2*axis+side+1)});
        }
    return result;
}

Scorer::Front frontFor(const mesh::Mesh& source){
    Scorer::Front front;
    for(const auto p:source.nodes())front.vertices.push_back({p,front.vertices.size()+1});
    for(Index i=0;i<source.cellCount();++i){const auto* v=source.cellNodes(i);
        front.boundary.push_back({Vertex(v[0]),Vertex(v[1]),Vertex(v[2])});}
    require(front.initialize(),"owned cube front is not oriented and closed");
    return front;
}

std::vector<std::uint64_t> carrierDigest(const Carrier& carrier){
    std::vector<std::uint64_t> values;
    const auto add=[&](auto v){values.push_back(std::uint64_t(v));};
    add(carrier.available);add(carrier.sourceNodes);add(carrier.liveCells);
    add(carrier.region.dim);add(carrier.region.tag);add(carrier.originalCoordinateRetryUsed);
    const auto& c=carrier.complex;add(c.vertices.size());add(c.cells.size());
    add(c.maximumStoredCells);add(c.walkConstraints);add(bool(c.chargeWork));
    for(const auto& v:c.vertices){for(auto x:bits(v.position))add(x);add(v.sequence);}
    for(const auto& cell:c.cells){for(auto v:cell.vertices)add(v);
        add(cell.sequence);add(cell.classification.dim);add(cell.classification.tag);}
    for(auto x:c.alive)add(x);
    for(const auto& row:c.neighbours)for(auto x:row)add(x);
    for(const auto& [f,owners]:c.faces){for(auto x:f)add(x);add(owners.size());for(auto x:owners)add(x);}
    for(const auto& [e,owners]:c.edges){for(auto x:e)add(x);add(owners.size());for(auto x:owners)add(x);}
    for(const auto& [key,directed]:carrier.sourceFaces){for(auto x:key)add(x);for(auto x:directed)add(x);}
    return values;
}

struct Fixture {mesh::Mesh source,placed;Carrier carrier;std::uint64_t constructionWork=0;};
Fixture fixture(){
    Fixture result;result.source=cube();result.placed=result.source;
    Scorer scorer;std::atomic_bool cancelled=false;
    const auto charge=[&](std::uint64_t n){require(n<=100'000'000-result.constructionWork,
        "fixture exceeded its construction work cap");result.constructionWork+=n;};
    require(scorer.initialize(result.source,region,4096,64*MiB,charge,inside,0,500'000,&cancelled,true),
        "owned initial topology did not bootstrap");
    auto front=frontFor(result.source);front.chargeWork=charge;
    const std::array<Vec3,4> points{{{.13,.21,.37},{-.27,.31,-.19},{.33,-.29,.11},{-.17,-.23,-.41}}};
    unsigned accepted=0;
    for(const auto p:points){
        const auto id=Vertex(front.vertices.size());front.vertices.push_back({p,front.vertices.size()+1});
        const auto plan=scorer.score(p,front);
        if(plan&&plan->valid&&scorer.commit(*plan,id)){result.placed.addNode(p,region);++accepted;}
        else front.vertices.pop_back();
    }
    require(accepted>=2,"owned fixture did not exercise admitted sites after exterior guards");
    const auto before=scorer.complexForDiagnostics().cells.size();
    require(scorer.releaseNumericalTopology(result.carrier,region)&&!scorer.initialized(),
        "certified owner did not release its numerical topology");
    require(result.carrier.available&&!result.carrier.complex.chargeWork&&
        result.carrier.complex.cells.size()==before&&result.carrier.sourceNodes==std::size_t(result.source.nodeCount()),
        "handoff retained an owner callback or lost source/cell identity");
    return result;
}

std::vector<Vec3> recoveryPoints(const mesh::Mesh& placed){
    std::vector<Vec3> points(placed.nodes().begin(),placed.nodes().end());
    // The known cube envelope remains unchanged by strictly interior sites.
    for(unsigned i=0;i<8;++i)points.push_back({(i&1)?3.:-3.,(i&2)?3.:-3.,(i&4)?3.:-3.});
    return points;
}

// Rebuild expected incidence independently from ordered live cell corners.
void coherent(const Complex& c){
    std::map<Complex::Face,std::set<std::size_t>> faces;
    std::map<Complex::Edge,std::set<std::size_t>> edges;
    require(c.cells.size()==c.alive.size()&&c.cells.size()==c.neighbours.size(),"incomplete adjacency arrays");
    for(std::size_t id=0;id<c.cells.size();++id)if(c.alive[id]){
        const auto& cell=c.cells[id];const auto& v=cell.vertices;
        for(auto x:v)require(x<c.vertices.size(),"live cell contains an invalid dense point ID");
        require(cell.classification==region&&predicates::orient3d(c.vertices[v[0]].position,c.vertices[v[1]].position,
            c.vertices[v[2]].position,c.vertices[v[3]].position)==predicates::Sign::Positive,
            "handoff changed body ownership or exact positive orientation");
        for(unsigned side=0;side<4;++side){
            std::array<Vertex,3> face{};unsigned next=0;
            for(unsigned j=0;j<4;++j)if(j!=side)face[next++]=v[j];
            std::sort(face.begin(),face.end());faces[face].insert(id);
        }
        for(unsigned a=0;a<4;++a)for(unsigned b=a+1;b<4;++b){
            auto edge=Complex::Edge{v[a],v[b]};std::sort(edge.begin(),edge.end());edges[edge].insert(id);
        }
    }
    require(c.faces==faces&&c.edges==edges,"imported incidence differs from live ordered geometry");
    for(const auto& [face,owners]:faces){
        require(!owners.empty()&&owners.size()<=2,"imported topology is nonmanifold");
        for(auto id:owners){
            const auto& vertices=c.cells[id].vertices;unsigned opposite=0;
            while(opposite<4&&std::ranges::find(face,vertices[opposite])!=face.end())++opposite;
            require(opposite<4,"independent opposite-corner lookup failed");
            const auto expected=owners.size()==1?Complex::noCell:(*owners.begin()==id?*owners.rbegin():*owners.begin());
            require(c.neighbours[id][opposite]==expected,"reciprocal adjacency changed at handoff");
        }
    }
}

void orderedImportAndLifetime(){
    auto f=fixture();const auto sourceHash=mesh::meshHash(f.placed);
    auto topology=std::move(f.carrier);coherent(topology.complex);
    const auto original=topology.sourceNodes,ns=std::size_t(f.placed.nodeCount());
    const auto remap=[&](Vertex id){return id<original?id:id<original+8?Vertex(ns+id-original):Vertex(id-8);};
    std::vector<Quartet> expected;
    for(std::size_t id=0;id<topology.complex.cells.size();++id)if(topology.complex.alive[id]){
        auto v=topology.complex.cells[id].vertices;for(auto& x:v)x=remap(x);expected.push_back(v);
    }
    const auto points=recoveryPoints(f.placed);const auto faces=topology.sourceFaces;
    const auto oldRetained=topology.retainedModeledBytes(),live=topology.liveCells;
    Complex imported;imported.walkConstraints=true;imported.maximumStoredCells=4096;
    bnd::detail::TetConstraintRecoveryReport report;
    imported.chargeWork=[&](std::uint64_t n){require(n<=100'000'000-report.work,"linear import exceeded work cap");report.work+=n;};
    const auto diagnostic=bnd::detail::constraint_detail::importNumericalTopology(topology,f.placed,region,
        faces,points,imported,128*MiB,32*MiB,report);
    require(!diagnostic&&report.initialTopologyUsed&&std::string_view(report.initialTopologyDiagnostic)=="preserved"&&
        report.initialTopologyCells==live&&!topology.available,"topology was not consumed after complete staged import");
    require(report.initialTopologyImportWork==report.work&&report.work>0&&
        report.initialTopologyPeakModeledBytes==32*MiB+oldRetained+4096*live,
        "import repeated bootstrap work or omitted old/new simultaneous storage");
    require(imported.vertices.size()==points.size()&&imported.cells.size()==expected.size(),"import changed point/cell count");
    for(std::size_t i=0;i<points.size();++i)require(bits(imported.vertices[i].position)==bits(points[i])&&
        imported.vertices[i].sequence==i+1,"bit-exact ordered point identity was not preserved");
    for(std::size_t i=0;i<expected.size();++i)require(imported.cells[i].vertices==expected[i],
        "handoff retriangulated or reordered an oriented live tetrahedron");
    coherent(imported);require(mesh::meshHash(f.placed)==sourceHash,"import mutated prescribed geometry");
    std::cout<<"ordered handoff: points="<<points.size()<<" live="<<live<<" import_work="<<report.work<<'\n';
}

mesh::Mesh held(){mesh::Mesh result;result.addNode({11,12,13},{0,91});return result;}
auto recover(const mesh::Mesh& source,Carrier& carrier,mesh::Mesh& output,std::uint64_t work=20'000'000,
             std::uint64_t bytes=128*MiB,std::uint64_t cells=100'000,const std::atomic_bool* stop=nullptr){
    return bnd::detail::recoverTetConstraints(source,region,work,bytes,cells,stop,output,
        0,false,true,true,true,&carrier);
}

void checkSeed(const mesh::Mesh& placed,const mesh::Mesh& seed){
    require(seed.nodeCount()==placed.nodeCount()&&seed.cellCount()>0,"private seed leaked guards or lost admitted points");
    for(Index i=0;i<placed.nodeCount();++i)require(bits(seed.node(i))==bits(placed.node(i))&&
        seed.nodeClass(i)==placed.nodeClass(i),"recovery changed a source or admitted-site identity");
    for(const auto p:seed.nodes())require(std::abs(p.x)<=1&&std::abs(p.y)<=1&&std::abs(p.z)<=1,
        "private exterior guard was published in the material seed");
    const auto quality=algo::analyzeMixedMeshQuality(seed);
    require(quality.valid&&std::abs(quality.totalVolume-8.)<1e-10,"handoff recovery changed cube material coverage");
    auto audited=seed;
    for(Index i=0;i<placed.cellCount();++i){const auto* v=placed.cellNodes(i);
        audited.addCell(mesh::CellType::Tri3,std::array<Index,3>{v[0],v[1],v[2]},placed.cell(i).onEntity);}
    require(mesh::auditLinearVolumeBoundary(audited).valid,"handoff seed does not conform to every original tagged facet");
}

void completeRecoveryAndGuards(){
    auto f=fixture();const auto sourceHash=mesh::meshHash(f.placed);mesh::Mesh output;
    struct Observed {unsigned topology=0,delaunay=0;} observed;
    const auto observer=[](const bnd::detail::TetConstructionObservation& row,void* context)noexcept{
        auto& count=*static_cast<Observed*>(context);
        if(std::string_view(row.stage)=="recovery_topology_import")++count.topology;
        if(std::string_view(row.stage)=="recovery_delaunay_import")++count.delaunay;
    };
    bnd::detail::ScopedTetConstructionObserver watch(observer,&observed);
    const auto result=recover(f.placed,f.carrier,output);
    require(std::string_view(result.diagnostic)=="complete_seed"&&result.initialTopologyUsed&&
        result.initialTopologyImportWork>0&&result.missingFacesAfter==0&&result.work<=20'000'000,
        "preserved topology did not reach complete constrained recovery");
    require(observed.topology==1&&observed.delaunay==0,
        "handoff recovery rebuilt the point-cloud triangulation");
    checkSeed(f.placed,output);require(mesh::meshHash(f.placed)==sourceHash,"recovery mutated the source");
}

void invalidCarrierRefusals(){
    const auto f=fixture();const auto originalSource=mesh::meshHash(f.placed);
    const auto check=[&](auto alter,const char* expected){
        auto carrier=f.carrier;alter(carrier);const auto before=carrierDigest(carrier);
        auto output=held();const auto old=mesh::meshHash(output);const auto result=recover(f.placed,carrier,output);
        require(std::string_view(result.diagnostic)==expected&&!result.initialTopologyUsed&&
            mesh::meshHash(output)==old&&carrierDigest(carrier)==before,
            "corrupt topology was accepted, reconstructed, consumed or partially published");
    };
    check([](Carrier& c){c.complex.vertices[0].position.x=std::nextafter(c.complex.vertices[0].position.x,0.);},
        "topology_point_mismatch");
    check([](Carrier& c){std::swap(c.sourceFaces.begin()->second[0],c.sourceFaces.begin()->second[1]);},
        "topology_source_mismatch");
    check([](Carrier& c){c.region.tag=8;},"topology_contract_mismatch");
    check([](Carrier& c){for(std::size_t i=0;i<c.complex.cells.size();++i)if(c.complex.alive[i]){
        c.complex.cells[i].classification.tag=8;break;}},"topology_region_mismatch");
    check([](Carrier& c){for(std::size_t i=0;i<c.complex.cells.size();++i)if(c.complex.alive[i]){
        c.complex.cells[i].vertices[0]=bnd::kInvalidTetFacetVertex;break;}},"topology_connectivity_mismatch");
    check([](Carrier& c){for(std::size_t i=0;i<c.complex.cells.size();++i)if(c.complex.alive[i]){
        std::swap(c.complex.cells[i].vertices[0],c.complex.cells[i].vertices[1]);break;}},"topology_orientation_mismatch");
    check([](Carrier& c){for(std::size_t i=0;i<c.complex.cells.size();++i)if(c.complex.alive[i]){
        for(auto& other:c.complex.neighbours[i])if(other!=Complex::noCell){other=Complex::noCell;return;}}},
        "topology_adjacency_mismatch");
    require(mesh::meshHash(f.placed)==originalSource,"failed carrier validation mutated its source");
}

void resourceAndCancellationRefusals(){
    const auto f=fixture();const auto check=[&](std::uint64_t work,std::uint64_t bytes,std::uint64_t cells,
        const std::atomic_bool* stop,const char* expected){
        auto carrier=f.carrier;const auto before=carrierDigest(carrier);auto output=held();const auto old=mesh::meshHash(output);
        const auto result=recover(f.placed,carrier,output,work,bytes,cells,stop);
        require(std::string_view(result.diagnostic)==expected&&result.work<=work&&mesh::meshHash(output)==old,
            "resource/cancellation refusal overcharged or published a partial seed");
        if(!result.initialTopologyUsed)require(carrierDigest(carrier)==before,"refused import consumed its old owner");
        return result;
    };
    const auto importCap=4*std::uint64_t(f.placed.cellCount())+std::uint64_t(f.placed.nodeCount())+
        1+f.carrier.sourceFaces.size()+f.carrier.complex.vertices.size()+1;
    const auto stopped=check(importCap,128*MiB,100'000,nullptr,"work_budget");
    require(stopped.initialTopologyImportWork>0&&!stopped.initialTopologyUsed,
        "small work allowance did not exercise a transactional stop during handoff");
    check(20'000'000,8*MiB,100'000,nullptr,"memory_budget");
    check(20'000'000,128*MiB,1,nullptr,"cell_budget");
    std::atomic_bool cancel=true;
    const auto cancelled=check(20'000'000,128*MiB,100'000,&cancel,"cancelled");
    require(cancelled.work==0,"pre-cancelled handoff performed charged work");

    // Exercise the overlap check itself, separately from the public64MiB floor.
    auto carrier=f.carrier;const auto before=carrierDigest(carrier);const auto points=recoveryPoints(f.placed);
    const auto faces=carrier.sourceFaces;Complex staged;staged.maximumStoredCells=4096;staged.walkConstraints=true;
    bnd::detail::TetConstraintRecoveryReport report;staged.chargeWork=[&](std::uint64_t n){report.work+=n;};
    const auto bytes=32*MiB+carrier.retainedModeledBytes()+4096*carrier.liveCells-1;
    const auto reason=bnd::detail::constraint_detail::importNumericalTopology(carrier,f.placed,region,faces,points,
        staged,bytes,32*MiB,report);
    require(reason&&std::string_view(reason)=="topology_memory_budget"&&staged.cells.empty()&&
        carrierDigest(carrier)==before,"simultaneous old/new import memory was not admitted transactionally");
    if constexpr(sizeof(std::size_t)==sizeof(std::uint64_t)){
        Carrier huge;huge.complex.maximumStoredCells=(std::numeric_limits<std::size_t>::max)();
        require(huge.retainedModeledBytes()==UINT64_MAX,"large retained capacity wrapped its modeled-byte counter");
    }
}

void exportRefusals(){
    const auto source=cube();Scorer scorer;std::atomic_bool stop=false;std::uint64_t work=0;
    require(scorer.initialize(source,region,4096,64*MiB,[&](std::uint64_t n){work+=n;},inside,
        0,500'000,&stop,true),"export-refusal owner did not initialize");
    const auto beforeWork=work,oldCells=scorer.storedCells();
    Carrier empty;stop=true;bool cancelled=false;
    try{(void)scorer.releaseNumericalTopology(empty,region);}catch(const OcctInitialConstructionCancelled&){cancelled=true;}
    require(cancelled&&!empty.available&&scorer.initialized()&&scorer.storedCells()==oldCells&&work==beforeWork,
        "cancelled export moved geometry or performed additional construction work");
    // Cancellation is intentionally sticky for an existing scorer lifetime.
    stop=false;cancelled=false;
    try{(void)scorer.releaseNumericalTopology(empty,region);}catch(const OcctInitialConstructionCancelled&){cancelled=true;}
    require(cancelled&&scorer.cancelled&&scorer.initialized()&&!empty.available&&work==beforeWork,
        "resetting a token revived an already cancelled numerical owner");
    Scorer fresh;std::uint64_t freshWork=0;
    require(fresh.initialize(source,region,4096,64*MiB,[&](std::uint64_t n){freshWork+=n;},inside,
        0,500'000,&stop,true),"fresh export owner did not initialize");
    const auto freshBefore=freshWork;Carrier occupied;occupied.available=true;occupied.sourceNodes=91;
    require(!fresh.releaseNumericalTopology(occupied,region)&&occupied.available&&occupied.sourceNodes==91&&
        fresh.initialized(),"export overwrote an occupied carrier");
    require(fresh.releaseNumericalTopology(empty,region)&&!fresh.initialized()&&!empty.complex.chargeWork&&
        freshWork==freshBefore,"ownership transfer rebuilt geometry or retained the destroyed owner's callback");
}

void budgetLimitedPlacementHandoff(){
    struct Field {double eval(Vec3)const{return .45;}std::uint64_t queryWork()const{return 1;}} field;
    const auto source=cube(4);const auto sourceHash=mesh::meshHash(source);
    OcctInitialPointFrontLimits limits;limits.sites=2;limits.cones=512;limits.frontStates=2048;
    limits.candidates=20'000;limits.queries=20'000;limits.work=16'000'000;
    limits.workingBytes=256*MiB;limits.progressiveResources=true;
    Carrier topology;OcctInitialPointFrontReport report;std::uint64_t charged=0;
    const auto placed=appendOcctInitialPointFrontSites(source,field,.45,region,limits,inside,
        [&](std::uint64_t n){require(n<=limits.work-charged,"front escaped its work allowance");charged+=n;},&report,&topology);
    require(report.initialized&&report.sites==2&&report.budgetLimited&&!report.frontClosed&&topology.available&&
        report.workUsed==charged&&mesh::meshHash(source)==sourceHash,
        "site-limited placement did not preserve a reusable valid numerical topology");
    require(placed.nodeCount()==source.nodeCount()+2&&topology.complex.vertices.size()==std::size_t(placed.nodeCount())+8&&
        !topology.complex.chargeWork,"partial placement exported temporary cones or a dangling owner callback");
    coherent(topology.complex);
    mesh::Mesh seed;const auto recovered=recover(placed,topology,seed,50'000'000,256*MiB);
    require(recovered.initialTopologyUsed&&std::string_view(recovered.diagnostic)=="complete_seed"&&
        recovered.initialTopologyPeakModeledBytes<=256*MiB,"bounded placement topology could not complete constrained recovery");
    checkSeed(placed,seed);
}
}

int main(){try{
    std::cout<<"orderedImportAndLifetime\n"<<std::flush;orderedImportAndLifetime();
    std::cout<<"completeRecoveryAndGuards\n"<<std::flush;completeRecoveryAndGuards();
    std::cout<<"invalidCarrierRefusals\n"<<std::flush;invalidCarrierRefusals();
    std::cout<<"resourceAndCancellationRefusals\n"<<std::flush;resourceAndCancellationRefusals();
    std::cout<<"exportRefusals\n"<<std::flush;exportRefusals();
    std::cout<<"budgetLimitedPlacementHandoff\n"<<std::flush;budgetLimitedPlacementHandoff();
    std::cout<<"initial topology identity, incidence, recovery and resource guards passed\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
catch(const bnd::detail::constraint_detail::WorkLimit&){std::cerr<<"uncaught recovery WorkLimit\n";return 1;}
catch(const bnd::detail::constraint_detail::Cancelled&){std::cerr<<"uncaught recovery Cancelled\n";return 1;}
catch(const bnd::detail::constraint_detail::CellLimit&){std::cerr<<"uncaught recovery CellLimit\n";return 1;}
catch(const OcctInitialConstructionCancelled&){std::cerr<<"uncaught initial construction cancellation\n";return 1;}
catch(...){std::cerr<<"unexpected initial-topology exception\n";return 1;}}
