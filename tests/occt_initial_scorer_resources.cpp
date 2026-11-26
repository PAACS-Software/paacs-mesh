// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctInitialDelaunayScorer.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <vector>

using namespace paacs;
using namespace paacs::geo::detail;
namespace {
using Scorer=OcctInitialDelaunayScorer;
constexpr std::uint64_t MiB=1024ULL*1024;
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
struct WorkStop {};
struct Work {
    std::uint64_t used=0,limit=100'000'000;
    bool refuse=false,stopOnCharge=false;
    std::atomic_bool cancellation=false;
    void operator()(std::uint64_t n){
        if(refuse)throw WorkStop{};
        require(used<=limit&&n<=limit-used,"scorer charged beyond the enclosing work allowance");
        used+=n;
        if(stopOnCharge)cancellation=true;
    }
};
bool inside(Vec3 p){return std::abs(p.x)<1&&std::abs(p.y)<1&&std::abs(p.z)<1;}
mesh::Mesh cube(){
    mesh::Mesh source;
    for(unsigned i=0;i<8;++i)
        source.addNode({(i&1)?1.:-1.,(i&2)?1.:-1.,(i&4)?1.:-1.},{0,int(i)+1});
    for(unsigned axis=0;axis<3;++axis)for(unsigned side=0;side<2;++side){
        const auto vertex=[&](unsigned u,unsigned v){
            return Index((side<<axis)|(u<<((axis+1)%3))|(v<<((axis+2)%3)));
        };
        std::array<Index,3> a{vertex(0,0),vertex(1,0),vertex(1,1)};
        std::array<Index,3> b{vertex(0,0),vertex(1,1),vertex(0,1)};
        if(!side){std::swap(a[1],a[2]);std::swap(b[1],b[2]);}
        source.addCell(mesh::CellType::Tri3,a,{2,int(2*axis+side+1)});
        source.addCell(mesh::CellType::Tri3,b,{2,int(2*axis+side+1)});
    }
    return source;
}
Scorer::Front makeFront(const mesh::Mesh& source,Work& work){
    Scorer::Front front;
    for(const auto p:source.nodes())front.vertices.push_back({p,front.vertices.size()+1});
    for(Index i=0;i<source.cellCount();++i){const auto* v=source.cellNodes(i);
        front.boundary.push_back({Scorer::Vertex(v[0]),Scorer::Vertex(v[1]),Scorer::Vertex(v[2])});}
    require(front.initialize(),"owned cube is not a closed oriented front");
    front.chargeWork=[&work](std::uint64_t n){work(n);};
    return front;
}
std::vector<std::uint64_t> geometry(const Scorer& scorer){
    const auto& c=scorer.complexForDiagnostics();std::vector<std::uint64_t> result;
    const auto add=[&](auto x){result.push_back(std::uint64_t(x));};
    add(c.vertices.size());for(const auto& v:c.vertices){
        add(std::bit_cast<std::uint64_t>(v.position.x));add(std::bit_cast<std::uint64_t>(v.position.y));
        add(std::bit_cast<std::uint64_t>(v.position.z));add(v.sequence);
    }
    add(c.cells.size());for(const auto& cell:c.cells){
        for(auto v:cell.vertices)add(v);add(cell.sequence);add(cell.classification.dim);add(cell.classification.tag);
    }
    for(auto v:c.alive)add(v);
    for(const auto& row:c.neighbours)for(auto v:row)add(v);
    for(const auto& [key,owners]:c.faces){for(auto v:key)add(v);add(owners.size());for(auto v:owners)add(v);}
    for(const auto& [key,owners]:c.edges){for(auto v:key)add(v);add(owners.size());for(auto v:owners)add(v);}
    for(auto v:scorer.auxiliaryToFront())add(v);
    for(auto v:scorer.relevanceForDiagnostics())add(v);
    return result;
}
using Quartet=std::array<Scorer::Vertex,4>;
std::vector<Quartet> liveQuartets(const Scorer& scorer){
    std::vector<Quartet> result;const auto& c=scorer.complexForDiagnostics();
    for(std::size_t i=0;i<c.cells.size();++i)if(c.alive[i]){
        const auto& cell=c.cells[i];const auto& v=cell.vertices;
        require(cell.classification==mesh::Classification{3,7},"resource operation changed body ownership");
        require(predicates::orient3dFiltered(c.vertices[v[0]].position,c.vertices[v[1]].position,
            c.vertices[v[2]].position,c.vertices[v[3]].position)==predicates::Sign::Positive,
            "resource operation left a nonpositive numerical tetrahedron");
        auto q=v;std::sort(q.begin(),q.end());result.push_back(q);
    }
    std::sort(result.begin(),result.end());return result;
}
bool initialize(Scorer& scorer,const mesh::Mesh& source,Work& work,std::size_t ceiling,
                std::uint64_t bytes=64*MiB,bool progressive=true,std::uint64_t bootstrap=500'000){
    return scorer.initialize(source,{3,7},ceiling,bytes,[&](std::uint64_t n){work(n);},
        inside,0,bootstrap,&work.cancellation,progressive);
}
double radical(unsigned n,unsigned base){
    double value=0.,weight=1.;while(n){weight/=base;value+=weight*(n%base);n/=base;}return value;
}
Vec3 candidate(unsigned n){return {1.6*radical(n,2)-.8,1.6*radical(n,3)-.8,1.6*radical(n,5)-.8};}
// Stop at a real scored patch requiring more admitted storage. The selected
// point is retained by the caller for a retry after resource preparation.
Vec3 reachStoragePressure(Scorer& scorer,Scorer::Front& front,unsigned& next,bool memoryStop=false){
    for(unsigned attempts=0;attempts<500;++attempts){
        const auto p=candidate(next++);require(inside(p),"trial point left owned material");
        const auto id=Scorer::Vertex(front.vertices.size());front.vertices.push_back({p,front.vertices.size()+1});
        const auto plan=scorer.score(p,front);
        if(!plan||!plan->valid){front.vertices.pop_back();
            if(scorer.resourceRetryRequested())return p;
            if(memoryStop&&scorer.budgetLimited&&
                std::string_view(scorer.resourceDiagnostic)=="auxiliary_memory_limit")return p;
            require(!scorer.budgetLimited,"resource trial hit an unrelated terminal allowance");
            continue;
        }
        if(!scorer.commit(*plan,id)){front.vertices.pop_back();continue;}
    }
    throw std::runtime_error("point sequence did not exercise actual capacity pressure");
}
void retryPoint(Scorer& scorer,Scorer::Front& front,Vec3 p){
    const auto id=Scorer::Vertex(front.vertices.size());front.vertices.push_back({p,front.vertices.size()+1});
    const auto plan=scorer.score(p,front);
    require(plan&&plan->valid&&scorer.commit(*plan,id),"admitted growth did not retry the same exact-positive point");
}
void exerciseActualBootstrap(){
    const auto source=cube();Work smallWork,largeWork,legacyWork;Scorer small,large,legacy;
    require(initialize(small,source,smallWork,4096,1024*MiB),"small-ceiling progressive bootstrap failed");
    require(initialize(large,source,largeWork,1'000'000,1024*MiB),"unused future ceiling blocked progressive bootstrap");
    require(!initialize(legacy,source,legacyWork,1'000'000,1024*MiB,false),
        "fixture did not expose historical future-capacity import reservation");
    require(legacy.budgetLimited&&!legacy.bootstrapWork,"historical refusal did not occur before construction");
    require(geometry(small)==geometry(large)&&smallWork.used==largeWork.used,
        "unused future capacity changed bootstrap geometry or charged work");
    require(large.initialStorageCapacity==small.initialStorageCapacity&&
        large.initialStorageCapacity<large.hardStorageLimit()/8,"bootstrap reserved the future ceiling");
    require(large.bootstrapImportedCells==large.storedCells()&&
        large.bootstrapImportWork==source.nodeCount()+8+12*large.storedCells(),
        "bootstrap import accounting did not use actual cells");
    require(largeWork.used==source.nodeCount()+source.cellCount()+large.bootstrapWork+large.bootstrapImportWork,
        "bootstrap accounting omitted or repeated work");
    require(large.maximumResourceModeledBytes<=large.workingBytes()&&
        large.bootstrapImportModeledBytes<=large.workingBytes(),"bootstrap exceeded admitted peak memory");
    const auto before=geometry(large);const auto epoch=large.planEpoch();const auto capacity=large.storageCapacity();
    const auto currentQuota=large.maximumResourceModeledBytes;
    require(large.setWorkingBudget(currentQuota)&&large.workingBytes()==currentQuota,
        "actual current quota could not replace the unused future budget");
    require(!large.setWorkingBudget(1)&&large.workingBytes()==currentQuota&&
        geometry(large)==before&&large.planEpoch()==epoch&&large.storageCapacity()==capacity,
        "refused working-budget shrink changed quota or geometry");
    std::cout<<"actual bootstrap: cells="<<large.bootstrapImportedCells<<" capacity="<<capacity
        <<" work="<<largeWork.used<<" modeled_peak="<<large.maximumResourceModeledBytes<<'\n';
}
void exerciseGrowthAndCompaction(){
    const auto source=cube();Work work;Scorer scorer;
    require(initialize(scorer,source,work,4096),"growth bootstrap failed");
    scorer.compactionEnabled=false;auto front=makeFront(source,work);unsigned next=1;
    auto pending=reachStoragePressure(scorer,front,next);
    const auto before=geometry(scorer);const auto capacity=scorer.storageCapacity();const auto epoch=scorer.planEpoch();
    const auto bootstrapWork=scorer.bootstrapWork,imported=scorer.bootstrapImportedCells;
    work.refuse=true;bool refused=false;
    try{(void)scorer.prepareForBatch();}catch(const WorkStop&){refused=true;}work.refuse=false;
    require(refused&&geometry(scorer)==before&&scorer.storageCapacity()==capacity&&scorer.planEpoch()==epoch,
        "work refusal partially changed live numerical storage");
    require(scorer.prepareForBatch()&&scorer.capacityGrowths==1&&scorer.storageCapacity()>capacity,
        "real storage pressure did not admit progressive growth");
    require(geometry(scorer)==before&&scorer.planEpoch()==epoch&&scorer.growthWork>0,
        "capacity-only growth changed cell IDs or mappings");
    require(scorer.bootstrapWork==bootstrapWork&&scorer.bootstrapImportedCells==imported,
        "capacity growth rebuilt the Delaunay seed");
    retryPoint(scorer,front,pending);
    for(unsigned pass=0;pass<2;++pass){
        pending=reachStoragePressure(scorer,front,next);
        const auto geometryBefore=geometry(scorer);const auto epochBefore=scorer.planEpoch();
        require(scorer.prepareForBatch()&&geometry(scorer)==geometryBefore&&scorer.planEpoch()==epochBefore,
            "later capacity admission changed exact geometry");
        retryPoint(scorer,front,pending);
    }
    const auto live=liveQuartets(scorer);const auto mapping=scorer.auxiliaryToFront();
    const auto oldEpoch=scorer.planEpoch();const auto oldCompactions=scorer.compactions;
    const auto oldSites=scorer.acceptedSites();
    scorer.compactionEnabled=true;
    require(scorer.prepareForBatch(true)&&scorer.compactions==oldCompactions+1&&
        scorer.planEpoch()==oldEpoch+1,"compaction did not invalidate old cell identities");
    require(liveQuartets(scorer)==live&&scorer.auxiliaryToFront()==mapping&&scorer.acceptedSites()==oldSites&&
        scorer.storedCells()==scorer.liveCells(),"compaction changed the live material point/cell set");
    require(scorer.maximumResourceModeledBytes<=scorer.workingBytes()&&
        scorer.maximumGrowthModeledBytes<=scorer.workingBytes(),"growth or compaction exceeded modeled memory");
    std::cout<<"progressive growth: sites="<<scorer.acceptedSites()<<" growths="<<scorer.capacityGrowths
        <<" copied_work="<<scorer.growthWork<<" reclaimed="<<scorer.deadCellsReclaimed<<'\n';
}
void exerciseRefusalsAndCancellation(){
    const auto source=cube();Work probeWork;Scorer probe;
    require(initialize(probe,source,probeWork,4096),"refusal probe bootstrap failed");
    Work hardWork;Scorer hard;
    require(initialize(hard,source,hardWork,probe.storedCells()),"exact seed ceiling was refused");
    auto hardFront=makeFront(source,hardWork);const auto before=geometry(hard);const auto epoch=hard.planEpoch();
    hardFront.vertices.push_back({{.13,.21,.37},hardFront.vertices.size()+1});
    require(!hard.score({.13,.21,.37},hardFront)&&hard.budgetLimited&&!hard.resourceRetryRequested()&&
        std::string_view(hard.resourceDiagnostic)=="auxiliary_cell_limit",
        "hard ceiling was treated as a retryable geometry failure");
    require(!hard.prepareForBatch()&&geometry(hard)==before&&hard.planEpoch()==epoch,
        "hard-cap refusal changed a certified seed");
    require(hard.setWorkingBudget(128*MiB)&&hard.budgetLimited&&!hard.resourceRetryRequested()&&
        std::string_view(hard.resourceDiagnostic)=="auxiliary_cell_limit",
        "increasing memory overrode the caller's hard cell ceiling");
    Work memoryWork;Scorer memory;
    require(!initialize(memory,source,memoryWork,4096,8*MiB)&&memory.budgetLimited&&
        !memory.initialized()&&memory.storedCells()==0,"insufficient fixed memory admitted a partial seed");
    Work tinyWork;Scorer tiny;
    require(!initialize(tiny,source,tinyWork,4096,64*MiB,true,1)&&tiny.budgetLimited&&tinyWork.used<=1&&
        tiny.storedCells()==0,"tiny explicit work limit admitted or overcharged a seed");
    Work cancelledWork;Scorer cancelled;cancelledWork.cancellation=true;bool stopped=false;
    try{(void)initialize(cancelled,source,cancelledWork,4096);}catch(const OcctInitialConstructionCancelled&){stopped=true;}
    require(stopped&&cancelledWork.used==0&&cancelled.storedCells()==0,"pre-cancelled bootstrap performed work");
    Work growthWork;Scorer growth;
    require(initialize(growth,source,growthWork,4096),"cancelled-growth bootstrap failed");
    growth.compactionEnabled=false;auto front=makeFront(source,growthWork);unsigned next=1;
    (void)reachStoragePressure(growth,front,next);
    const auto old=geometry(growth);const auto oldCapacity=growth.storageCapacity();const auto oldEpoch=growth.planEpoch();
    growthWork.stopOnCharge=true;stopped=false;
    try{(void)growth.prepareForBatch();}catch(const OcctInitialConstructionCancelled&){stopped=true;}
    require(stopped&&geometry(growth)==old&&growth.storageCapacity()==oldCapacity&&growth.planEpoch()==oldEpoch,
        "cancellation installed staged growth arrays");
}
void exerciseMemoryShareResume(){
    const auto source=cube();Work work;Scorer scorer;
    require(initialize(scorer,source,work,4096),"memory-share bootstrap failed");
    scorer.compactionEnabled=false;auto front=makeFront(source,work);unsigned next=1;
    require(scorer.setWorkingBudget(scorer.currentModeledBytes()),"retained scratch envelope could not be admitted");
    const auto pending=reachStoragePressure(scorer,front,next,true);
    require(scorer.budgetLimited&&!scorer.resourceRetryRequested()&&
        std::string_view(scorer.resourceDiagnostic)=="auxiliary_memory_limit",
        "tight shared memory did not produce its distinct resource stop");
    const auto before=geometry(scorer);const auto epoch=scorer.planEpoch();
    const auto bootstrap=scorer.bootstrapWork;
    require(scorer.setWorkingBudget(64*MiB)&&!scorer.budgetLimited&&scorer.growthRequested()&&
        scorer.resourceRetryRequested()&&geometry(scorer)==before&&scorer.planEpoch()==epoch,
        "restoring shared memory did not resume the pending resource request");
    require(scorer.prepareForBatch(),"restored memory did not admit staged growth");
    retryPoint(scorer,front,pending);
    require(scorer.bootstrapWork==bootstrap&&scorer.maximumResourceModeledBytes<=scorer.workingBytes(),
        "memory-share recovery rebuilt the seed or exceeded admitted memory");
}
}
int main(){try{
    exerciseActualBootstrap();exerciseGrowthAndCompaction();exerciseRefusalsAndCancellation();exerciseMemoryShareResume();
    std::cout<<"initial scorer progressive resources passed\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
catch(...){std::cerr<<"unexpected scorer resource exception\n";return 1;}}
