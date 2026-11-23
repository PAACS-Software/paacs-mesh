// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/bnd/TetConstraintRecovery.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <new>
#include <stdexcept>
#include <string_view>

namespace {
std::atomic<std::size_t> allocations{0},cancelAt{0},failAt{0};
std::atomic<std::atomic_bool*> cancellationTarget{nullptr};
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
}

void* operator new(std::size_t bytes){
    const auto ordinal=++allocations;
    if(ordinal==cancelAt.load())if(auto* target=cancellationTarget.load())target->store(true);
    if(ordinal==failAt.load())throw std::bad_alloc{};
    if(auto* result=std::malloc(bytes?bytes:1))return result;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t bytes){return ::operator new(bytes);}
void operator delete(void* value)noexcept{std::free(value);}
void operator delete[](void* value)noexcept{std::free(value);}
void operator delete(void* value,std::size_t)noexcept{std::free(value);}
void operator delete[](void* value,std::size_t)noexcept{std::free(value);}

using namespace paacs;
using bnd::detail::recoverTetConstraints;
namespace {
constexpr auto mib=1024ULL*1024;
constexpr mesh::Classification region{3,7};

double radicalInverse(unsigned ordinal,unsigned base){
    double value=0.,factor=1.;
    while(ordinal){factor/=base;value+=factor*double(ordinal%base);ordinal/=base;}
    return value;
}

mesh::Mesh cube(unsigned divisions=1,unsigned interior=0){
    mesh::Mesh source;
    std::map<std::array<unsigned,3>,Index> ids;
    const auto node=[&](std::array<unsigned,3> key){
        if(const auto found=ids.find(key);found!=ids.end())return found->second;
        const auto id=source.addNode({double(key[0])/divisions,double(key[1])/divisions,
                                      double(key[2])/divisions},{2,int(ids.size()+1)});
        ids.emplace(key,id);return id;
    };
    for(unsigned axis=0;axis<3;++axis)for(unsigned side=0;side<2;++side)
        for(unsigned u=0;u<divisions;++u)for(unsigned v=0;v<divisions;++v){
            std::array<Index,4> vertices{};
            const std::array du{0U,1U,1U,0U},dv{0U,0U,1U,1U};
            for(unsigned j=0;j<4;++j){
                std::array<unsigned,3> key{};key[axis]=side*divisions;
                key[(axis+1)%3]=u+du[j];key[(axis+2)%3]=v+dv[j];vertices[j]=node(key);
            }
            std::array<Index,3> a{vertices[0],vertices[1],vertices[2]},b{vertices[0],vertices[2],vertices[3]};
            if(!side){std::swap(a[1],a[2]);std::swap(b[1],b[2]);}
            source.addCell(mesh::CellType::Tri3,a,{2,int(axis*2+side+1)});
            source.addCell(mesh::CellType::Tri3,b,{2,int(axis*2+side+1)});
        }
    for(unsigned i=1;i<=interior;++i)
        source.addNode({.05+.9*radicalInverse(i,2),.05+.9*radicalInverse(i,3),
                        .05+.9*radicalInverse(i,5)},region);
    return source;
}

mesh::Mesh tetrahedron(unsigned interior=0){
    mesh::Mesh source;
    const std::array<Vec3,4> corners{{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}};
    for(auto point:corners)source.addNode(point,{0,int(source.nodeCount()+1)});
    bnd::TetFacetCell cell{{0,1,2,3},1,region};
    if(predicates::orient3d(corners[0],corners[1],corners[2],corners[3])!=predicates::Sign::Positive)
        std::swap(cell.vertices[0],cell.vertices[1]);
    for(unsigned side=0;side<4;++side){
        const auto face=bnd::tet_facet_detail::outwardFace(cell,side);
        const std::array<Index,3> ids{Index(face[0]),Index(face[1]),Index(face[2])};
        source.addCell(mesh::CellType::Tri3,ids,{2,int(side+1)});
    }
    for(unsigned i=1;i<=interior;++i){
        const double a=.2+radicalInverse(i,2),b=.2+radicalInverse(i,3),
            c=.2+radicalInverse(i,5),d=.2+radicalInverse(i,7),sum=a+b+c+d;
        source.addNode({b/sum,c/sum,d/sum},region);
    }
    return source;
}

mesh::Mesh heldDestination(){mesh::Mesh held;held.addNode({11,12,13},{0,91});return held;}

void sourceCountAndRealBudgetStops(){
    const auto source=cube(1,12001);const auto before=mesh::meshHash(source);
    auto held=heldDestination();const auto original=mesh::meshHash(held);
    const auto legacy=recoverTetConstraints(source,region,2'000'000,128*mib,100'000,nullptr,held);
    require(std::string_view(legacy.diagnostic)=="admission_budget"&&legacy.work==0,
            "historical source-count admission changed");
    const auto normal=recoverTetConstraints(source,region,100'000,128*mib,100'000,nullptr,held,
                                            0,false,true,false,true);
    require((std::string_view(normal.diagnostic)=="construction_budget"||
             std::string_view(normal.diagnostic)=="work_budget")&&normal.work>0&&normal.work<=100'000,
            "large Normal input hit a fixed source-count cutoff instead of real charged work");
    require(mesh::meshHash(held)==original&&mesh::meshHash(source)==before,
            "resource refusal changed the source or partially replaced the held destination");
    const auto memory=recoverTetConstraints(source,region,2'000'000,8*mib,100'000,nullptr,held,
                                            0,false,true,false,true);
    require(std::string_view(memory.diagnostic)=="memory_budget"&&mesh::meshHash(held)==original,
            "large Normal input ignored modeled-memory admission");
    const auto cells=recoverTetConstraints(source,region,2'000'000,128*mib,0,nullptr,held,
                                           0,false,true,false,true);
    require(std::string_view(cells.diagnostic)=="cell_budget"&&mesh::meshHash(held)==original,
            "zero cell allowance was not a distinct resource refusal");
    const auto facets=cube(46);
    require(facets.cellCount()>24000,"fixture does not exceed the old source-face ceiling");
    const auto faceStop=recoverTetConstraints(facets,region,1,128*mib,100'000,nullptr,held,
                                              0,false,true,false,true);
    require(std::string_view(faceStop.diagnostic)=="work_budget"&&faceStop.work<=1&&
            mesh::meshHash(held)==original,"large face source bypassed its explicit work cap");
    std::cout<<"large source admitted to bounded construction: nodes="<<source.nodeCount()
             <<" diagnostic="<<normal.diagnostic<<" work="<<normal.work<<'\n';
}

void invalidSourceAndCancellation(){
    auto source=tetrahedron();auto held=heldDestination();const auto original=mesh::meshHash(held);
    const auto duplicate=source.cellNodes(0);
    source.addCell(mesh::CellType::Tri3,std::array<Index,3>{duplicate[0],duplicate[1],duplicate[2]},{2,99});
    const auto invalid=recoverTetConstraints(source,region,2'000'000,128*mib,100'000,nullptr,held,
                                             0,false,true,false,true);
    require(std::string_view(invalid.diagnostic)=="invalid_source"&&mesh::meshHash(held)==original,
            "invalid source constraints were treated as an admissible Normal seed");
    source=tetrahedron();source.setNode(1,{std::numeric_limits<double>::quiet_NaN(),0,0});
    const auto nonfinite=recoverTetConstraints(source,region,2'000'000,128*mib,100'000,nullptr,held,
                                               0,false,true,false,true);
    require(std::string_view(nonfinite.diagnostic)=="invalid_source"&&mesh::meshHash(held)==original,
            "nonfinite source escaped Normal geometry admission");
    source=tetrahedron();std::atomic_bool stop{true};
    const auto cancelled=recoverTetConstraints(source,region,2'000'000,128*mib,100'000,&stop,held,
                                               0,false,true,false,true);
    require(std::string_view(cancelled.diagnostic)=="cancelled"&&cancelled.work==0&&
            mesh::meshHash(held)==original,"cancelled construction published a seed");
    const auto aliasHash=mesh::meshHash(source);
    const auto alias=recoverTetConstraints(source,region,2'000'000,128*mib,100'000,nullptr,source,
                                           0,false,true,false,true);
    require(std::string_view(alias.diagnostic)=="invalid_source"&&mesh::meshHash(source)==aliasHash,
            "source/destination aliasing changed input geometry");
}

void completeSeedAndCellLimit(){
    const auto source=tetrahedron();mesh::Mesh output;
    const auto result=recoverTetConstraints(source,region,2'000'000,128*mib,100'000,nullptr,output,
                                            0,false,true,false,true);
    const auto quality=algo::analyzeMixedMeshQuality(output);
    require(std::string_view(result.diagnostic)=="complete_seed"&&result.seedCells>0&&quality.valid&&
            std::abs(quality.totalVolume-1./6)<1e-12,"Normal changed valid constraint construction");
    auto held=heldDestination();const auto original=mesh::meshHash(held);
    const auto cubeSource=cube();
    const auto capped=recoverTetConstraints(cubeSource,region,5'000'000,128*mib,1,nullptr,held,
                                            0,false,true,false,true);
    require(std::string_view(capped.diagnostic)=="cell_budget"&&mesh::meshHash(held)==original,
            "the final private seed exceeded its cell ceiling or replaced output partially");
    failAt.store(allocations.load()+1);
    const auto failed=recoverTetConstraints(source,region,2'000'000,128*mib,100'000,nullptr,held,
                                            0,false,true,false,true);
    failAt.store(0);
    require(std::string_view(failed.diagnostic)=="allocation_failed"&&mesh::meshHash(held)==original,
            "allocation failure escaped or published a partial seed");
}

void cancelledBuilderRetainsWork(){
    const auto source=cube(1,12001);auto held=heldDestination();const auto original=mesh::meshHash(held);
    std::atomic_bool stop{false};
    cancellationTarget.store(&stop);cancelAt.store(allocations.load()+200);
    const auto cancelled=recoverTetConstraints(source,region,100'000'000,512*mib,100'000,&stop,held,
                                               0,false,true,false,true);
    cancelAt.store(0);cancellationTarget.store(nullptr);
    require(stop.load()&&std::string_view(cancelled.diagnostic)=="cancelled"&&
            cancelled.candidateCells==0&&
            cancelled.work>std::uint64_t(source.nodeCount()+4*source.cellCount())&&
            cancelled.work<=100'000'000&&mesh::meshHash(held)==original,
            "cancelled bounded builder forgot completed work or partially published output");
}

void largeCompleteSeed(){
    auto source=cube();
    // Keep this resource-admission test independent of difficult nonlocal
    // facet recovery: central points leave empty supporting spheres between
    // each prescribed cube face and the eight exterior guard corners.
    for(unsigned i=1;i<=11997;++i)
        source.addNode({.25+.5*radicalInverse(i,2),.25+.5*radicalInverse(i,3),
                        .25+.5*radicalInverse(i,5)},region);
    const auto original=mesh::meshHash(source);mesh::Mesh output;
    const auto result=recoverTetConstraints(source,region,500'000'000,1024*mib,200'000,nullptr,output,
                                            0,false,true,false,true);
    const auto quality=algo::analyzeMixedMeshQuality(output);
    std::cout<<"large complete seed: nodes="<<source.nodeCount()<<" cells="<<output.cellCount()
             <<" diagnostic="<<result.diagnostic<<" work="<<result.work<<" seconds="<<result.seconds<<'\n';
    require(std::string_view(result.diagnostic)=="complete_seed"&&result.seedCells>0&&quality.valid&&
            std::abs(quality.totalVolume-1.)<1e-10&&result.work<=500'000'000&&
            mesh::meshHash(source)==original,"large Normal source failed real constrained construction");
}
}

int main(int argc,char** argv){try{
    sourceCountAndRealBudgetStops();invalidSourceAndCancellation();completeSeedAndCellLimit();
    cancelledBuilderRetainsWork();
    if(argc==2&&std::string_view(argv[1])=="--large")largeCompleteSeed();
    else require(argc==1,"only optional --large is accepted");
    std::cout<<"constraint resource admission, refusal and cancellation checks passed\n";return 0;
}catch(const std::exception& error){failAt.store(0);cancelAt.store(0);cancellationTarget.store(nullptr);
    std::cerr<<error.what()<<'\n';return 1;}}
