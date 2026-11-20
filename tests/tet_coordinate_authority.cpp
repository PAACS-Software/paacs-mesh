// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/TetCandidate.hpp"
#include "paacs/predicates/Predicates.hpp"
#include <atomic>
#include <bit>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <stdexcept>

static std::atomic<std::size_t> allocationCalls{0},failOrdinal{0};
void* operator new(std::size_t n){const auto ordinal=++allocationCalls;if(ordinal==failOrdinal.load())throw std::bad_alloc{};if(auto* p=std::malloc(n?n:1))return p;throw std::bad_alloc{};}
void* operator new[](std::size_t n){return ::operator new(n);}
void operator delete(void* p)noexcept{std::free(p);}
void operator delete[](void* p)noexcept{std::free(p);}
void operator delete(void* p,std::size_t)noexcept{std::free(p);}
void operator delete[](void* p,std::size_t)noexcept{std::free(p);}
using namespace paacs;
void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
int main(){try{
    // Six ordinary finite binary doubles form a nearly coplanar curved-wall
    // patch. Rounded centering/scaling previously emitted a physically
    // negative Tet and a same-winding interior face (AdjacencyMismatch).
    const std::vector<Vec3> points{
        {10.856053358834119,16.79720528749192,10.0},
        {15.937782728263846,10.343455984662688,10.0},
        {15.984496655211151,11.002740097382754,10.91419796248158},
        {15.76057402256841,10.848605589300725,10.991058967623033},
        {15.57124299398417,11.580190611268433,10.91419796248158},
        {15.353109523788468,11.417966749978774,10.991058967623033}};
    const auto before=points;
    algo::Delaunay3DStats stats;const auto raw=algo::delaunay3D(points,&stats);
    require(raw.size()==4&&stats.originalCoordinateRetryUsed&&stats.discardedConstructionWork>0,
            "bounded exact-coordinate retry must replace rejected physical output and retain discarded work");
    for(const auto& cell:raw)
        require(predicates::orient3d(points[cell.v[0]],points[cell.v[1]],points[cell.v[2]],points[cell.v[3]])==predicates::Sign::Positive,
                "public Delaunay returned a nonpositive original-coordinate tetrahedron");
    const auto candidate=algo::buildTetCandidate(points);
    require(candidate.built()&&algo::validateTetCandidate(*candidate.candidate).valid&&
            candidate.metrics.originalCoordinateRetryUsed&&candidate.metrics.construction.discardedConstructionWork>0,
            "six-point candidate must pass original-space topology and exact hull audits");
    for(std::size_t i=0;i<points.size();++i){
        require(std::bit_cast<std::uint64_t>(points[i].x)==std::bit_cast<std::uint64_t>(before[i].x)&&
                std::bit_cast<std::uint64_t>(points[i].y)==std::bit_cast<std::uint64_t>(before[i].y)&&
                std::bit_cast<std::uint64_t>(points[i].z)==std::bit_cast<std::uint64_t>(before[i].z),
                "candidate retry changed caller-owned coordinates");
    }
    const auto replay=algo::buildTetCandidate(points);require(replay.built()&&
        replay.candidate->cells.size()==candidate.candidate->cells.size(),"retry replay counts changed");
    for(std::size_t i=0;i<replay.candidate->cells.size();++i)
        require(replay.candidate->cells[i].vertices==candidate.candidate->cells[i].vertices&&
                replay.candidate->cells[i].neighbours==candidate.candidate->cells[i].neighbours,
                "exact-coordinate retry replay changed connectivity");
    const std::vector<Vec3> ordinary{{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1},{0,0,0},{.125,.25,-.375},{-.375,.5,.625},{.75,-.625,.375}};
    constexpr std::array<Index,4> golden[]={
{{0,3,9,2}},
{{0,8,9,3}},
{{0,1,9,8}},
{{0,2,9,1}},
{{0,8,11,1}},
{{0,4,11,8}},
{{1,8,11,9}},
{{1,9,11,2}},
{{0,4,5,11}},
{{0,11,5,1}},
{{11,5,10,4}},
{{8,11,10,4}},
{{2,9,10,3}},
{{8,3,10,9}},
{{0,8,10,4}},
{{0,3,10,8}},
{{0,3,7,10}},
{{0,10,7,4}},
{{4,10,6,5}},
{{4,7,6,10}},
{{5,10,6,11}},
{{1,11,6,2}},
{{1,5,6,11}},
{{2,9,6,10}},
{{2,11,6,9}},
{{11,10,6,8}},
{{9,8,6,10}},
{{9,11,6,8}},
{{10,7,6,3}},
{{2,10,6,3}},
};
    algo::Delaunay3DStats ordinaryStats;const auto ordinaryRaw=algo::delaunay3D(ordinary,&ordinaryStats);
    require(ordinaryRaw.size()==std::size(golden)&&!ordinaryStats.originalCoordinateRetryUsed&&
            ordinaryStats.discardedConstructionWork==0,"already-valid output unnecessarily retried");
    for(std::size_t i=0;i<ordinaryRaw.size();++i)
        require(ordinaryRaw[i].v==golden[i],"already-valid baseline tetrahedron order or winding changed");
    std::size_t failures=0;bool failedDuringRetry=false,finished=false;
    for(std::size_t ordinal=1;ordinal<=1024;++ordinal){
        allocationCalls=0;failOrdinal=ordinal;
        const auto attempted=algo::buildTetCandidate(points);
        const auto calls=allocationCalls.load();failOrdinal=0;
        if(attempted.built()){
            require(ordinal>calls,"allocation failure was swallowed and partially accepted");finished=true;break;
        }
        ++failures;require(attempted.reason==algo::TetCandidateReason::AllocationFailed&&!attempted.candidate,
                           "allocation replay lost its distinct failure reason or exposed a partial candidate");
        failedDuringRetry|=attempted.metrics.construction.originalCoordinateRetryUsed&&
                          attempted.metrics.construction.discardedConstructionWork>0;
    }
    require(finished&&failures>0&&failedDuringRetry,"allocation replay did not reach exact-retry failures and completion");
    const std::vector<Vec3> nonfinite{{0,0,0},{1,0,0},{0,1,0},{0,0,(std::numeric_limits<double>::quiet_NaN)()}};
    const auto rejected=algo::buildTetCandidate(nonfinite);
    require(rejected.reason==algo::TetCandidateReason::NonfiniteInput&&!rejected.metrics.originalCoordinateRetryUsed,
            "invalid input must retain its reason without an exact-coordinate retry");
    std::puts("Original-coordinate positivity, six-point topology, baseline bytes, retry work and allocation checks passed");return 0;
}catch(const std::exception& e){failOrdinal=0;std::fprintf(stderr,"%s\n",e.what());return 1;}}
