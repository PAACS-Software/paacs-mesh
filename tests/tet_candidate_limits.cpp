// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0
#include "paacs/algo/TetCandidate.hpp"

#include <atomic>
#include <bit>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <stdexcept>

namespace {
std::atomic<std::size_t> allocations{0}, cancelAt{0}, failAt{0};
std::atomic_bool cancellation{false};
void require(bool value,const char* message) {
    if (!value) throw std::runtime_error(message);
}
}
void* operator new(std::size_t bytes) {
    const auto ordinal=++allocations;
    if (ordinal==cancelAt.load()) cancellation.store(true);
    if (ordinal==failAt.load()) throw std::bad_alloc{};
    if (auto* result=std::malloc(bytes ? bytes : 1)) return result;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p,std::size_t) noexcept { std::free(p); }
void operator delete[](void* p,std::size_t) noexcept { std::free(p); }

using namespace paacs;
void same(const algo::TetCandidate& a,const algo::TetCandidate& b) {
    require(a.points.size()==b.points.size() &&
            std::memcmp(a.points.data(),b.points.data(),a.points.size()*sizeof(Vec3))==0,
            "bounded construction changed coordinate bytes");
    require(a.inputToPoint==b.inputToPoint && a.pointToInput==b.pointToInput,
            "bounded construction changed identity mappings");
    require(a.cells.size()==b.cells.size() && a.boundaryFaces.size()==b.boundaryFaces.size(),
            "bounded construction changed topology counts");
    for (std::size_t i=0;i<a.cells.size();++i)
        require(a.cells[i].vertices==b.cells[i].vertices &&
                a.cells[i].neighbours==b.cells[i].neighbours &&
                a.cells[i].sourceOrdinal==b.cells[i].sourceOrdinal,
                "bounded construction changed canonical cells");
    for (std::size_t i=0;i<a.boundaryFaces.size();++i)
        require(a.boundaryFaces[i].vertices==b.boundaryFaces[i].vertices &&
                a.boundaryFaces[i].owner==b.boundaryFaces[i].owner &&
                a.boundaryFaces[i].ownerLocalFace==b.boundaryFaces[i].ownerLocalFace,
                "bounded construction changed directed boundary");
}
void stopped(const algo::TetCandidateBuildOutcome& outcome,algo::TetCandidateReason reason) {
    require(!outcome.built() && !outcome.candidate && outcome.reason==reason,
            "bounded stop published a partial candidate or lost its reason");
}
int main() { try {
    const std::vector<Vec3> tetra{{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
    const auto original=tetra;
    algo::TetCandidateBuildLimits limits;
    limits.cancellation=&cancellation;
    limits.maxCells=1; // Four symbolic hull guards do not consume this cap.
    const auto unlimited=algo::buildTetCandidate(tetra);
    const auto bounded=algo::buildTetCandidate(tetra,limits);
    require(unlimited.built() && bounded.built(),"one finite tetrahedron must fit maxCells=1");
    same(*unlimited.candidate,*bounded.candidate);
    require(bounded.metrics.boundedWork>0 && bounded.metrics.boundedPeakTemporaryBytes>0 &&
            unlimited.metrics.boundedWork==0 && unlimited.metrics.boundedPeakTemporaryBytes==0,
            "bounded accounting contaminated the ordinary constructor");

    cancellation.store(true);
    auto before=allocations.load();
    stopped(algo::buildTetCandidate(tetra,limits),algo::TetCandidateReason::Cancelled);
    require(allocations.load()==before,"pre-cancelled construction allocated payload");
    cancellation.store(false);

    auto small=limits;
    small.maxWork=1;
    const auto workStop=algo::buildTetCandidate(tetra,small);
    stopped(workStop,algo::TetCandidateReason::BudgetExceeded);
    require(workStop.metrics.boundedWork<=small.maxWork,"work quota overshot");
    small=limits; small.maxCells=0;
    stopped(algo::buildTetCandidate(tetra,small),algo::TetCandidateReason::BudgetExceeded);
    small=limits; small.maxTemporaryBytes=0;
    before=allocations.load();
    stopped(algo::buildTetCandidate(tetra,small),algo::TetCandidateReason::BudgetExceeded);
    require(allocations.load()==before,"zero payload quota allocated before refusal");

    // A quota just below measured total work reaches canonicalization or its
    // final topology audit; admission must still publish no partial result.
    small=limits; small.maxWork=bounded.metrics.boundedWork-1;
    const auto late=algo::buildTetCandidate(tetra,small);
    stopped(late,algo::TetCandidateReason::BudgetExceeded);
    require(late.metrics.construction.runComplete,"late quota never reached completed raw construction");

    // Cancellation is flipped by a real constructor allocation, without a
    // racing timer or a public callback that could escape a noexcept API.
    cancelAt.store(allocations.load()+4);
    const auto interrupted=algo::buildTetCandidate(tetra,limits);
    cancelAt.store(0); cancellation.store(false);
    stopped(interrupted,algo::TetCandidateReason::Cancelled);
    require(interrupted.metrics.boundedWork>0,"in-flight cancellation was only checked before work");

    failAt.store(allocations.load()+4);
    const auto allocationStop=algo::buildTetCandidate(tetra,limits);
    failAt.store(0);
    stopped(allocationStop,algo::TetCandidateReason::AllocationFailed);
    const auto afterFailure=algo::buildTetCandidate(tetra,limits);
    require(afterFailure.built(),"failed reserve left budget state attached to later construction");
    same(*bounded.candidate,*afterFailure.candidate);

    std::vector<Vec3> plane;
    for (unsigned i=0;i<32;++i) plane.push_back({double(i%8),double(i/8),0});
    small=limits; small.maxCells=1000; small.maxWork=3000;
    const auto coplanarStop=algo::buildTetCandidate(plane,small);
    stopped(coplanarStop,algo::TetCandidateReason::BudgetExceeded);
    require(!coplanarStop.metrics.construction.runComplete &&
            coplanarStop.metrics.construction.bootstrapOrientationTests>0,
            "coplanar bootstrap quota did not stop its nested search");
    auto planeLimits=limits; planeLimits.maxCells=1000;
    const auto planeFull=algo::buildTetCandidate(plane,planeLimits);
    stopped(planeFull,algo::TetCandidateReason::NonVolumetricInput);
    require(planeFull.metrics.originalCoordinateRetryUsed &&
            planeFull.metrics.boundedWork%2==0,
            "coplanar retry must perform two identical complete searches");
    planeLimits.maxWork=planeFull.metrics.boundedWork/2;
    const auto planeRetryStop=algo::buildTetCandidate(plane,planeLimits);
    stopped(planeRetryStop,algo::TetCandidateReason::BudgetExceeded);
    require(planeRetryStop.metrics.originalCoordinateRetryUsed &&
            planeRetryStop.metrics.initialCoordinateFailure==algo::TetCandidateReason::NonVolumetricInput &&
            planeRetryStop.metrics.boundedWork==planeLimits.maxWork,
            "second constructor attempt received a fresh work allowance");

    const std::vector<Vec3> retry{
        {10.856053358834119,16.79720528749192,10.0},
        {15.937782728263846,10.343455984662688,10.0},
        {15.984496655211151,11.002740097382754,10.91419796248158},
        {15.76057402256841,10.848605589300725,10.991058967623033},
        {15.57124299398417,11.580190611268433,10.91419796248158},
        {15.353109523788468,11.417966749978774,10.991058967623033}};
    limits.maxCells=1000;
    const auto retryFull=algo::buildTetCandidate(retry,limits);
    const auto retryUnlimited=algo::buildTetCandidate(retry);
    require(retryFull.built() && retryUnlimited.built() &&
            retryFull.metrics.originalCoordinateRetryUsed &&
            retryFull.metrics.construction.discardedConstructionWork>0,
            "exact-coordinate retry fixture did not retry");
    same(*retryFull.candidate,*retryUnlimited.candidate);
    const std::vector<Vec3> cloud{
        {-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},
        {-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1},{0,0,0},
        {.125,.25,-.375},{-.375,.5,.625},{.75,-.625,.375},{0,0,0}};
    const auto cloudBounded=algo::buildTetCandidate(cloud,limits);
    const auto cloudUnlimited=algo::buildTetCandidate(cloud);
    require(cloudBounded.built() && cloudUnlimited.built(),"ordinary cloud construction failed");
    same(*cloudBounded.candidate,*cloudUnlimited.candidate);
    small=limits; small.maxWork=retryFull.metrics.boundedWork-1;
    const auto retryStop=algo::buildTetCandidate(retry,small);
    stopped(retryStop,algo::TetCandidateReason::BudgetExceeded);
    require(retryStop.metrics.construction.originalCoordinateRetryUsed &&
            retryStop.metrics.boundedWork<=small.maxWork,
            "coordinate retry reset or lost the shared work budget");
    small=limits; small.maxTemporaryBytes=retryFull.metrics.boundedPeakTemporaryBytes-1;
    const auto memoryStop=algo::buildTetCandidate(retry,small);
    stopped(memoryStop,algo::TetCandidateReason::BudgetExceeded);
    require(memoryStop.metrics.boundedPeakTemporaryBytes<=small.maxTemporaryBytes,
            "payload quota was admitted after growth");
    require(std::memcmp(tetra.data(),original.data(),tetra.size()*sizeof(Vec3))==0,
            "bounded construction mutated caller input");
    std::vector<Vec3> representative;
    std::uint32_t state=0x6d2b79f5U;
    auto coordinate=[&]() {
        state=1664525U*state+1013904223U;
        return double(state)/2147483648.0-1.0;
    };
    for (unsigned i=0;i<256;++i)
        representative.push_back({coordinate(),coordinate(),coordinate()});
    limits.maxCells=10000;
    const auto representativeBuild=algo::buildTetCandidate(representative,limits);
    require(representativeBuild.built(),"representative bounded cloud failed");
    std::printf("representative points=%zu cells=%zu bounded_work=%llu peak_payload_bytes=%llu\n",
                representative.size(),representativeBuild.candidate->cells.size(),
                static_cast<unsigned long long>(representativeBuild.metrics.boundedWork),
                static_cast<unsigned long long>(representativeBuild.metrics.boundedPeakTemporaryBytes));
    std::puts("tet candidate limits passed");
    return 0;
} catch (const std::exception& error) {
    std::fprintf(stderr,"tet candidate limits: %s\n",error.what()); return 1;
} }
