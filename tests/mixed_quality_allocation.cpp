// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MixedMeshQuality.hpp"

#include <array>
#include <bit>
#include <cstdlib>
#include <iostream>
#include <new>
#include <stdexcept>

namespace {
// Armed only around single-threaded operations under test. Once the selected
// allocation fails, subsequent allocations fail too until the call returns.
long allocationsBeforeFailure = -1;
std::size_t allocationCount = 0;
bool countAllocations = false;
}
void* operator new(std::size_t size) {
    if (countAllocations) ++allocationCount;
    if (allocationsBeforeFailure == 0) throw std::bad_alloc{};
    if (allocationsBeforeFailure > 0) --allocationsBeforeFailure;
    if (void* p = std::malloc(size ? size : 1)) return p;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }

using namespace paacs;
namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
mesh::Mesh fixture(bool inverted = false) {
    mesh::Mesh value;
    for (const auto p : std::array<Vec3, 4>{{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}})
        (void)value.addNode(p, {2,7});
    const std::array<Index, 4> ids = inverted ? std::array<Index,4>{0,2,1,3}
                                             : std::array<Index,4>{0,1,2,3};
    (void)value.addCell(mesh::CellType::Tet4, ids, {3,1});
    return value;
}
bool unchanged(const mesh::Mesh& a, const mesh::Mesh& b) {
    if (a.nodeCount()!=b.nodeCount() || a.cellCount()!=b.cellCount() ||
        a.connectivityCount()!=b.connectivityCount()) return false;
    for (Index i=0;i<a.nodeCount();++i) {
        const auto x=a.node(i), y=b.node(i);
        if (std::bit_cast<std::uint64_t>(x.x)!=std::bit_cast<std::uint64_t>(y.x) ||
            std::bit_cast<std::uint64_t>(x.y)!=std::bit_cast<std::uint64_t>(y.y) ||
            std::bit_cast<std::uint64_t>(x.z)!=std::bit_cast<std::uint64_t>(y.z) ||
            a.nodeClass(i)!=b.nodeClass(i)) return false;
    }
    for (Index i=0;i<a.cellCount();++i) {
        if (a.cell(i).type!=b.cell(i).type || a.cell(i).onEntity!=b.cell(i).onEntity)
            return false;
        for (int j=0;j<mesh::nodesPerCell(a.cell(i).type);++j)
            if (a.cellNodes(i)[j]!=b.cellNodes(i)[j]) return false;
    }
    return true;
}
}
int main() {
    try {
        auto value=fixture(); const auto saved=value;
        const auto healthy=algo::analyzeMixedMeshQuality(value);
        require(healthy.complete && healthy.valid && !healthy.allocationFailed,
                "healthy quality report is complete without allocation failure");
        const auto invalid=algo::analyzeMixedMeshQuality(fixture(true));
        require(invalid.complete && !invalid.valid && !invalid.allocationFailed,
                "inverted geometry is distinguished from allocation failure");
        allocationsBeforeFailure=0;
        const auto failed=algo::analyzeMixedMeshQuality(value);
        allocationsBeforeFailure=-1;
        require(!failed.complete && !failed.valid && failed.allocationFailed && unchanged(value,saved),
                "quality allocation failure is precise and leaves source unchanged");

        // Boundary-classified nodes cannot move. Exhaust every allocation in
        // both quality analyses and optimizer adjacency preparation, including
        // the final report after an otherwise successful optimization.
        countAllocations=true; allocationCount=0;
        const auto healthyOptimization=algo::optimizeMixedMeshQuality(value);
        countAllocations=false;
        const auto count=allocationCount;
        require(healthyOptimization.completed && count>2 && unchanged(value,saved),
                "optimizer fixture succeeds and retains protected nodes");
        bool sawBefore=false, sawAfter=false;
        for (std::size_t i=0;i<count;++i) {
            allocationsBeforeFailure=static_cast<long>(i);
            const auto out=algo::optimizeMixedMeshQuality(value);
            allocationsBeforeFailure=-1;
            require(!out.completed && out.reason==algo::MixedMeshOptimizationReason::AllocationFailed &&
                    unchanged(value,saved), "every injected optimizer allocation is typed and protected source retained");
            sawBefore=sawBefore||out.before.allocationFailed;
            sawAfter=sawAfter||out.after.allocationFailed;
        }
        require(sawBefore && sawAfter, "initial and final swallowed quality allocation failures both propagated");
        std::cout<<"quality reports: valid/inverted/persistent-OOM distinct; optimizer "<<count
                 <<" allocation sites including both reports preserve protected source\n";
        return 0;
    } catch (const std::exception& e) {
        allocationsBeforeFailure=-1; countAllocations=false;
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
