// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/mesh/Mesh.hpp"
#include <cstddef>
#include <span>

namespace paacs::bnd::detail {
class BistellarComplex;

// Private synchronous views for bounded diagnostic executables. The observer
// owns its work/output limits and must not retain pointers, mutate construction,
// perform CAD queries, or throw. No observer means no scans or allocations.
// These numerical complexes may include guards/exterior cells; this is not a
// geometry, material, quality or publication certificate.
struct TetConstructionObservation {
    const char* stage;
    int bodyTag;
    const BistellarComplex* complex=nullptr;
    const mesh::Mesh* mesh=nullptr;
    std::span<const unsigned char> relevance{};
    std::size_t guardBegin=0,guardCount=0;
    std::size_t firstCell=0,cellCount=0;
    std::span<const std::size_t> removedCells{};
    bool completeSnapshot=true;
    bool originalCoordinateRetryUsed=false;
};
using TetConstructionObserver=void (*)(const TetConstructionObservation&,void*) noexcept;
struct TetConstructionObserverState {
    TetConstructionObserver observer=nullptr;
    void* context=nullptr;
};
inline thread_local TetConstructionObserverState tetConstructionObserver;

class ScopedTetConstructionObserver {
public:
    ScopedTetConstructionObserver(TetConstructionObserver observer,void* context) noexcept
        : previous_(tetConstructionObserver) { tetConstructionObserver={observer,context}; }
    ~ScopedTetConstructionObserver() { tetConstructionObserver=previous_; }
    ScopedTetConstructionObserver(const ScopedTetConstructionObserver&)=delete;
    ScopedTetConstructionObserver& operator=(const ScopedTetConstructionObserver&)=delete;
private:
    TetConstructionObserverState previous_;
};
inline void observeTetConstruction(const TetConstructionObservation& observation) noexcept {
    if(tetConstructionObserver.observer)
        tetConstructionObserver.observer(observation,tetConstructionObserver.context);
}
inline void observeTetMesh(const char* stage,int bodyTag,const mesh::Mesh& mesh) noexcept {
    if(tetConstructionObserver.observer)
        observeTetConstruction({stage,bodyTag,nullptr,&mesh});
}
} // namespace paacs::bnd::detail
