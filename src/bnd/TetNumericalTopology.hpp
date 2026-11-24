// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "TetBistellarRecovery.hpp"

namespace paacs::bnd::detail {

// Private in-process handoff from the certified initial-point owner. This is
// the complete numerical complex, including exterior cells and eight guards;
// it is neither a constrained material mesh nor an external import format.
// No callbacks capturing the former owner may survive its move into this box.
struct TetNumericalTopology {
    BistellarComplex complex;
    std::map<BistellarComplex::Face,BistellarComplex::Face> sourceFaces;
    std::size_t sourceNodes=0,liveCells=0;
    mesh::Classification region{};
    bool available=false,originalCoordinateRetryUsed=false;

    std::uint64_t retainedModeledBytes()const noexcept {
        // Same 4096/cell unprotected-complex model as initial compaction,
        // plus retained vertex/source maps and a fixed scratch envelope.
        std::uint64_t bytes=8ULL*1024*1024;
        for(const auto [count,unit]:std::array<std::pair<std::uint64_t,std::uint64_t>,3>{{
            {complex.maximumStoredCells,4096},{complex.vertices.size(),512},{sourceFaces.size(),512}}}){
            if(count>(UINT64_MAX-bytes)/unit)return UINT64_MAX;
            bytes+=count*unit;
        }
        return bytes;
    }
};

} // namespace paacs::bnd::detail
