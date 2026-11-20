// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <paacs/bnd/TetPlcCavityRecovery.hpp>
#include <paacs/mesh/LinearVolumeAudit.hpp>
#include <iostream>
#include <array>
int main(){using namespace paacs;mesh::Mesh source,partial,completed;
 for(auto p:std::array<Vec3,4>{{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}})source.addNode(p,{0,static_cast<int>(source.nodeCount())+1});
 for(auto f:std::array<std::array<Index,3>,4>{{{1,2,3},{0,3,2},{0,1,3},{0,2,1}}})source.addCell(mesh::CellType::Tri3,f,{2,static_cast<int>(source.cellCount())+1});
 bnd::TetPlcCavityOptions options;options.region={3,1};options.maximumRandomCandidatePoints=8;options.randomSeed=17;options.maximumMergedComponents=1;const auto result=bnd::recoverTetPlcCavities(source,partial,options,completed);
 if(!result.published||!mesh::auditLinearVolumeBoundary(completed).valid||std::abs(result.meshVolume-1./6)>1e-14)return 1;
 std::cout<<"installed public API: complete volume, source classes retained\n";return 0;
}
