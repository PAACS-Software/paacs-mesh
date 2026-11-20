// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/bnd/TetPlcCavityRecovery.hpp"
#include "paacs/bnd/TetFacetGeometry.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <iostream>
#include <stdexcept>
using namespace paacs;
void check(bool value,const char* why){if(!value)throw std::runtime_error(why);}
std::vector<std::uint8_t> bytes(const mesh::Mesh& m){std::vector<std::uint8_t> out;check(sdk::serializeMeshArchive(m,out).completed,"encode");return out;}
mesh::Mesh source(bool pair=false){
    mesh::Mesh m;const auto add=[&](Vec3 offset){const Index begin=m.nodeCount();
        for(auto p:std::array<Vec3,4>{{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}})m.addNode(p+offset,{0,int(m.nodeCount()+1)});
        bnd::TetFacetCell c{{unsigned(begin),unsigned(begin+1),unsigned(begin+2),unsigned(begin+3)},1,{3,9}};
        for(unsigned side=0;side<4;++side){const auto f=bnd::tet_facet_detail::outwardFace(c,side);m.addCell(mesh::CellType::Tri3,std::array<Index,3>{f[0],f[1],f[2]},{2,int(m.cellCount()+1)});}
    };add({});if(pair)add({4,0,0});return m;
}
mesh::Mesh prefix(const mesh::Mesh& s){mesh::Mesh m;for(Index i=0;i<s.nodeCount();++i)m.addNode(s.node(i),s.nodeClass(i));return m;}
int main(){try{
    const auto s=source();auto partial=prefix(s);bnd::TetPlcCavityOptions options;options.region={3,9};options.maximumWork=2000000;
    options.patchSearch.maximumStates=64;options.useVisibilityExpansion=false;
    options.interiorAdmission=[](Vec3 p){return p.x<.24;};mesh::Mesh output=s;const auto sentinel=bytes(output);
    auto result=bnd::recoverTetPlcCavities(s,partial,options,output);
    check(!result.published&&bytes(output)==sentinel,"strict default must reject the excluded centroid transactionally");
    options.centroidAdmission=[](Vec3 p){return p.x<=.26;};result=bnd::recoverTetPlcCavities(s,partial,options,output);
    check(result.published&&mesh::auditLinearVolumeBoundary(output).valid&&std::abs(result.meshVolume-1./6)<1e-12,"explicit cell policy must complete the exact source tetrahedron");
    check(output.nodeCount()==4,"cell policy must not require an excluded interior node");auto complete=bytes(output);
    auto retained=partial;retained.addCell(mesh::CellType::Tet4,std::array<Index,4>{0,1,2,3},{3,9});
    result=bnd::recoverTetPlcCavities(s,retained,options,output);
    check(result.published&&mesh::auditLinearVolumeBoundary(output).valid&&std::abs(result.meshVolume-1./6)<1e-12,"retained partial cells must use the same centroid policy");complete=bytes(output);
    auto invalid=partial;invalid.addNode({.25,.25,.25},{3,9});result=bnd::recoverTetPlcCavities(s,invalid,options,output);
    check(!result.published&&result.reason==bnd::TetPlcCavityReason::InvalidPartialMesh&&bytes(output)==complete,"cell policy must not admit an excluded interior node");
    options.centroidAdmission=[](Vec3 p){return p.x<=.2;};result=bnd::recoverTetPlcCavities(s,partial,options,output);
    check(!result.published&&bytes(output)==complete,"explicit centroid bound must remain enforced");
    options.centroidAdmission=[](Vec3)->bool{throw std::runtime_error("query failure");};result=bnd::recoverTetPlcCavities(s,partial,options,output);
    check(!result.published&&bytes(output)==complete,"throwing cell policy must not publish");
    options.interiorAdmission={};options.centroidAdmission=[](Vec3){return true;};const auto separate=source(true);auto p2=prefix(separate);
    result=bnd::recoverTetPlcCavities(separate,p2,options,output);
    check(result.published&&mesh::auditLinearVolumeBoundary(output).valid&&std::abs(result.meshVolume-1./3)<1e-12,"cell policy must retain disconnected volumes and their physical gap");
    const auto saved=bytes(output);std::stop_source stop;stop.request_stop();options.cancellation=stop.get_token();
    result=bnd::recoverTetPlcCavities(separate,p2,options,output);
    check(!result.published&&result.reason==bnd::TetPlcCavityReason::Cancelled&&bytes(output)==saved,"cancelled policy recovery must retain output");
    std::cout<<"Separated centroid/node admission, retained cells, bounds, exceptions, physical gap and cancellation checks passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
