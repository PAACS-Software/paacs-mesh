// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctSurfaceAliasRepair.hpp"
#include "OcctStarMesher.hpp"
#include <stdexcept>
#include <cstdio>
using namespace paacs;using namespace paacs::geo::detail;
namespace {void need(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}}
int main(){try{
    // Two charts share a bent native curve A-C-B. Both triangulations chose
    // the shortcut A-B and produced opposite copies of triangle A-C-B.
    // Face10 represents a shared cut; only source face20 may change.
    const std::vector<Vec3> p{{0,0,0},{0,1,0},{0,.5,.1},{1,.5,.5},{-1,.5,.5}};
    const std::vector<OcctAliasTriangle> source{{1,0,2},{3,0,1},{1,2,0},{4,1,0},{0,3,4},{3,1,4}};
    const std::array tags{10,10,20,20,30,31};const std::array allowed{20};
    const std::set<OcctAliasEdge> wires{{0,2},{1,2},{0,3},{1,3},{0,4},{1,4},{3,4}};
    const auto position=[&](int i){return p[i];};std::uint64_t work=0;
    const auto charge=[&](std::uint64_t n){work+=n;};
    const auto admit=[](int,auto,auto){return true;};
    const auto meshOf=[&](const auto& faces){mesh::Mesh m;for(unsigned i=0;i<p.size();++i)m.addNode(p[i],{0,int(i)+1});for(unsigned i=0;i<faces.size();++i){const auto f=faces[i];m.addCell(mesh::CellType::Tri3,std::array<Index,3>{f[0],f[1],f[2]},{2,tags[i]});}return m;};
    OcctStarOptions limits;const auto broken=meshOf(source);need(!certifyOcctBoundarySurface(broken,limits).certified,"fixture must reproduce invalid boundary");
    auto repaired=source;need(repairOcctSurfaceAliases(repaired,tags,wires,allowed,1,position,charge,admit)==1,"one bounded diagonal repair");
    need(repaired[0]==source[0]&&repaired[1]==source[1]&&repaired[4]==source[4]&&repaired[5]==source[5],"shared cut and other faces changed");
    const auto edgePresent=[&](OcctAliasEdge edge){for(auto f:repaired)for(unsigned j=0;j<3;++j)if(OcctAliasEdge{std::min(f[j],f[(j+1)%3]),std::max(f[j],f[(j+1)%3])}==edge)return true;return false;};
    for(auto edge:wires)need(edgePresent(edge),"native curve segment lost");
    auto fixed=meshOf(repaired);const auto certificate=certifyOcctBoundarySurface(fixed,limits);need(certificate.certified&&certificate.eulerCharacteristic==2&&certificate.connectedShells==1,"repair must produce an embedded closed sphere");
    mesh::Mesh volume;const auto filled=buildOcctStarVolume(fixed,1,limits,volume);need(filled.published&&filled.minimumQuality>0,"repaired region must fill with positive valid tetrahedra");
    need(work>0&&repairOcctSurfaceAliases(repaired,tags,wires,allowed,1,position,charge,admit)==0,"completed repair must be idempotent");
    auto unchanged=source;need(repairOcctSurfaceAliases(unchanged,tags,wires,allowed,0,position,charge,admit)==0&&unchanged==source,"zero disables changes");
    auto locked=wires;locked.insert({0,1});need(repairOcctSurfaceAliases(unchanged,tags,locked,allowed,1,position,charge,admit)==0&&unchanged==source,"protected native diagonal must stay fixed");
    const std::array foreign{30};need(repairOcctSurfaceAliases(unchanged,tags,wires,foreign,1,position,charge,admit)==0&&unchanged==source,"unselected CAD face must stay fixed");
    need(repairOcctSurfaceAliases(unchanged,tags,wires,allowed,1,position,charge,[](int,auto,auto){return false;})==0&&unchanged==source,"failed geometry admission must preserve triangles");
    bool stopped=false;try{repairOcctSurfaceAliases(unchanged,tags,wires,allowed,1,position,[](std::uint64_t){throw 1;},admit);}catch(int){stopped=true;}
    need(stopped&&unchanged==source,"work/cancellation stop before mutation");
    std::puts("Opposite-facet diagonal repair, fixed curves/cuts, geometry refusal, and embedded positive refill passed");return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
