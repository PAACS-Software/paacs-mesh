// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0
#include "paacs/geometry/CadBoundaryExport.hpp"
#include "paacs/geometry/OcctBoundary.hpp"
#include "../src/geometry/OcctModelInternal.hpp"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

using namespace paacs;
namespace {
void check(bool value,const char* message) { if(!value)throw std::runtime_error(message); }
auto bits(Vec3 p) { return std::array{std::bit_cast<std::uint64_t>(p.x),
    std::bit_cast<std::uint64_t>(p.y),std::bit_cast<std::uint64_t>(p.z)}; }
auto edgeKey(Index a,Index b) { return std::array{std::min(a,b),std::max(a,b)}; }
mesh::Mesh sentinel() {
    mesh::Mesh result;result.addNode({37,19,-11},{0,91});return result;
}
void unchanged(const mesh::Mesh& value) {
    check(value.nodeCount()==1&&value.cellCount()==0&&bits(value.node(0))==bits({37,19,-11})&&
          value.nodeClass(0)==mesh::Classification{0,91},"failed extraction changed destination");
}
mesh::Mesh boundaryOf(const geo::OcctModel& model,int body) {
    mesh::Mesh result;geo::OcctBoundaryOptions options;
    options.maximumWork=50'000'000;options.maximumChordError=.06;
    const auto built=geo::buildOcctBoundary(model,body,options,result);
    if(!built.published)throw std::runtime_error(std::string("boundary: ")+geo::occtBoundaryReasonName(built.reason));
    return result;
}
void prove(const geo::OcctModel& model,int body,const mesh::Mesh& source,bool expectClosed,bool expectSeam,bool expectDegenerate) {
    geo::CadBoundaryEdgeOptions options;options.maximumChordError=.06;
    mesh::Mesh lines;const auto before=source;
    const auto outcome=geo::exportCadBoundaryEdges(model,body,source,options,lines);
    if(!outcome.complete)throw std::runtime_error(std::string("CAD edges: ")+geo::cadBoundaryEdgeReasonName(outcome.reason));
    check(!outcome.destinationUnchanged&&outcome.reason==geo::CadBoundaryEdgeReason::None,"completion state inconsistent");
    check(lines.nodeCount()==source.nodeCount()&&std::uint64_t(lines.cellCount())==outcome.segments,"node space/count changed");
    for(Index i=0;i<source.nodeCount();++i)
        check(bits(lines.node(i))==bits(source.node(i))&&lines.nodeClass(i)==source.nodeClass(i)&&
              bits(source.node(i))==bits(before.node(i)),"edge extraction moved or reclassified a source node");
    std::set<std::array<Index,2>> triangleEdges,usedEdges;
    for(Index i=0;i<source.cellCount();++i){const auto* v=source.cellNodes(i);
        for(unsigned j=0;j<3;++j)triangleEdges.insert(edgeKey(v[j],v[(j+1)%3]));}
    bool closed=false,seam=false,degenerated=false;std::set<Index> usedNodes;
    Index next=0;
    for(const auto& curve:outcome.edges) {
        check(std::binary_search(curve.bodyTags.begin(),curve.bodyTags.end(),body),"CAD body identity missing");
        check(!curve.faceTags.empty()&&!curve.uses.empty(),"CAD face/wire provenance missing");
        for(const auto& use:curve.uses)seam|=use.seam;
        if(curve.degenerated){degenerated=true;check(curve.naturalLineOrder.empty(),"degenerate edge became a segment");continue;}
        check(curve.naturalParameters.size()==curve.naturalNodeOrder.size()&&
              curve.naturalLineOrder.size()+1==curve.naturalNodeOrder.size(),"ordered chain schema mismatch");
        if(curve.closed){closed=true;check(curve.naturalNodeOrder.front()==curve.naturalNodeOrder.back(),"closed edge lost closure");}
        for(std::size_t j=1;j<curve.naturalNodeOrder.size();++j){
            check(curve.naturalParameters[j]>curve.naturalParameters[j-1],"native parameter order lost");
            check(curve.naturalLineOrder[j-1]==next,"line IDs not preserved in deterministic edge order");
            const auto c=lines.cell(next);const auto* v=lines.cellNodes(next++);
            check(c.type==mesh::CellType::Line2&&c.onEntity==mesh::Classification{1,curve.edgeTag},"not a tagged CAD Line2");
            check(v[0]==curve.naturalNodeOrder[j-1]&&v[1]==curve.naturalNodeOrder[j],"line connectivity differs from canonical record");
            const auto key=edgeKey(v[0],v[1]);
            check(triangleEdges.contains(key)&&usedEdges.insert(key).second,"invented or duplicate CAD segment");
            usedNodes.insert(v[0]);usedNodes.insert(v[1]);
        }
    }
    check(next==lines.cellCount()&&usedNodes.size()+outcome.unusedNodes==std::size_t(lines.nodeCount()),"unused-node evidence differs");
    check(!expectClosed||closed,"closed periodic curve not exercised");
    check(!expectSeam||seam,"seam use not retained");check(!expectDegenerate||degenerated,"degenerate curve not retained");
    check(usedEdges.size()<triangleEdges.size(),"CAD edges were confused with all triangulation edges");
}
void refusals(const geo::OcctModel& model,const mesh::Mesh& source) {
    geo::CadBoundaryEdgeOptions options;options.maximumChordError=.06;auto output=sentinel();
    options.maximumWork=1;auto result=geo::exportCadBoundaryEdges(model,1,source,options,output);
    check(!result.complete&&result.reason==geo::CadBoundaryEdgeReason::BudgetExceeded,"work allowance ignored");unchanged(output);
    options.maximumWork=20'000'000;options.maximumOutputBytes=40;
    result=geo::exportCadBoundaryEdges(model,1,source,options,output);
    check(!result.complete&&result.reason==geo::CadBoundaryEdgeReason::BudgetExceeded,"output allowance ignored");unchanged(output);
    options.maximumOutputBytes=256ULL*1024*1024;
    std::atomic<bool> cancel{true};options.cancellation=&cancel;
    result=geo::exportCadBoundaryEdges(model,1,source,options,output);
    check(!result.complete&&result.reason==geo::CadBoundaryEdgeReason::Cancelled,"cancellation ignored");unchanged(output);
    options.cancellation=nullptr;
    auto alias=source;result=geo::exportCadBoundaryEdges(model,1,alias,options,alias);
    check(!result.complete&&result.reason==geo::CadBoundaryEdgeReason::InvalidRequest&&alias.cellCount()==source.cellCount(),"input alias accepted");
    mesh::Mesh corrupt;
    bool changed=false;
    for(Index i=0;i<source.nodeCount();++i){auto cls=source.nodeClass(i);auto p=source.node(i);
        if(!changed&&cls.dim==1){p.x+=1000;changed=true;}corrupt.addNode(p,cls);}
    check(changed,"projection refusal did not exercise an edge-interior anchor");
    for(Index i=0;i<source.cellCount();++i){const auto* ids=source.cellNodes(i);
        corrupt.addCell(mesh::CellType::Tri3,{ids[0],ids[1],ids[2]},source.cell(i).onEntity);}
    result=geo::exportCadBoundaryEdges(model,1,corrupt,options,output);
    check(!result.complete,"off-curve CAD anchor accepted");unchanged(output);
}
}
int main() {
    try {
        auto box=geo::detail::OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(1.,1.,1.).Shape());
        auto cylinder=geo::detail::OcctModelAccess::fromShape(BRepPrimAPI_MakeCylinder(1.,2.).Shape());
        auto sphere=geo::detail::OcctModelAccess::fromShape(BRepPrimAPI_MakeSphere(1.).Shape());
        check(box.imported()&&cylinder.imported()&&sphere.imported(),"owned CAD import failed");
        const auto boxBoundary=boundaryOf(*box.model,1),cylinderBoundary=boundaryOf(*cylinder.model,1);
        prove(*box.model,1,boxBoundary,false,false,false);
        prove(*cylinder.model,1,cylinderBoundary,true,true,false);
        prove(*sphere.model,1,boundaryOf(*sphere.model,1),false,true,true);
        refusals(*cylinder.model,cylinderBoundary);
        BRep_Builder builder;TopoDS_Compound assembly;builder.MakeCompound(assembly);
        builder.Add(assembly,BRepPrimAPI_MakeBox(1.,1.,1.).Shape());
        builder.Add(assembly,BRepPrimAPI_MakeBox(gp_Pnt(3,0,0),1.,1.,1.).Shape());
        auto parts=geo::detail::OcctModelAccess::fromShape(assembly);
        check(parts.imported()&&parts.solids==2,"independent body fixture failed");
        prove(*parts.model,2,boundaryOf(*parts.model,2),false,false,false);
        std::cout<<"CAD edge export: actual surface segments, periodic/seam/degenerate uses, body tags and atomic refusals passed\n";
        return 0;
    } catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
