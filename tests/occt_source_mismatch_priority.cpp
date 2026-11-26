// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/geometry/OcctSourceMismatchPriority.hpp"
#include <iostream>
#include <numeric>
#include <stdexcept>

using namespace paacs;
using namespace paacs::geo::detail;
namespace {
unsigned checks=0;
void check(bool value,const char* message){++checks;if(!value)throw std::runtime_error(message);}
mesh::Mesh source(unsigned components){
    mesh::Mesh mesh;
    for(unsigned component=0;component<components;++component){
        const auto first=mesh.nodeCount();
        for(auto p:std::array{Vec3{0,0,0},Vec3{1,0,0},Vec3{0,1,0},Vec3{0,0,1}}){
            p.x+=2*component;mesh.addNode(p,{0,int(mesh.nodeCount())+1});
        }
        bnd::TetFacetCell tet;
        for(unsigned j=0;j<4;++j)tet.vertices[j]=static_cast<bnd::TetFacetVertexId>(first+j);
        if(predicates::orient3d(mesh.node(tet.vertices[0]),mesh.node(tet.vertices[1]),mesh.node(tet.vertices[2]),mesh.node(tet.vertices[3]))==predicates::Sign::Negative)
            std::swap(tet.vertices[1],tet.vertices[2]);
        check(predicates::orient3d(mesh.node(tet.vertices[0]),mesh.node(tet.vertices[1]),mesh.node(tet.vertices[2]),mesh.node(tet.vertices[3]))==predicates::Sign::Positive,
            "owned tetrahedron orientation");
        for(unsigned side=0;side<4;++side){const auto f=bnd::tet_facet_detail::outwardFace(tet,side);
            mesh.addCell(mesh::CellType::Tri3,{f[0],f[1],f[2]},{2,int(mesh.cellCount())+1});}
    }return mesh;
}
std::vector<Index> candidateIds(const mesh::Mesh& source){
    std::vector<Index> result(source.cellCount());std::iota(result.begin(),result.end(),0);return result;
}
bool sameRanks(const OcctSourceMismatchOutcome& a,const OcctSourceMismatchOutcome& b){
    if(a.facets.size()!=b.facets.size())return false;
    for(std::size_t i=0;i<a.facets.size();++i)if(a.facets[i].triangle!=b.facets[i].triangle||
        a.facets[i].insideSourceWitnesses!=b.facets[i].insideSourceWitnesses||
        a.facets[i].nearestDistanceSquared!=b.facets[i].nearestDistanceSquared)return false;
    return true;
}
void exercise(){
    auto boundary=source(1);const auto hash=mesh::meshHash(boundary);auto ids=candidateIds(boundary);
    const OcctSourceMismatchStamp identity{7,3,hash};OcctSourceMismatchLimits limits;
    const std::array points{Vec3{.1,.1,.1},Vec3{2,2,2},Vec3{.25,.25,0},Vec3{0,0,0},Vec3{.5,0,0}};
    const auto result=rankOcctSourceMismatchFacets(boundary,identity,identity,points,ids,limits);
    check(result.complete&&result.reason==OcctSourceMismatchReason::None,"owned source mismatch diagnosis completed");
    check(result.insideSourceOutsideCad==1&&result.outsideSource==1&&result.ambiguousSource==3&&result.nativeRejectedPoints==points.size(),
        "strict source inclusion distinguishes exterior and face/edge/vertex contacts");
    check(result.facets.size()==ids.size()&&result.prioritizedFacets==3&&result.insideWitnessesNearMissingFacets==1,
        "only inside-source Native-rejected points prioritize nearby missing facets");
    check(result.facets.front().insideSourceWitnesses==1&&result.facets.back().insideSourceWitnesses==0,
        "supported triangles precede unsupported fallback without implying acceptance");
    const auto repeated=rankOcctSourceMismatchFacets(boundary,identity,identity,points,ids,limits);
    check(repeated.complete&&sameRanks(result,repeated)&&result.workUsed==repeated.workUsed,
        "source mismatch ranking and logical work are deterministic");
    std::reverse(ids.begin(),ids.end());
    const auto reversed=rankOcctSourceMismatchFacets(boundary,identity,identity,points,ids,limits);
    check(reversed.complete&&sameRanks(result,reversed),"input candidate order cannot override geometric witness ranking");
    const std::array exterior{Vec3{2,2,2}};
    const auto outside=rankOcctSourceMismatchFacets(boundary,identity,identity,exterior,ids,limits);
    check(outside.complete&&outside.outsideSource==1&&!outside.prioritizedFacets&&
        std::is_sorted(outside.facets.begin(),outside.facets.end(),[](const auto& a,const auto& b){return a.triangle<b.triangle;}),
        "outside-source bridges receive the original deterministic fallback order");
    auto many=source(12);auto manyIds=candidateIds(many);
    const OcctSourceMismatchStamp manyIdentity{19,0,mesh::meshHash(many)};
    const std::array later{Vec3{18.1,.1,.1},Vec3{18.15,.1,.1}};
    const auto localized=rankOcctSourceMismatchFacets(many,manyIdentity,manyIdentity,later,manyIds,limits);
    check(localized.complete&&localized.insideSourceOutsideCad==2&&localized.prioritizedFacets>0,
        "disconnected owned source components retain exact inclusion");
    check(localized.facets.front().triangle>=36&&localized.facets.front().triangle<=39&&
        localized.facets.front().insideSourceWitnesses==2,"geometric mismatch beyond the first32 facets takes priority");
    auto sorted=localized.facets;std::sort(sorted.begin(),sorted.end(),[](const auto& a,const auto& b){return a.triangle<b.triangle;});
    check(sorted.size()==manyIds.size(),"ranked output preserves all supplied candidates");
    for(std::size_t i=0;i<sorted.size();++i)check(sorted[i].triangle==manyIds[i],"ranked output neither duplicates nor loses source IDs");
    const auto refuse=[&](const mesh::Mesh& input,OcctSourceMismatchStamp current,OcctSourceMismatchStamp expected,
                          std::span<const Vec3> query,std::span<const Index> candidates,OcctSourceMismatchLimits cap,
                          OcctSourceMismatchReason reason,const char* diagnostic){
        const auto refused=rankOcctSourceMismatchFacets(input,current,expected,query,candidates,cap);
        check(!refused.complete&&refused.reason==reason&&refused.facets.empty()&&!refused.prioritizedFacets,
            "refused optional diagnosis returns no partial priorities");
        check(std::string_view(refused.diagnostic)==diagnostic,"source mismatch refusal retains a precise diagnostic");
    };
    auto stale=identity;++stale.generation;
    refuse(boundary,identity,stale,points,ids,limits,OcctSourceMismatchReason::InvalidInput,"source-mismatch-stale-identity");
    stale=identity;++stale.bodyTag;
    refuse(boundary,identity,stale,points,ids,limits,OcctSourceMismatchReason::InvalidInput,"source-mismatch-stale-identity");
    auto changed=boundary;changed.setNode(0,{.001,0,0});
    refuse(changed,identity,identity,points,ids,limits,OcctSourceMismatchReason::InvalidInput,"source-mismatch-source-changed");
    const std::array badIds{Index{-1}};
    refuse(boundary,identity,identity,points,badIds,limits,OcctSourceMismatchReason::InvalidInput,"source-mismatch-candidate-index");
    const std::array duplicates{Index{0},Index{0}};
    refuse(boundary,identity,identity,points,duplicates,limits,OcctSourceMismatchReason::InvalidInput,"source-mismatch-candidate-index");
    const std::array nonfinite{Vec3{std::numeric_limits<double>::quiet_NaN(),0,0}};
    refuse(boundary,identity,identity,nonfinite,ids,limits,OcctSourceMismatchReason::InvalidInput,"source-mismatch-point");
    auto cap=limits;cap.maximumWork=reversed.workUsed-1;
    refuse(boundary,identity,identity,points,ids,cap,OcctSourceMismatchReason::BudgetExceeded,"source-mismatch-work-budget");
    cap=limits;cap.maximumWork=reversed.workUsed;
    check(rankOcctSourceMismatchFacets(boundary,identity,identity,points,ids,cap).complete,"exact logical work allowance remains sufficient");
    cap=limits;cap.maximumWorkingBytes=result.modeledBytes-1;
    refuse(boundary,identity,identity,points,ids,cap,OcctSourceMismatchReason::BudgetExceeded,"source-mismatch-working-budget");
    cap=limits;cap.maximumRejectedPoints=points.size()-1;
    refuse(boundary,identity,identity,points,ids,cap,OcctSourceMismatchReason::BudgetExceeded,"source-mismatch-count-budget");
    cap=limits;cap.maximumCandidates=ids.size()-1;
    refuse(boundary,identity,identity,points,ids,cap,OcctSourceMismatchReason::BudgetExceeded,"source-mismatch-count-budget");
    std::atomic<bool> cancelled{true};cap=limits;cap.cancellation=&cancelled;
    refuse(boundary,identity,identity,points,ids,cap,OcctSourceMismatchReason::Cancelled,"cancelled");
    check(mesh::meshHash(boundary)==hash&&mesh::meshHash(many)==manyIdentity.sourceMeshHash,
        "all diagnostics preserve immutable source coordinates, classifications and connectivity");
}
}
int main(){try{exercise();std::cout<<checks<<" source mismatch priority checks passed\n";return 0;}
catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
