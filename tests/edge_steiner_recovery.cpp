// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Retained V0.10 non-star witness, exercised through the current indexed route.
#include "paacs/bnd/TetFacetSetRecovery.hpp"
#include "paacs/bnd/TetRecoveryIndex.hpp"
#include "paacs/topo/TetRefinement.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include "allocation_nothrow.hpp"

namespace {
thread_local std::int64_t failOrdinal=-1;
thread_local std::uint64_t allocations=0;
}
void* operator new(std::size_t size) {
    if (failOrdinal>=0) {
        ++allocations;
        if (failOrdinal--==0) { failOrdinal=-1;throw std::bad_alloc(); }
    }
    if (auto* p=std::malloc(size?size:1)) return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p,std::size_t) noexcept { std::free(p); }
void operator delete[](void* p,std::size_t) noexcept { std::free(p); }

using namespace paacs;
namespace {
void require(bool yes,const char* why) { if(!yes) throw std::runtime_error(why); }
struct Witness {
    std::vector<bnd::TetFacetVertex> vertices{
        {{1,-2.6000000000000001,-2},24},{{2,-4,-4},44},
        {{2,-4,-1.3333333333333335},45},{{2,-1.9000000000000004,-1},48},
        {{3,-4,-4},59},{{3,-4,-1.3333333333333335},60},
        {{3,-1.9000000000000004,-1},63}};
    std::vector<bnd::TetFacetCell> cells{
        {{{0,4,5,1}},260,{3,0}},{{{0,2,6,3}},262,{3,0}},
        {{{0,5,6,2}},263,{3,0}},{{{4,5,6,0}},266,{3,0}}};
    bnd::TetEdgeGeneralCavityInput input() const {
        return {vertices,cells,{1,3},{},{},{3,0}};
    }
};
constexpr bnd::TetEdgeGeneralCavityLimits limits{16,16,256,4000000,40000000,8};
struct Fixture {
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::unique_ptr<bnd::TetRecoveryContext> context;
    bnd::TetFacetSetState state;
    bnd::TetEdgeRecoveryRequest request;
    Fixture() {
        Witness witness;topo::NodePool np;std::vector<topo::NodeHandle> handles;
        for(const auto& vertex:witness.vertices) handles.push_back(np.add(vertex.position,{3,0}));
        topo::CellPool cp;
        for(const auto& source:witness.cells) {
            std::array<topo::NodeHandle,4> nodes;
            for(unsigned i=0;i<4;++i)nodes[i]=handles[source.vertices[i]];
            (void)cp.add(topo::ElementType::Tet4,nodes,{3,0},np);
        }
        topo::EdgePool ep;for(const auto c:cp.activeBySequence())(void)ep.insertCellEdges(np,cp,c);
        topo::FacePool fp;for(const auto c:cp.activeBySequence())(void)fp.insertCellFaces(np,cp,ep,c);
        workspace=std::make_unique<topo::TopologyWorkspace>(std::move(np),std::move(cp),std::move(ep),std::move(fp));
        auto built=bnd::buildTetRecoveryContext(*workspace,{3,0});
        require(built.built(),"retained source fixture failed its current context audit");
        context=std::move(built.context);state.workspaceStamp=topo::WorkspaceStamp::of(*workspace);
        state.artifactStamp={27,1};request.stamp=state.workspaceStamp;request.artifactStamp=state.artifactStamp;
        request.vertices={handles[1],handles[3]};request.key=bnd::VolumeArtifactEdgeKey::of(1,3);
        request.region={3,0};request.classification={1,91};
    }
};
double volume(const topo::TopologyWorkspace& workspace) {
    double v=0;
    for(const auto c:workspace.cells().activeBySequence()) {
        const auto cell=workspace.cells().get(c);
        const auto p=[&](unsigned i){return workspace.nodes().get(cell.nodes[i]).coord;};
        v+=std::abs((p(1)-p(0)).dot((p(2)-p(0)).cross(p(3)-p(0))))/6.;
    }
    return v;
}
void plainValueContracts() {
    Witness source;const auto input=source.input();
    require(bnd::solveTetEdgeGeneralCavity(input,limits).status==bnd::TetEdgeGeneralCavityStatus::Unavailable,
            "existing-node no-solution witness changed");
    const auto solved=bnd::solveTetEdgeGeneralCavityOneSteiner(input,limits,8);
    require(solved.solved()&&solved.cavity.cells.size()==12,"restored one-point solve failed");
    auto vertices=source.vertices;vertices.push_back({solved.steinerPosition,64});
    bnd::TetFacetComplexView view{vertices,solved.cavity.cells,{},{},{},{3,0}};
    std::vector<std::size_t> all(solved.cavity.cells.size());
    for(std::size_t i=0;i<all.size();++i)all[i]=i;
    require(bnd::auditTetFacetBall(view,all).valid,"replacement is not a positive manifold ball");
    for(std::size_t i=0;i<all.size();++i)for(std::size_t j=i+1;j<all.size();++j)
        require(!bnd::tet_facet_detail::tetImproperlyIntersects(view,view.cells[i],view.cells[j]),
                "replacement tetrahedra intersect");
    const auto shortRun=bnd::solveTetEdgeGeneralCavityOneSteiner(input,limits,5);
    require(shortRun.status==bnd::TetEdgeGeneralCavityStatus::BudgetReached,"candidate cap lost");
    const auto work=bnd::solveTetEdgeGeneralCavityOneSteiner(input,limits,8,1);
    require(work.status==bnd::TetEdgeGeneralCavityStatus::BudgetReached&&work.cumulativeWorkUnits<=1,"cumulative work cap lost");
    auto small=limits;small.maxCavityVertices=7;
    require(bnd::solveTetEdgeGeneralCavityOneSteiner(input,small,8).status==bnd::TetEdgeGeneralCavityStatus::BudgetReached,
            "staged point bypassed vertex cap");
    std::reverse(source.cells.begin(),source.cells.end());
    const auto reversed=bnd::solveTetEdgeGeneralCavityOneSteiner(source.input(),limits,8);
    require(reversed.solved()&&reversed.cavity.orderedTetKeys==solved.cavity.orderedTetKeys&&
            reversed.steinerPosition.x==solved.steinerPosition.x&&
            reversed.steinerPosition.y==solved.steinerPosition.y&&
            reversed.steinerPosition.z==solved.steinerPosition.z,"source presentation changes recovery");
    source.cells[0].vertices[0]=100;
    require(bnd::solveTetEdgeGeneralCavityOneSteiner(source.input(),limits,8).status==bnd::TetEdgeGeneralCavityStatus::InvalidInput,
            "invalid source connectivity was dereferenced");
}
void indexedContracts(const char* artifact) {
    Fixture f;const auto before=f.workspace->strictWorkspaceHash();const auto v=volume(*f.workspace);
    auto old=limits;old.maximumSteinerCandidates=0;
    require(!bnd::recoverTetEdgeIndexedDetailed(f.request,f.state,*f.workspace,*f.context,old).recovery.succeeded(),
            "disabled route unexpectedly adds a point");
    require(f.workspace->strictWorkspaceHash()==before,"disabled refusal changed topology");
    bnd::TetEdgeGeneralRecoveryDiagnostics evidence;
    const auto assessed=bnd::assessTetEdgeIndexedDetailed(f.request,f.state,*f.workspace,*f.context,limits,&evidence);
    require(assessed.recovery.succeeded()&&evidence.stagedCandidateAvailable,"assessment missed staged recovery");
    require(f.workspace->strictWorkspaceHash()==before&&f.context->currentFor(*f.workspace)&&f.state.edgeBindings.empty(),
            "assessment published state");
    const auto recovered=bnd::recoverTetEdgeIndexedDetailed(f.request,f.state,*f.workspace,*f.context,limits,&evidence);
    require(recovered.recovery.succeeded()&&recovered.recovery.stagedNodeCreated,"indexed route did not commit staged recovery");
    require(f.workspace->nodes().activeBySequence().size()==8&&recovered.createdCells.size()==12&&recovered.retiredCells.size()==4,
            "staged commit lost nodes/cells");
    require(std::abs(volume(*f.workspace)-v)<1e-12,"staged replacement changed the cavity volume");
    require(f.context->currentFor(*f.workspace)&&bnd::buildTetRecoveryContext(*f.workspace,{3,0}).built(),
            "staged context update failed");
    const auto edge=f.workspace->edges().get(recovered.recovery.edge);
    require(edge.requiredConstraint&&edge.cls==f.request.classification&&
            f.state.edgeBindings.at(f.request.key).edge==recovered.recovery.edge,"target edge metadata was not published atomically");
    if(artifact) {
        std::ofstream file(artifact);require(bool(file),"cannot create witness artifact");
        file<<std::setprecision(17);
        const auto nodes=f.workspace->nodes().activeBySequence();std::map<topo::NodeHandle,std::size_t> indices;
        file<<"{\"nodes\":[";
        for(std::size_t i=0;i<nodes.size();++i){indices.emplace(nodes[i],i);const auto p=f.workspace->nodes().get(nodes[i]).coord;
            file<<(i?",":"")<<'['<<p.x<<','<<p.y<<','<<p.z<<']';}
        file<<"],\"cells\":[";bool first=true;
        for(const auto c:f.workspace->cells().activeBySequence()){const auto cell=f.workspace->cells().get(c);
            file<<(first?"":",")<<'['<<indices.at(cell.nodes[0])<<','<<indices.at(cell.nodes[1])<<','<<indices.at(cell.nodes[2])<<','<<indices.at(cell.nodes[3])<<']';first=false;}
        file<<"],\"source_cells\":[[0,4,5,1],[0,2,6,3],[0,5,6,2],[4,5,6,0]],\"target\":[1,3],\"volume\":"<<v<<"}\n";
    }
}
#include "legacy_edge_recovery_witnesses.hpp"
void allocationContracts() {
    Fixture baseline;allocations=0;failOrdinal=(std::numeric_limits<std::int64_t>::max)();
    const auto success=bnd::recoverTetEdgeIndexedDetailed(baseline.request,baseline.state,*baseline.workspace,*baseline.context,limits);
    const auto count=allocations;failOrdinal=-1;require(success.recovery.succeeded(),"allocation baseline failed");
    std::printf("one-point allocation baseline=%llu\n",static_cast<unsigned long long>(count));
    std::uint64_t unchanged=0,poisoned=0;
    const char* selected=std::getenv("PAACS_STEINER_FAULT_ORDINAL");
    const auto begin=selected?std::strtoull(selected,nullptr,10):0;
    const char* selectedEnd=std::getenv("PAACS_STEINER_FAULT_END");
    const auto end=selectedEnd?std::strtoull(selectedEnd,nullptr,10):(selected?begin+1:count);
    for(std::uint64_t ordinal=begin;ordinal<end;++ordinal) {
        if(ordinal%32==0)std::printf("fault ordinal=%llu\n",static_cast<unsigned long long>(ordinal));
        Fixture f;const auto before=f.workspace->strictWorkspaceHash();allocations=0;failOrdinal=static_cast<std::int64_t>(ordinal);
        const auto result=bnd::recoverTetEdgeIndexedDetailed(f.request,f.state,*f.workspace,*f.context,limits);failOrdinal=-1;
        if(result.recovery.succeeded())std::printf("unexpected success at allocation ordinal=%llu\n",static_cast<unsigned long long>(ordinal));
        require(!result.recovery.succeeded(),"allocation failure was silently ignored");
        if(f.workspace->strictWorkspaceHash()==before) {
            require(f.state.edgeBindings.empty()&&f.state.workspaceStamp==f.request.stamp&&f.context->currentFor(*f.workspace),
                    "pre-commit allocation changed state/context");++unchanged;
        } else {
            require(f.context->poisoned&&!f.context->currentFor(*f.workspace),"post-commit failure left a usable stale index");++poisoned;
        }
    }
    std::printf("one-point allocation ordinals=%llu rollback=%llu post-commit-poisoned=%llu\n",static_cast<unsigned long long>(count),static_cast<unsigned long long>(unchanged),static_cast<unsigned long long>(poisoned));
}
}
int main(int argc,char** argv) {
    std::setvbuf(stdout,nullptr,_IONBF,0);
    try {std::puts("plain-value witnesses");plainValueContracts();std::puts("indexed staged recovery");indexedContracts(argc>1?argv[1]:nullptr);
        abcBoundaryEdgeCavityWitness();abcNonStarOneSteinerSearchWitness();
        abcProtectedCandidateOneSteinerWitness();abcFailedEdgeOneSteinerSearchWitness();
        abcCompletedPathOneSteinerWitness();allocationContracts();
        std::puts("retained one-point recovery contracts passed");return 0;
    } catch(const std::exception& e) {failOrdinal=-1;std::fprintf(stderr,"%s\n",e.what());return 1;}
}
