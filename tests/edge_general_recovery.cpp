// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Small exact fixture captured from the first refused ABC5454 edge cavity.
// No downloaded STEP data or OCCT dependency is needed for this regression.
#include "paacs/bnd/TetFacetSetRecovery.hpp"
#include "../src/bnd/TetEdgeManifoldClosure.hpp"
#include "../src/bnd/TetEdgeGeneralCavityCadPolicy.hpp"
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <stdexcept>
#include "allocation_nothrow.hpp"

static std::atomic_bool failNextAllocation{false};
void* operator new(std::size_t n) {
    if(failNextAllocation.exchange(false))throw std::bad_alloc{};
    if(void* p=std::malloc(n?n:1))return p;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t n){return ::operator new(n);}
void operator delete(void* p) noexcept {std::free(p);}
void operator delete[](void* p) noexcept {std::free(p);}
void operator delete(void* p,std::size_t) noexcept {std::free(p);}
void operator delete[](void* p,std::size_t) noexcept {std::free(p);}

using namespace paacs;
namespace {
void require(bool condition,const char* text) {
    if(!condition)throw std::runtime_error(text);
}
constexpr bnd::TetEdgeGeneralCavityLimits complete{64,16,5000,10000,500000};
struct Fixture {
    std::vector<topo::NodeHandle> nodes;
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::unique_ptr<bnd::TetRecoveryContext> context;
    bnd::TetFacetSetState state;
    bnd::TetEdgeRecoveryRequest request;
    topo::CellHandle outside{};
    topo::EdgeHandle protectedEdge{};
    topo::FaceHandle protectedFace{};
    Fixture(bool protectBoundary=false,bool protectInterior=false) {
        topo::NodePool nodePool;
        for(const Vec3& p:std::array<Vec3,12>{{
            {-6,-13.899999999999999,3}, {-5,-10.300000000000001,3},
            {-2.2000000000000002,-16.899999999999999,1.8},
            {-2.2000000000000002,-16.899999999999999,6.2000000000000002},
            {2.2000000000000002,-16.899999999999999,1.8},
            {2.2000000000000002,-16.899999999999999,6.2000000000000002},
            {5,-10.300000000000001,3},{6,-13.899999999999999,3},
            {50,50,50},{51,50,50},{50,51,50},{50,50,51}}})
            nodes.push_back(nodePool.add(p,{3,9}));
        topo::CellPool cellPool;
        for(auto ids:std::array<std::array<std::size_t,4>,6>{{
            {{2,0,3,1}},{{2,3,5,1}},{{2,1,5,6}},{{2,6,5,4}},{{4,6,5,7}},{{8,9,10,11}}}}) {
            std::array<topo::NodeHandle,4> handles{};
            for(std::size_t i=0;i<4;++i)handles[i]=nodes[ids[i]];
            if(predicates::orient3d(nodePool.get(handles[0]).coord,nodePool.get(handles[1]).coord,
                nodePool.get(handles[2]).coord,nodePool.get(handles[3]).coord)==predicates::Sign::Negative)
                std::swap(handles[0],handles[1]);
            outside=cellPool.add(topo::ElementType::Tet4,handles,{3,9},nodePool);
        }
        topo::EdgePool edgePool;
        for(const auto cell:cellPool.activeBySequence())
            (void)edgePool.insertCellEdges(nodePool,cellPool,cell);
        if(protectBoundary) {
            const auto edge=edgePool.find(nodePool,nodes[0],nodes[1]);
            require(edge.has_value(),"missing boundary protection fixture");
            protectedEdge=*edge;
            edgePool.setClassification(*edge,{1,77});
            edgePool.setRequiredConstraint(*edge,true);
        }
        topo::FacePool facePool;
        for(const auto cell:cellPool.activeBySequence())
            (void)facePool.insertCellFaces(nodePool,cellPool,edgePool,cell);
        if(protectBoundary||protectInterior) {
            const auto ids=protectInterior?std::array<std::size_t,3>{2,3,1}:
                                           std::array<std::size_t,3>{0,3,1};
            const auto face=facePool.find(nodePool,
                std::array<topo::NodeHandle,3>{nodes[ids[0]],nodes[ids[1]],nodes[ids[2]]});
            require(face.has_value(),"missing protected face fixture");
            protectedFace=*face;facePool.setClassification(*face,{2,78});
        }
        workspace=std::make_unique<topo::TopologyWorkspace>(
            std::move(nodePool),std::move(cellPool),std::move(edgePool),std::move(facePool));
        auto built=bnd::buildTetRecoveryContext(*workspace,{3,9});
        require(built.built(),"fixture context audit failed");context=std::move(built.context);
        state.workspaceStamp=topo::WorkspaceStamp::of(*workspace);state.artifactStamp={5454,1};
        request.stamp=state.workspaceStamp;request.artifactStamp=state.artifactStamp;
        request.key=bnd::VolumeArtifactEdgeKey::of(0,7);
        request.vertices={nodes[0],nodes[7]};request.classification={1,193};request.region={3,9};
    }
};
template<class F> double volume(const F& f) {
    double sum=0;
    for(const auto cell:f.workspace->cells().activeBySequence()) {
        const auto c=f.workspace->cells().get(cell);
        const auto p=[&](std::size_t i){return f.workspace->nodes().get(c.nodes[i]).coord;};
        sum+=std::abs((p(1)-p(0)).dot((p(2)-p(0)).cross(p(3)-p(0))))/6;
    }
    return sum;
}
template<class F> void unchanged(const F& f,std::uint64_t before) {
    require(f.workspace->strictWorkspaceHash()==before,"failed recovery changed topology identity/content");
    require(f.state.workspaceStamp==f.request.stamp&&f.state.edgeBindings.empty(),"failed recovery changed bindings");
    require(f.context->currentFor(*f.workspace),"failed recovery poisoned current context");
}
struct ExpansionFixture {
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::unique_ptr<bnd::TetRecoveryContext> context;
    bnd::TetFacetSetState state;
    bnd::TetEdgeRecoveryRequest request;
    ExpansionFixture(bool pair,bool protectBridge=false) {
        const std::vector<Vec3> points=pair?std::vector<Vec3>{
            {5,-10.300000000000001,3},{5,-1.7763568394002505e-15,3},
            {7.6384650723747658,1.0387160589922075,0},
            {12.574999999999989,-10.300000000000001,.40000000000000002},
            {12.575000000000014,-10.300000000000001,3},
            {14.38440328417316,3.7309177010787886,.40000000000000002},
            {15.765596715826838,-3.6534855830943727,0}}:std::vector<Vec3>{
            {-25.149999999999999,-7.2999999999999998,10.800000000000001},
            {-23.150000000000002,-10.300000000000001,8.4000000000000004},
            {-23.150000000000002,-10.300000000000001,10.200000000000001},
            {-23.150000000000002,-7.2999999999999998,8.4000000000000004},
            {-23.150000000000002,-7.2999999999999998,10.200000000000001},
            {-20.150000000000002,-10.300000000000001,8.4000000000000004},
            {-20.150000000000002,-10.300000000000001,10.200000000000001},
            {-20.150000000000002,-7.2999999999999998,8.4000000000000004}};
        const std::vector<std::array<std::size_t,4>> connectivity=pair?
            std::vector<std::array<std::size_t,4>>{{4,1,2,0},{3,4,2,0},{1,2,6,4},{2,6,5,1},{3,4,6,2}}:
            std::vector<std::array<std::size_t,4>>{{2,0,4,1},{3,1,4,0},{1,6,4,2},{1,3,4,7},{1,5,4,6},{1,7,4,5},{7,6,4,5}};
        topo::NodePool np;std::vector<topo::NodeHandle> nodes;
        for(const auto& p:points)nodes.push_back(np.add(p,{3,9}));
        topo::CellPool cp;
        for(const auto& ids:connectivity) {
            std::array<topo::NodeHandle,4> handles{};
            for(std::size_t i=0;i<4;++i)handles[i]=nodes[ids[i]];
            if(predicates::orient3d(np.get(handles[0]).coord,np.get(handles[1]).coord,
                np.get(handles[2]).coord,np.get(handles[3]).coord)==predicates::Sign::Negative)
                std::swap(handles[0],handles[1]);
            (void)cp.add(topo::ElementType::Tet4,handles,{3,9},np);
        }
        topo::EdgePool ep;for(const auto cell:cp.activeBySequence())(void)ep.insertCellEdges(np,cp,cell);
        topo::FacePool fp;for(const auto cell:cp.activeBySequence())(void)fp.insertCellFaces(np,cp,ep,cell);
        if(protectBridge) {
            const auto ids=pair?std::array<std::size_t,3>{1,2,4}:std::array<std::size_t,3>{4,5,7};
            const auto face=fp.find(np,std::array<topo::NodeHandle,3>{nodes[ids[0]],nodes[ids[1]],nodes[ids[2]]});
            require(face.has_value(),"missing exact adjacency bridge");fp.setClassification(*face,{2,78});
        }
        workspace=std::make_unique<topo::TopologyWorkspace>(std::move(np),std::move(cp),std::move(ep),std::move(fp));
        auto built=bnd::buildTetRecoveryContext(*workspace,{3,9});require(built.built(),"expansion context audit failed");
        context=std::move(built.context);state.workspaceStamp=topo::WorkspaceStamp::of(*workspace);state.artifactStamp={5454,1};
        request.stamp=state.workspaceStamp;request.artifactStamp=state.artifactStamp;request.key=bnd::VolumeArtifactEdgeKey::of(1,6);
        request.vertices=pair?std::array<topo::NodeHandle,2>{nodes[3],nodes[5]}:
                              std::array<topo::NodeHandle,2>{nodes[2],nodes[3]};
        request.classification={1,193};request.region={3,9};
    }
};
void exactAdjacentExpansion() {
    for(bool pair:{false,true}) {
        ExpansionFixture f(pair);const auto before=f.workspace->strictWorkspaceHash();const auto v=volume(f);
        const auto noExpansion=bnd::recoverTetEdgeIndexedDetailed(f.request,f.state,*f.workspace,*f.context,complete);
        require(noExpansion.recovery.reason==bnd::TetEdgeRecoveryReason::NonCycleLink,"captured no-solution cavity changed");
        unchanged(f,before);
        if(pair) {
            const auto shortRun=bnd::recoverTetEdgeIndexedDetailed(f.request,f.state,*f.workspace,*f.context,complete,nullptr,nullptr,1);
            require(shortRun.recovery.reason==bnd::TetEdgeRecoveryReason::ResourceBudgetReached,"adjacent attempt cap not distinct");
            unchanged(f,before);
        }
        bnd::TetEdgeGeneralRecoveryDiagnostics evidence;
        std::array<topo::CellHandle,16> frontier{};std::size_t used=0;bool witness=true;
        const auto assessed=bnd::assessTetEdgeIndexedFrontier(f.request,f.state,*f.workspace,*f.context,
            frontier,used,witness,complete,&evidence,nullptr,80);
        require(assessed.succeeded()&&witness&&used==(pair?5:7),"expanded frontier omitted retired source cells");
        unchanged(f,before);
        const auto result=bnd::recoverTetEdgeIndexedDetailed(f.request,f.state,*f.workspace,*f.context,complete,&evidence,nullptr,80);
        require(result.recovery.succeeded()&&evidence.adjacentCellsAdded==(pair?2:1),"captured bounded adjacent recovery failed");
        require(result.retiredCells.size()==(pair?5:7)&&result.createdCells.size()==(pair?6:7),"expanded exact replacement counts changed");
        require(evidence.metrics.workUnits<=complete.maxWorkUnits&&evidence.adjacentExpansionAttempts<=80,"expansion budgets exceeded");
        require(std::abs(volume(f)-v)<=1e-12*std::max(1.,v),"expanded cavity changed volume");
        require(bnd::buildTetRecoveryContext(*f.workspace,{3,9}).built(),"expanded global topology audit failed");
        std::printf("ABC5454 %s-adjacent cavity: attempts=%llu work=%llu, exact global audit passed\n",pair?"two":"one",
            static_cast<unsigned long long>(evidence.adjacentExpansionAttempts),static_cast<unsigned long long>(evidence.metrics.workUnits));
    }
    ExpansionFixture blocked(false,true);const auto before=blocked.workspace->strictWorkspaceHash();
    const auto refusal=bnd::recoverTetEdgeIndexedDetailed(blocked.request,blocked.state,*blocked.workspace,*blocked.context,complete,nullptr,nullptr,80);
    require(!refusal.recovery.succeeded(),"expansion crossed a protected dual face");unchanged(blocked,before);
}
void successAndDefault() {
    Fixture f(true);const auto before=f.workspace->strictWorkspaceHash();const auto v=volume(f);
    const auto old=bnd::recoverTetEdgeIndexedDetailed(f.request,f.state,*f.workspace,*f.context);
    require(old.recovery.reason==bnd::TetEdgeRecoveryReason::NonCycleLink,"default route changed");unchanged(f,before);
    bnd::TetEdgeGeneralRecoveryDiagnostics evidence;
    const auto assessed=bnd::assessTetEdgeIndexedDetailed(
        f.request,f.state,*f.workspace,*f.context,complete,&evidence);
    require(assessed.recovery.succeeded()&&evidence.attempted&&evidence.status==bnd::TetEdgeGeneralCavityStatus::Solved,
            "bounded assessment did not solve actual strip");
    require(assessed.retiredCells.size()==5,"assessment did not witness the exact selected cells");unchanged(f,before);
    const auto outside=f.workspace->cells().get(f.outside).nodes;
    const auto result=bnd::recoverTetEdgeIndexedDetailed(
        f.request,f.state,*f.workspace,*f.context,complete,&evidence);
    require(result.recovery.status==bnd::TetEdgeRecoveryStatus::Recovered,"bounded actual strip not recovered");
    require(result.retiredCells.size()==5&&result.createdCells.size()==9,"unexpected exact replacement counts");
    require(f.workspace->nodes().activeBySequence().size()==12,"existing-node recovery changed nodes");
    require(f.workspace->cells().isCurrent(f.outside)&&f.workspace->cells().get(f.outside).nodes==outside,"outside cell changed");
    require(f.workspace->edges().isCurrent(f.protectedEdge)&&f.workspace->edges().get(f.protectedEdge).requiredConstraint,
            "protected boundary edge identity/required bit lost");
    require(f.workspace->faces().isCurrent(f.protectedFace)&&f.workspace->faces().get(f.protectedFace).cls.tag==78,
            "protected boundary face identity/classification lost");
    require(f.context->currentFor(*f.workspace),"committed context stale");
    require(bnd::buildTetRecoveryContext(*f.workspace,{3,9}).built(),"committed global topology audit failed");
    require(std::abs(volume(f)-v)<=1e-12*std::max(1.,v),"replacement volume differs");
    const auto binding=f.state.edgeBindings.find(f.request.key);
    require(binding!=f.state.edgeBindings.end()&&binding->second.edge==result.recovery.edge,"missing exact persistent binding");
    const auto edge=f.workspace->edges().get(result.recovery.edge);
    require(edge.requiredConstraint&&edge.cls.dim==1&&edge.cls.tag==193,"target metadata not atomic");
    std::printf("ABC5454 strip: 5 -> 9 Tets; states=%llu work=%llu; outside/protected identities retained\n",
        static_cast<unsigned long long>(evidence.metrics.searchStates),
        static_cast<unsigned long long>(evidence.metrics.workUnits));
}
void budgetsAndProtection() {
    for(int mode=0;mode<5;++mode) {
        Fixture f;auto limits=complete;const auto before=f.workspace->strictWorkspaceHash();
        if(mode==0)limits.maxCavityCells=4;
        if(mode==1)limits.maxCavityVertices=7;
        if(mode==2)limits.maxEnumeratedTets=69;
        if(mode==3)limits.maxSearchStates=1;
        if(mode==4)limits.maxWorkUnits=1;
        bnd::TetEdgeGeneralRecoveryDiagnostics evidence;
        const auto result=bnd::recoverTetEdgeIndexedDetailed(f.request,f.state,*f.workspace,*f.context,limits,&evidence);
        require(result.recovery.reason==bnd::TetEdgeRecoveryReason::ResourceBudgetReached,"hard limit did not refuse distinctly");
        require(evidence.metrics.workUnits<=limits.maxWorkUnits,"solver exceeded checked work budget");unchanged(f,before);
    }
    Fixture invalid;auto limits=complete;limits.maxWorkUnits=0;
    const auto before=invalid.workspace->strictWorkspaceHash();
    require(bnd::recoverTetEdgeIndexedDetailed(invalid.request,invalid.state,*invalid.workspace,*invalid.context,limits)
        .recovery.status==bnd::TetEdgeRecoveryStatus::Rejected,"nonpositive limits accepted");unchanged(invalid,before);
    Fixture blocked(false,true);const auto blockedBefore=blocked.workspace->strictWorkspaceHash();
    const auto result=bnd::recoverTetEdgeIndexedDetailed(blocked.request,blocked.state,*blocked.workspace,*blocked.context,complete);
    require(!result.recovery.succeeded(),"crossed protected interior face accepted");unchanged(blocked,blockedBefore);
}
void manifoldPreparation() {
  for(unsigned depth:{1U,2U}) {
    topo::NodePool np;std::vector<topo::NodeHandle> nodes;
    const auto index=[](unsigned x,unsigned y,unsigned z){return z*16+y*4+x;};
    for(unsigned z=0;z<=depth;++z)for(unsigned y=0;y<4;++y)for(unsigned x=0;x<4;++x)
        nodes.push_back(np.add({double(x),double(y),double(z)},{3,9}));
    topo::CellPool cp;std::set<topo::CellHandle> selected;
    const std::set<std::array<unsigned,3>> columns=depth==1?
        std::set<std::array<unsigned,3>>{{1,0,0},{0,0,0},{0,1,0},{0,2,0},{1,2,0},{2,2,0},{2,1,0}}:
        std::set<std::array<unsigned,3>>{{1,0,0},{0,0,0},{0,1,0},{0,2,0},{1,2,0},{2,2,0},{2,2,1},{2,1,1}};
    for(unsigned z=0;z<depth;++z)for(unsigned y=0;y<3;++y)for(unsigned x=0;x<3;++x){
        const std::array<unsigned,8> cube{index(x,y,z),index(x+1,y,z),index(x+1,y+1,z),index(x,y+1,z),
            index(x,y,z+1),index(x+1,y,z+1),index(x+1,y+1,z+1),index(x,y+1,z+1)};
        for(const auto& t:std::array<std::array<unsigned,4>,6>{{{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}}}){
            std::array<topo::NodeHandle,4> handles{nodes[cube[t[0]]],nodes[cube[t[1]]],nodes[cube[t[2]]],nodes[cube[t[3]]]};
            if(predicates::orient3d(np.get(handles[0]).coord,np.get(handles[1]).coord,np.get(handles[2]).coord,np.get(handles[3]).coord)==predicates::Sign::Negative)std::swap(handles[0],handles[1]);
            const auto cell=cp.add(topo::ElementType::Tet4,handles,{3,9},np);
            if(columns.contains({x,y,z}))selected.insert(cell);
        }
    }
    topo::EdgePool ep;for(const auto cell:cp.activeBySequence())(void)ep.insertCellEdges(np,cp,cell);
    topo::FacePool fp;for(const auto cell:cp.activeBySequence())(void)fp.insertCellFaces(np,cp,ep,cell);
    topo::TopologyWorkspace workspace(std::move(np),std::move(cp),std::move(ep),std::move(fp));
    auto built=bnd::buildTetRecoveryContext(workspace,{3,9});require(built.built(),"grid context");
    auto& context=*built.context;const auto before=workspace.strictWorkspaceHash();
    const auto audit=[&](const std::set<topo::CellHandle>& patch){
        std::vector<bnd::TetFacetCell> cells;
        for(const auto cell:patch){const auto& row=context.cells.at(cell);cells.push_back({row.nodes,row.sequence,row.classification});}
        const bnd::TetFacetComplexView view{context.plainVertices,cells,{},{},{},{3,9}};
        std::vector<std::size_t> all(cells.size());for(std::size_t i=0;i<all.size();++i)all[i]=i;
        return bnd::auditTetFacetBall(view,all);
    };
    const auto originalAudit=audit(selected);
    require(!originalAudit.valid&&originalAudit.reason==(depth==1?bnd::TetFacetGeometryReason::BoundaryNotClosed:bnd::TetFacetGeometryReason::InvalidVertexLink),"analytic edge/vertex pinch differs");
    const auto closed=bnd::detail::prepareTetEdgeManifoldClosure(selected,context,{3,9},128,64,100000);
    require(closed.status==bnd::detail::ManifoldClosureStatus::Closed&&closed.passes>0&&closed.cells.size()>selected.size(),"pinched edge/vertex boundary not closed");
    require(audit(closed.cells).valid,"prepared boundary did not pass independent exact ball audit");
    const auto replay=bnd::detail::prepareTetEdgeManifoldClosure(selected,context,{3,9},128,64,100000);
    require(replay.cells==closed.cells&&replay.work==closed.work&&replay.passes==closed.passes,"closure replay differs");
    for(unsigned mode=0;mode<3;++mode){
        const auto stopped=bnd::detail::prepareTetEdgeManifoldClosure(selected,context,{3,9},mode==0?selected.size():128,64,mode==1?1:100000,mode==2?0:32);
        require(stopped.status==bnd::detail::ManifoldClosureStatus::BudgetReached&&stopped.cells.empty(),"bounded closure published partial proposal");
        require(stopped.work<=(mode==1?1U:100000U),"closure exceeded work budget");
    }
    require(workspace.strictWorkspaceHash()==before,"pure preparation changed topology");
    // Supply exact existing-face identities to the pure preparation input.
    // Metadata identity preservation during a commit is tested separately.
    for(const auto& [key,row]:context.faces){if(row.owners.size()!=2)continue;
        if(selected.contains(row.owners[0].cell)!=selected.contains(row.owners[1].cell))
            context.protectedFaces.push_back({bnd::canonicalTetFacetFace(key.nodes),static_cast<std::uint64_t>(row.handle.slot)+1,{2,78}});}
    const auto stopped=bnd::detail::prepareTetEdgeManifoldClosure(selected,context,{3,9},128,64,100000);
    require(stopped.status==bnd::detail::ManifoldClosureStatus::Unavailable&&stopped.cells.empty(),"closure crossed protected boundary");
    require(workspace.strictWorkspaceHash()==before,"protected refusal changed topology");
    std::printf("analytic %s closure: %zu -> %zu cells, %llu passes; exact ball, protection and budgets passed\n",depth==1?"edge":"vertex",selected.size(),closed.cells.size(),static_cast<unsigned long long>(closed.passes));
  }
}
void cadRequiredEdgeSearch() {
    require(!bnd::detail::cadGeneralPruningEnabled(),"CAD policy leaked before test");
    Fixture f(true);const auto v=volume(f);
    {
        bnd::detail::ScopedCadGeneralPruning policy(true);
        bnd::TetEdgeGeneralRecoveryDiagnostics evidence;
        const auto result=bnd::recoverTetEdgeIndexedDetailed(f.request,f.state,*f.workspace,*f.context,complete,&evidence);
        require(result.recovery.succeeded(),"required-edge-first search refused valid strip");
        require(evidence.metrics.workUnits<=complete.maxWorkUnits,"CAD search work exceeded");
        require(f.workspace->cells().isCurrent(f.outside)&&f.workspace->faces().isCurrent(f.protectedFace)&&f.workspace->edges().isCurrent(f.protectedEdge),"CAD search lost outside or protected identity");
        require(std::abs(volume(f)-v)<=1e-12*std::max(1.,v),"CAD search changed volume");
        require(bnd::buildTetRecoveryContext(*f.workspace,{3,9}).built(),"CAD search global audit failed");
    }
    require(!bnd::detail::cadGeneralPruningEnabled(),"CAD policy scope did not restore");
}
#if defined(PAACS_TOPO_TESTING)
bnd::TetEdgePreparedGuardDecisionForTesting defer(
    const bnd::TetEdgeRecoveryRequest&,topo::WorkspaceStamp,topo::EdgeHandle,
    const topo::TetEdgeCavityPreview& preview,const bnd::TetRecoveryContext&,void*) noexcept {
    return preview.addedConnectivity.size()==9?bnd::TetEdgePreparedGuardDecisionForTesting::Defer:
                                              bnd::TetEdgePreparedGuardDecisionForTesting::Allow;
}
bnd::TetEdgePreparedGuardDecisionForTesting allocationAfterPrepare(
    const bnd::TetEdgeRecoveryRequest&,topo::WorkspaceStamp,topo::EdgeHandle,
    const topo::TetEdgeCavityPreview&,const bnd::TetRecoveryContext&,void*) noexcept {
    failNextAllocation.store(true);return bnd::TetEdgePreparedGuardDecisionForTesting::Allow;
}
void preparedRollback() {
    for(bool allocation:{false,true}) {
        Fixture f;const auto before=f.workspace->strictWorkspaceHash();
        const auto result=bnd::recoverTetEdgeIndexedPreparedGuardedForTesting(
            f.request,f.state,*f.workspace,*f.context,allocation?allocationAfterPrepare:defer,nullptr,complete);
        failNextAllocation.store(false);
        require(result.recovery.reason==(allocation?bnd::TetEdgeRecoveryReason::AllocationFailed:
                                                   bnd::TetEdgeRecoveryReason::PendingConstraint),"prepared stop taxonomy changed");
        unchanged(f,before);
    }
}
#endif
}
int main() {
    try {
        successAndDefault();budgetsAndProtection();exactAdjacentExpansion();manifoldPreparation();cadRequiredEdgeSearch();
#if defined(PAACS_TOPO_TESTING)
        preparedRollback();
#endif
        std::puts("bounded general edge recovery contracts passed");return 0;
    } catch(const std::exception& e) {failNextAllocation.store(false);std::fprintf(stderr,"%s\n",e.what());return 1;}
}
