// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "TetPlcRetainedRecovery.hpp"
#include "TetPlcCavityExpansion.hpp"
#include "TetPlcLocalFront.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include <bit>
#include <numeric>
namespace paacs::bnd {
namespace {
struct Refusal {TetPlcCavityReason reason;const char* diagnostic;};
void require(bool ok,TetPlcCavityReason reason,const char* diagnostic){if(!ok)throw Refusal{reason,diagnostic};}
bool same(Vec3 a,Vec3 b){return std::bit_cast<std::uint64_t>(a.x)==std::bit_cast<std::uint64_t>(b.x)&&std::bit_cast<std::uint64_t>(a.y)==std::bit_cast<std::uint64_t>(b.y)&&std::bit_cast<std::uint64_t>(a.z)==std::bit_cast<std::uint64_t>(b.z);}
std::uint64_t checkedSum(std::uint64_t a,std::uint64_t b){require(b<=UINT64_MAX-a,TetPlcCavityReason::BudgetExceeded,"storage-overflow");return a+b;}
std::uint64_t checkedMul(std::uint64_t a,std::uint64_t b){require(!b||a<=UINT64_MAX/b,TetPlcCavityReason::BudgetExceeded,"storage-overflow");return a*b;}
}
TetPlcCavityOutcome recoverTetPlcCavitiesRetaining(const mesh::Mesh& source,const mesh::Mesh& partial,
    const TetPlcCavityOptions& options,mesh::Mesh& destination,mesh::Mesh* retainedProgress) noexcept {
    TetPlcCavityOutcome out;using Reason=TetPlcCavityReason;detail::PlcFrontExpansionMetrics metrics;
    using Clock=std::chrono::steady_clock;auto phaseBegan=Clock::now();double* phase=&out.preparationSeconds;
    const auto enterPhase=[&](double& next){const auto now=Clock::now();*phase+=std::chrono::duration<double>(now-phaseBegan).count();phaseBegan=now;phase=&next;};
    try {
        require(&source!=&partial&&&destination!=&source&&&destination!=&partial,Reason::InvalidRequest,"input-destination-alias");
        require(!retainedProgress||(retainedProgress!=&source&&retainedProgress!=&partial&&retainedProgress!=&destination),Reason::InvalidRequest,"progress-input-alias");
        require(options.region.dim==3&&options.region.tag>=0&&options.maximumNodes&&options.maximumCells&&options.maximumSourceTriangles&&options.maximumWork&&options.maximumWorkingBytes&&options.maximumPatchCells>=4&&options.maximumComponentFaces>=4&&options.maximumCandidatePoints&&options.maximumPatches&&std::isfinite(options.warningQuality)&&options.warningQuality>=0&&options.warningQuality<=1,Reason::InvalidRequest,"invalid-options");
        require(options.maximumRandomCandidatePoints<=options.maximumCandidatePoints&&options.maximumRandomCandidatePoints<=4096,Reason::InvalidRequest,"random-candidate-budget");
        require(!options.patchSearch.maximumStates||(options.patchSearch.maximumStates<=65536&&options.patchSearch.maximumCandidateTets&&options.patchSearch.maximumCandidateTets<=32768&&options.patchSearch.maximumVertices>=4&&options.patchSearch.maximumVertices<=512&&options.patchSearch.maximumShells<=8),Reason::InvalidRequest,"invalid-region-search-limits");
        require(options.maximumLocalFrontTets<=4096,Reason::InvalidRequest,"local-front-limit");
        require(options.useVisibilityExpansion||options.patchSearch.maximumStates||options.maximumLocalFrontTets,Reason::InvalidRequest,"no-recovery-route");
        detail::PlcFrontExpansionLimits limits;limits.maximumRemovedCells=options.maximumExpansionCells;limits.maximumCandidatePoints=options.maximumCandidatePoints;limits.maximumRandomCandidatePoints=options.maximumRandomCandidatePoints;limits.randomSeed=options.randomSeed;limits.maximumMergedComponents=options.maximumMergedComponents;limits.maximumNewCells=options.maximumPatchCells;limits.maximumComponentFaces=options.maximumComponentFaces;limits.maximumWork=options.maximumWork;limits.work=&out.workUsed;limits.cancellation=options.cancellation;limits.tick();
        limits.patchSearch=options.patchSearch;limits.maximumNodes=options.maximumNodes;
        const auto ns=source.nodeCount(),fs=source.cellCount();const bool empty=partial.nodeCount()==0&&partial.cellCount()==0;const auto np=empty?ns:partial.nodeCount();
        require(ns>=4&&fs>=4&&np>=ns,Reason::InvalidSource,"source-or-prefix-count");
        require(std::uint64_t(np)<=options.maximumNodes&&std::uint64_t(fs)<=options.maximumSourceTriangles&&std::uint64_t(partial.cellCount())<=options.maximumCells&&std::uint64_t(np)<kInvalidTetFacetVertex,Reason::BudgetExceeded,"input-count-budget");
        const auto workspace=checkedSum(checkedMul(np,512),checkedSum(checkedMul(fs,2048),checkedMul(partial.cellCount(),4096)));
        out.workingPayloadEstimate=checkedSum(workspace,checkedSum(checkedMul(options.maximumPatchCells,4096),checkedMul(options.maximumExpansionCells,4096)));
        const bool localFront=options.maximumLocalFrontTets!=0;
        const auto searchStorage=checkedSum(options.patchSearch.maximumStates?checkedSum(65536ULL*96,checkedSum(checkedMul(options.patchSearch.maximumCandidateTets,2048),checkedMul(options.patchSearch.maximumStates,256))):0,(localFront?24ULL:16ULL)*1024*1024);
        out.workingPayloadEstimate=checkedSum(out.workingPayloadEstimate,searchStorage);
        require(out.workingPayloadEstimate<=options.maximumWorkingBytes,Reason::BudgetExceeded,"working-payload-budget");
        detail::PlcGiftWrapping g;g.useLocalCandidates=options.patchSearch.maximumStates!=0;g.useExactContactShortcut=localFront;g.cancellation=options.cancellation;g.interiorAdmission=options.interiorAdmission;g.centroidAdmission=options.centroidAdmission;g.vertices.reserve(static_cast<std::size_t>(np));
        for(Index i=0;i<np;++i){limits.tick();const auto p=empty?source.node(i):partial.node(i);require(tet_facet_detail::finite(p),Reason::InvalidPartialMesh,"nonfinite-coordinate");
            if(i<ns){require(source.nodeClass(i).dim>=0&&source.nodeClass(i).dim<=2&&source.nodeClass(i).tag>=0,Reason::InvalidSource,"source-node-classification");if(!empty)require(same(source.node(i),p)&&source.nodeClass(i)==partial.nodeClass(i),Reason::InvalidPartialMesh,"source-node-prefix-changed");}
            else {require(partial.nodeClass(i)==options.region,Reason::InvalidPartialMesh,"interior-region-classification");require(!options.interiorAdmission||options.interiorAdmission(p),Reason::InvalidPartialMesh,"interior-admission-refused-node");}
            g.vertices.push_back({p,static_cast<std::uint64_t>(i)+1});}
        std::vector<unsigned char> sourceUsed(static_cast<std::size_t>(ns),0);g.boundary.reserve(static_cast<std::size_t>(fs));
        for(Index i=0;i<fs;++i){limits.tick();const auto& c=source.cell(i);require(c.type==mesh::CellType::Tri3&&c.onEntity.dim==2&&c.onEntity.tag>=0,Reason::InvalidSource,"non-triangular-or-unclassified-source");std::array<TetFacetVertexId,3> f;const auto* ids=source.cellNodes(i);for(unsigned j=0;j<3;++j){require(ids[j]>=0&&ids[j]<ns,Reason::InvalidSource,"source-index");f[j]=static_cast<TetFacetVertexId>(ids[j]);sourceUsed[ids[j]]=1;}g.boundary.push_back(f);}
        require(std::find(sourceUsed.begin(),sourceUsed.end(),0)==sourceUsed.end(),Reason::InvalidSource,"unused-source-node");
        require(g.initialize(),Reason::InvalidSource,"source-is-not-closed-oriented");g.indexManifoldFront();
        for(const auto& [vertex,faces]:g.frontVertexFaces){limits.tick();std::map<unsigned,std::set<unsigned>> link;for(auto f:faces){std::array<unsigned,2> other;unsigned n=0;for(auto v:f)if(v!=vertex)other[n++]=v;require(n==2,Reason::InvalidSource,"source-vertex-link");link[other[0]].insert(other[1]);link[other[1]].insert(other[0]);}
            for(const auto& [v,next]:link){(void)v;require(next.size()==2,Reason::InvalidSource,"source-vertex-link-degree");}std::set<unsigned> seen;std::vector<unsigned> queue{link.begin()->first};seen.insert(queue.front());for(std::size_t j=0;j<queue.size();++j)for(auto v:link.at(queue[j]))if(seen.insert(v).second)queue.push_back(v);require(seen.size()==link.size(),Reason::InvalidSource,"disconnected-source-vertex-link");}
        // Reuse the same exact narrow predicates behind a three-dimensional
        // broad phase. A long x projection must not compare unrelated faces
        // across every y/z layer of a dense CAD surface.
        detail::PlcGrowingAabbIndex sourcePairs;
        const auto sourceBox=[&](std::size_t i){const auto b=g.boundaryBoxes[i];return detail::PlcGrowingAabbIndex::Box{b.lo,b.hi};};
        const auto chargeSource=[&](std::uint64_t n){limits.tick(n);};
        TetFacetComplexView view{g.vertices,{},{},{},{},{3,0}};
        for(std::size_t i=0;i<g.boundary.size();++i){
            const auto compatible=sourcePairs.all(g.boundary.size(),sourceBox(i),sourceBox,chargeSource,[&](std::size_t j){
                if(j<=i)return true;
                limits.tick();
                const auto a=g.boundary[i],b=g.boundary[j];
                bool bad=false;
                for(unsigned e=0;e<3;++e){bad|=tet_facet_detail::segmentCrossesOpenTriangle(g.vertices[a[e]].position,g.vertices[a[(e+1)%3]].position,g.vertices[b[0]].position,g.vertices[b[1]].position,g.vertices[b[2]].position);bad|=tet_facet_detail::segmentCrossesOpenTriangle(g.vertices[b[e]].position,g.vertices[b[(e+1)%3]].position,g.vertices[a[0]].position,g.vertices[a[1]].position,g.vertices[a[2]].position);}
                return !bad;
            });
            require(compatible,Reason::InvalidSource,"intersecting-source-triangles");
        }
        const auto origin=g.vertices.front().position;long double sourceVolume=0;for(auto f:g.boundary){limits.tick();const auto a=g.vertices[f[0]].position-origin,b=g.vertices[f[1]].position-origin,c=g.vertices[f[2]].position-origin;sourceVolume+=static_cast<long double>(a.dot(b.cross(c)))/6;}
        out.sourceVolume=static_cast<double>(sourceVolume);require(std::isfinite(out.sourceVolume)&&out.sourceVolume>0,Reason::InvalidSource,"source-volume-or-orientation");
        out.inputCells=partial.cellCount();
        for(Index i=0;i<partial.cellCount();++i){limits.tick();const auto& row=partial.cell(i);require(row.type==mesh::CellType::Tet4&&row.onEntity==options.region,Reason::InvalidPartialMesh,"partial-cell-type-or-region");TetFacetCell cell;const auto* ids=partial.cellNodes(i);std::set<Index> unique;for(unsigned j=0;j<4;++j){require(ids[j]>=0&&ids[j]<np&&unique.insert(ids[j]).second,Reason::InvalidPartialMesh,"partial-cell-indices");cell.vertices[j]=static_cast<TetFacetVertexId>(ids[j]);}cell.classification={3,0};cell.sequence=g.cells.size()+1;
            const auto v=cell.vertices;const std::array<Vec3,4> corners{g.vertices[v[0]].position,g.vertices[v[1]].position,g.vertices[v[2]].position,g.vertices[v[3]].position};
            require(predicates::orient3d(corners[0],corners[1],corners[2],corners[3])==predicates::Sign::Positive,Reason::InvalidPartialMesh,"nonpositive-partial-tet");
            for(unsigned side=0;side<4;++side){limits.tick();const auto directed=tet_facet_detail::outwardFace(cell,side);const auto at=g.front.find(detail::PlcGiftWrapping::key(directed));if(at!=g.front.end())require(tet_facet_detail::faceSense(view,at->second.directed)==tet_facet_detail::faceSense(view,directed),Reason::InvalidPartialMesh,"partial-source-or-tet-intersection");}
            // Every source vertex belongs to a source facet. If a vertex lies
            // strictly inside this tet, that facet's box overlaps the tet and
            // the exact contact test rejects it. Its supporting-plane shortcut
            // cannot separate an interior vertex. No second point scan is needed.
            const auto box=g.box(v);const auto before=g.metrics.boundaryTests;const bool sourceValid=g.boundaryAccepts(cell,box,0,true);limits.tick(g.metrics.boundaryTests-before);
            require(sourceValid,options.cancellation.stop_requested()?Reason::Cancelled:Reason::InvalidPartialMesh,"partial-source-or-tet-intersection");
            require(g.admitCentroid((corners[0]+corners[1]+corners[2]+corners[3])*.25),Reason::InvalidPartialMesh,"partial-source-or-tet-intersection");g.insert(cell);
        }
        const TetFacetComplexView partialView{g.vertices,g.cells,{},{},{},{3,0}};
        if(localFront){
            const auto cellBounds=[&](std::size_t i){const auto b=g.cellBoxes[i];return detail::PlcGrowingAabbIndex::Box{b.lo,b.hi};};
            for(std::size_t i=0;i<g.cells.size();++i)
                require(g.occupiedIndex.all(g.cells.size(),cellBounds(i),cellBounds,chargeSource,[&](std::size_t j){
                    return j<=i||!g.cellsConflict(partialView,g.cells[i],g.cells[j]);
                }),Reason::InvalidPartialMesh,"partial-source-or-tet-intersection");
        }else{
            std::vector<std::size_t> cellOrder(g.cells.size());std::iota(cellOrder.begin(),cellOrder.end(),0);
            std::sort(cellOrder.begin(),cellOrder.end(),[&](auto a,auto b){return std::pair(g.cellBoxes[a].lo.x,a)<std::pair(g.cellBoxes[b].lo.x,b);});
            for(std::size_t ai=0;ai<cellOrder.size();++ai){const auto i=cellOrder[ai];for(std::size_t aj=ai+1;aj<cellOrder.size();++aj){limits.tick();const auto j=cellOrder[aj];if(g.cellBoxes[j].lo.x>g.cellBoxes[i].hi.x)break;if(!detail::PlcGiftWrapping::overlap(g.cellBoxes[i],g.cellBoxes[j]))continue;require(!tet_facet_detail::tetImproperlyIntersects(partialView,g.cells[i],g.cells[j]),Reason::InvalidPartialMesh,"partial-source-or-tet-intersection");}}
        }
        out.initialFrontFaces=g.front.size();
        const auto retainProgress=[&](){
            if(!retainedProgress||(!metrics.closed&&!out.localFrontTets))return;
            std::vector<unsigned char> used(g.vertices.size(),0);for(auto c:g.cells)for(auto v:c.vertices)used[v]=1;
            mesh::Mesh staged;std::vector<Index> map(g.vertices.size(),kInvalidIndex);
            for(Index i=0;i<ns;++i){limits.tick();map[i]=staged.addNode(g.vertices[i].position,source.nodeClass(i));}
            for(std::size_t i=std::size_t(ns);i<g.vertices.size();++i)if(used[i]){limits.tick();map[i]=staged.addNode(g.vertices[i].position,options.region);}
            for(auto c:g.cells){limits.tick();std::array<Index,4> ids;for(unsigned j=0;j<4;++j)ids[j]=map[c.vertices[j]];staged.addCell(mesh::CellType::Tet4,ids,options.region);}
            retainedProgress->swap(staged);
        };
        if(localFront&&!g.front.empty()){
            enterPhase(out.localFrontSeconds);
            g.chargeWork=chargeSource;
            detail::advancePlcLocalFront(g,static_cast<unsigned>(std::min<std::uint64_t>(
                options.maximumLocalFrontTets,options.maximumCells-g.cells.size())),chargeSource,
                [&](bool closed){++out.localFrontTets;++metrics.added;out.localFrontClosed+=closed;});
            g.chargeWork={};
            metrics.closed+=static_cast<unsigned>(out.localFrontClosed);
            if(!g.front.empty()&&out.localFrontTets)retainProgress();
        }
        using CachedFace=detail::PlcExpansionFace;
        std::map<std::vector<CachedFace>,std::set<CachedFace>> failedRegions;
        std::size_t cachedDependencies=0,cachedFaces=0;
        const bool cacheFailures=true;
        for(unsigned pass=0;!g.front.empty()&&pass<options.maximumPatches;++pass){
            limits.tick();out.workingPayloadEstimate=std::max(out.workingPayloadEstimate,checkedSum(searchStorage,checkedSum(checkedSum(checkedMul(g.vertices.size(),512),checkedMul(fs,2048)),checkedSum(checkedMul(g.cells.size(),4096),checkedSum(checkedMul(options.maximumPatchCells,4096),checkedMul(options.maximumExpansionCells,4096))))));
            require(out.workingPayloadEstimate<=options.maximumWorkingBytes,Reason::BudgetExceeded,"expanded-working-payload-budget");
            require(g.vertices.size()<options.maximumNodes||options.patchSearch.maximumStates,Reason::BudgetExceeded,"node-budget-before-cavity");bool progress=false;
            const auto components=detail::plcFrontComponents(g);
            std::optional<detail::PlcCellOwners> owners;
            if(localFront&&(options.useVisibilityExpansion||options.patchSearch.maximumStates))
                owners=detail::plcCellOwners(g,limits);
            const auto sharedOwners=owners?&*owners:nullptr;
            if(options.useVisibilityExpansion)for(const auto& component:components){limits.tick(component.size());enterPhase(out.visibilitySeconds);
                if(g.vertices.size()<options.maximumNodes&&detail::closePlcFrontComponent(g,component,metrics,limits,sharedOwners)){progress=true;failedRegions.clear();cachedDependencies=0;cachedFaces=0;break;}}
            if(!progress&&options.patchSearch.maximumStates){enterPhase(out.regionSearchSeconds);
                for(const auto& component:components){
                    limits.tick(component.size());
                    if(cacheFailures&&failedRegions.contains(component))continue;
                    if(detail::closePlcFrontRegion(g,component,metrics,limits,sharedOwners)){
                        progress=true;
                        if(cacheFailures&&metrics.lastRegionTouched)for(auto at=failedRegions.begin();at!=failedRegions.end();){bool affected=false;for(auto v:*metrics.lastRegionTouched)if(at->second.contains(v)){affected=true;break;}if(affected){cachedDependencies-=at->second.size();cachedFaces-=at->first.size();at=failedRegions.erase(at);}else ++at;}
                        break;
                    }
                    if(cacheFailures&&metrics.lastRegionDependencies&&failedRegions.size()<4096&&cachedDependencies+metrics.lastRegionDependencies->size()<=65536&&cachedFaces+component.size()<=65536){cachedDependencies+=metrics.lastRegionDependencies->size();cachedFaces+=component.size();failedRegions.emplace(component,*metrics.lastRegionDependencies);}
                }
            }
            if(!progress){out.remainingFrontFaces=g.front.size();out.cavitiesClosed=metrics.closed;out.removedCells=metrics.removed;out.addedCells=metrics.added;
                for(Index i=0;i<fs;++i){limits.tick();if(g.front.contains(detail::PlcGiftWrapping::key(g.boundary[i])))out.unresolvedSourceTriangles.push_back(i);}
                retainProgress();
                throw Refusal{metrics.expansionBudget||metrics.regionBoundHits?Reason::BudgetExceeded:Reason::Unavailable,options.patchSearch.maximumStates?"no-admissible-region-refill":"no-admissible-expanded-kernel"};}
            require(g.cells.size()<=options.maximumCells&&g.vertices.size()<=options.maximumNodes,Reason::BudgetExceeded,"recovered-output-count-budget");
        }
        out.remainingFrontFaces=g.front.size();if(!g.front.empty())retainProgress();require(g.front.empty(),Reason::BudgetExceeded,"patch-count-budget");out.cavitiesClosed=metrics.closed;out.removedCells=metrics.removed;out.addedCells=metrics.added;
        enterPhase(out.publicationSeconds);std::vector<unsigned char> used(g.vertices.size(),0);for(auto c:g.cells)for(auto v:c.vertices)used[v]=1;
        mesh::Mesh result=source;std::vector<Index> remap(g.vertices.size(),kInvalidIndex);for(Index i=0;i<ns;++i)remap[i]=i;for(std::size_t i=static_cast<std::size_t>(ns);i<g.vertices.size();++i)if(used[i]){limits.tick();remap[i]=result.addNode(g.vertices[i].position,options.region);}
        for(auto c:g.cells){limits.tick();std::array<Index,4> ids;for(unsigned j=0;j<4;++j)ids[j]=remap[c.vertices[j]];result.addCell(mesh::CellType::Tet4,ids,options.region);}
        enterPhase(out.validationSeconds);limits.tick(result.cellCount()*4);mesh::LinearVolumeAuditLimits auditLimits;auditLimits.maximumTemporaryBytes=options.maximumWorkingBytes-out.workingPayloadEstimate;std::atomic_bool cancelled{options.cancellation.stop_requested()};std::stop_callback callback(options.cancellation,[&]{cancelled.store(true,std::memory_order_relaxed);});auditLimits.cancellation=&cancelled;const auto audit=mesh::auditLinearVolumeBoundary(result,auditLimits);require(!audit.cancelled,Reason::Cancelled,"completed-audit-cancelled");require(!audit.budgetExceeded,Reason::BudgetExceeded,"completed-audit-budget");require(!audit.allocationFailed,Reason::AllocationFailed,"completed-audit-allocation");require(audit.valid,Reason::InvalidPartialMesh,"completed-boundary-audit-failed");
        limits.tick(result.cellCount());const auto quality=algo::analyzeMixedMeshQuality(result);require(!quality.allocationFailed,Reason::AllocationFailed,"quality-allocation-failed");require(quality.complete&&quality.valid,Reason::InvalidPartialMesh,"completed-positive-volume-audit-failed");out.meshVolume=quality.totalVolume;out.minimumQuality=quality.minimumQuality;out.qualityWarning=quality.minimumQuality<options.warningQuality;
        require(std::abs(out.meshVolume-out.sourceVolume)<=std::max(1.,out.sourceVolume)*1e-10,Reason::InvalidPartialMesh,"completed-plc-volume-mismatch");limits.tick();out.outputCells=g.cells.size();destination.swap(result);out.published=true;out.destinationUnchanged=false;out.reason=Reason::None;out.diagnostic="complete";
    }catch(const detail::PlcExpansionStop& stop){out.reason=stop.budget?Reason::BudgetExceeded:Reason::Cancelled;out.diagnostic=stop.budget?"work-or-patch-budget":"cancelled";}
     catch(const Refusal& failure){out.reason=failure.reason;out.diagnostic=failure.diagnostic;}
     catch(const std::bad_alloc&){out.reason=Reason::AllocationFailed;out.diagnostic="allocation-failed";}
     catch(...){out.reason=Reason::InternalError;out.diagnostic="internal-error";}
    *phase+=std::chrono::duration<double>(Clock::now()-phaseBegan).count();
    out.regionSearchAttempts=metrics.regionAttempts;out.regionSearchStates=metrics.regionStates;out.regionSearchBacktracks=metrics.regionBacktracks;out.regionSearchCandidates=metrics.regionCandidates;out.regionSearchPairCacheHits=metrics.regionPairCacheHits;out.regionSearchPatches=metrics.regionPatches;out.regionSearchBoundHits=metrics.regionBoundHits;
    out.regionVertexBoundHits=metrics.regionVertexBoundHits;out.regionStateBoundHits=metrics.regionStateBoundHits;out.largestRegionVertices=metrics.largestRegionVertices;out.largestRegionFaces=metrics.largestRegionFaces;
    out.candidatePointsTested=metrics.points;out.randomCandidatePointsProposed=metrics.randomProposed;out.patchesCommitted=metrics.closed;out.cavitiesClosed=metrics.closed+metrics.merged;out.mergedComponents=metrics.merged;out.sourceBlockedCandidates=metrics.sourceBlocked;out.otherComponentBlockedCandidates=metrics.otherHoleBlocked;
    try{if(!out.published&&metrics.blockedSources&&!metrics.blockedSources->empty()&&!options.cancellation.stop_requested()){std::vector<std::pair<std::uint64_t,Index>> ranked;for(Index i=0;i<source.cellCount()&&out.workUsed<options.maximumWork;++i){++out.workUsed;const auto* p=source.cellNodes(i);const auto found=metrics.blockedSources->find(detail::PlcGiftWrapping::key({static_cast<unsigned>(p[0]),static_cast<unsigned>(p[1]),static_cast<unsigned>(p[2])}));if(found!=metrics.blockedSources->end())ranked.emplace_back(found->second,i);}std::sort(ranked.begin(),ranked.end(),[](auto a,auto b){return a.first!=b.first?a.first>b.first:a.second<b.second;});for(auto [count,index]:ranked){(void)count;out.blockedSourceTriangles.push_back(index);}}}
    catch(const std::bad_alloc&){out.reason=Reason::AllocationFailed;out.diagnostic="diagnostic-allocation-failed";}catch(...){out.reason=Reason::InternalError;out.diagnostic="diagnostic-internal-error";}return out;
}
TetPlcCavityOutcome recoverTetPlcCavities(const mesh::Mesh& source,const mesh::Mesh& partial,
    const TetPlcCavityOptions& options,mesh::Mesh& destination) noexcept {
    return recoverTetPlcCavitiesRetaining(source,partial,options,destination,nullptr);
}
} // namespace paacs::bnd
