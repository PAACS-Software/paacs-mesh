// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/bnd/TetEdgeGeneralCavity.hpp"
#include "TetEdgeGeneralCavityCadPolicy.hpp"
#include "PlcSpatialIndex.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <map>
#include <new>
#include <numeric>
#include <optional>
#include <queue>
#include <set>

namespace paacs::bnd {
namespace {
thread_local bool cadPruning=false;

using tet_facet_detail::FaceUse;
using tet_facet_detail::FaceUseMap;
using ContactBox=detail::PlcGrowingAabbIndex::Box;

template<class NodeIds>
ContactBox contactBounds(std::span<const TetFacetVertex> vertices,const NodeIds& ids) {
    const auto first=vertices[ids[0]].position;
    ContactBox box{first,first};
    for(auto id:ids){
        const auto p=vertices[id].position;
        box.lo={std::min(box.lo.x,p.x),std::min(box.lo.y,p.y),std::min(box.lo.z,p.z)};
        box.hi={std::max(box.hi.x,p.x),std::max(box.hi.y,p.y),std::max(box.hi.z,p.z)};
    }
    return box;
}

struct CavityContactFace {
    std::array<TetFacetVertexId,3> vertices;
    ContactBox box;
    bool valid=false;
};

bool cavityCellsConflict(const TetFacetComplexView& view,const TetFacetCell& a,const TetFacetCell& b) {
    if(tet_facet_detail::separatedByTetFace(view,a,b.vertices,[]{})||
       tet_facet_detail::separatedByTetFace(view,b,a.vertices,[]{}))return false;
    return tet_facet_detail::tetImproperlyIntersects(view,a,b);
}

struct LedgerRow {
    std::uint8_t count = 0;
    int firstSense = 0;
};

struct UndoRow {
    std::size_t faceIndex = 0;
    LedgerRow old{};
};

struct SearchFrame {
    const std::vector<std::size_t>* options = nullptr;
    std::size_t next = 0;
    std::size_t appliedCandidate = (std::numeric_limits<std::size_t>::max)();
    std::array<UndoRow,4> undo{};
    std::size_t undoCount = 0;
};

bool checkedChoose4(std::uint64_t n,std::uint64_t& out) noexcept {
    if(n<4) { out=0; return true; }
    std::array<std::uint64_t,4> terms{{n,n-1,n-2,n-3}};
    std::uint64_t divisor=24;
    for(auto& term:terms) {
        const auto g=std::gcd(term,divisor);
        term/=g;
        divisor/=g;
    }
    out=1;
    for(const auto term:terms) {
        if(term!=0&&out>(std::numeric_limits<std::uint64_t>::max)()/term)
            return false;
        out*=term;
    }
    return divisor==1;
}

bool consume(TetEdgeGeneralCavityMetrics& metrics,
             const TetEdgeGeneralCavityLimits& limits,
             std::uint64_t amount=1) noexcept {
    if(amount>limits.maxWorkUnits-metrics.workUnits)return false;
    metrics.workUnits+=amount;
    return true;
}

std::array<TetFacetVertexId,4> keyOf(const TetFacetCell& cell) noexcept {
    auto key=cell.vertices;
    std::sort(key.begin(),key.end());
    return key;
}

bool sameClass(const mesh::Classification& a,
               const mesh::Classification& b) noexcept {
    return a.dim==b.dim&&a.tag==b.tag;
}

bool canonicalInput(const TetEdgeGeneralCavityInput& input) noexcept {
    if(input.vertices.size()<4||input.removedCells.empty()||
       input.vertices.size()>kInvalidTetFacetVertex||
       input.target.lo>=input.vertices.size()||
       input.target.hi>=input.vertices.size()||
       input.target.lo>=input.target.hi||input.region.dim!=3||
       input.region.tag<0)return false;
    for(std::size_t i=0;i<input.vertices.size();++i) {
        if(!tet_facet_detail::finite(input.vertices[i].position))return false;
        if(i!=0&&input.vertices[i-1].sequence>=input.vertices[i].sequence)
            return false;
    }
    return true;
}

TetFacetComplexView viewOf(const TetEdgeGeneralCavityInput& input,
                           std::span<const TetFacetCell> cells) noexcept {
    return {input.vertices,cells,{},input.protectedEdges,input.protectedFaces,
            input.region};
}

bool targetCycle(const TetEdgeGeneralCavityInput& input,
                 std::span<const TetFacetCell> cells,
                 TetEdgeGeneralCavityMetrics& metrics,
                 const TetEdgeGeneralCavityLimits& limits) {
    std::map<TetFacetVertexId,std::set<TetFacetVertexId>> link;
    std::size_t star=0;
    for(const auto& cell:cells) {
        if(!consume(metrics,limits))return false;
        ++metrics.validatorVisits;
        if(std::find(cell.vertices.begin(),cell.vertices.end(),input.target.lo)==
               cell.vertices.end()||
           std::find(cell.vertices.begin(),cell.vertices.end(),input.target.hi)==
               cell.vertices.end())continue;
        std::array<TetFacetVertexId,2> opposite{};
        std::size_t count=0;
        for(const auto vertex:cell.vertices)
            if(vertex!=input.target.lo&&vertex!=input.target.hi)
                opposite[count++]=vertex;
        if(count!=2||opposite[0]==opposite[1])return false;
        link[opposite[0]].insert(opposite[1]);
        link[opposite[1]].insert(opposite[0]);
        ++star;
    }
    if(star<3||link.size()<3)return false;
    for(const auto& [vertex,neighbours]:link) {
        (void)vertex;
        if(neighbours.size()!=2)return false;
    }
    std::set<TetFacetVertexId> seen;
    std::queue<TetFacetVertexId> pending;
    pending.push(link.begin()->first);
    seen.insert(link.begin()->first);
    while(!pending.empty()) {
        const auto vertex=pending.front();
        pending.pop();
        for(const auto neighbour:link[vertex])
            if(seen.insert(neighbour).second)pending.push(neighbour);
    }
    return seen.size()==link.size();
}

bool completeCandidate(
    const TetEdgeGeneralCavityInput& input,
    const std::map<TetFacetFaceKey,int>& expectedBoundary,
    std::span<const TetFacetCell> candidate,
    TetEdgeGeneralCavityMetrics& metrics,
    const TetEdgeGeneralCavityLimits& limits,
    TetFacetGeometryReason& reason) {
    if(candidate.empty()) { reason=TetFacetGeometryReason::InvalidInput; return false; }
    const auto view=viewOf(input,candidate);
    std::vector<std::size_t> all(candidate.size());
    for(std::size_t i=0;i<all.size();++i)all[i]=i;
    if(!consume(metrics,limits,candidate.size()))return false;
    metrics.validatorVisits+=candidate.size();
    const auto ball=auditTetFacetBall(view,all);
    if(!ball.valid) { reason=ball.reason; return false; }
    const auto faces=tet_facet_detail::candidateFaceMap(view,candidate);
    if(tet_facet_detail::boundarySenses(faces)!=expectedBoundary) {
        reason=TetFacetGeometryReason::BoundaryMismatch;
        return false;
    }
    std::vector<unsigned char> used(input.vertices.size(),0);
    std::set<TetFacetEdgeKey> edges;
    for(const auto& cell:candidate) {
        if(!sameClass(cell.classification,input.region)) {
            reason=TetFacetGeometryReason::RegionMismatch;
            return false;
        }
        for(const auto vertex:cell.vertices)used[vertex]=1;
        for(const auto& local:tet_facet_detail::kTetEdges)
            edges.insert(tet_facet_detail::edgeKey(
                cell.vertices[static_cast<std::size_t>(local[0])],
                cell.vertices[static_cast<std::size_t>(local[1])]));
    }
    if(std::find(used.begin(),used.end(),0)!=used.end()) {
        reason=TetFacetGeometryReason::OrphanNode;
        return false;
    }
    for(const auto edge:input.protectedEdges)
        if(edges.find(edge)==edges.end()) {
            reason=TetFacetGeometryReason::CrossedProtectedEdge;
            return false;
        }
    for(const auto& face:input.protectedFaces)
        if(faces.find(face.key)==faces.end()) {
            reason=TetFacetGeometryReason::ProtectedFaceLost;
            return false;
        }
    for(std::size_t i=0;i<candidate.size();++i)
        for(std::size_t j=i+1;j<candidate.size();++j) {
            if(!consume(metrics,limits))return false;
            ++metrics.pairContactTests;
            if(cavityCellsConflict(
                   view,candidate[i],candidate[j])) {
                reason=TetFacetGeometryReason::CandidateIntersection;
                return false;
            }
        }
    reason=TetFacetGeometryReason::None;
    return true;
}

}  // namespace

TetEdgeGeneralCavityOutcome solveTetEdgeGeneralCavity(
    const TetEdgeGeneralCavityInput& input,
    const TetEdgeGeneralCavityLimits& limits) noexcept {
    TetEdgeGeneralCavityOutcome out;
    try {
        out.metrics.cavityCells=input.removedCells.size();
        out.metrics.cavityVertices=input.vertices.size();
        if(!canonicalInput(input)||!limits.positive())return out;
        if(input.removedCells.size()>limits.maxCavityCells||
           input.vertices.size()>limits.maxCavityVertices) {
            out.status=TetEdgeGeneralCavityStatus::BudgetReached;
            out.geometryReason=TetFacetGeometryReason::None;
            return out;
        }
        if(!checkedChoose4(input.vertices.size(),out.metrics.combinationCount))
            return out;
        if(out.metrics.combinationCount>limits.maxEnumeratedTets) {
            out.status=TetEdgeGeneralCavityStatus::BudgetReached;
            out.geometryReason=TetFacetGeometryReason::None;
            return out;
        }
        const auto sourceView=viewOf(input,input.removedCells);
        std::vector<std::size_t> removed(input.removedCells.size());
        for(std::size_t i=0;i<removed.size();++i)removed[i]=i;
        const auto sourceBall=auditTetFacetBall(sourceView,removed);
        if(!sourceBall.valid) {
            out.geometryReason=sourceBall.reason;
            if(sourceBall.reason==TetFacetGeometryReason::AllocationFailed)
                out.status=TetEdgeGeneralCavityStatus::AllocationFailed;
            return out;
        }
        const auto expectedBoundary=tet_facet_detail::boundarySenses(
            tet_facet_detail::selectedFaceMap(sourceView,removed));
        if(expectedBoundary.empty())return out;

        std::vector<CavityContactFace> fixedFaces;
        if(cadPruning){
            fixedFaces.reserve(expectedBoundary.size()+input.protectedFaces.size());
            auto contactView=sourceView;
            const auto registerFace=[&](std::array<TetFacetVertexId,3> ids){
                contactView.target.vertices=ids;
                const bool valid=tet_facet_detail::validTriangle(contactView);
                fixedFaces.push_back({ids,valid?contactBounds(input.vertices,ids):ContactBox{},valid});
            };
            for(const auto& [face,sense]:expectedBoundary){(void)sense;registerFace(face.vertices);}
            for(const auto& face:input.protectedFaces)registerFace(face.key.vertices);
        }

        std::vector<TetFacetCell> candidates;
        candidates.reserve(static_cast<std::size_t>(out.metrics.combinationCount));
        std::uint64_t sequence=1;
        for(TetFacetVertexId a=0;a<input.vertices.size();++a)
            for(TetFacetVertexId b=a+1;b<input.vertices.size();++b)
                for(TetFacetVertexId c=b+1;c<input.vertices.size();++c)
                    for(TetFacetVertexId d=c+1;d<input.vertices.size();++d) {
                        if(!consume(out.metrics,limits)) {
                            out.status=TetEdgeGeneralCavityStatus::BudgetReached;
                            out.geometryReason=TetFacetGeometryReason::None;
                            return out;
                        }
                        ++out.metrics.enumeratedTets;
                        TetFacetCell cell{{a,b,c,d},sequence++,input.region};
                        const auto sign=predicates::orient3d(
                            input.vertices[a].position,input.vertices[b].position,
                            input.vertices[c].position,input.vertices[d].position);
                        if(sign==predicates::Sign::Zero)continue;
                        if(sign==predicates::Sign::Negative)
                            std::swap(cell.vertices[0],cell.vertices[1]);
                        if(cadPruning) {
                        bool enclosedOther=false;
                        const std::array<Vec3,4> corners{{input.vertices[cell.vertices[0]].position,input.vertices[cell.vertices[1]].position,input.vertices[cell.vertices[2]].position,input.vertices[cell.vertices[3]].position}};
                        for(TetFacetVertexId v=0;v<input.vertices.size();++v) {
                            if(std::find(cell.vertices.begin(),cell.vertices.end(),v)!=cell.vertices.end())continue;
                            if(!consume(out.metrics,limits,4)) {out.status=TetEdgeGeneralCavityStatus::BudgetReached;out.geometryReason=TetFacetGeometryReason::None;return out;}
                            out.metrics.validatorVisits+=4;
                            if(tet_facet_detail::pointStrictlyInsideTet(input.vertices[v].position,corners)) {enclosedOther=true;break;}
                        }
                        if(enclosedOther)continue;
                        bool crossesRequired=false;
                        auto contactView=sourceView;
                        const auto box=contactBounds(input.vertices,cell.vertices);
                        for(const auto& face:fixedFaces) {
                            if(!consume(out.metrics,limits)) {out.status=TetEdgeGeneralCavityStatus::BudgetReached;out.geometryReason=TetFacetGeometryReason::None;return out;}
                            ++out.metrics.validatorVisits;
                            if(!face.valid){crossesRequired=true;break;}
                            if(box.hi.x<face.box.lo.x||face.box.hi.x<box.lo.x||
                               box.hi.y<face.box.lo.y||face.box.hi.y<box.lo.y||
                               box.hi.z<face.box.lo.z||face.box.hi.z<box.lo.z)continue;
                            if(tet_facet_detail::separatedByTetFace(contactView,cell,face.vertices,[]{}))continue;
                            contactView.target.vertices=face.vertices;
                            const auto contact=classifyTetFacetContact(contactView,cell);
                            if(contact!=TetFacetContact::Disjoint&&contact!=TetFacetContact::CoplanarExactFace&&contact!=TetFacetContact::DeclaredBoundaryTouch) {crossesRequired=true;break;}
                        }
                        if(crossesRequired)continue;
                        }
                        candidates.push_back(cell);
                        ++out.metrics.positiveTets;
                    }
        std::sort(candidates.begin(),candidates.end(),[](const auto& a,const auto& b) {
            return keyOf(a)<keyOf(b);
        });
        for(std::size_t i=0;i<candidates.size();++i)candidates[i].sequence=i+1;
        out.enumeratedTetKeys.reserve(candidates.size());
        for(const auto& cell:candidates)
            out.enumeratedTetKeys.push_back(keyOf(cell));

        std::vector<std::size_t> targetCandidates;
        if(cadPruning)for(std::size_t ci=0;ci<candidates.size();++ci) {
            if(!consume(out.metrics,limits)) {out.status=TetEdgeGeneralCavityStatus::BudgetReached;out.geometryReason=TetFacetGeometryReason::None;return out;}
            const auto& ids=candidates[ci].vertices;
            if(std::find(ids.begin(),ids.end(),input.target.lo)!=ids.end()&&
               std::find(ids.begin(),ids.end(),input.target.hi)!=ids.end())targetCandidates.push_back(ci);
        }
        const auto candidateView=viewOf(input,candidates);
        std::map<TetFacetFaceKey,std::vector<std::size_t>> incident;
        std::vector<std::array<std::pair<TetFacetFaceKey,int>,4>> cellFaces;
        cellFaces.resize(candidates.size());
        for(std::size_t ci=0;ci<candidates.size();++ci)
            for(std::size_t f=0;f<4;++f) {
                const auto directed=tet_facet_detail::outwardFace(candidates[ci],f);
                const auto key=tet_facet_detail::faceKey(directed);
                const int sense=tet_facet_detail::faceSense(candidateView,directed);
                cellFaces[ci][f]={key,sense};
                incident[key].push_back(ci);
            }

        std::vector<TetFacetFaceKey> faceKeys;
        faceKeys.reserve(incident.size());
        std::map<TetFacetFaceKey,std::size_t> faceIndex;
        for(const auto& [key,rows]:incident) {
            (void)rows;
            faceIndex.emplace(key,faceKeys.size());
            faceKeys.push_back(key);
        }
        std::vector<std::vector<std::size_t>> incidentByFace(faceKeys.size());
        for(auto& [key,rows]:incident)
            incidentByFace[faceIndex.at(key)]=std::move(rows);
        std::vector<std::array<std::pair<std::size_t,int>,4>> denseCellFaces(
            candidates.size());
        for(std::size_t ci=0;ci<candidates.size();++ci)
            for(std::size_t f=0;f<4;++f)
                denseCellFaces[ci][f]={
                    faceIndex.at(cellFaces[ci][f].first),
                    cellFaces[ci][f].second};
        std::vector<unsigned char> boundaryFace(faceKeys.size(),0);
        std::vector<int> boundarySense(faceKeys.size(),0);
        std::vector<std::size_t> orderedBoundary;
        orderedBoundary.reserve(expectedBoundary.size());
        for(const auto& [key,sense]:expectedBoundary) {
            const auto index=faceIndex.find(key);
            if(index==faceIndex.end())return out;
            boundaryFace[index->second]=1;
            boundarySense[index->second]=sense;
            orderedBoundary.push_back(index->second);
        }
        // Synthetic search constraint: select a target-edge incident cell
        // first. No ledger face or acceptance condition is changed.
        incidentByFace.push_back(std::move(targetCandidates));
        std::vector<LedgerRow> ledger(faceKeys.size());
        std::vector<unsigned char> selected(candidates.size(),0);
        std::vector<TetFacetCell> selectedCells;
        selectedCells.reserve(candidates.size());
        std::vector<SearchFrame> stack;
        stack.reserve(candidates.size()+1);

        auto unresolved=[&]() -> std::optional<std::size_t> {
            if(cadPruning) {
                bool targetPresent=false;
                for(const auto& cell:selectedCells)
                    targetPresent=targetPresent||
                        (std::find(cell.vertices.begin(),cell.vertices.end(),input.target.lo)!=cell.vertices.end()&&
                         std::find(cell.vertices.begin(),cell.vertices.end(),input.target.hi)!=cell.vertices.end());
                if(!targetPresent)return faceKeys.size();
            }
            for(const auto index:orderedBoundary)
                if(ledger[index].count==0)return index;
            for(std::size_t index=0;index<ledger.size();++index)
                if(!boundaryFace[index]&&ledger[index].count==1)return index;
            return std::nullopt;
        };
        auto apply=[&](std::size_t ci,std::array<UndoRow,4>& undo,
                       std::size_t& undoCount) {
            if(selected[ci])return false;
            if(!cadPruning) {
            for(std::size_t other=0;other<candidates.size();++other)
                if(selected[other]) {
                    if(!consume(out.metrics,limits))return false;
                    ++out.metrics.pairContactTests;
                    if(cavityCellsConflict(
                           candidateView,candidates[ci],candidates[other]))
                        return false;
                }
            }
            undoCount=0;
            const auto restore=[&]() {
                while(undoCount!=0) {
                    --undoCount;
                    ledger[undo[undoCount].faceIndex]=undo[undoCount].old;
                }
            };
            for(const auto& [index,sense]:denseCellFaces[ci]) {
                if(!consume(out.metrics,limits)) { restore(); return false; }
                ++out.metrics.faceLedgerVisits;
                const LedgerRow old=ledger[index];
                undo[undoCount++]={index,old};
                if(boundaryFace[index]) {
                    if(old.count!=0||sense!=boundarySense[index]) {
                        restore();
                        return false;
                    }
                } else if(old.count>=2||
                          (old.count==1&&old.firstSense==sense)) {
                    restore();
                    return false;
                }
                auto& row=ledger[index];
                if(row.count==0)row.firstSense=sense;
                ++row.count;
            }
            if(cadPruning) {
            for(std::size_t other=0;other<candidates.size();++other)
                if(selected[other]) {
                    if(!consume(out.metrics,limits)){restore();return false;}
                    ++out.metrics.pairContactTests;
                    if(cavityCellsConflict(
                           candidateView,candidates[ci],candidates[other]))
                        {restore();return false;}
                }
            }
            selected[ci]=1;
            selectedCells.push_back(candidates[ci]);
            return true;
        };
        auto undo=[&](SearchFrame& frame) {
            if(frame.appliedCandidate==(std::numeric_limits<std::size_t>::max)())
                return;
            selected[frame.appliedCandidate]=0;
            selectedCells.pop_back();
            while(frame.undoCount!=0) {
                --frame.undoCount;
                const auto& row=frame.undo[frame.undoCount];
                ledger[row.faceIndex]=row.old;
            }
        };

        const auto first=unresolved();
        if(!first) {
            out.status=TetEdgeGeneralCavityStatus::Unavailable;
            out.metrics.runComplete=true;
            return out;
        }
        stack.push_back({&incidentByFace[*first],0});
        out.metrics.stackHighWater=1;
        bool finalValidationFailed=false;
        TetFacetGeometryReason finalReason=TetFacetGeometryReason::None;
        while(!stack.empty()) {
            auto& frame=stack.back();
            if(frame.options==nullptr||frame.next>=frame.options->size()) {
                SearchFrame finished=std::move(frame);
                stack.pop_back();
                undo(finished);
                continue;
            }
            if(out.metrics.searchStates>=limits.maxSearchStates||
               !consume(out.metrics,limits)) {
                out.status=TetEdgeGeneralCavityStatus::BudgetReached;
                out.geometryReason=TetFacetGeometryReason::None;
                return out;
            }
            ++out.metrics.searchStates;
            const auto candidate=(*frame.options)[frame.next++];
            std::array<UndoRow,4> appliedUndo{};
            std::size_t appliedUndoCount=0;
            if(!apply(candidate,appliedUndo,appliedUndoCount)) {
                if(out.metrics.workUnits>=limits.maxWorkUnits) {
                    out.status=TetEdgeGeneralCavityStatus::BudgetReached;
                    out.geometryReason=TetFacetGeometryReason::None;
                    return out;
                }
                continue;
            }
            const auto nextFace=unresolved();
            if(!nextFace) {
                TetFacetGeometryReason reason=TetFacetGeometryReason::None;
                const bool cycle=targetCycle(
                    input,selectedCells,out.metrics,limits);
                if(cycle&&completeCandidate(input,expectedBoundary,selectedCells,
                                     out.metrics,limits,reason)) {
                    out.cells=selectedCells;
                    out.orderedTetKeys.reserve(out.cells.size());
                    for(const auto& cell:out.cells)
                        out.orderedTetKeys.push_back(keyOf(cell));
                    std::sort(out.orderedTetKeys.begin(),out.orderedTetKeys.end());
                    out.status=TetEdgeGeneralCavityStatus::Solved;
                    out.geometryReason=TetFacetGeometryReason::None;
                    out.metrics.runComplete=true;
                    return out;
                }
                // An allocation failure in the exact validator is not a
                // geometric rejection. Do not hide it by trying another cover.
                if(cycle&&reason==TetFacetGeometryReason::AllocationFailed) {
                    out.status=TetEdgeGeneralCavityStatus::AllocationFailed;
                    out.geometryReason=reason;
                    return out;
                }
                if(out.metrics.workUnits>=limits.maxWorkUnits) {
                    out.status=TetEdgeGeneralCavityStatus::BudgetReached;
                    out.geometryReason=TetFacetGeometryReason::None;
                    return out;
                }
                if(cycle) {
                    finalValidationFailed=true;
                    finalReason=reason;
                }
                selected[candidate]=0;
                selectedCells.pop_back();
                while(appliedUndoCount!=0) {
                    --appliedUndoCount;
                    const auto& row=appliedUndo[appliedUndoCount];
                    ledger[row.faceIndex]=row.old;
                }
                continue;
            }
            const auto* options=&incidentByFace[*nextFace];
            if(options->empty()) {
                selected[candidate]=0;
                selectedCells.pop_back();
                while(appliedUndoCount!=0) {
                    --appliedUndoCount;
                    const auto& row=appliedUndo[appliedUndoCount];
                    ledger[row.faceIndex]=row.old;
                }
                continue;
            }
            SearchFrame child;
            child.options=options;
            child.appliedCandidate=candidate;
            child.undo=appliedUndo;
            child.undoCount=appliedUndoCount;
            stack.push_back(std::move(child));
            out.metrics.stackHighWater=(std::max)(
                out.metrics.stackHighWater,
                static_cast<std::uint64_t>(stack.size()));
        }
        out.status=finalValidationFailed
            ? TetEdgeGeneralCavityStatus::InvalidCandidate
            : TetEdgeGeneralCavityStatus::Unavailable;
        out.geometryReason=finalValidationFailed?finalReason:
            TetFacetGeometryReason::None;
        out.metrics.runComplete=true;
        return out;
    } catch(const std::bad_alloc&) {
        out.status=TetEdgeGeneralCavityStatus::AllocationFailed;
        out.geometryReason=TetFacetGeometryReason::AllocationFailed;
        return out;
    } catch(...) {
        out.status=TetEdgeGeneralCavityStatus::InvalidInput;
        out.geometryReason=TetFacetGeometryReason::InvalidInput;
        return out;
    }
}

TetEdgeOneSteinerOutcome solveTetEdgeGeneralCavityOneSteiner(
    const TetEdgeGeneralCavityInput& input,
    const TetEdgeGeneralCavityLimits& perCandidateLimits,
    std::uint64_t maximumCandidateAttempts,
    std::uint64_t maximumCumulativeWorkUnits) noexcept {
    TetEdgeOneSteinerOutcome out;
    try {
        if (!canonicalInput(input) || !perCandidateLimits.positive() ||
            maximumCandidateAttempts == 0 ||
            maximumCumulativeWorkUnits == 0 ||
            input.vertices.size() == kInvalidTetFacetVertex ||
            maximumCandidateAttempts > 4096)
            return out;

        // Admission precedes proposal allocation and every source-cell access.
        if (input.vertices.size()+1 > perCandidateLimits.maxCavityVertices ||
            input.removedCells.size() > perCandidateLimits.maxCavityCells) {
            out.status = TetEdgeGeneralCavityStatus::BudgetReached;
            out.geometryReason = TetFacetGeometryReason::None;
            return out;
        }
        for (const auto& cell : input.removedCells)
            for (const auto vertex : cell.vertices)
                if (vertex >= input.vertices.size()) return out;
        Vec3 mean{};
        for (const auto& vertex : input.vertices)
            mean = mean + vertex.position;
        mean = mean / static_cast<double>(input.vertices.size());
        const Vec3 midpoint =
            (input.vertices[input.target.lo].position +
             input.vertices[input.target.hi].position) * 0.5;
        std::vector<Vec3> candidates{
            mean, (mean + midpoint) * 0.5,
            mean * 0.75 + midpoint * 0.25,
            mean * 0.25 + midpoint * 0.75};

        std::vector<const TetFacetCell*> orderedCells;
        orderedCells.reserve(input.removedCells.size());
        for (const auto& cell : input.removedCells)
            orderedCells.push_back(&cell);
        std::sort(orderedCells.begin(), orderedCells.end(),
                  [](const auto* a, const auto* b) {
                      const auto ak = keyOf(*a);
                      const auto bk = keyOf(*b);
                      return ak != bk ? ak < bk : a->sequence < b->sequence;
                  });
        for (const auto* cell : orderedCells) {
            Vec3 center{};
            for (const auto vertex : cell->vertices) {
                if (vertex >= input.vertices.size()) return out;
                center = center + input.vertices[vertex].position;
            }
            candidates.push_back(center * 0.25);
        }

        // For a larger recovered ring, cover every source Tet before the
        // global vertex-blend family consumes the checked candidate cap.
        // The previous category-major ordering reached only the first half
        // of an eleven-cell ABC cavity at 96 attempts.  Small accepted
        // witnesses retain their historical ordering exactly.
        constexpr std::size_t kInterleavedSourceTetThreshold=8;
        if(orderedCells.size()>=kInterleavedSourceTetThreshold) {
            for(const auto* cell:orderedCells) {
                Vec3 sum{};
                for(const auto vertex:cell->vertices)
                    sum=sum+input.vertices[vertex].position;
                const Vec3 center=sum*0.25;
                candidates.push_back((center+midpoint)*0.5);
                candidates.push_back((center+mean)*0.5);
                for(const auto vertex:cell->vertices)
                    candidates.push_back(
                        (sum+input.vertices[vertex].position)*0.2);
            }
        }

        // Keep the accepted original sequence above byte-for-byte first, then
        // broaden the bounded deterministic search for CAD cavities whose
        // feasible kernel is missed by means and cell centroids. Every point
        // remains only a numerical proposal: solveTetEdgeGeneralCavity exact-
        // certifies all replacement Tets, protected features and closure.
        std::vector<const TetFacetVertex*> orderedVertices;
        orderedVertices.reserve(input.vertices.size());
        for (const auto& vertex : input.vertices)
            orderedVertices.push_back(&vertex);
        std::sort(orderedVertices.begin(), orderedVertices.end(),
                  [](const auto* a, const auto* b) {
                      if (a->sequence != b->sequence)
                          return a->sequence < b->sequence;
                      if (a->position.x != b->position.x)
                          return a->position.x < b->position.x;
                      if (a->position.y != b->position.y)
                          return a->position.y < b->position.y;
                      return a->position.z < b->position.z;
                  });
        for (const auto* vertex : orderedVertices) {
            candidates.push_back((vertex->position + midpoint) * 0.5);
            candidates.push_back(vertex->position * 0.25 + midpoint * 0.75);
            candidates.push_back((vertex->position + mean) * 0.5);
        }
        for (const auto* cell : orderedCells) {
            Vec3 center{};
            for (const auto vertex : cell->vertices)
                center = center + input.vertices[vertex].position;
            center = center * 0.25;
            candidates.push_back((center + midpoint) * 0.5);
            candidates.push_back((center + mean) * 0.5);
        }

        // A centroid alone can sit on the wrong side of the visibility
        // kernel in a short non-star strip. Try the four deterministic
        // 1/5-interior offsets of each source Tet before lower-authority
        // boundary-face averages: barycentric weights (2,1,1,1) and their
        // permutations. These remain strictly interior numerical proposals;
        // the exact solver below still proves every publication invariant.
        for (const auto* cell : orderedCells) {
            Vec3 sum{};
            for (const auto vertex : cell->vertices)
                sum = sum + input.vertices[vertex].position;
            for (const auto vertex : cell->vertices)
                candidates.push_back(
                    (sum + input.vertices[vertex].position) * 0.2);
        }

        std::map<std::array<TetFacetVertexId, 3>, std::size_t> faceUses;
        constexpr std::array<std::array<std::size_t, 3>, 4> localFaces{{
            {{1, 2, 3}}, {{0, 3, 2}}, {{0, 1, 3}}, {{0, 2, 1}}}};
        for (const auto* cell : orderedCells)
            for (const auto& local : localFaces) {
                std::array<TetFacetVertexId, 3> face{
                    cell->vertices[local[0]], cell->vertices[local[1]],
                    cell->vertices[local[2]]};
                std::sort(face.begin(), face.end());
                ++faceUses[face];
            }
        for (const auto& [face, uses] : faceUses) {
            if (uses != 1U) continue;
            const Vec3 center =
                (input.vertices[face[0]].position +
                 input.vertices[face[1]].position +
                 input.vertices[face[2]].position) / 3.0;
            candidates.push_back((center + midpoint) * 0.5);
            candidates.push_back((center + mean) * 0.5);
        }

        std::set<std::array<std::uint64_t, 3>> exactCandidates;
        std::vector<TetFacetVertex> vertices(input.vertices.begin(),
                                             input.vertices.end());
        const std::uint64_t nextSequence =
            input.vertices.empty() ? 1 : input.vertices.back().sequence + 1;
        if (nextSequence == 0) return out;
        vertices.push_back({{}, nextSequence});
        const std::uint64_t available = static_cast<std::uint64_t>(
            candidates.size());
        const std::uint64_t attemptedLimit =
            (std::min)(maximumCandidateAttempts, available);
        bool sawInvalidCandidate = false;
        bool sawCandidateBudget = false;
        TetFacetGeometryReason lastInvalidReason =
            TetFacetGeometryReason::None;
        for (std::uint64_t ordinal = 0; ordinal < attemptedLimit; ++ordinal) {
            if (out.cumulativeWorkUnits >= maximumCumulativeWorkUnits) {
                out.status = TetEdgeGeneralCavityStatus::BudgetReached;
                out.geometryReason = TetFacetGeometryReason::None;
                return out;
            }
            const Vec3 point = candidates[static_cast<std::size_t>(ordinal)];
            if (!tet_facet_detail::finite(point)) continue;
            const std::array<std::uint64_t, 3> bits{
                std::bit_cast<std::uint64_t>(point.x == 0.0 ? 0.0 : point.x),
                std::bit_cast<std::uint64_t>(point.y == 0.0 ? 0.0 : point.y),
                std::bit_cast<std::uint64_t>(point.z == 0.0 ? 0.0 : point.z)};
            if (!exactCandidates.insert(bits).second) continue;
            vertices.back().position = point;
            auto expanded = input;
            expanded.vertices = vertices;
            auto candidateLimits = perCandidateLimits;
            candidateLimits.maxWorkUnits = (std::min)(
                candidateLimits.maxWorkUnits,
                maximumCumulativeWorkUnits - out.cumulativeWorkUnits);
            auto candidate = solveTetEdgeGeneralCavity(
                expanded, candidateLimits);
            ++out.candidateAttempts;
            if (out.cumulativeWorkUnits >
                (std::numeric_limits<std::uint64_t>::max)() -
                    candidate.metrics.workUnits)
                out.cumulativeWorkUnits =
                    (std::numeric_limits<std::uint64_t>::max)();
            else
                out.cumulativeWorkUnits += candidate.metrics.workUnits;
            out.status = candidate.status;
            out.geometryReason = candidate.geometryReason;
            out.cavity = std::move(candidate);
            if (out.status == TetEdgeGeneralCavityStatus::Solved) {
                out.steinerPosition = point;
                out.selectedCandidateOrdinal = ordinal;
                return out;
            }
            if (out.status ==
                TetEdgeGeneralCavityStatus::InvalidCandidate) {
                sawInvalidCandidate = true;
                lastInvalidReason = out.geometryReason;
                continue;
            }
            if(out.status==TetEdgeGeneralCavityStatus::BudgetReached) {
                sawCandidateBudget=true;
                if(out.cumulativeWorkUnits>=maximumCumulativeWorkUnits)
                    return out;
                continue;
            }
            if (out.status != TetEdgeGeneralCavityStatus::Unavailable)
                return out;
        }
        if (maximumCandidateAttempts < available) {
            out.status = TetEdgeGeneralCavityStatus::BudgetReached;
            out.geometryReason = TetFacetGeometryReason::None;
        } else if(sawCandidateBudget) {
            out.status=TetEdgeGeneralCavityStatus::BudgetReached;
            out.geometryReason=TetFacetGeometryReason::None;
        } else if (sawInvalidCandidate) {
            out.status = TetEdgeGeneralCavityStatus::InvalidCandidate;
            out.geometryReason = lastInvalidReason;
        } else {
            out.status = TetEdgeGeneralCavityStatus::Unavailable;
            out.geometryReason = TetFacetGeometryReason::None;
        }
        return out;
    } catch (const std::bad_alloc&) {
        out.status = TetEdgeGeneralCavityStatus::AllocationFailed;
        out.geometryReason = TetFacetGeometryReason::AllocationFailed;
        return out;
    } catch (...) {
        out.status = TetEdgeGeneralCavityStatus::InvalidInput;
        out.geometryReason = TetFacetGeometryReason::InvalidInput;
        return out;
    }
}

}  // namespace paacs::bnd

namespace paacs::bnd::detail {
bool cadGeneralPruningEnabled() noexcept { return cadPruning; }
bool exchangeCadGeneralPruning(bool enabled) noexcept {
    const bool previous=cadPruning;cadPruning=enabled;return previous;
}
}
