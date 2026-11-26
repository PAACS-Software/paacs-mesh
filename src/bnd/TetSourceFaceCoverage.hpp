// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/bnd/TetFacetGeometry.hpp"
#include <atomic>
#include <bit>
#include <map>
#include <new>

namespace paacs::bnd::detail {
struct SourceFaceCoverageLimits {
    std::uint64_t maximumWork=10'000'000,maximumWorkingBytes=256ULL*1024*1024;
    std::uint64_t maximumNodes=1'000'000,maximumSourceTriangles=1'000'000,maximumPartialCells=5'000'000;
    const std::atomic<bool>* cancellation=nullptr;
};
enum class SourceFaceCoverageReason {None,InvalidInput,BudgetExceeded,Cancelled,AllocationFailed,InternalError};
struct SourceFaceCoverageOutcome {
    bool complete=false;
    SourceFaceCoverageReason reason=SourceFaceCoverageReason::InternalError;
    const char* diagnostic="not-attempted";
    std::uint64_t workUsed=0,modeledBytes=0,sourceTriangles=0,partialCells=0,sourceFacesCovered=0,tetFacesVisited=0;
    // Original body-local source triangle IDs, sorted. This receipt describes
    // connectivity only, never a complete or geometrically valid volume mesh.
    std::vector<Index> missingSourceTriangles;
};
namespace source_face_coverage_detail {
struct Stop {SourceFaceCoverageReason reason;const char* diagnostic;};
using Face=std::array<Index,3>;
struct Entry {Face directed{};Index source=-1;unsigned count=0;};
inline bool same(Vec3 a,Vec3 b){
    return std::bit_cast<std::uint64_t>(a.x)==std::bit_cast<std::uint64_t>(b.x)&&
           std::bit_cast<std::uint64_t>(a.y)==std::bit_cast<std::uint64_t>(b.y)&&
           std::bit_cast<std::uint64_t>(a.z)==std::bit_cast<std::uint64_t>(b.z);
}
inline bool sameDirected(Face a,Face b){
    for(unsigned shift=0;shift<3;++shift)
        if(a[0]==b[shift]&&a[1]==b[(shift+1)%3]&&a[2]==b[(shift+2)%3])return true;
    return false;
}
}

// Called only on an independently certified source and retained partial mesh.
// It maps source face identities only, not every interior face. The unchanged
// source prefix and source-face orientation/owner counts are checked here;
// complete source/contact/centroid/pair and Native gates remain mandatory.
inline SourceFaceCoverageOutcome collectMissingSourceTriangles(const mesh::Mesh& source,
    const mesh::Mesh& partial,mesh::Classification region,const SourceFaceCoverageLimits& limits) noexcept {
    using CoverageStop=source_face_coverage_detail::Stop;
    using Face=source_face_coverage_detail::Face;
    using Entry=source_face_coverage_detail::Entry;
    SourceFaceCoverageOutcome out;
    const auto need=[](bool ok,SourceFaceCoverageReason reason,const char* diagnostic){if(!ok)throw CoverageStop{reason,diagnostic};};
    const auto tick=[&](std::uint64_t n=1){
        if(limits.cancellation&&limits.cancellation->load())throw CoverageStop{SourceFaceCoverageReason::Cancelled,"cancelled"};
        need(out.workUsed<=limits.maximumWork&&n<=limits.maximumWork-out.workUsed,
             SourceFaceCoverageReason::BudgetExceeded,"source-coverage-work-budget");out.workUsed+=n;
    };
    try {
        tick(0);
        const auto ns=source.nodeCount(),np=partial.nodeCount(),fs=source.cellCount(),cells=partial.cellCount();
        need(&source!=&partial&&ns>=3&&np>=ns&&fs>0&&region.dim==3&&region.tag>=0,
             SourceFaceCoverageReason::InvalidInput,"source-coverage-input");
        out.sourceTriangles=static_cast<std::uint64_t>(fs);out.partialCells=static_cast<std::uint64_t>(cells);
        need(std::uint64_t(np)<kInvalidTetFacetVertex&&std::uint64_t(np)<=limits.maximumNodes&&
             std::uint64_t(fs)<=limits.maximumSourceTriangles&&std::uint64_t(cells)<=limits.maximumPartialCells,
             SourceFaceCoverageReason::BudgetExceeded,"source-coverage-count-budget");
        // Conservative map-node, flag and returned-ID storage plus fixed
        // predicate scratch. There is no full tetrahedron-face map/copy.
        need(std::uint64_t(fs)<=(UINT64_MAX-32*1024)/256,
             SourceFaceCoverageReason::BudgetExceeded,"source-coverage-byte-overflow");
        out.modeledBytes=std::uint64_t(fs)*256+32*1024;
        need(out.modeledBytes<=limits.maximumWorkingBytes,SourceFaceCoverageReason::BudgetExceeded,"source-coverage-working-budget");
        const auto lookupWork=std::uint64_t(2*std::bit_width(std::uint64_t(fs))+1);
        for(Index i=0;i<ns;++i){tick();
            const auto entity=source.nodeClass(i);
            need(tet_facet_detail::finite(source.node(i))&&entity.dim>=0&&entity.dim<=2&&entity.tag>=0&&
                 source_face_coverage_detail::same(source.node(i),partial.node(i))&&entity==partial.nodeClass(i),
                 SourceFaceCoverageReason::InvalidInput,"source-coverage-prefix");
        }
        std::map<Face,Entry> faces;
        std::vector<unsigned char> missing(static_cast<std::size_t>(fs),1);
        for(Index i=0;i<fs;++i){tick(lookupWork+4);const auto& cell=source.cell(i);
            need(cell.type==mesh::CellType::Tri3&&cell.onEntity.dim==2&&cell.onEntity.tag>=0,
                 SourceFaceCoverageReason::InvalidInput,"source-coverage-source-row");
            const auto* ids=source.cellNodes(i);const Face directed{ids[0],ids[1],ids[2]};
            for(auto id:directed)need(id>=0&&id<ns,SourceFaceCoverageReason::InvalidInput,"source-coverage-source-index");
            auto key=directed;std::sort(key.begin(),key.end());
            need(key[0]!=key[1]&&key[1]!=key[2],SourceFaceCoverageReason::InvalidInput,"source-coverage-source-degenerate");
            need(faces.emplace(key,Entry{directed,i,0}).second,SourceFaceCoverageReason::InvalidInput,"source-coverage-duplicate-source-face");
        }
        for(Index ci=0;ci<cells;++ci){tick(6);const auto& row=partial.cell(ci);
            need(row.type==mesh::CellType::Tet4&&row.onEntity==region,
                 SourceFaceCoverageReason::InvalidInput,"source-coverage-partial-row");
            const auto* ids=partial.cellNodes(ci);TetFacetCell cell;
            std::array<Index,4> unique;
            for(unsigned j=0;j<4;++j){
                need(ids[j]>=0&&ids[j]<np,SourceFaceCoverageReason::InvalidInput,"source-coverage-tet-index");
                need(tet_facet_detail::finite(partial.node(ids[j])),SourceFaceCoverageReason::InvalidInput,"source-coverage-nonfinite-tet");
                need(ids[j]<ns||partial.nodeClass(ids[j])==region,SourceFaceCoverageReason::InvalidInput,"source-coverage-interior-region");
                unique[j]=ids[j];cell.vertices[j]=static_cast<TetFacetVertexId>(ids[j]);
            }
            std::sort(unique.begin(),unique.end());
            need(std::adjacent_find(unique.begin(),unique.end())==unique.end(),SourceFaceCoverageReason::InvalidInput,"source-coverage-repeated-tet-node");
            need(predicates::orient3d(partial.node(ids[0]),partial.node(ids[1]),partial.node(ids[2]),partial.node(ids[3]))==predicates::Sign::Positive,
                 SourceFaceCoverageReason::InvalidInput,"source-coverage-nonpositive-tet");
            for(unsigned side=0;side<4;++side){tick(lookupWork);++out.tetFacesVisited;
                const auto directedIds=tet_facet_detail::outwardFace(cell,side);
                const Face directed{directedIds[0],directedIds[1],directedIds[2]};
                auto key=directed;std::sort(key.begin(),key.end());const auto found=faces.find(key);
                if(found==faces.end())continue;
                auto& face=found->second;
                need(face.count==0,SourceFaceCoverageReason::InvalidInput,"source-coverage-multiple-source-owners");
                need(source_face_coverage_detail::sameDirected(face.directed,directed),SourceFaceCoverageReason::InvalidInput,"source-coverage-source-winding");
                ++face.count;++out.sourceFacesCovered;missing[static_cast<std::size_t>(face.source)]=0;
            }
        }
        std::vector<Index> result;result.reserve(static_cast<std::size_t>(fs)-static_cast<std::size_t>(out.sourceFacesCovered));
        for(Index i=0;i<fs;++i){tick();if(missing[static_cast<std::size_t>(i)])result.push_back(i);}
        tick(0);out.missingSourceTriangles.swap(result);out.complete=true;out.reason=SourceFaceCoverageReason::None;
        out.diagnostic="source-face-coverage-complete";
    }catch(const CoverageStop& stop){out.reason=stop.reason;out.diagnostic=stop.diagnostic;}
    catch(const std::bad_alloc&){out.reason=SourceFaceCoverageReason::AllocationFailed;out.diagnostic="allocation-failed";}
    catch(...){out.reason=SourceFaceCoverageReason::InternalError;out.diagnostic="internal-error";}
    return out;
}
}
