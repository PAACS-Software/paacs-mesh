// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/ProtectedFeatureRefinement.hpp"

#include "CheckedSizing.hpp"
#include "paacs/algo/FacetBarycentre.hpp"
#include "paacs/algo/FacetConflictResolver.hpp"
#include "T2CSurfaceCache.hpp"
#include "paacs/algo/Meshers.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <new>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <tuple>

namespace paacs::algo {
namespace {

bool finite(const Vec3& point) noexcept {
    return std::isfinite(point.x) && std::isfinite(point.y) &&
           std::isfinite(point.z);
}

bool finite(const Vec2& point) noexcept {
    return std::isfinite(point.u) && std::isfinite(point.v);
}

double canonicalZero(double value) noexcept {
    return value==0.0?0.0:value;
}

std::uint64_t parameterBits(double value) noexcept {
    return std::bit_cast<std::uint64_t>(canonicalZero(value));
}

bool exactPoint(const Vec3& a,const Vec3& b) noexcept {
    return a.x==b.x&&a.y==b.y&&a.z==b.z;
}

bool checkedMetricAdd(std::uint64_t& total,std::uint64_t value) noexcept {
    if(value>(std::numeric_limits<std::uint64_t>::max)()-total)return false;
    total+=value;
    return true;
}

bool checkedSizeAdd(std::size_t& total,std::size_t value) noexcept {
    if(value>(std::numeric_limits<std::size_t>::max)()-total)return false;
    total+=value;
    return true;
}

bool checkedSizeBytes(std::size_t& total,std::size_t count,
                      std::size_t width) noexcept {
    if(count!=0&&width>(std::numeric_limits<std::size_t>::max)()/count)
        return false;
    return checkedSizeAdd(total,count*width);
}

bool meshAccounted(const mesh::MeshCapacitySnapshot& snapshot,
                   std::size_t& total) noexcept {
    if(!snapshot.complete)return false;
    total=0;
    constexpr std::size_t entryWidth=sizeof(std::uint64_t)+
        sizeof(std::vector<Index>)+2*sizeof(void*);
    return checkedSizeBytes(total,snapshot.nodeCoordinateCapacity,sizeof(Vec3))&&
           checkedSizeBytes(total,snapshot.nodeClassificationCapacity,
                            sizeof(mesh::Classification))&&
           checkedSizeBytes(total,snapshot.cellCapacity,sizeof(mesh::Cell))&&
           checkedSizeBytes(total,snapshot.connectivityCapacity,sizeof(Index))&&
           checkedSizeBytes(total,snapshot.dedupBucketCount,sizeof(void*))&&
           checkedSizeBytes(total,snapshot.dedupEntryCount,entryWidth)&&
           checkedSizeBytes(total,snapshot.dedupIndexCapacity,sizeof(Index));
}

std::vector<double> edgeParameters(const geo::GEdge& edge,
                                   const sizing::MeshSpec& spec) {
    const double length=edge.length();
    if(!std::isfinite(length)||length<0.0)
        throw detail::SizingEvaluationError{};
    if(length<=0.0)return {};
    constexpr int samples=64;
    std::array<double,samples+1> cumulative{};
    double integral=0.0;
    const double sampleLength=length/samples;
    if(!detail::positiveFinite(sampleLength))
        throw detail::SizingEvaluationError{};
    for(int i=0;i<samples;++i) {
        const double t=(i+0.5)/samples;
        const double h=(std::max)(detail::evaluatedGeometrySize(
            spec,edge.pointAt(t),edge.curvatureAt(t)),1.0e-12);
        const double contribution=sampleLength/h;
        const double next=integral+contribution;
        if(!detail::positiveFinite(contribution)||!std::isfinite(next)||
           !(next>integral))throw detail::SizingEvaluationError{};
        integral=next;
        cumulative[static_cast<std::size_t>(i)+1]=integral;
    }
    const double rounded=std::round(integral);
    if(!std::isfinite(rounded)||rounded>
       static_cast<double>((std::numeric_limits<int>::max)()))
        throw detail::SizingEvaluationError{};
    const int count=(std::max)(1,static_cast<int>(rounded));
    std::vector<double> parameters;
    parameters.reserve(static_cast<std::size_t>(count)+1);
    parameters.push_back(0.0);
    for(int k=1;k<count;++k) {
        const double target=integral*k/count;
        int i=0;
        while(i<samples&&cumulative[static_cast<std::size_t>(i)+1]<target)++i;
        const double denominator=(std::max)(
            cumulative[static_cast<std::size_t>(i)+1]-
            cumulative[static_cast<std::size_t>(i)],1.0e-30);
        const double fraction=(target-cumulative[static_cast<std::size_t>(i)])/
            denominator;
        const double parameter=(i+fraction)/samples;
        if(!std::isfinite(parameter)||parameter<0.0||parameter>1.0)
            throw detail::SizingEvaluationError{};
        parameters.push_back(parameter);
    }
    parameters.push_back(1.0);
    return parameters;
}

ArtifactStamp entityStamp(int tag) noexcept {
    return {static_cast<std::uint64_t>(static_cast<std::uint32_t>(tag))+1U,1};
}

SurfaceNodeIdentity curveIdentity(ArtifactStamp owningCurve,
                                  std::uint64_t baseSegment,
                                  std::uint64_t numerator,
                                  std::uint32_t level,double parameter) noexcept {
    SurfaceNodeIdentity result;
    result.domain=SurfaceNodeIdentityDomain::RefinedCurve;
    result.parentArtifact=owningCurve;
    result.parentOrdinal=baseSegment;
    result.parameterUBits=parameterBits(parameter);
    result.dyadicNumerator=numerator;
    result.dyadicLevel=level;
    return result;
}

SurfaceNodeIdentity faceIdentity(SurfaceNodeIdentityDomain domain,
                                 ArtifactStamp faceStamp,
                                 std::uint64_t cellOrdinal,
                                 std::uint64_t subordinal,
                                 const Vec2& parameter) noexcept {
    SurfaceNodeIdentity result;
    result.domain=domain;
    result.parentArtifact=faceStamp;
    result.parentOrdinal=cellOrdinal;
    result.parentSubordinal=subordinal;
    result.parameterUBits=parameterBits(parameter.u);
    result.parameterVBits=parameterBits(parameter.v);
    return result;
}

struct SelectedSurfaceEdge {
    ArtifactStamp faceStamp{};
    FacePersistentId facePersistentId = 0;
    int faceTag = -1;
    std::uint64_t localCellOrdinal = 0;
    std::uint64_t localEdgeOrdinal = 0;
    std::array<Index,3> parentFacetGlobalNodes{
        kInvalidIndex,kInvalidIndex,kInvalidIndex};
    std::array<Index,2> globalNodes{kInvalidIndex,kInvalidIndex};
    Vec2 firstParameter{};
    Vec2 secondParameter{};
    std::size_t incidentChartCount = 0;
};

struct SurfaceSupportIndex {
    std::map<std::array<Index,2>,SelectedSurfaceEdge> edges;
    std::map<int,ArtifactStamp> owningCurveByTag;
    std::uint64_t workUnits = 0;
};

std::optional<SurfaceSupportIndex> buildSurfaceSupportIndex(
    const BoundaryConformingV3Outcome& source) {
    if(!source.surfaceOutcome||!source.surfaceOutcome->retainedArtifact)
        return std::nullopt;
    SurfaceSupportIndex result;
    std::map<std::array<Index,2>,std::set<ArtifactStamp>> chartsByEdge;
    for(const auto& face:source.surfaceOutcome->retainedArtifact->faces) {
        if(!face||face->publishedGlobalByLocalNode.size()!=
                     face->localMappingParams.size())continue;
        for(const auto& chain:face->normalizedChains) {
            if(!checkedMetricAdd(result.workUnits,1))return std::nullopt;
            if(chain.classification.dim!=1||!chain.owningCurve)continue;
            const auto [it,inserted]=result.owningCurveByTag.emplace(
                chain.classification.tag,*chain.owningCurve);
            if(!inserted&&it->second!=*chain.owningCurve)return std::nullopt;
        }
        for(std::size_t ci=0;ci<face->localCells.size();++ci) {
            const auto& cell=face->localCells[ci];
            constexpr std::array<std::array<std::size_t,2>,3> edges{{
                {{0,1}},{{1,2}},{{2,0}}}};
            for(std::size_t ei=0;ei<edges.size();++ei) {
                if(!checkedMetricAdd(result.workUnits,1))return std::nullopt;
                const auto a=cell.nodes[edges[ei][0]];
                const auto b=cell.nodes[edges[ei][1]];
                if(a>=face->publishedGlobalByLocalNode.size()||
                   b>=face->publishedGlobalByLocalNode.size())continue;
                std::array<Index,2> key{{
                    face->publishedGlobalByLocalNode[a],
                    face->publishedGlobalByLocalNode[b]}};
                if(key[1]<key[0])std::swap(key[0],key[1]);
                chartsByEdge[key].insert(face->faceStamp);
                SelectedSurfaceEdge row;
                row.faceStamp=face->faceStamp;
                row.facePersistentId=face->faceStamp.persistentId;
                row.faceTag=face->faceClass.tag;
                row.localCellOrdinal=ci;
                row.localEdgeOrdinal=ei;
                for(std::size_t vi=0;vi<3;++vi) {
                    if(cell.nodes[vi]>=face->publishedGlobalByLocalNode.size())
                        return std::nullopt;
                    row.parentFacetGlobalNodes[vi]=
                        face->publishedGlobalByLocalNode[cell.nodes[vi]];
                }
                std::sort(row.parentFacetGlobalNodes.begin(),
                          row.parentFacetGlobalNodes.end());
                row.globalNodes=key;
                row.firstParameter=face->localMappingParams[a];
                row.secondParameter=face->localMappingParams[b];
                const auto found=result.edges.find(key);
                if(found==result.edges.end())result.edges.emplace(key,row);
                else {
                    if(std::tie(row.faceStamp,row.parentFacetGlobalNodes,
                                row.localEdgeOrdinal)<
                       std::tie(found->second.faceStamp,
                                found->second.parentFacetGlobalNodes,
                                found->second.localEdgeOrdinal))
                        found->second=row;
                }
            }
        }
    }
    for(auto& [key,row]:result.edges)
        row.incidentChartCount=chartsByEdge.at(key).size();
    return result;
}

struct PlannedFeaturePlacement {
    std::size_t frontierOrdinal = 0;
    bool curve = false;
    ProtectedFeatureEntityKind kind = ProtectedFeatureEntityKind::FaceEdge;
    int entityTag = -1;
    FacePersistentId facePersistentId = 0;
    ArtifactStamp entityStamp{};
    SurfaceNodeIdentity identity{};
    double curveParameter = 0.0;
    Vec2 faceParameter{};
    t2c_detail::RefinedSurfaceSegment segment{};
    SelectedSurfaceEdge surfaceEdge{};
    std::uint64_t witness = 0;
    const geo::GEdge* geometryEdge = nullptr;
    const geo::GFace* geometryFace = nullptr;
};

struct PreparedFeaturePlacement {
    PlannedFeaturePlacement plan{};
    ProtectedFeaturePlacement placement{};
};

struct SourceOutcomeRetention {
    std::unique_ptr<BoundaryConformingV3Outcome>& pending;
    BoundaryConformingV4Outcome& outcome;
    ~SourceOutcomeRetention() noexcept {
        if(pending&&!outcome.lastSourceOutcome)
            outcome.lastSourceOutcome=std::move(pending);
    }
};

bool sameIdentity(const PlannedFeaturePlacement& row,
                  const SurfaceNodeIdentity& identity) noexcept {
    return row.identity==identity;
}

bool supportsOverlap(const bnd::TetFacetSetFrontierSink& frontier,
                     const bnd::TetFacetSetFrontierRow& a,
                     const bnd::TetFacetSetFrontierRow& b) noexcept {
    if(a.supportOffset>frontier.supportCellsUsed||
       a.supportCount>frontier.supportCellsUsed-a.supportOffset||
       b.supportOffset>frontier.supportCellsUsed||
       b.supportCount>frontier.supportCellsUsed-b.supportOffset)return true;
    const auto first=frontier.supportCells.subspan(a.supportOffset,
                                                   a.supportCount);
    const auto second=frontier.supportCells.subspan(b.supportOffset,
                                                    b.supportCount);
    for(const auto cell:first)
        if(std::find(second.begin(),second.end(),cell)!=second.end())
            return true;
    return false;
}

std::uint64_t frontierWitness(
    const bnd::TetFacetSetFrontierRow& row) noexcept {
    auto mix=[](std::uint64_t hash,std::uint64_t value) noexcept {
        for(unsigned byte=0;byte<8;++byte) {
            hash^=(value>>(byte*8U))&0xffU;
            hash*=1099511628211ULL;
        }
        return hash;
    };
    std::uint64_t digest=1469598103934665603ULL;
    digest=mix(digest,static_cast<std::uint64_t>(
        row.kind==bnd::TetFacetSetWorkKind::Edge
            ? VolumeWorkKind::ArtifactEdge : VolumeWorkKind::Facet));
    if(row.kind==bnd::TetFacetSetWorkKind::Edge) {
        digest=mix(digest,static_cast<std::uint64_t>(
            row.edgeKey.canonicalGlobalNodes[0]));
        digest=mix(digest,static_cast<std::uint64_t>(
            row.edgeKey.canonicalGlobalNodes[1]));
        digest=mix(digest,static_cast<std::uint64_t>(
            row.edgeOutcome.reason));
    } else {
        digest=mix(digest,row.facetKey.facePersistentId);
        digest=mix(digest,row.facetKey.normalizedLocalCellOrdinal);
        digest=mix(digest,static_cast<std::uint64_t>(
            row.facetOutcome.reason));
    }
    return digest;
}

std::uint64_t facetConflictWitness(
    const bnd::TetFacetConflictCycleSink& sink) noexcept {
    auto mix=[](std::uint64_t hash,std::uint64_t value) noexcept {
        for(unsigned byte=0;byte<8;++byte) {
            hash^=(value>>(byte*8U))&0xffU;
            hash*=1099511628211ULL;
        }
        return hash;
    };
    std::uint64_t digest=1469598103934665603ULL;
    digest=mix(digest,sink.artifactStamp.persistentId);
    digest=mix(digest,sink.artifactStamp.revision);
    digest=mix(digest,static_cast<std::uint64_t>(sink.membersUsed));
    for(std::size_t i=0;i<sink.membersUsed;++i) {
        digest=mix(digest,sink.members[i].facePersistentId);
        digest=mix(digest,static_cast<std::uint64_t>(
            sink.members[i].normalizedLocalCellOrdinal));
    }
    digest=mix(digest,static_cast<std::uint64_t>(sink.dependenciesUsed));
    for(std::size_t i=0;i<sink.dependenciesUsed;++i)
        for(const auto& key:sink.dependencies[i]) {
            digest=mix(digest,key.facePersistentId);
            digest=mix(digest,static_cast<std::uint64_t>(
                key.normalizedLocalCellOrdinal));
        }
    return digest;
}

#if defined(PAACS_TOPO_TESTING)
thread_local std::optional<SurfaceNodeIdentity>
    facetConflictRepeatedIdentityForTesting;
#endif

BoundaryConformingV4Status sourceStatus(BoundaryConformingV3Status status) {
    switch(status) {
        case BoundaryConformingV3Status::Published:
            return BoundaryConformingV4Status::Published;
        case BoundaryConformingV3Status::Rejected:
            return BoundaryConformingV4Status::Rejected;
        case BoundaryConformingV3Status::Unsupported:
            return BoundaryConformingV4Status::Unsupported;
        case BoundaryConformingV3Status::Refused:
            return BoundaryConformingV4Status::Refused;
        case BoundaryConformingV3Status::Failed:
            return BoundaryConformingV4Status::Failed;
    }
    return BoundaryConformingV4Status::Failed;
}

BoundaryConformingV4Reason sourceReason(BoundaryConformingV3Status status) {
    switch(status) {
        case BoundaryConformingV3Status::Rejected:
            return BoundaryConformingV4Reason::SourceRejected;
        case BoundaryConformingV3Status::Unsupported:
            return BoundaryConformingV4Reason::SourceUnsupported;
        case BoundaryConformingV3Status::Refused:
            return BoundaryConformingV4Reason::SourceRefused;
        case BoundaryConformingV3Status::Failed:
            return BoundaryConformingV4Reason::SourceFailed;
        case BoundaryConformingV3Status::Published:
            return BoundaryConformingV4Reason::None;
    }
    return BoundaryConformingV4Reason::InternalError;
}

}  // namespace

namespace t2c_detail {

#if defined(PAACS_TOPO_TESTING)
void setFacetConflictRepeatedIdentityForTesting(
    const SurfaceNodeIdentity& identity) noexcept {
    facetConflictRepeatedIdentityForTesting=identity;
}

void resetFacetConflictRepeatedIdentityForTesting() noexcept {
    facetConflictRepeatedIdentityForTesting.reset();
}
#endif

void meshEdgesWithProtectedSplits(
    const geo::Model& model,const sizing::MeshSpec& spec,
    std::span<const DyadicEdgeSplit> splits,mesh::Mesh& output,
    PrivateBuildTrace& trace,
    std::span<const DyadicEdgeSplit> proposedSplits) {
    if(proposedSplits.size()>
       (std::numeric_limits<std::size_t>::max)()-splits.size())
        throw std::length_error("protected split overlay overflow");
    const std::size_t splitRows=splits.size()+proposedSplits.size();
    mesh::Mesh base;
    mesh::Mesh* source=&base;
    if(splitRows==0) {
        // Fast trace-only case used by single-pass face prefill. The former
        // route meshed every edge twice and copied the complete first mesh.
        // Meshing directly into the requested output preserves the ordinary
        // edge order and identities while eliminating both costs.
        output=mesh::Mesh{};
        meshEdges(model,spec,output);
        source=&output;
    } else {
        meshEdges(model,spec,base);
        output=mesh::Mesh{};
        output.reserveAppend(static_cast<std::size_t>(base.nodeCount())+
                                 splitRows,
                             static_cast<std::size_t>(base.cellCount())+
                                 splitRows,
                             2*(static_cast<std::size_t>(base.cellCount())+
                                splitRows));
        for(Index node=0;node<base.nodeCount();++node)
            (void)output.addNode(base.node(node),base.nodeClass(node));
    }

    std::map<int,std::vector<double>> parameters;
    for(const auto* edge:model.edges()) {
        if(!edge||edge->tag()<0)throw std::invalid_argument("invalid edge");
        parameters.emplace(edge->tag(),edgeParameters(*edge,spec));
    }
    std::map<int,std::size_t> ordinal;
    trace.segments.clear();
    trace.curveNodeIdentities.clear();
    trace.segments.reserve(static_cast<std::size_t>(source->cellCount())+
                           splitRows);
    trace.curveNodeIdentities.reserve(splitRows);
    for(Index cell=0;cell<source->cellCount();++cell) {
        const auto& row=source->cell(cell);
        if(row.type!=mesh::CellType::Line2||row.onEntity.dim!=1)
            throw std::logic_error("edge mesh contains non-Line2 cell");
        const int tag=row.onEntity.tag;
        const auto found=parameters.find(tag);
        if(found==parameters.end())throw std::logic_error("unknown edge tag");
        const std::size_t oi=ordinal[tag]++;
        if(oi+1>=found->second.size())
            throw std::logic_error("edge ordinal exceeds parameters");
        const double lo=found->second[oi],hi=found->second[oi+1];
        std::vector<const DyadicEdgeSplit*> inside;
        const auto collect=[&](std::span<const DyadicEdgeSplit> rows) {
            for(const auto& split:rows)
                if(split.edgeTag==tag&&split.baseSegmentOrdinal==oi&&
                   split.parameter>lo&&split.parameter<hi)
                    inside.push_back(&split);
        };
        collect(splits); collect(proposedSplits);
        std::sort(inside.begin(),inside.end(),[](const auto* a,const auto* b) {
            return std::tie(a->parameter,a->identity)<
                   std::tie(b->parameter,b->identity);
        });
        for(std::size_t i=1;i<inside.size();++i)
            if(inside[i-1]->parameter==inside[i]->parameter)
                throw std::invalid_argument("duplicate edge split");
        struct Leaf {
            std::uint64_t path = 0;
            std::uint32_t level = 0;
            double lo = 0.0;
            double hi = 0.0;
        };
        std::vector<Leaf> leaves{{0,0,lo,hi}};
        std::vector<const DyadicEdgeSplit*> byLevel=inside;
        std::sort(byLevel.begin(),byLevel.end(),[](const auto* a,const auto* b) {
            return std::tie(a->identity.dyadicLevel,
                            a->identity.dyadicNumerator)<
                   std::tie(b->identity.dyadicLevel,
                            b->identity.dyadicNumerator);
        });
        for(const auto* split:byLevel) {
            const auto& id=split->identity;
            if(id.domain!=SurfaceNodeIdentityDomain::RefinedCurve||
               id.parentArtifact!=entityStamp(tag)||
               id.parentOrdinal!=oi||id.dyadicLevel==0||
               id.dyadicLevel>63||(id.dyadicNumerator&1U)==0)
                throw std::invalid_argument("invalid curve split identity");
            const std::uint32_t parentLevel=id.dyadicLevel-1U;
            const std::uint64_t parentPath=id.dyadicNumerator>>1U;
            const auto leaf=std::find_if(leaves.begin(),leaves.end(),
                [&](const Leaf& row) {
                    return row.level==parentLevel&&row.path==parentPath;
                });
            if(leaf==leaves.end()||std::midpoint(leaf->lo,leaf->hi)!=
                                      split->parameter)
                throw std::invalid_argument("non-canonical curve split path");
            const auto position=static_cast<std::size_t>(leaf-leaves.begin());
            const Leaf left{parentPath*2U,id.dyadicLevel,
                            leaf->lo,split->parameter};
            const Leaf right{parentPath*2U+1U,id.dyadicLevel,
                             split->parameter,leaf->hi};
            *leaf=left;
            leaves.insert(leaves.begin()+static_cast<std::ptrdiff_t>(position+1),
                          right);
        }
        std::sort(leaves.begin(),leaves.end(),[](const Leaf& a,const Leaf& b) {
            return a.lo<b.lo;
        });
        if(leaves.size()!=inside.size()+1)
            throw std::logic_error("curve split leaf mismatch");

        const Index* sourceNodes=source->cellNodes(cell);
        Index previous=sourceNodes[0];
        for(std::size_t i=0;i<inside.size();++i) {
            const auto* split=inside[i];
            const Index middle=output.addNode(split->coordinate,{1,tag});
            trace.curveNodeIdentities.push_back({middle,split->identity});
            (void)output.addCell(mesh::CellType::Line2,{previous,middle},{1,tag});
            RefinedSurfaceSegment segment;
            segment.globalNodes={previous,middle};
            segment.edgeTag=tag;
            segment.owningCurve=entityStamp(tag);
            segment.baseSegmentOrdinal=oi;
            segment.dyadicPath=leaves[i].path;
            segment.dyadicLevel=leaves[i].level;
            segment.parameterLo=leaves[i].lo;
            segment.parameterHi=leaves[i].hi;
            trace.segments.push_back(segment);
            previous=middle;
        }
        if(splitRows!=0)
            (void)output.addCell(mesh::CellType::Line2,
                                 {previous,sourceNodes[1]},{1,tag});
        RefinedSurfaceSegment tail;
        tail.globalNodes={previous,sourceNodes[1]};
        tail.edgeTag=tag;
        tail.owningCurve=entityStamp(tag);
        tail.baseSegmentOrdinal=oi;
        tail.dyadicPath=leaves.back().path;
        tail.dyadicLevel=leaves.back().level;
        tail.parameterLo=leaves.back().lo;
        tail.parameterHi=leaves.back().hi;
        trace.segments.push_back(tail);
    }
    for(const auto& [tag,count]:ordinal) {
        const auto& values=parameters.at(tag);
        if(count+1!=values.size())
            throw std::logic_error("edge parameter count mismatch");
    }
    trace.surfaceNodes=static_cast<std::uint64_t>(output.nodeCount());
    trace.surfaceSegments=static_cast<std::uint64_t>(output.cellCount());
}

}  // namespace t2c_detail

BoundaryConformingV4Outcome meshModelBoundaryConformingV4(
    const geo::Model& model,const sizing::MeshSpec& spec,
    const BoundaryConformingV4Request& request,mesh::Mesh& destination,
    const VolumeBoundaryOwnershipProvider& ownershipProvider,
    const ProtectedFeaturePlacementProvider& placementProvider) noexcept {
    BoundaryConformingV4Outcome out;
    std::unique_ptr<BoundaryConformingV4Result> result;
    std::vector<SurfaceCacheFaceReportRow> cacheProcessedStorage;
    std::vector<FacePersistentId> cacheUnattemptedStorage;
    std::vector<FacePersistentId> cacheCanonicalFaces;
    bool cacheEnvelopePrepared=false;
    BoundaryConformingV4Phase activePhase=BoundaryConformingV4Phase::Preflight;
    try {
        const auto sizeMax=static_cast<std::uint64_t>(
            (std::numeric_limits<std::size_t>::max)());
        const bool facetConflictEnabled=request.facetConflictMode==
            FacetConflictRefinementMode::BoundedCanonicalSccSubdivision;
        const bool facetConflictModeValid=facetConflictEnabled||
            request.facetConflictMode==FacetConflictRefinementMode::Refuse;
        const bool mixedPresentEnabled=request.recoverySchedulingMode==
            TetRecoverySchedulingMode::MixedExactPresentBatching;
        const bool recoverySchedulingModeValid=mixedPresentEnabled||
            request.recoverySchedulingMode==
                TetRecoverySchedulingMode::AcceptedScheduler;
        const bool directSinglePass=request.executionStrategy==
            ProtectedFeatureExecutionStrategy::
                SinglePassInterfacePrefillDirect;
        const bool preparedConstraintSinglePass=request.executionStrategy==
            ProtectedFeatureExecutionStrategy::
                SinglePassPreparedConstraintPrefill;
        const bool bypassCacheSinglePass=directSinglePass||
            preparedConstraintSinglePass;
        const bool interfacePrefillEnabled=directSinglePass||
            preparedConstraintSinglePass||
            request.executionStrategy==
                ProtectedFeatureExecutionStrategy::SinglePassInterfacePrefill;
        const bool executionStrategyValid=interfacePrefillEnabled||
            request.executionStrategy==
                ProtectedFeatureExecutionStrategy::IterativeExact;
        // Single-pass prefill is deliberately independent of the iterative
        // SCC facet subdivision route.  Combining both would make the
        // one-source-attempt contract ambiguous.
        const bool executionStrategyCompatible=
            !interfacePrefillEnabled||!facetConflictEnabled;
        const std::uint32_t maximumPrefillU=
            request.interfacePrefill.maximumPlannedU();
        const std::uint32_t maximumPrefillV=
            request.interfacePrefill.maximumPlannedV();
        const bool prefillGridFits=!interfacePrefillEnabled||
            (maximumPrefillV!=0&&maximumPrefillU<=
                (std::numeric_limits<std::uint64_t>::max)()/
                    maximumPrefillV);
        const std::uint64_t prefillPointsPerFace=
            interfacePrefillEnabled&&prefillGridFits
                ? static_cast<std::uint64_t>(maximumPrefillU)*
                    maximumPrefillV:0U;
        const bool prefillLimitsFit=!interfacePrefillEnabled||
            (request.interfacePrefill.positive()&&prefillGridFits&&
             prefillPointsPerFace<=request.interfacePrefill.maxPoints&&
             request.interfacePrefill.maxInterfaceFaces<=sizeMax&&
             request.interfacePrefill.maxPoints<=sizeMax&&
             request.interfacePrefill.maxPoints<=
                 request.limits.maxAddedNodes);
        const bool facetLimitFits=facetConflictEnabled&&
            request.facetConflictLimits.maxSccMembers!=0&&
            request.facetConflictLimits.maxSccMembers<=sizeMax&&
            request.facetConflictLimits.maxFacetInteriorPoints<=sizeMax&&
            request.facetConflictLimits.maxSccMembers<=
                sizeMax/request.facetConflictLimits.maxSccMembers;
        const std::uint64_t facetRows=facetConflictEnabled
            ? request.facetConflictLimits.maxFacetInteriorPoints:0U;
        const std::uint64_t prefillRows=interfacePrefillEnabled
            ? request.interfacePrefill.maxPoints:0U;
        const auto metricMax=(std::numeric_limits<std::uint64_t>::max)();
        const bool combinedLimitFits=
            request.limits.maxEdgeSplits<=metricMax-facetRows&&
            request.limits.maxEdgeSplits+facetRows<=metricMax-prefillRows;
        if(!request.limits.positive()||
           !request.surfaceCacheLimits.positive()||
           !facetConflictModeValid||
           !recoverySchedulingModeValid||
           !executionStrategyValid||!executionStrategyCompatible||
           !prefillLimitsFit||
           (facetConflictEnabled&&
                (!request.facetConflictLimits.positive()||
                 !facetLimitFits))||
           !combinedLimitFits||
           request.expectedPlacementProviderStamp.persistentId==0||
           request.limits.maxEdgeSplits>sizeMax||
           request.limits.maxAddedNodes>sizeMax||
           request.limits.maxBatchSize>sizeMax||
           request.limits.maxFrontierSupportRows>sizeMax||
           request.limits.maxFrontierDependencyRows>sizeMax||
           request.limits.maxFrontierRows>sizeMax||
           request.surfaceCacheLimits.maxFaces>sizeMax||
           request.surfaceCacheLimits.maxCachedNodeRows>sizeMax||
           request.surfaceCacheLimits.maxCachedCellRows>sizeMax||
           request.surfaceCacheLimits.maxCachedConstraintRows>sizeMax||
           request.surfaceCacheLimits.maxCachedMappingRows>sizeMax||
           request.surfaceCacheLimits.maxDirtyFacesPerBatch>sizeMax||
           request.surfaceCacheLimits.maxSurfaceCacheAccountedBytes>sizeMax) {
            out.status=BoundaryConformingV4Status::Rejected;
            out.reason=BoundaryConformingV4Reason::InvalidRequest;
            out.phase=BoundaryConformingV4Phase::Preflight;
            return out;
        }
        const std::uint64_t maximumRefinementRows=
            request.limits.maxEdgeSplits+facetRows+prefillRows;
        if(maximumRefinementRows>sizeMax) {
            out.status=BoundaryConformingV4Status::Rejected;
            out.reason=BoundaryConformingV4Reason::InvalidRequest;
            out.phase=BoundaryConformingV4Phase::Preflight;
            return out;
        }
        {
#if defined(PAACS_TOPO_TESTING)
            const t2c_detail::SurfaceCacheAllocationZoneScope allocationZone(
                t2c_detail::SurfaceCacheAllocationZone::DiagnosticEnvelope);
#endif
            result=std::make_unique<BoundaryConformingV4Result>();
            result->executionStrategy=request.executionStrategy;
        }
        std::size_t faceCount=0;
        {
#if defined(PAACS_TOPO_TESTING)
            const t2c_detail::SurfaceCacheAllocationZoneScope allocationZone(
                t2c_detail::SurfaceCacheAllocationZone::DiagnosticEnvelope);
#endif
            auto reportFaces=model.faces();
            faceCount=reportFaces.size();
            std::sort(reportFaces.begin(),reportFaces.end(),
                [](const geo::GFace* a,const geo::GFace* b) {
                    return (a?a->tag():-1)<(b?b->tag():-1);
                });
            if(faceCount>request.surfaceCacheLimits.maxFaces) {
                result->surfaceCacheReport.processed.reserve(faceCount);
                result->surfaceCacheReport.reportComplete=true;
                for(const auto* face:reportFaces) {
                    SurfaceCacheFaceReportRow row;
                    row.facePersistentId=face&&face->tag()>=0
                        ? static_cast<FacePersistentId>(face->tag())+1U : 0;
                    row.disposition=SurfaceCacheFaceDisposition::Stopped;
                    row.surfaceStatus=SurfaceFaceStatus::Failed;
                    row.originalSurfaceFailure.reason=
                        SurfaceReason::InternalError;
                    result->surfaceCacheReport.processed.push_back(row);
                }
                out.status=BoundaryConformingV4Status::Incomplete;
                out.reason=BoundaryConformingV4Reason::
                    SurfaceCacheBudgetReached;
                out.phase=BoundaryConformingV4Phase::SurfaceCachePreparation;
                out.retainedResult=std::move(result);
                return out;
            }
            cacheProcessedStorage.resize(faceCount);
            cacheUnattemptedStorage.resize(faceCount);
            cacheCanonicalFaces.resize(faceCount);
            result->surfaceCacheReport.processed.reserve(faceCount);
            result->surfaceCacheReport.unattempted.reserve(faceCount);
            for(std::size_t i=0;i<faceCount;++i) {
                const auto* face=reportFaces[i];
                cacheUnattemptedStorage[i]=face&&face->tag()>=0
                    ? static_cast<FacePersistentId>(face->tag())+1U : 0;
                cacheCanonicalFaces[i]=cacheUnattemptedStorage[i];
            }
            cacheEnvelopePrepared=true;
        }
        t2c_detail::SurfaceCacheDiagnosticChannel cacheChannel{
            cacheProcessedStorage,cacheUnattemptedStorage};
        const auto retainWholeCacheStop=[&](SurfaceFaceStatus status,
                const SurfaceFailure& failure) noexcept {
            auto& report=result->surfaceCacheReport;
            report.processed.clear();
            report.unattempted.clear();
            report.hasFailedFace=false;
            report.failedFacePersistentId=0;
            report.reportComplete=true;
            for(const auto face:cacheUnattemptedStorage)
                report.processed.push_back({face,
                    SurfaceCacheFaceDisposition::Stopped,status,failure});
        };
        const auto consumeWork=[&](std::uint64_t amount) noexcept {
            if(amount>request.limits.maxCoordinatorWorkUnits-
                      result->metrics.coordinatorWorkUnits)return false;
            result->metrics.coordinatorWorkUnits+=amount;
            std::uint64_t measured=result->metrics.nestedConstructionWorkUnits;
            if(!checkedMetricAdd(measured,
                    result->metrics.nestedRecoveryWorkUnits)||
               !checkedMetricAdd(measured,
                                 result->metrics.coordinatorWorkUnits)||
               !checkedMetricAdd(measured,
                    result->metrics.surfaceCachePreparationWorkUnits)||
               !checkedMetricAdd(measured,
                    result->metrics.surfaceCacheReplayWorkUnits)||
               !checkedMetricAdd(measured,
                    result->metrics.surfaceCacheSetAuditWorkUnits))
                return false;
            result->metrics.totalMeasuredWorkUnits=measured;
            return true;
        };
        const auto consumeFacetConflictWork=[&](std::uint64_t amount) noexcept {
            if(!facetConflictEnabled)return false;
            if(result->metrics.facetConflictCoordinatorWorkUnits>
                   request.facetConflictLimits.maxCoordinatorWorkUnits||
               amount>request.facetConflictLimits.maxCoordinatorWorkUnits-
                    result->metrics.facetConflictCoordinatorWorkUnits)
                return false;
            if(!consumeWork(amount))return false;
            result->metrics.facetConflictCoordinatorWorkUnits+=amount;
            return true;
        };
        const auto consumeInterfacePrefillWork=
            [&](std::uint64_t amount) noexcept {
                if(!interfacePrefillEnabled)return false;
                if(!consumeWork(amount))return false;
                return checkedMetricAdd(
                    result->metrics.interfacePrefillWorkUnits,amount);
            };
        if(!consumeWork(1)) {
            SurfaceFailure failure;
            failure.reason=SurfaceReason::InternalError;
            retainWholeCacheStop(SurfaceFaceStatus::Failed,failure);
            out.status=BoundaryConformingV4Status::Incomplete;
            out.reason=BoundaryConformingV4Reason::SubdivisionBudgetReached;
            out.phase=BoundaryConformingV4Phase::Preflight;
            out.retainedResult=std::move(result);
            return out;
        }
        if(facetConflictEnabled) {
            const auto extraRows=static_cast<std::size_t>(
                request.facetConflictLimits.maxFacetInteriorPoints);
            const auto memberRows=static_cast<std::size_t>(
                request.facetConflictLimits.maxSccMembers);
            const auto dependencyRows=memberRows*memberRows;
            const auto edgeBatch=static_cast<std::size_t>((std::min)(
                request.limits.maxBatchSize,request.limits.maxEdgeSplits));
            const auto facetBatch=(std::max)(edgeBatch,memberRows);
            const auto extraBatch=facetBatch-edgeBatch;
            std::size_t lowerBound=0;
            using SemanticStateValue=
                std::pair<SurfaceNodeIdentity,std::uint64_t>;
            if(!checkedSizeBytes(lowerBound,extraRows,
                                 sizeof(t2c_detail::PrivateFacePoint))||
               !checkedSizeBytes(lowerBound,extraRows,
                                 sizeof(ProtectedFeatureChildIdentity))||
               !checkedSizeBytes(lowerBound,extraRows,
                                 sizeof(SemanticStateValue))||
               !checkedSizeBytes(lowerBound,memberRows,
                                 sizeof(bnd::VolumeFacetKey))||
               !checkedSizeBytes(lowerBound,dependencyRows,
                    sizeof(std::array<bnd::VolumeFacetKey,2>))||
               !checkedSizeBytes(lowerBound,extraBatch,
                                 sizeof(std::size_t))||
               !checkedSizeBytes(lowerBound,extraBatch,
                                 sizeof(PlannedFeaturePlacement))||
               !checkedSizeBytes(lowerBound,extraBatch,
                                 sizeof(PreparedFeaturePlacement))||
               !checkedSizeBytes(lowerBound,extraBatch,
                                 sizeof(t2c_detail::DyadicEdgeSplit))||
               !checkedSizeBytes(lowerBound,extraBatch,
                                 sizeof(t2c_detail::PrivateFacePoint))||
               !checkedSizeBytes(lowerBound,extraBatch,
                                 sizeof(ProtectedFeatureChildIdentity))||
               !checkedSizeBytes(lowerBound,extraBatch,
                                 sizeof(SemanticStateValue))||
               lowerBound>request.limits.maxAccountedBytes) {
                SurfaceFailure failure;
                failure.reason=SurfaceReason::InternalError;
                retainWholeCacheStop(SurfaceFaceStatus::Failed,failure);
                out.status=BoundaryConformingV4Status::Incomplete;
                out.reason=BoundaryConformingV4Reason::
                    FacetConflictBudgetReached;
                out.phase=BoundaryConformingV4Phase::Preflight;
                out.retainedResult=std::move(result);
                return out;
            }
        }
        if(interfacePrefillEnabled) {
            const auto rows=static_cast<std::size_t>(
                request.interfacePrefill.maxPoints);
            std::size_t lowerBound=0;
            using SemanticStateValue=
                std::pair<SurfaceNodeIdentity,std::uint64_t>;
            if(!checkedSizeBytes(lowerBound,rows,
                                 sizeof(t2c_detail::PrivateFacePoint))||
               !checkedSizeBytes(lowerBound,rows,
                                 sizeof(ProtectedFeatureChildIdentity))||
               !checkedSizeBytes(lowerBound,rows,
                                 sizeof(SemanticStateValue))||
               lowerBound>request.limits.maxAccountedBytes) {
                SurfaceFailure failure;
                failure.reason=SurfaceReason::InternalError;
                retainWholeCacheStop(SurfaceFaceStatus::Failed,failure);
                out.status=BoundaryConformingV4Status::Incomplete;
                out.reason=BoundaryConformingV4Reason::
                    InterfacePrefillBudgetReached;
                out.phase=BoundaryConformingV4Phase::InterfacePrefill;
                out.retainedResult=std::move(result);
                return out;
            }
        }
        std::vector<t2c_detail::DyadicEdgeSplit> splits;
        std::vector<t2c_detail::PrivateFacePoint> facePoints;
#if defined(PAACS_TOPO_TESTING)
        {
            const t2c_detail::SurfaceCacheAllocationZoneScope allocationZone(
                facetConflictEnabled
                    ? t2c_detail::SurfaceCacheAllocationZone::
                        FacetConflictPreparedStateStorage
                    : t2c_detail::SurfaceCacheAllocationZone::
                        ProposedStateIntegration);
#endif
            splits.reserve(static_cast<std::size_t>(
                request.limits.maxEdgeSplits));
            result->subdivisionLedger.reserve(
                static_cast<std::size_t>(maximumRefinementRows));
            facePoints.reserve(static_cast<std::size_t>(
                maximumRefinementRows));
#if defined(PAACS_TOPO_TESTING)
        }
#endif
        const auto frontierRowCapacity=static_cast<std::size_t>(
            request.limits.maxFrontierRows);
        const auto frontierSupportCapacity=static_cast<std::size_t>(
            request.limits.maxFrontierSupportRows);
        const auto frontierDependencyCapacity=static_cast<std::size_t>(
            request.limits.maxFrontierDependencyRows);
        std::vector<bnd::TetFacetSetFrontierRow> frontierRows(
            frontierRowCapacity);
        std::vector<topo::CellHandle> frontierSupport(
            frontierSupportCapacity);
        std::vector<std::size_t> frontierDependencies(
            frontierDependencyCapacity);
        bnd::TetFacetSetFrontierSink frontier{
            frontierRows,frontierSupport,frontierDependencies};
        const auto facetMemberCapacity=facetConflictEnabled
            ? static_cast<std::size_t>(
                request.facetConflictLimits.maxSccMembers):0U;
        const auto facetDependencyCapacity=facetConflictEnabled
            ? facetMemberCapacity*facetMemberCapacity:0U;
        std::vector<bnd::VolumeFacetKey> facetConflictMembers;
        std::vector<std::array<bnd::VolumeFacetKey,2>>
            facetConflictDependencies;
        std::vector<FacetConflictRankedFacetScratch> facetResolverScratch;
        std::size_t facetStorageBytes=0;
        if(facetConflictEnabled&&
           (!checkedSizeBytes(facetStorageBytes,facetMemberCapacity,
                              sizeof(bnd::VolumeFacetKey))||
            !checkedSizeBytes(facetStorageBytes,facetDependencyCapacity,
                              sizeof(std::array<bnd::VolumeFacetKey,2>))||
            facetMemberCapacity>facetConflictMembers.max_size()||
            facetDependencyCapacity>facetConflictDependencies.max_size()||
            facetStorageBytes>request.limits.maxAccountedBytes)) {
            SurfaceFailure failure;
            failure.reason=SurfaceReason::InternalError;
            retainWholeCacheStop(SurfaceFaceStatus::Failed,failure);
            out.status=BoundaryConformingV4Status::Incomplete;
            out.reason=BoundaryConformingV4Reason::FacetConflictBudgetReached;
            out.phase=BoundaryConformingV4Phase::Preflight;
            out.retainedResult=std::move(result);
            return out;
        }
#if defined(PAACS_TOPO_TESTING)
        {
            const t2c_detail::SurfaceCacheAllocationZoneScope allocationZone(
                t2c_detail::SurfaceCacheAllocationZone::FacetConflictStorage);
#endif
            facetConflictMembers.resize(facetMemberCapacity);
            facetConflictDependencies.resize(facetDependencyCapacity);
#if defined(PAACS_TOPO_TESTING)
        }
#endif
        bnd::TetFacetConflictCycleSink facetConflictCycle{
            facetConflictMembers,facetConflictDependencies};
        bnd::TetFacetSetMixedPresentBatchSink mixedPresentBatch;
        using SemanticState=std::pair<SurfaceNodeIdentity,std::uint64_t>;
        std::vector<SemanticState> semanticStates;
        semanticStates.reserve(static_cast<std::size_t>(
            maximumRefinementRows));
#if defined(PAACS_TOPO_TESTING)
        if(facetConflictRepeatedIdentityForTesting)
            semanticStates.push_back(
                {*facetConflictRepeatedIdentityForTesting,0});
#endif
        const auto edgeBatchCapacity=static_cast<std::size_t>((std::min)(
            request.limits.maxBatchSize,request.limits.maxEdgeSplits));
        const auto batchCapacity=(std::max)(edgeBatchCapacity,
                                            facetMemberCapacity);
        std::vector<std::size_t> selectedFrontierRows;
        std::vector<PlannedFeaturePlacement> planned;
        std::vector<PreparedFeaturePlacement> prepared;
        std::vector<t2c_detail::DyadicEdgeSplit> pendingSplits;
        std::vector<t2c_detail::PrivateFacePoint> pendingFacePoints;
        std::vector<ProtectedFeatureChildIdentity> pendingLedger;
        std::vector<SemanticState> pendingSemanticStates;
#if defined(PAACS_TOPO_TESTING)
        {
            const t2c_detail::SurfaceCacheAllocationZoneScope allocationZone(
                facetConflictEnabled
                    ? t2c_detail::SurfaceCacheAllocationZone::
                        FacetConflictPreparedStateStorage
                    : t2c_detail::SurfaceCacheAllocationZone::
                        ProposedStateIntegration);
#endif
            selectedFrontierRows.reserve(batchCapacity);
            planned.reserve(batchCapacity);
            prepared.reserve(batchCapacity);
            pendingSplits.reserve(batchCapacity);
            pendingFacePoints.reserve(batchCapacity);
            pendingLedger.reserve(batchCapacity);
            pendingSemanticStates.reserve(batchCapacity);
#if defined(PAACS_TOPO_TESTING)
        }
#endif
        t2c_detail::SurfaceFaceCacheSession surfaceCache(
            request.surfaceCacheLimits);
        surfaceCache.attachDiagnosticChannel(&cacheChannel);

        // The complete cache/report envelope and every fixed-capacity v4
        // owner are secured before any placement-provider call.  Authority is
        // checked by the per-placement pre/place/result/post protocol.  The
        // facet-conflict opt-in makes no provider call before a fully prepared
        // SCC proposal.  The non-opt-in checkpoint also secures the fast
        // prefill provider before any deterministic interface placement.

        const auto retainCacheReport=[&]() {
#if defined(PAACS_TOPO_TESTING)
            const t2c_detail::SurfaceCacheAllocationZoneScope allocationZone(
                t2c_detail::SurfaceCacheAllocationZone::DiagnosticEnvelope);
#endif
            auto& retained=result->surfaceCacheReport;
            retained.reportComplete=cacheChannel.reportComplete;
            retained.hasFailedFace=cacheChannel.hasFailedFace;
            retained.failedFacePersistentId=
                cacheChannel.failedFacePersistentId;
            retained.processed.assign(
                cacheProcessedStorage.begin(),
                cacheProcessedStorage.begin()+
                    static_cast<std::ptrdiff_t>(cacheChannel.processedUsed));
            retained.unattempted.assign(
                cacheUnattemptedStorage.begin(),
                cacheUnattemptedStorage.begin()+
                    static_cast<std::ptrdiff_t>(cacheChannel.unattemptedUsed));
        };
        const auto retainDirectSurfaceReport=[&](
                const BoundaryConformingV3Outcome& source) noexcept {
            auto& retained=result->surfaceCacheReport;
            retained.processed.clear();
            retained.unattempted.clear();
            retained.reportComplete=false;
            retained.hasFailedFace=false;
            retained.failedFacePersistentId=0;
            if(!source.surfaceOutcome) {
                SurfaceFailure failure;
                failure.reason=SurfaceReason::InternalError;
                for(const auto face:cacheCanonicalFaces)
                    retained.processed.push_back({face,
                        SurfaceCacheFaceDisposition::Stopped,
                        SurfaceFaceStatus::Failed,failure});
                retained.reportComplete=true;
                return;
            }
            const auto& surface=*source.surfaceOutcome;
            retained.reportComplete=surface.reportComplete;
            for(const auto& row:surface.processed) {
                retained.processed.push_back({row.facePersistentId,
                    SurfaceCacheFaceDisposition::DirtyBuild,
                    row.status,row.failure});
                if(row.status!=SurfaceFaceStatus::Ready&&
                   !retained.hasFailedFace) {
                    retained.hasFailedFace=true;
                    retained.failedFacePersistentId=row.facePersistentId;
                }
            }
            for(const auto face:surface.unattempted)
                retained.unattempted.push_back(face);
        };
        const auto accountAttempt=[&](const mesh::Mesh& privateMesh,
                const BoundaryConformingV3Outcome& source,
                std::size_t& accounted) noexcept {
            bool ok=meshAccounted(privateMesh.capacitySnapshot(),accounted)&&
                checkedSizeBytes(accounted,splits.capacity(),
                                 sizeof(t2c_detail::DyadicEdgeSplit))&&
                checkedSizeBytes(accounted,facePoints.capacity(),
                                 sizeof(t2c_detail::PrivateFacePoint))&&
                checkedSizeBytes(accounted,
                    result->subdivisionLedger.capacity(),
                    sizeof(ProtectedFeatureChildIdentity))&&
                checkedSizeBytes(accounted,semanticStates.capacity(),
                                 sizeof(SemanticState))&&
                checkedSizeBytes(accounted,frontierRows.capacity(),
                                 sizeof(bnd::TetFacetSetFrontierRow))&&
                checkedSizeBytes(accounted,frontierSupport.capacity(),
                                 sizeof(topo::CellHandle))&&
                checkedSizeBytes(accounted,frontierDependencies.capacity(),
                                 sizeof(std::size_t))&&
                checkedSizeBytes(accounted,facetConflictMembers.capacity(),
                                 sizeof(bnd::VolumeFacetKey))&&
                checkedSizeBytes(accounted,
                    facetConflictDependencies.capacity(),
                    sizeof(std::array<bnd::VolumeFacetKey,2>))&&
                checkedSizeBytes(accounted,facetResolverScratch.capacity(),
                                 sizeof(FacetConflictRankedFacetScratch))&&
                checkedSizeBytes(accounted,selectedFrontierRows.capacity(),
                                 sizeof(std::size_t))&&
                checkedSizeBytes(accounted,planned.capacity(),
                                 sizeof(PlannedFeaturePlacement))&&
                checkedSizeBytes(accounted,prepared.capacity(),
                                 sizeof(PreparedFeaturePlacement))&&
                checkedSizeBytes(accounted,pendingSplits.capacity(),
                                 sizeof(t2c_detail::DyadicEdgeSplit))&&
                checkedSizeBytes(accounted,pendingFacePoints.capacity(),
                                 sizeof(t2c_detail::PrivateFacePoint))&&
                checkedSizeBytes(accounted,pendingLedger.capacity(),
                                 sizeof(ProtectedFeatureChildIdentity))&&
                checkedSizeBytes(accounted,pendingSemanticStates.capacity(),
                                 sizeof(SemanticState))&&
                checkedSizeBytes(accounted,cacheProcessedStorage.capacity(),
                                 sizeof(SurfaceCacheFaceReportRow))&&
                checkedSizeBytes(accounted,cacheUnattemptedStorage.capacity(),
                                 sizeof(FacePersistentId))&&
                checkedSizeBytes(accounted,cacheCanonicalFaces.capacity(),
                                 sizeof(FacePersistentId))&&
                checkedSizeBytes(accounted,
                    result->surfaceCacheReport.processed.capacity(),
                    sizeof(SurfaceCacheFaceReportRow))&&
                checkedSizeBytes(accounted,
                    result->surfaceCacheReport.unattempted.capacity(),
                    sizeof(FacePersistentId))&&
                checkedSizeAdd(accounted,
                    result->metrics.surfaceCacheAccountedBytes);
            if(source.retainedResult)
                ok=ok&&checkedSizeAdd(accounted,source.retainedResult->
                    memoryMetrics.finalRetainedOverlapBytes);
            result->metrics.peakAccountedBytes=(std::max)(
                result->metrics.peakAccountedBytes,accounted);
            return ok;
        };

        // Preserve the already-established non-opt-in v4 authority
        // checkpoint exactly.  The SCC opt-in instead secures and validates
        // its complete bounded proposal before its first facet-provider
        // checkpoint.
        if(!facetConflictEnabled) {
            if(!consumeWork(1)) {
                SurfaceFailure failure;
                failure.reason=SurfaceReason::InternalError;
                retainWholeCacheStop(SurfaceFaceStatus::Failed,failure);
                out.status=BoundaryConformingV4Status::Incomplete;
                out.reason=BoundaryConformingV4Reason::SubdivisionBudgetReached;
                out.phase=BoundaryConformingV4Phase::Preflight;
                out.retainedResult=std::move(result);
                return out;
            }
            ProtectedFeatureProviderStamp initialStamp{};
            try {
                ++result->metrics.providerStampCalls;
                initialStamp=placementProvider.stamp();
            } catch(const std::bad_alloc&) {
                SurfaceFailure failure;
                failure.reason=SurfaceReason::AllocationFailed;
                retainWholeCacheStop(SurfaceFaceStatus::Failed,failure);
                out.status=BoundaryConformingV4Status::Failed;
                out.reason=BoundaryConformingV4Reason::AllocationFailed;
                out.phase=BoundaryConformingV4Phase::Preflight;
                out.retainedResult=std::move(result);
                return out;
            } catch(...) {
                SurfaceFailure failure;
                failure.reason=SurfaceReason::InternalError;
                retainWholeCacheStop(SurfaceFaceStatus::Failed,failure);
                out.status=BoundaryConformingV4Status::Failed;
                out.reason=BoundaryConformingV4Reason::ProviderFailed;
                out.phase=BoundaryConformingV4Phase::Preflight;
                out.retainedResult=std::move(result);
                return out;
            }
            if(initialStamp!=request.expectedPlacementProviderStamp) {
                SurfaceFailure failure;
                failure.reason=SurfaceReason::InternalError;
                retainWholeCacheStop(SurfaceFaceStatus::Failed,failure);
                out.status=BoundaryConformingV4Status::Failed;
                out.reason=BoundaryConformingV4Reason::ProviderRevisionChanged;
                out.phase=BoundaryConformingV4Phase::Preflight;
                out.retainedResult=std::move(result);
                return out;
            }
        }

        if(interfacePrefillEnabled) {
            activePhase=BoundaryConformingV4Phase::InterfacePrefill;
            const auto stopPrefill=[&](BoundaryConformingV4Status status,
                    BoundaryConformingV4Reason reason,
                    SurfaceReason surfaceReason) noexcept {
                SurfaceFailure failure;
                failure.reason=surfaceReason;
                retainWholeCacheStop(SurfaceFaceStatus::Failed,failure);
                out.status=status;
                out.reason=reason;
                out.phase=BoundaryConformingV4Phase::InterfacePrefill;
                out.retainedResult=std::move(result);
            };
            if(!consumeInterfacePrefillWork(1)) {
                stopPrefill(BoundaryConformingV4Status::Incomplete,
                    BoundaryConformingV4Reason::
                        InterfacePrefillBudgetReached,
                    SurfaceReason::InternalError);
                return out;
            }
            const auto expectedOwnership=
                request.sourceRequest.expectedOwnershipProviderStamp;
            if(ownershipProvider.stamp()!=expectedOwnership) {
                stopPrefill(BoundaryConformingV4Status::Failed,
                    BoundaryConformingV4Reason::
                        InterfacePrefillAuthorityChanged,
                    SurfaceReason::InternalError);
                return out;
            }
            std::vector<const geo::GFace*> orderedFaces;
            try {
                orderedFaces=model.faces();
                std::sort(orderedFaces.begin(),orderedFaces.end(),
                    [](const geo::GFace* a,const geo::GFace* b) {
                        return (a?a->tag():-1)<(b?b->tag():-1);
                    });
            } catch(const std::bad_alloc&) {
                stopPrefill(BoundaryConformingV4Status::Failed,
                    BoundaryConformingV4Reason::AllocationFailed,
                    SurfaceReason::AllocationFailed);
                return out;
            } catch(...) {
                stopPrefill(BoundaryConformingV4Status::Failed,
                    BoundaryConformingV4Reason::InterfacePrefillUnsupported,
                    SurfaceReason::InternalError);
                return out;
            }
            std::uint64_t interfaceFaces=0;
            for(const auto* face:orderedFaces) {
                if(!consumeInterfacePrefillWork(1)) {
                    stopPrefill(BoundaryConformingV4Status::Incomplete,
                        BoundaryConformingV4Reason::
                            InterfacePrefillBudgetReached,
                        SurfaceReason::InternalError);
                    return out;
                }
                if(!face||face->tag()<0)continue;
                const auto faceId=static_cast<FacePersistentId>(
                    static_cast<std::uint32_t>(face->tag()))+1U;
                VolumeFaceOwnership ownership;
                try {
                    ++result->metrics.interfacePrefillOwnershipCalls;
                    ownership=ownershipProvider.ownership(faceId);
                } catch(const std::bad_alloc&) {
                    stopPrefill(BoundaryConformingV4Status::Failed,
                        BoundaryConformingV4Reason::AllocationFailed,
                        SurfaceReason::AllocationFailed);
                    return out;
                } catch(...) {
                    stopPrefill(BoundaryConformingV4Status::Failed,
                        BoundaryConformingV4Reason::ProviderFailed,
                        SurfaceReason::InternalError);
                    return out;
                }
                if(ownershipProvider.stamp()!=expectedOwnership||
                   ownership.providerStamp!=expectedOwnership) {
                    stopPrefill(BoundaryConformingV4Status::Failed,
                        BoundaryConformingV4Reason::
                            InterfacePrefillAuthorityChanged,
                        SurfaceReason::InternalError);
                    return out;
                }
                if(ownership.facePersistentId!=faceId) {
                    stopPrefill(BoundaryConformingV4Status::Failed,
                        BoundaryConformingV4Reason::
                            PlacementContractViolation,
                        SurfaceReason::InternalError);
                    return out;
                }
                if(ownership.availability!=
                       VolumeProviderAvailability::Available||
                   ownership.owners.size()!=2)
                    continue;
                ++interfaceFaces;
                result->metrics.interfacePrefillCandidateFaces=interfaceFaces;
                if(interfaceFaces>
                   request.interfacePrefill.maxInterfaceFaces) {
                    stopPrefill(BoundaryConformingV4Status::Incomplete,
                        BoundaryConformingV4Reason::
                            InterfacePrefillBudgetReached,
                        SurfaceReason::InternalError);
                    return out;
                }
                Vec2 lo{},hi{},metric{};
                double estimatedU=0.0,estimatedV=0.0;
                bool controlRequestsFiner=false;
                try {
                    face->paramBounds(lo,hi);
                    metric=face->metricScale();
                    const Vec2 centre=(lo+hi)*0.5;
                    const auto point=face->pointAt(centre);
                    const auto curvature=face->curvatureAt(centre);
                    double h=detail::evaluatedGeometrySize(
                        spec,point,curvature);
                    const double baseH=h;
                    if(request.meshControls) {
                        h=request.meshControls->effectiveInterfaceSize(
                            face->tag(),ownership.owners[0].regionTag,
                            ownership.owners[1].regionTag,h);
                        controlRequestsFiner=h<baseH;
                    }
                    if(!finite(lo)||!finite(hi)||!finite(metric)||
                       !finite(point)||!std::isfinite(h)||!(h>0.0))
                        throw detail::SizingEvaluationError{};
                    estimatedU=std::abs(hi.u-lo.u)*
                        std::abs(metric.u)/h;
                    estimatedV=std::abs(hi.v-lo.v)*
                        std::abs(metric.v)/h;
                } catch(const std::bad_alloc&) {
                    stopPrefill(BoundaryConformingV4Status::Failed,
                        BoundaryConformingV4Reason::AllocationFailed,
                        SurfaceReason::AllocationFailed);
                    return out;
                } catch(...) {
                    stopPrefill(BoundaryConformingV4Status::Unsupported,
                        BoundaryConformingV4Reason::
                            InterfacePrefillUnsupported,
                        SurfaceReason::InternalError);
                    return out;
                }
                if(!std::isfinite(estimatedU)||
                   !std::isfinite(estimatedV)) {
                    stopPrefill(BoundaryConformingV4Status::Unsupported,
                        BoundaryConformingV4Reason::
                            InterfacePrefillUnsupported,
                        SurfaceReason::InternalError);
                    return out;
                }
                // The base surface builder has already consumed spec.h(), so
                // its density can justify the historical skip.  A finer
                // entity control has not been consumed by that builder yet;
                // force one-pass planning instead of incorrectly declaring
                // the face dense from a target count alone.
                if(!controlRequestsFiner&&
                   estimatedU>=
                       request.interfacePrefill.minimumEstimatedIntervals&&
                   estimatedV>=
                       request.interfacePrefill.minimumEstimatedIntervals) {
                    ++result->metrics.interfacePrefillDenseFacesSkipped;
                    continue;
                }
                if(controlRequestsFiner)
                    ++result->metrics.interfacePrefillControlForcedFaces;
                std::uint32_t plannedU=request.interfacePrefill.interiorU;
                std::uint32_t plannedV=request.interfacePrefill.interiorV;
                if(request.interfacePrefill.densityMode==
                       InterfacePrefillDensityMode::SizingDriven) {
                    const auto derived=[&](double estimate,
                            std::uint32_t minimum,
                            std::uint32_t maximum) noexcept {
                        const double scaled=std::ceil(
                            estimate*request.interfacePrefill.sizingScale)-1.0;
                        const double bounded=std::clamp(
                            scaled,static_cast<double>(minimum),
                            static_cast<double>(maximum));
                        return static_cast<std::uint32_t>(bounded);
                    };
                    plannedU=derived(estimatedU,
                        request.interfacePrefill.minimumInteriorU,
                        request.interfacePrefill.maximumInteriorU);
                    plannedV=derived(estimatedV,
                        request.interfacePrefill.minimumInteriorV,
                        request.interfacePrefill.maximumInteriorV);
                    ++result->metrics.interfacePrefillSizingDrivenFaces;
                }
                result->metrics.interfacePrefillMaximumPlannedU=std::max(
                    result->metrics.interfacePrefillMaximumPlannedU,
                    static_cast<std::uint64_t>(plannedU));
                result->metrics.interfacePrefillMaximumPlannedV=std::max(
                    result->metrics.interfacePrefillMaximumPlannedV,
                    static_cast<std::uint64_t>(plannedV));
                ++result->metrics.interfacePrefillFaces;
                for(std::uint32_t iu=1;iu<=plannedU;++iu) {
                    for(std::uint32_t iv=1;iv<=plannedV;++iv) {
                        if(result->metrics.interfacePrefillPoints>=
                               request.interfacePrefill.maxPoints||
                           result->metrics.uniquePlacements>=
                               request.limits.maxAddedNodes) {
                            stopPrefill(BoundaryConformingV4Status::Incomplete,
                                BoundaryConformingV4Reason::
                                    InterfacePrefillBudgetReached,
                                SurfaceReason::InternalError);
                            return out;
                        }
                        const double uFraction=static_cast<double>(iu)/
                            static_cast<double>(
                                plannedU+1U);
                        const double vFraction=static_cast<double>(iv)/
                            static_cast<double>(
                                plannedV+1U);
                        const Vec2 uv{
                            lo.u+(hi.u-lo.u)*uFraction,
                            lo.v+(hi.v-lo.v)*vFraction};
                        bool contained=false;
                        try {
                            contained=finite(uv)&&face->containsParam(uv);
                        } catch(const std::bad_alloc&) {
                            stopPrefill(BoundaryConformingV4Status::Failed,
                                BoundaryConformingV4Reason::AllocationFailed,
                                SurfaceReason::AllocationFailed);
                            return out;
                        } catch(...) {
                            stopPrefill(BoundaryConformingV4Status::Unsupported,
                                BoundaryConformingV4Reason::
                                    InterfacePrefillUnsupported,
                                SurfaceReason::InternalError);
                            return out;
                        }
                        if(!contained) {
                            ++result->metrics.
                                interfacePrefillSkippedParameters;
                            continue;
                        }
                        const auto identity=faceIdentity(
                            SurfaceNodeIdentityDomain::InterfacePrefill,
                            entityStamp(face->tag()),iu,iv,uv);
                        if(!validSurfaceNodeIdentity(identity)) {
                            stopPrefill(BoundaryConformingV4Status::Failed,
                                BoundaryConformingV4Reason::InternalError,
                                SurfaceReason::InternalError);
                            return out;
                        }
                        ProtectedFeaturePlacementRequest placeRequest;
                        placeRequest.kind=
                            ProtectedFeatureEntityKind::FacetInterior;
                        placeRequest.entityStamp=entityStamp(face->tag());
                        placeRequest.requestedIdentity=identity;
                        placeRequest.entityTag=face->tag();
                        placeRequest.faceParameter=uv;
                        if(!consumeInterfacePrefillWork(3)) {
                            stopPrefill(BoundaryConformingV4Status::Incomplete,
                                BoundaryConformingV4Reason::
                                    InterfacePrefillBudgetReached,
                                SurfaceReason::InternalError);
                            return out;
                        }
                        ProtectedFeatureProviderStamp preStamp{};
                        try {
                            ++result->metrics.providerStampCalls;
                            preStamp=placementProvider.stamp();
                        } catch(...) {
                            stopPrefill(BoundaryConformingV4Status::Failed,
                                BoundaryConformingV4Reason::ProviderFailed,
                                SurfaceReason::InternalError);
                            return out;
                        }
                        if(preStamp!=request.expectedPlacementProviderStamp) {
                            stopPrefill(BoundaryConformingV4Status::Failed,
                                BoundaryConformingV4Reason::
                                    ProviderRevisionChanged,
                                SurfaceReason::InternalError);
                            return out;
                        }
                        ProtectedFeaturePlacement placement;
                        try {
                            ++result->metrics.providerPlacementCalls;
                            ++result->metrics.
                                interfacePrefillProviderCalls;
                            placement=placementProvider.placeOnFace(
                                *face,placeRequest);
                        } catch(const std::bad_alloc&) {
                            stopPrefill(BoundaryConformingV4Status::Failed,
                                BoundaryConformingV4Reason::AllocationFailed,
                                SurfaceReason::AllocationFailed);
                            return out;
                        } catch(...) {
                            stopPrefill(BoundaryConformingV4Status::Failed,
                                BoundaryConformingV4Reason::ProviderFailed,
                                SurfaceReason::InternalError);
                            return out;
                        }
                        if(placement.providerStamp!=
                               request.expectedPlacementProviderStamp) {
                            stopPrefill(BoundaryConformingV4Status::Failed,
                                BoundaryConformingV4Reason::
                                    ProviderRevisionChanged,
                                SurfaceReason::InternalError);
                            return out;
                        }
                        if(placement.entityStamp!=placeRequest.entityStamp||
                           placement.requestedIdentity!=identity) {
                            stopPrefill(BoundaryConformingV4Status::Failed,
                                BoundaryConformingV4Reason::
                                    PlacementContractViolation,
                                SurfaceReason::InternalError);
                            return out;
                        }
                        const auto postStamp=placementProvider.stamp();
                        ++result->metrics.providerStampCalls;
                        if(postStamp!=
                           request.expectedPlacementProviderStamp) {
                            stopPrefill(BoundaryConformingV4Status::Failed,
                                BoundaryConformingV4Reason::
                                    ProviderRevisionChanged,
                                SurfaceReason::InternalError);
                            return out;
                        }
                        if(placement.status!=
                           ProtectedFeaturePlacementStatus::Placed) {
                            const auto status=placement.status==
                                    ProtectedFeaturePlacementStatus::Rejected
                                ? BoundaryConformingV4Status::Rejected
                                : placement.status==
                                      ProtectedFeaturePlacementStatus::
                                          Unsupported
                                    ? BoundaryConformingV4Status::Unsupported
                                    : BoundaryConformingV4Status::Failed;
                            const auto reason=placement.status==
                                    ProtectedFeaturePlacementStatus::Rejected
                                ? BoundaryConformingV4Reason::ProviderRejected
                                : placement.status==
                                      ProtectedFeaturePlacementStatus::
                                          Unsupported
                                    ? BoundaryConformingV4Reason::
                                          ProviderUnsupported
                                    : BoundaryConformingV4Reason::
                                          ProviderFailed;
                            stopPrefill(status,reason,
                                SurfaceReason::InternalError);
                            return out;
                        }
                        Vec3 authoritative{};
                        try {
                            authoritative=face->pointAt(uv);
                        } catch(const std::bad_alloc&) {
                            stopPrefill(BoundaryConformingV4Status::Failed,
                                BoundaryConformingV4Reason::AllocationFailed,
                                SurfaceReason::AllocationFailed);
                            return out;
                        } catch(...) {
                            stopPrefill(BoundaryConformingV4Status::Failed,
                                BoundaryConformingV4Reason::ProviderFailed,
                                SurfaceReason::InternalError);
                            return out;
                        }
                        if(!finite(placement.coordinate)||
                           !finite(authoritative)||
                           !exactPoint(placement.coordinate,authoritative)) {
                            stopPrefill(BoundaryConformingV4Status::Failed,
                                BoundaryConformingV4Reason::
                                    PlacementContractViolation,
                                SurfaceReason::InternalError);
                            return out;
                        }
                        const std::uint64_t ordinal=
                            static_cast<std::uint64_t>(
                                result->subdivisionLedger.size()+1U);
                        facePoints.push_back({faceId,uv,
                            placement.coordinate,identity});
                        semanticStates.push_back({identity,ordinal});
                        ProtectedFeatureChildIdentity ledger;
                        ledger.parentEntityPersistentId=
                            identity.parentArtifact.persistentId;
                        ledger.childIdentity=identity;
                        ledger.splitOrdinal=ordinal;
                        ledger.faceTag=face->tag();
                        ledger.faceParameter=uv;
                        ledger.coordinate=placement.coordinate;
                        ledger.classification={2,face->tag()};
                        ledger.providerStamp=placement.providerStamp;
                        result->subdivisionLedger.push_back(ledger);
                        ++result->metrics.facetSplits;
                        ++result->metrics.semanticStates;
                        ++result->metrics.uniquePlacements;
                        ++result->metrics.interfacePrefillPoints;
                    }
                }
            }
            if(ownershipProvider.stamp()!=expectedOwnership) {
                stopPrefill(BoundaryConformingV4Status::Failed,
                    BoundaryConformingV4Reason::
                        InterfacePrefillAuthorityChanged,
                    SurfaceReason::InternalError);
                return out;
            }
        }

        for(std::uint64_t round=0;round<request.limits.maxRounds;++round) {
            activePhase=BoundaryConformingV4Phase::SourceAttempt;
            ++result->metrics.rounds;
            if(mixedPresentEnabled)mixedPresentBatch.reset();
            if(!consumeWork(1)) {
                out.status=BoundaryConformingV4Status::Incomplete;
                out.reason=BoundaryConformingV4Reason::SubdivisionBudgetReached;
                out.phase=BoundaryConformingV4Phase::SourceAttempt;
                out.retainedResult=std::move(result);
                return out;
            }
            mesh::Mesh privateDestination;
            frontier.reset();
            if(facetConflictEnabled)facetConflictCycle.reset();
            t2c_detail::PrivateBuildTrace trace;
            trace.frontierSink=&frontier;
            trace.facetConflictCycleSink=facetConflictEnabled
                ? &facetConflictCycle:nullptr;
            trace.mixedPresentBatchSink=mixedPresentEnabled
                ? &mixedPresentBatch:nullptr;
            trace.surfaceCacheSession=bypassCacheSinglePass?nullptr:&surfaceCache;
            trace.surfaceCacheChannel=bypassCacheSinglePass?nullptr:&cacheChannel;
            trace.preferDirectConstraintInitialization=preparedConstraintSinglePass;
            if(bypassCacheSinglePass)
                ++result->metrics.directSurfaceBuildAttempts;
            trace.proposedEdgeSplits=pendingSplits;
            trace.proposedFacePoints=pendingFacePoints;
            auto source=std::make_unique<BoundaryConformingV3Outcome>(
                t2c_detail::meshModelBoundaryConformingWithProtectedSplits(
                    model,spec,request.sourceRequest,privateDestination,
                    ownershipProvider,splits,facePoints,&trace));
            if(source->surfaceOutcome&&
               source->surfaceOutcome->retainedArtifact) {
                for(const auto& face:
                    source->surfaceOutcome->retainedArtifact->faces)
                    if(face&&face->preparedRecoveryConsumed)
                        ++result->metrics.preparedSurfaceRecoveryFaces;
            }
            if(mixedPresentEnabled) {
                auto& total=result->metrics.mixedPresentScheduling;
                const auto& current=mixedPresentBatch.metrics;
                result->metrics.latestSourceMixedPresentScheduling=current;
                bool mixedComplete=current.complete;
                const auto addGeometryCounters=[&](
                    auto& destination,const auto& source) {
                    return checkedMetricAdd(destination.cellsScanned,
                                            source.cellsScanned)&&
                        checkedMetricAdd(destination.crossedEdgeStars,
                                         source.crossedEdgeStars)&&
                        checkedMetricAdd(destination.closureInsertions,
                                         source.closureInsertions)&&
                        checkedMetricAdd(destination.fixedPointPasses,
                                         source.fixedPointPasses)&&
                        checkedMetricAdd(destination.apexAttempts,
                                         source.apexAttempts)&&
                        checkedMetricAdd(destination.tetTetIntersectionTests,
                                         source.tetTetIntersectionTests);
                };
                const auto addBindingWork=[&](auto& destination,
                                              const auto& source) {
                    bool complete=destination.complete&&source.complete&&
                        checkedMetricAdd(
                            destination.edgeBindingAuditCalls,
                            source.edgeBindingAuditCalls)&&
                        checkedMetricAdd(
                            destination.edgeBindingRowsVisited,
                            source.edgeBindingRowsVisited)&&
                        checkedMetricAdd(
                            destination.edgeAliasRowsVisited,
                            source.edgeAliasRowsVisited)&&
                        checkedMetricAdd(
                            destination.facetBindingAuditCalls,
                            source.facetBindingAuditCalls)&&
                        checkedMetricAdd(
                            destination.facetBindingRowsVisited,
                            source.facetBindingRowsVisited)&&
                        checkedMetricAdd(
                            destination.facetBindingRequestRowsVisited,
                            source.facetBindingRequestRowsVisited)&&
                        checkedMetricAdd(
                            destination.facetBindingDuplicateRowsCompared,
                            source.facetBindingDuplicateRowsCompared)&&
                        checkedMetricAdd(
                            destination.facetBindingReceiptBuilds,
                            source.facetBindingReceiptBuilds)&&
                        checkedMetricAdd(
                            destination.facetBindingReceiptReuses,
                            source.facetBindingReceiptReuses)&&
                        checkedMetricAdd(
                            destination.facetBindingReceiptRowsCopied,
                            source.facetBindingReceiptRowsCopied)&&
                        checkedMetricAdd(
                            destination.facetBindingReceiptLookups,
                            source.facetBindingReceiptLookups)&&
                        checkedMetricAdd(
                            destination.facetBindingReceiptReserveCalls,
                            source.facetBindingReceiptReserveCalls)&&
                        checkedMetricAdd(destination.edgeDataCopyRows,
                                         source.edgeDataCopyRows)&&
                        checkedMetricAdd(destination.faceDataCopyRows,
                                         source.faceDataCopyRows)&&
                        checkedMetricAdd(
                            destination.faceIncidenceRowsCopied,
                            source.faceIncidenceRowsCopied)&&
                        checkedMetricAdd(destination.stateCopies,
                                         source.stateCopies)&&
                        checkedMetricAdd(destination.stateCopyRows,
                                         source.stateCopyRows);
                    destination.complete=complete;
                    return complete;
                };
                for(std::size_t i=0;
                    i<total.edgeAssessmentEventStatusCounts.size();++i)
                    mixedComplete=mixedComplete&&checkedMetricAdd(
                        total.edgeAssessmentEventStatusCounts[i],
                        current.edgeAssessmentEventStatusCounts[i]);
                for(std::size_t i=0;
                    i<total.facetAssessmentEventStatusCounts.size();++i)
                    mixedComplete=mixedComplete&&checkedMetricAdd(
                        total.facetAssessmentEventStatusCounts[i],
                        current.facetAssessmentEventStatusCounts[i]);
                mixedComplete=mixedComplete&&
                    checkedMetricAdd(total.presentBindingBatches,
                        current.presentBindingBatches)&&
                    checkedMetricAdd(total.presentBindingEntities,
                        current.presentBindingEntities)&&
                    checkedMetricAdd(total.stateCopies,current.stateCopies)&&
                    checkedMetricAdd(total.copiedRows,current.copiedRows)&&
                    checkedMetricAdd(total.copiedLogicalBytes,
                        current.copiedLogicalBytes)&&
                    checkedMetricAdd(total.schedulerBindingAuditRows,
                        current.schedulerBindingAuditRows)&&
                    checkedMetricAdd(total.pendingRowsVisited,
                        current.pendingRowsVisited)&&
                    checkedMetricAdd(total.dependencyPairsTested,
                        current.dependencyPairsTested)&&
                    checkedMetricAdd(total.exactPresentSelectedRebuilds,
                        current.exactPresentSelectedRebuilds)&&
                    checkedMetricAdd(total.exactPresentDependencyTests,
                        current.exactPresentDependencyTests)&&
                    checkedMetricAdd(total.nontrivialRecoveries,
                        current.nontrivialRecoveries)&&
                    checkedMetricAdd(
                        total.modeledBatchPreparationAllocationRequests,
                        current.modeledBatchPreparationAllocationRequests)&&
                    checkedMetricAdd(
                        total.modeledBatchPreparationLogicalBytes,
                        current.modeledBatchPreparationLogicalBytes)&&
                    addGeometryCounters(total.assessmentCounters,
                                        current.assessmentCounters)&&
                    addGeometryCounters(total.recoveryCounters,
                                        current.recoveryCounters)&&
                    addBindingWork(total.assessmentBindingWork,
                                   current.assessmentBindingWork)&&
                    addBindingWork(total.recoveryBindingWork,
                                   current.recoveryBindingWork);
                total.complete=total.complete&&mixedComplete;
                if(!mixedComplete)trace.workComplete=false;
            }
            SourceOutcomeRetention sourceRetention{source,out};
            if(bypassCacheSinglePass) {
                retainDirectSurfaceReport(*source);
                result->metrics.directSurfaceBuildFaces=
                    static_cast<std::uint64_t>(
                        result->surfaceCacheReport.processed.size());
            } else {
                retainCacheReport();
            }
            result->metrics.surfaceCacheFacesBuilt=
                trace.surfaceCacheMetrics.facesBuilt;
            result->metrics.surfaceCacheFacesReused=
                trace.surfaceCacheMetrics.facesReused;
            result->metrics.surfaceCacheFacesReplayed=
                trace.surfaceCacheMetrics.facesReplayed;
            result->metrics.surfaceCacheFacesRejected=
                trace.surfaceCacheMetrics.facesRejected;
            result->metrics.surfaceCacheKeyRowsVisited=
                trace.surfaceCacheMetrics.keyRowsVisited;
            result->metrics.surfaceCacheDirtyPropagationVisits=
                trace.surfaceCacheMetrics.dirtyPropagationVisits;
            result->metrics.surfaceCacheSizingReplayCalls=
                trace.surfaceCacheMetrics.sizingReplayCalls;
            result->metrics.surfaceCachePreparationWorkUnits=
                trace.surfaceCacheMetrics.preparationWorkUnits;
            result->metrics.surfaceCacheReplayWorkUnits=
                trace.surfaceCacheMetrics.replayWorkUnits;
            result->metrics.surfaceCacheSetAuditWorkUnits=
                trace.surfaceCacheMetrics.setAuditWorkUnits;
            result->metrics.surfaceCacheCommits=
                trace.surfaceCacheMetrics.cacheCommits;
            result->metrics.surfaceCacheCleanNodeRows=
                trace.surfaceCacheMetrics.cleanNodeRows;
            result->metrics.surfaceCacheDirtyNodeRows=
                trace.surfaceCacheMetrics.dirtyNodeRows;
            result->metrics.surfaceCacheCleanCellRows=
                trace.surfaceCacheMetrics.cleanCellRows;
            result->metrics.surfaceCacheDirtyCellRows=
                trace.surfaceCacheMetrics.dirtyCellRows;
            result->metrics.surfaceCacheCleanConstraintRows=
                trace.surfaceCacheMetrics.cleanConstraintRows;
            result->metrics.surfaceCacheDirtyConstraintRows=
                trace.surfaceCacheMetrics.dirtyConstraintRows;
            result->metrics.surfaceCacheCleanMappingRows=
                trace.surfaceCacheMetrics.cleanMappingRows;
            result->metrics.surfaceCacheDirtyMappingRows=
                trace.surfaceCacheMetrics.dirtyMappingRows;
            result->metrics.surfaceCacheCopiedBytes=
                trace.surfaceCacheMetrics.copiedBytes;
            result->metrics.surfaceCacheHighWaterSlots=
                trace.surfaceCacheMetrics.highWaterSlots;
            result->metrics.surfaceCacheHighWaterRows=
                trace.surfaceCacheMetrics.highWaterRows;
            result->metrics.surfaceCacheAllFaceScansPerDirtyFace=
                trace.surfaceCacheMetrics.allFaceScansPerDirtyFace;
            result->metrics.surfaceCacheAccountedBytes=
                trace.surfaceCacheMetrics.accountedByteHighWater;
            if(!bypassCacheSinglePass&&
               (trace.surfaceCacheDiagnostic.status==
                    t2c_detail::SurfaceFaceCacheStatus::Incomplete||
                (trace.surfaceCacheDiagnostic.status==
                     t2c_detail::SurfaceFaceCacheStatus::Failed&&
                 trace.surfaceCacheDiagnostic.reason!=
                     t2c_detail::SurfaceFaceCacheReason::FaceBuildStopped))) {
                out.status=trace.surfaceCacheDiagnostic.status==
                        t2c_detail::SurfaceFaceCacheStatus::Incomplete
                    ? BoundaryConformingV4Status::Incomplete
                    : BoundaryConformingV4Status::Failed;
                switch(trace.surfaceCacheDiagnostic.reason) {
                    case t2c_detail::SurfaceFaceCacheReason::AuthorityChanged:
                        out.reason=BoundaryConformingV4Reason::
                            SurfaceCacheAuthorityChanged;
                        break;
                    case t2c_detail::SurfaceFaceCacheReason::BudgetReached:
                        out.reason=BoundaryConformingV4Reason::
                            SurfaceCacheBudgetReached;
                        break;
                    case t2c_detail::SurfaceFaceCacheReason::AllocationFailed:
                        out.reason=BoundaryConformingV4Reason::AllocationFailed;
                        break;
                    default:
                        out.reason=BoundaryConformingV4Reason::
                            SurfaceCacheInvalid;
                        break;
                }
                out.phase=trace.surfaceCacheDiagnostic.phase==
                        t2c_detail::SurfaceFaceCachePhase::Preparation||
                    trace.surfaceCacheDiagnostic.phase==
                        t2c_detail::SurfaceFaceCachePhase::Envelope
                    ? BoundaryConformingV4Phase::SurfaceCachePreparation
                    : BoundaryConformingV4Phase::SurfaceCacheAssembly;
                out.lastSourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }
            const bool proposedSurfaceAccepted=bypassCacheSinglePass
                ? source->surfaceOutcome&&
                    source->surfaceOutcome->setStatus==
                        SurfaceSetStatus::Published&&
                    source->surfaceOutcome->reportComplete
                : trace.surfaceCacheDiagnostic.status==
                    t2c_detail::SurfaceFaceCacheStatus::Ready;
            if(source->published()&&!proposedSurfaceAccepted) {
                out.status=BoundaryConformingV4Status::Failed;
                out.reason=BoundaryConformingV4Reason::SurfaceCacheInvalid;
                out.phase=BoundaryConformingV4Phase::SurfaceCacheAssembly;
                out.lastSourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }
            if(proposedSurfaceAccepted) {
                // The proposed rows have passed one complete surface rebuild,
                // cache/set audit, and the downstream v3 immutable preflight.
                // All owning vectors were reserved to their hard limits before
                // the first provider call, so this private commit is
                // nonthrowing.  A stopped/incomplete cache attempt leaves the
                // live split, identity and ledger state untouched.
#if defined(PAACS_TOPO_TESTING)
                const t2c_detail::SurfaceCacheAllocationZoneScope allocationZone(
                    t2c_detail::SurfaceCacheAllocationZone::
                        ProposedStateIntegration);
#endif
                for(const auto& row:pendingSplits)splits.push_back(row);
                for(const auto& row:pendingFacePoints)
                    facePoints.push_back(row);
                for(const auto& row:pendingSemanticStates)
                    semanticStates.push_back(row);
                for(const auto& row:pendingLedger) {
                    result->subdivisionLedger.push_back(row);
                    if(row.childIdentity.domain==
                       SurfaceNodeIdentityDomain::RefinedFacetInterior)
                        ++result->metrics.facetSplits;
                    else
                        ++result->metrics.edgeSplits;
                    ++result->metrics.semanticStates;
                    ++result->metrics.uniquePlacements;
                }
                pendingSplits.clear(); pendingFacePoints.clear();
                pendingSemanticStates.clear(); pendingLedger.clear();
            }
            result->metrics.rebuiltSurfaceNodes=trace.surfaceNodes;
            result->metrics.rebuiltSurfaceSegments=trace.surfaceSegments;
            result->metrics.candidateTetCells=trace.candidateTetCells;
#if defined(PAACS_TOPO_TESTING)
            for(std::size_t i=0;i<trace.phaseNanoseconds.size();++i)
                if(!checkedMetricAdd(
                        result->metrics.nestedPhaseNanoseconds[i],
                        trace.phaseNanoseconds[i])) {
                    trace.workComplete=false;
                    break;
                }
#endif
            if(!trace.workComplete||
               !checkedMetricAdd(result->metrics.nestedConstructionWorkUnits,
                                 trace.constructionWorkUnits)||
               !checkedMetricAdd(result->metrics.nestedRecoveryWorkUnits,
                                 trace.recoveryWorkUnits)) {
                out.status=BoundaryConformingV4Status::Failed;
                out.reason=BoundaryConformingV4Reason::InternalError;
                out.phase=BoundaryConformingV4Phase::SourceAttempt;
                out.lastSourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }
            if(!consumeWork(0)) {
                out.status=BoundaryConformingV4Status::Failed;
                out.reason=BoundaryConformingV4Reason::InternalError;
                out.phase=BoundaryConformingV4Phase::SourceAttempt;
                out.lastSourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }
            if(trace.candidateTetCells>request.limits.maxTetCells) {
                out.status=BoundaryConformingV4Status::Incomplete;
                out.reason=BoundaryConformingV4Reason::SubdivisionBudgetReached;
                out.phase=BoundaryConformingV4Phase::SourceAttempt;
                out.lastSourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }
            if(source->published()) {
                const auto tets=static_cast<std::uint64_t>(
                    privateDestination.countCells(mesh::CellType::Tet4));
                std::size_t accounted=0;
                const bool accountedOk=accountAttempt(
                    privateDestination,*source,accounted);
                if(!accountedOk||tets>request.limits.maxTetCells||
                   accounted>request.limits.maxAccountedBytes) {
                    out.status=BoundaryConformingV4Status::Incomplete;
                    out.reason=BoundaryConformingV4Reason::SubdivisionBudgetReached;
                    out.phase=BoundaryConformingV4Phase::BoundaryAudit;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                const bool needsFacetFinalCheckpoint=
                    result->metrics.facetConflictRounds!=0;
                const bool finalCallCapacityReady=!needsFacetFinalCheckpoint||
                    (result->metrics.facetConflictProviderCalls<
                         request.facetConflictLimits.maxProviderCalls);
                if(!finalCallCapacityReady||
                   !(needsFacetFinalCheckpoint
                        ? consumeFacetConflictWork(1):consumeWork(1))) {
                    out.status=BoundaryConformingV4Status::Incomplete;
                    out.reason=needsFacetFinalCheckpoint
                        ? BoundaryConformingV4Reason::
                            FacetConflictBudgetReached
                        : BoundaryConformingV4Reason::
                            SubdivisionBudgetReached;
                    out.phase=needsFacetFinalCheckpoint
                        ? BoundaryConformingV4Phase::ConstraintRecovery
                        : BoundaryConformingV4Phase::Publication;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                ProtectedFeatureProviderStamp finalPlacementStamp{};
                try {
                    ++result->metrics.providerStampCalls;
                    if(needsFacetFinalCheckpoint)
                        ++result->metrics.facetConflictProviderCalls;
                    finalPlacementStamp=placementProvider.stamp();
                } catch(const std::bad_alloc&) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::AllocationFailed;
                    out.phase=BoundaryConformingV4Phase::Publication;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                } catch(...) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::ProviderFailed;
                    out.phase=BoundaryConformingV4Phase::Publication;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                if(finalPlacementStamp!=
                   request.expectedPlacementProviderStamp) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::
                        ProviderRevisionChanged;
                    out.phase=BoundaryConformingV4Phase::Publication;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                result->metrics.runComplete=true;
                result->sourceOutcome=std::move(source);
#if defined(PAACS_TOPO_TESTING)
                {
                    const t2c_detail::SurfaceCacheAllocationZoneScope
                        allocationZone(
                            t2c_detail::SurfaceCacheAllocationZone::CommitTail);
#endif
                    destination.swap(privateDestination);
                    out.status=BoundaryConformingV4Status::Published;
                    out.reason=BoundaryConformingV4Reason::None;
                    out.phase=BoundaryConformingV4Phase::Publication;
                    out.destinationUnchanged=false;
                    out.retainedResult=std::move(result);
#if defined(PAACS_TOPO_TESTING)
                }
#endif
                return out;
            }

            if(interfacePrefillEnabled) {
                out.status=BoundaryConformingV4Status::Incomplete;
                out.reason=BoundaryConformingV4Reason::
                    SinglePassRecoveryIncomplete;
                out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                result->sourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }

            const bool facetConflictCaptured=
                facetConflictEnabled&&facetConflictCycle.captured;
            if(facetConflictCaptured) {
                activePhase=BoundaryConformingV4Phase::ConstraintRecovery;
                if(!checkedMetricAdd(result->metrics.facetConflictSccMembers,
                                     facetConflictCycle.membersRequired)||
                   !checkedMetricAdd(
                       result->metrics.facetConflictDependencyRows,
                       facetConflictCycle.dependenciesRequired)) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::InternalError;
                    out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                if(!facetConflictCycle.complete) {
                    out.status=BoundaryConformingV4Status::Incomplete;
                    out.reason=BoundaryConformingV4Reason::
                        FacetConflictBudgetReached;
                    out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
            }
            if(!frontier.complete) {
                out.status=BoundaryConformingV4Status::Incomplete;
                out.reason=BoundaryConformingV4Reason::FrontierCapacityReached;
                out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                result->sourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }
            if(frontier.rowsUsed==0) {
                out.status=sourceStatus(source->status);
                out.reason=sourceReason(source->status);
                out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                result->sourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }
            if(!checkedMetricAdd(result->metrics.frontierRows,
                                 frontier.rowsUsed)) {
                out.status=BoundaryConformingV4Status::Failed;
                out.reason=BoundaryConformingV4Reason::InternalError;
                out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                result->sourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }
            selectedFrontierRows.clear();
            planned.clear();
            prepared.clear();
            if(!facetConflictCaptured) {
            std::vector<std::size_t> canonicalRows(frontier.rowsUsed);
            std::iota(canonicalRows.begin(),canonicalRows.end(),0);
            std::sort(canonicalRows.begin(),canonicalRows.end(),
                [&](std::size_t a,std::size_t b) {
                    const auto& left=frontier.rows[a];
                    const auto& right=frontier.rows[b];
                    if(left.kind!=right.kind)return left.kind<right.kind;
                    return left.kind==bnd::TetFacetSetWorkKind::Edge
                        ? left.edgeKey<right.edgeKey
                        : left.facetKey<right.facetKey;
                });
            const auto refinable=[](const bnd::TetFacetSetFrontierRow& row) {
                return row.kind==bnd::TetFacetSetWorkKind::Edge&&
                       row.edgeOutcome.status==
                           bnd::TetEdgeRecoveryStatus::Refused&&
                       row.edgeOutcome.reason==
                           bnd::TetEdgeRecoveryReason::NonCycleLink&&
                       row.edgeOutcome.phase==
                           bnd::TetEdgeRecoveryPhase::Discovery;
            };
            const auto firstNonrefinable=std::find_if(
                canonicalRows.begin(),canonicalRows.end(),
                [&](std::size_t i) { return !refinable(frontier.rows[i]); });
            if(firstNonrefinable!=canonicalRows.end()) {
                const auto& failure=source->failure;
                out.status=sourceStatus(source->status);
                out.reason=(!splits.empty()||!facePoints.empty())
                    ? (failure.workKind==VolumeWorkKind::ArtifactEdge
                        ? BoundaryConformingV4Reason::ChildEdgeRecoveryRefused
                        : failure.workKind==VolumeWorkKind::Facet
                            ? BoundaryConformingV4Reason::ChildFacetRecoveryRefused
                            : sourceReason(source->status))
                    : sourceReason(source->status);
                out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                result->sourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }

            std::map<std::array<Index,2>,const t2c_detail::RefinedSurfaceSegment*>
                curveByEdge;
            for(const auto& segment:trace.segments) {
                auto key=segment.globalNodes;
                if(key[1]<key[0])std::swap(key[0],key[1]);
                const auto [it,inserted]=curveByEdge.emplace(key,&segment);
                if(!inserted&&it->second!=&segment) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::InternalError;
                    out.phase=BoundaryConformingV4Phase::EdgeSelection;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
            }
            const auto surfaceSupport=buildSurfaceSupportIndex(*source);
            if(!surfaceSupport) {
                out.status=BoundaryConformingV4Status::Failed;
                out.reason=BoundaryConformingV4Reason::InternalError;
                out.phase=BoundaryConformingV4Phase::EdgeSelection;
                result->sourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }
            std::uint64_t selectionWork=surfaceSupport->workUnits;
            if(!checkedMetricAdd(selectionWork,trace.segments.size())||
               !checkedMetricAdd(selectionWork,canonicalRows.size())||
               !consumeWork(selectionWork)) {
                out.status=BoundaryConformingV4Status::Incomplete;
                out.reason=BoundaryConformingV4Reason::SubdivisionBudgetReached;
                out.phase=BoundaryConformingV4Phase::EdgeSelection;
                result->sourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }

            for(const auto rowOrdinal:canonicalRows) {
                if(selectedFrontierRows.size()>=batchCapacity) {
                    ++result->metrics.deferredFrontierRows;
                    continue;
                }
                const auto& frontierRow=frontier.rows[rowOrdinal];
                bool overlap=false;
                for(const auto selected:selectedFrontierRows)
                    if(supportsOverlap(frontier,frontierRow,
                                       frontier.rows[selected])) {
                        overlap=true;
                        break;
                    }
                if(overlap) {
                    ++result->metrics.deferredFrontierRows;
                    continue;
                }
                PlannedFeaturePlacement plan;
                plan.frontierOrdinal=rowOrdinal;
                plan.witness=frontierWitness(frontierRow);
                const auto edge=frontierRow.edgeKey.canonicalGlobalNodes;
                const auto curveIt=curveByEdge.find(edge);
                const auto surfaceIt=surfaceSupport->edges.find(edge);
                if(curveIt==curveByEdge.end()&&
                   surfaceIt==surfaceSupport->edges.end()) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::InternalError;
                    out.phase=BoundaryConformingV4Phase::EdgeSelection;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                plan.curve=curveIt!=curveByEdge.end();
                if(plan.curve) {
                    plan.kind=ProtectedFeatureEntityKind::Curve;
                    plan.segment=*curveIt->second;
                    plan.entityTag=plan.segment.edgeTag;
                    const auto stampIt=surfaceSupport->owningCurveByTag.find(
                        plan.entityTag);
                    if(stampIt==surfaceSupport->owningCurveByTag.end()||
                       plan.segment.dyadicLevel>=63U||
                       plan.segment.dyadicPath>
                         ((std::numeric_limits<std::uint64_t>::max)()-1U)/2U) {
                        out.status=BoundaryConformingV4Status::Unsupported;
                        out.reason=stampIt==surfaceSupport->owningCurveByTag.end()
                            ? BoundaryConformingV4Reason::
                                SharedChartReconciliationUnsupported
                            : BoundaryConformingV4Reason::
                                SubdivisionDepthUnsupported;
                        out.phase=BoundaryConformingV4Phase::EdgeSelection;
                        result->sourceOutcome=std::move(source);
                        out.retainedResult=std::move(result);
                        return out;
                    }
                    plan.entityStamp=stampIt->second;
                    plan.curveParameter=std::midpoint(
                        plan.segment.parameterLo,plan.segment.parameterHi);
                    plan.identity=curveIdentity(
                        plan.entityStamp,plan.segment.baseSegmentOrdinal,
                        plan.segment.dyadicPath*2U+1U,
                        plan.segment.dyadicLevel+1U,plan.curveParameter);
                    plan.geometryEdge=model.edge(plan.entityTag);
                    if(!(plan.curveParameter>plan.segment.parameterLo&&
                         plan.curveParameter<plan.segment.parameterHi)) {
                        out.status=BoundaryConformingV4Status::Incomplete;
                        out.reason=BoundaryConformingV4Reason::NoProgressCycle;
                        out.phase=BoundaryConformingV4Phase::EdgeSelection;
                        result->sourceOutcome=std::move(source);
                        out.retainedResult=std::move(result);
                        return out;
                    }
                } else {
                    plan.kind=ProtectedFeatureEntityKind::FaceEdge;
                    plan.surfaceEdge=surfaceIt->second;
                    plan.facePersistentId=plan.surfaceEdge.facePersistentId;
                    if(plan.surfaceEdge.incidentChartCount!=1) {
                        out.status=BoundaryConformingV4Status::Unsupported;
                        out.reason=BoundaryConformingV4Reason::
                            SharedChartReconciliationUnsupported;
                        out.phase=BoundaryConformingV4Phase::EdgeSelection;
                        result->sourceOutcome=std::move(source);
                        out.retainedResult=std::move(result);
                        return out;
                    }
                    plan.entityTag=plan.surfaceEdge.faceTag;
                    plan.entityStamp=plan.surfaceEdge.faceStamp;
                    plan.faceParameter={
                        std::midpoint(plan.surfaceEdge.firstParameter.u,
                                      plan.surfaceEdge.secondParameter.u),
                        std::midpoint(plan.surfaceEdge.firstParameter.v,
                                      plan.surfaceEdge.secondParameter.v)};
                    plan.identity=faceIdentity(
                        SurfaceNodeIdentityDomain::RefinedFaceEdge,
                        plan.entityStamp,0,0,plan.faceParameter);
                    plan.identity.parentFacetGlobalNodes=
                        plan.surfaceEdge.parentFacetGlobalNodes;
                    plan.geometryFace=model.face(plan.entityTag);
                }
                if((plan.curve&&!plan.geometryEdge)||
                   (!plan.curve&&!plan.geometryFace)) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::InternalError;
                    out.phase=BoundaryConformingV4Phase::EdgeSelection;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                const auto prior=std::find_if(
                    semanticStates.begin(),semanticStates.end(),
                    [&](const SemanticState& state) {
                        return state.first==plan.identity;
                    });
                if(prior!=semanticStates.end()) {
                    out.status=BoundaryConformingV4Status::Incomplete;
                    out.reason=BoundaryConformingV4Reason::NoProgressCycle;
                    out.phase=BoundaryConformingV4Phase::EdgeSelection;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                const auto duplicate=std::find_if(
                    planned.begin(),planned.end(),
                    [&](const PlannedFeaturePlacement& current) {
                        return sameIdentity(current,plan.identity);
                    });
                selectedFrontierRows.push_back(rowOrdinal);
                if(duplicate==planned.end())planned.push_back(plan);
            }
            } else {
                // The canonical SCC branch is prepared below from the
                // production-safe cycle sink. It intentionally bypasses the
                // single terminal frontier row, which is not a complete
                // representation of the blocker graph.
#if defined(PAACS_TOPO_TESTING)
                const t2c_detail::SurfaceCacheAllocationZoneScope
                    facetPreparationZone(
                        t2c_detail::SurfaceCacheAllocationZone::
                            FacetConflictPreparation);
#endif
                if(facetConflictCycle.membersUsed==0||
                   facetConflictCycle.membersUsed!=
                       facetConflictCycle.membersRequired||
                   facetConflictCycle.dependenciesUsed!=
                       facetConflictCycle.dependenciesRequired||
                   !source->surfaceOutcome||
                   !source->surfaceOutcome->retainedArtifact||
                   source->surfaceOutcome->retainedArtifact->setStamp!=
                       facetConflictCycle.artifactStamp) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::InternalError;
                    out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                if(result->metrics.facetSplits>
                       request.facetConflictLimits.maxFacetInteriorPoints||
                   facetConflictCycle.membersUsed>
                       request.facetConflictLimits.maxFacetInteriorPoints-
                           result->metrics.facetSplits||
                   facePoints.size()>request.limits.maxAddedNodes||
                   facetConflictCycle.membersUsed>
                       request.limits.maxAddedNodes-facePoints.size()) {
                    out.status=BoundaryConformingV4Status::Incomplete;
                    out.reason=BoundaryConformingV4Reason::
                        FacetConflictBudgetReached;
                    out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                if(frontier.rowsUsed!=0&&
                   (frontier.rows[0].workspaceStamp!=
                        facetConflictCycle.workspaceStamp||
                    frontier.rows[0].artifactStamp!=
                        facetConflictCycle.artifactStamp)) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::InternalError;
                    out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                const auto memberBegin=facetConflictCycle.members.begin();
                const auto memberEnd=memberBegin+
                    static_cast<std::ptrdiff_t>(facetConflictCycle.membersUsed);
                if(!std::is_sorted(memberBegin,memberEnd)||
                   std::adjacent_find(memberBegin,memberEnd)!=memberEnd) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::InternalError;
                    out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                for(std::size_t i=0;
                    i<facetConflictCycle.dependenciesUsed;++i) {
                    const auto& dependency=facetConflictCycle.dependencies[i];
                    if(!std::binary_search(memberBegin,memberEnd,
                                           dependency[0])||
                       !std::binary_search(memberBegin,memberEnd,
                                           dependency[1])) {
                        out.status=BoundaryConformingV4Status::Failed;
                        out.reason=BoundaryConformingV4Reason::InternalError;
                        out.phase=BoundaryConformingV4Phase::
                            ConstraintRecovery;
                        result->sourceOutcome=std::move(source);
                        out.retainedResult=std::move(result);
                        return out;
                    }
                }
                std::size_t requiredResolverRows=0;
                for(const auto& faceOwner:
                    source->surfaceOutcome->retainedArtifact->faces) {
                    bool relevant=false;
                    for(std::size_t i=0;i<facetConflictCycle.membersUsed;++i) {
                        if(!consumeFacetConflictWork(1)) {
                            out.status=BoundaryConformingV4Status::Incomplete;
                            out.reason=BoundaryConformingV4Reason::
                                FacetConflictBudgetReached;
                            out.phase=BoundaryConformingV4Phase::
                                ConstraintRecovery;
                            result->sourceOutcome=std::move(source);
                            out.retainedResult=std::move(result);
                            return out;
                        }
                        if(faceOwner&&faceOwner->faceStamp.persistentId==
                               facetConflictCycle.members[i].facePersistentId)
                            relevant=true;
                    }
                    if(relevant&&faceOwner)
                        requiredResolverRows=(std::max)(requiredResolverRows,
                                                       faceOwner->localCells.size());
                }
                std::size_t currentAccounted=0;
                std::size_t projectedAccounted=0;
                const std::size_t additionalResolverRows=
                    requiredResolverRows>facetResolverScratch.capacity()
                        ? requiredResolverRows-
                            facetResolverScratch.capacity():0U;
                if(requiredResolverRows>facetResolverScratch.max_size()||
                   !accountAttempt(privateDestination,*source,
                                   currentAccounted)||
                   !checkedSizeBytes(projectedAccounted,
                                     additionalResolverRows,
                                     sizeof(FacetConflictRankedFacetScratch))||
                   !checkedSizeAdd(projectedAccounted,currentAccounted)||
                   projectedAccounted>request.limits.maxAccountedBytes) {
                    out.status=BoundaryConformingV4Status::Incomplete;
                    out.reason=BoundaryConformingV4Reason::
                        FacetConflictBudgetReached;
                    out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                facetResolverScratch.resize(requiredResolverRows);
                if(!accountAttempt(privateDestination,*source,
                                   currentAccounted)||
                   currentAccounted>request.limits.maxAccountedBytes) {
                    out.status=BoundaryConformingV4Status::Incomplete;
                    out.reason=BoundaryConformingV4Reason::
                        FacetConflictBudgetReached;
                    out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                const std::uint64_t cycleWitness=
                    facetConflictWitness(facetConflictCycle);
                std::uint64_t preparationWork=1;
                if(!checkedMetricAdd(preparationWork,
                        facetConflictCycle.membersUsed)||
                   !checkedMetricAdd(preparationWork,
                        facetConflictCycle.dependenciesUsed)||
                   !consumeFacetConflictWork(preparationWork)) {
                    out.status=BoundaryConformingV4Status::Incomplete;
                    out.reason=BoundaryConformingV4Reason::
                        FacetConflictBudgetReached;
                    out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                for(std::size_t i=0;i<facetConflictCycle.membersUsed;++i) {
                    const auto& member=facetConflictCycle.members[i];
                    const auto facetRemaining=
                        request.facetConflictLimits.maxCoordinatorWorkUnits-
                            result->metrics.facetConflictCoordinatorWorkUnits;
                    const auto globalRemaining=
                        request.limits.maxCoordinatorWorkUnits-
                            result->metrics.coordinatorWorkUnits;
                    const auto resolved=resolveFacetConflict(
                        *source->surfaceOutcome->retainedArtifact,member,
                        facetResolverScratch,
                        (std::min)(facetRemaining,globalRemaining));
                    if(!consumeFacetConflictWork(resolved.workUnits)) {
                        out.status=BoundaryConformingV4Status::Incomplete;
                        out.reason=BoundaryConformingV4Reason::
                            FacetConflictBudgetReached;
                        out.phase=BoundaryConformingV4Phase::
                            ConstraintRecovery;
                        result->sourceOutcome=std::move(source);
                        out.retainedResult=std::move(result);
                        return out;
                    }
                    if(resolved.reason==
                           FacetConflictResolveReason::WorkBudgetReached||
                       resolved.reason==FacetConflictResolveReason::
                           WorkAccountingOverflow) {
                        out.status=BoundaryConformingV4Status::Incomplete;
                        out.reason=BoundaryConformingV4Reason::
                            FacetConflictBudgetReached;
                        out.phase=BoundaryConformingV4Phase::
                            ConstraintRecovery;
                        result->sourceOutcome=std::move(source);
                        out.retainedResult=std::move(result);
                        return out;
                    }
                    if(!resolved.resolved()) {
                        out.status=resolved.reason==
                                FacetConflictResolveReason::AllocationFailed
                            ? BoundaryConformingV4Status::Failed
                            : resolved.status==
                                  FacetConflictResolveStatus::Unsupported
                                ? BoundaryConformingV4Status::Unsupported
                                : BoundaryConformingV4Status::Failed;
                        out.reason=resolved.reason==
                                FacetConflictResolveReason::
                                    InvalidChartCertificate
                            ? BoundaryConformingV4Reason::
                                SurfaceCacheAuthorityChanged
                            : resolved.reason==
                                FacetConflictResolveReason::AllocationFailed
                            ? BoundaryConformingV4Reason::AllocationFailed
                            : resolved.status==
                                  FacetConflictResolveStatus::Unsupported
                                ? BoundaryConformingV4Reason::
                                    FacetConflictSubdivisionUnsupported
                                : BoundaryConformingV4Reason::InternalError;
                        out.phase=BoundaryConformingV4Phase::
                            ConstraintRecovery;
                        result->sourceOutcome=std::move(source);
                        out.retainedResult=std::move(result);
                        return out;
                    }
                    FacetBarycentreRequest barycentreRequest;
                    barycentreRequest.faceStamp=resolved.faceStamp;
                    barycentreRequest.chart=resolved.chart;
                    barycentreRequest.expectedProviderStamp=
                        resolved.chartCertificate.providerStamp;
                    barycentreRequest.expectedCertificateStamp=
                        resolved.chartCertificate.certificateStamp;
                    barycentreRequest.vertices=resolved.vertices;
                    const auto baryFacetRemaining=
                        request.facetConflictLimits.maxCoordinatorWorkUnits-
                            result->metrics.facetConflictCoordinatorWorkUnits;
                    const auto baryGlobalRemaining=
                        request.limits.maxCoordinatorWorkUnits-
                            result->metrics.coordinatorWorkUnits;
                    const auto barycentre=computeFacetBarycentre(
                        barycentreRequest,resolved.chartCertificate,
                        (std::min)(baryFacetRemaining,baryGlobalRemaining));
                    if(!consumeFacetConflictWork(barycentre.workUnits)) {
                        out.status=BoundaryConformingV4Status::Incomplete;
                        out.reason=BoundaryConformingV4Reason::
                            FacetConflictBudgetReached;
                        out.phase=BoundaryConformingV4Phase::
                            ConstraintRecovery;
                        result->sourceOutcome=std::move(source);
                        out.retainedResult=std::move(result);
                        return out;
                    }
                    if(barycentre.status==
                           FacetBarycentreStatus::WorkBudgetReached||
                       barycentre.status==
                           FacetBarycentreStatus::WorkAccountingOverflow) {
                        out.status=BoundaryConformingV4Status::Incomplete;
                        out.reason=BoundaryConformingV4Reason::
                            FacetConflictBudgetReached;
                        out.phase=BoundaryConformingV4Phase::
                            ConstraintRecovery;
                        result->sourceOutcome=std::move(source);
                        out.retainedResult=std::move(result);
                        return out;
                    }
                    if(!barycentre.ready()) {
                        const bool authorityFailure=barycentre.status==
                            FacetBarycentreStatus::InvalidCertificateAuthority;
                        out.status=authorityFailure
                            ? BoundaryConformingV4Status::Failed
                            : BoundaryConformingV4Status::Unsupported;
                        out.reason=authorityFailure
                            ? BoundaryConformingV4Reason::
                                SurfaceCacheAuthorityChanged
                            : BoundaryConformingV4Reason::
                                FacetConflictSubdivisionUnsupported;
                        out.phase=BoundaryConformingV4Phase::
                            ConstraintRecovery;
                        result->sourceOutcome=std::move(source);
                        out.retainedResult=std::move(result);
                        return out;
                    }
                    const auto* geometryFace=model.face(resolved.faceTag);
                    if(!consumeFacetConflictWork(1)) {
                        out.status=BoundaryConformingV4Status::Incomplete;
                        out.reason=BoundaryConformingV4Reason::
                            FacetConflictBudgetReached;
                        out.phase=BoundaryConformingV4Phase::
                            ConstraintRecovery;
                        result->sourceOutcome=std::move(source);
                        out.retainedResult=std::move(result);
                        return out;
                    }
                    if(!geometryFace||
                       !geometryFace->containsParam(barycentre.parameter)) {
                        out.status=BoundaryConformingV4Status::Unsupported;
                        out.reason=BoundaryConformingV4Reason::
                            FacetConflictSubdivisionUnsupported;
                        out.phase=BoundaryConformingV4Phase::
                            ConstraintRecovery;
                        result->sourceOutcome=std::move(source);
                        out.retainedResult=std::move(result);
                        return out;
                    }
                    bool repeated=false;
                    for(const auto& state:semanticStates) {
                        if(!consumeFacetConflictWork(1)) {
                            out.status=BoundaryConformingV4Status::Incomplete;
                            out.reason=BoundaryConformingV4Reason::
                                FacetConflictBudgetReached;
                            out.phase=BoundaryConformingV4Phase::
                                ConstraintRecovery;
                            result->sourceOutcome=std::move(source);
                            out.retainedResult=std::move(result);
                            return out;
                        }
                        if(state.first==barycentre.identity) {
                            repeated=true;
                            break;
                        }
                    }
                    for(const auto& row:planned) {
                        if(repeated)break;
                        if(!consumeFacetConflictWork(1)) {
                            out.status=BoundaryConformingV4Status::Incomplete;
                            out.reason=BoundaryConformingV4Reason::
                                FacetConflictBudgetReached;
                            out.phase=BoundaryConformingV4Phase::
                                ConstraintRecovery;
                            result->sourceOutcome=std::move(source);
                            out.retainedResult=std::move(result);
                            return out;
                        }
                        if(row.identity==barycentre.identity)repeated=true;
                    }
                    if(repeated) {
                        out.status=BoundaryConformingV4Status::Incomplete;
                        out.reason=BoundaryConformingV4Reason::NoProgressCycle;
                        out.phase=BoundaryConformingV4Phase::
                            ConstraintRecovery;
                        result->sourceOutcome=std::move(source);
                        out.retainedResult=std::move(result);
                        return out;
                    }
                    PlannedFeaturePlacement plan;
                    plan.frontierOrdinal=i;
                    plan.kind=ProtectedFeatureEntityKind::FacetInterior;
                    plan.curve=false;
                    plan.entityTag=resolved.faceTag;
                    plan.facePersistentId=member.facePersistentId;
                    plan.entityStamp=resolved.faceStamp;
                    plan.identity=barycentre.identity;
                    plan.faceParameter=barycentre.parameter;
                    plan.witness=cycleWitness;
                    plan.geometryFace=geometryFace;
                    planned.push_back(plan);
                }
            }
            if(planned.empty()) {
                out.status=BoundaryConformingV4Status::Incomplete;
                out.reason=BoundaryConformingV4Reason::NoProgressCycle;
                out.phase=BoundaryConformingV4Phase::EdgeSelection;
                result->sourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }
            // A repeated structured identity takes precedence over the hard
            // round stop.  At this point every planned identity is new; now
            // secure both the facet-round and outer source-attempt capacity
            // before the first provider call.
            if((facetConflictCaptured&&
                result->metrics.facetConflictRounds>=
                    request.facetConflictLimits.maxFacetConflictRounds)||
               round+1U==request.limits.maxRounds) {
                out.status=BoundaryConformingV4Status::Incomplete;
                out.reason=facetConflictCaptured
                    ? BoundaryConformingV4Reason::FacetConflictBudgetReached
                    : BoundaryConformingV4Reason::SubdivisionBudgetReached;
                out.phase=facetConflictCaptured
                    ? BoundaryConformingV4Phase::ConstraintRecovery
                    : BoundaryConformingV4Phase::SourceAttempt;
                result->sourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }
            if(facetConflictCaptured) {
                t2c_detail::SurfaceFacetCacheProposalCapacities capacities;
                capacities.committedFacePointRows=facePoints.size();
                capacities.committedFacePointCapacity=facePoints.capacity();
                capacities.preparedFacePointRows=pendingFacePoints.size();
                capacities.preparedFacePointCapacity=pendingFacePoints.capacity();
                capacities.committedLedgerRows=
                    result->subdivisionLedger.size();
                capacities.committedLedgerCapacity=
                    result->subdivisionLedger.capacity();
                capacities.preparedLedgerRows=pendingLedger.size();
                capacities.preparedLedgerCapacity=pendingLedger.capacity();
                capacities.committedSemanticStateRows=semanticStates.size();
                capacities.committedSemanticStateCapacity=
                    semanticStates.capacity();
                capacities.preparedSemanticStateRows=
                    pendingSemanticStates.size();
                capacities.preparedSemanticStateCapacity=
                    pendingSemanticStates.capacity();
                t2c_detail::SurfaceFacetCacheProposal cacheProposal;
                const auto plannedFacets=facetConflictCycle.members.first(
                    facetConflictCycle.membersUsed);
                if(!surfaceCache.prepareFacetProposal(
                       cacheCanonicalFaces,plannedFacets,capacities,
                       cacheProposal)) {
                    const auto& diagnostic=surfaceCache.diagnostic();
                    result->metrics.surfaceCachePreparationWorkUnits=
                        surfaceCache.metrics().preparationWorkUnits;
                    result->metrics.surfaceCacheAccountedBytes=
                        surfaceCache.metrics().accountedByteHighWater;
                    (void)consumeWork(0);
                    retainCacheReport();
                    out.status=diagnostic.status==
                            t2c_detail::SurfaceFaceCacheStatus::Incomplete
                        ? BoundaryConformingV4Status::Incomplete
                        : BoundaryConformingV4Status::Failed;
                    out.reason=diagnostic.reason==
                            t2c_detail::SurfaceFaceCacheReason::AllocationFailed
                        ? BoundaryConformingV4Reason::AllocationFailed
                        : diagnostic.status==
                              t2c_detail::SurfaceFaceCacheStatus::Incomplete
                            ? BoundaryConformingV4Reason::
                                SurfaceCacheBudgetReached
                            : BoundaryConformingV4Reason::SurfaceCacheInvalid;
                    out.phase=BoundaryConformingV4Phase::
                        SurfaceCachePreparation;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                result->metrics.surfaceCachePreparationWorkUnits=
                    surfaceCache.metrics().preparationWorkUnits;
                result->metrics.surfaceCacheAccountedBytes=
                    surfaceCache.metrics().accountedByteHighWater;
                if(!consumeWork(0)) {
                    surfaceCache.abandonFacetProposal();
                    out.status=BoundaryConformingV4Status::Incomplete;
                    out.reason=BoundaryConformingV4Reason::
                        FacetConflictBudgetReached;
                    out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
            }
            if(facetConflictCaptured) {
                std::size_t accounted=0;
                if(!accountAttempt(privateDestination,*source,accounted)||
                   accounted>request.limits.maxAccountedBytes) {
                    surfaceCache.abandonFacetProposal();
                    out.status=BoundaryConformingV4Status::Incomplete;
                    out.reason=BoundaryConformingV4Reason::
                        FacetConflictBudgetReached;
                    out.phase=BoundaryConformingV4Phase::ConstraintRecovery;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                ++result->metrics.facetConflictRounds;
            }
            const std::size_t splitCount=splits.size()+facePoints.size();
            const std::uint64_t applicableRefinementLimit=
                facetConflictCaptured?maximumRefinementRows:
                    request.limits.maxEdgeSplits;
            if(splitCount>applicableRefinementLimit||
               splitCount>request.limits.maxAddedNodes||
               planned.size()>applicableRefinementLimit-splitCount||
               planned.size()>request.limits.maxAddedNodes-splitCount) {
                out.status=BoundaryConformingV4Status::Incomplete;
                out.reason=facetConflictCaptured
                    ? BoundaryConformingV4Reason::FacetConflictBudgetReached
                    : BoundaryConformingV4Reason::SubdivisionBudgetReached;
                out.phase=facetConflictCaptured
                    ? BoundaryConformingV4Phase::ConstraintRecovery
                    : BoundaryConformingV4Phase::EdgeSelection;
                result->sourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }
            if(planned.size()>
               (std::numeric_limits<std::uint64_t>::max)()/3U) {
                out.status=BoundaryConformingV4Status::Incomplete;
                out.reason=BoundaryConformingV4Reason::SubdivisionBudgetReached;
                out.phase=BoundaryConformingV4Phase::EdgePlacement;
                result->sourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }
            const std::uint64_t providerUnits=
                static_cast<std::uint64_t>(planned.size())*3U;
            const bool facetProviderReservationOverflow=
                facetConflictCaptured&&providerUnits==
                    (std::numeric_limits<std::uint64_t>::max)();
            const std::uint64_t facetProviderReservation=
                facetConflictCaptured&&!facetProviderReservationOverflow
                    ? providerUnits+1U:providerUnits;
            if((facetConflictCaptured&&
                (facetProviderReservationOverflow||
                 result->metrics.facetConflictProviderCalls>
                    request.facetConflictLimits.maxProviderCalls||
                 facetProviderReservation>
                    request.facetConflictLimits.maxProviderCalls-
                        result->metrics.facetConflictProviderCalls))||
               providerUnits>request.limits.maxCoordinatorWorkUnits-
                    result->metrics.coordinatorWorkUnits) {
                out.status=BoundaryConformingV4Status::Incomplete;
                out.reason=facetConflictCaptured
                    ? BoundaryConformingV4Reason::FacetConflictBudgetReached
                    : BoundaryConformingV4Reason::SubdivisionBudgetReached;
                out.phase=facetConflictCaptured
                    ? BoundaryConformingV4Phase::ConstraintRecovery
                    : BoundaryConformingV4Phase::EdgePlacement;
                result->sourceOutcome=std::move(source);
                out.retainedResult=std::move(result);
                return out;
            }

            activePhase=BoundaryConformingV4Phase::EdgePlacement;
#if defined(PAACS_TOPO_TESTING)
            const t2c_detail::SurfaceCacheAllocationZoneScope
                facetProviderBatchZone(facetConflictCaptured
                    ? t2c_detail::SurfaceCacheAllocationZone::
                        FacetConflictProviderBatch
                    : t2c_detail::SurfaceCacheAllocationZone::None);
#endif
            for(const auto& plan:planned) {
                ProtectedFeaturePlacementRequest placementRequest;
                placementRequest.kind=plan.kind;
                placementRequest.entityStamp=plan.entityStamp;
                placementRequest.requestedIdentity=plan.identity;
                placementRequest.entityTag=plan.entityTag;
                placementRequest.curveParameter=plan.curveParameter;
                placementRequest.faceParameter=plan.faceParameter;
                // Secure the fixed pre/place/post call budget before entering
                // the provider protocol.  A placement exception never causes
                // the post-call checkpoint to execute.
                if(!(facetConflictCaptured
                        ? consumeFacetConflictWork(3):consumeWork(3))) {
                    out.status=BoundaryConformingV4Status::Incomplete;
                    out.reason=facetConflictCaptured
                        ? BoundaryConformingV4Reason::
                            FacetConflictBudgetReached
                        : BoundaryConformingV4Reason::
                            SubdivisionBudgetReached;
                    out.phase=facetConflictCaptured
                        ? BoundaryConformingV4Phase::ConstraintRecovery
                        : BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                ProtectedFeatureProviderStamp preStamp{};
                try {
                    ++result->metrics.providerStampCalls;
                    if(facetConflictCaptured)
                        ++result->metrics.facetConflictProviderCalls;
                    preStamp=placementProvider.stamp();
                } catch(const std::bad_alloc&) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::AllocationFailed;
                    out.phase=BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                } catch(...) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::ProviderFailed;
                    out.phase=BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                if(preStamp!=request.expectedPlacementProviderStamp) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::
                        ProviderRevisionChanged;
                    out.phase=BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                ProtectedFeaturePlacement placement;
                try {
                    ++result->metrics.providerPlacementCalls;
                    if(facetConflictCaptured)
                        ++result->metrics.facetConflictProviderCalls;
                    placement=plan.curve
                        ? placementProvider.placeOnEdge(*plan.geometryEdge,
                                                       placementRequest)
                        : placementProvider.placeOnFace(*plan.geometryFace,
                                                       placementRequest);
                } catch(const std::bad_alloc&) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::AllocationFailed;
                    out.phase=BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                  catch(...) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::ProviderFailed;
                    out.phase=BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                // The returned certificate is checked before the post-call
                // stamp exactly as the placement-provider contract requires.
                if(placement.providerStamp!=
                        request.expectedPlacementProviderStamp) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::
                        ProviderRevisionChanged;
                    out.phase=BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                if(placement.entityStamp!=plan.entityStamp||
                   placement.requestedIdentity!=plan.identity) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::
                        PlacementContractViolation;
                    out.phase=BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                ProtectedFeatureProviderStamp postStamp{};
                try {
                    ++result->metrics.providerStampCalls;
                    if(facetConflictCaptured)
                        ++result->metrics.facetConflictProviderCalls;
                    postStamp=placementProvider.stamp();
                } catch(const std::bad_alloc&) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::AllocationFailed;
                    out.phase=BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                } catch(...) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::ProviderFailed;
                    out.phase=BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                if(postStamp!=request.expectedPlacementProviderStamp) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::
                        ProviderRevisionChanged;
                    out.phase=BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                if(placement.status!=ProtectedFeaturePlacementStatus::Placed) {
                    out.status=placement.status==
                            ProtectedFeaturePlacementStatus::Rejected
                        ? BoundaryConformingV4Status::Rejected
                        : placement.status==
                              ProtectedFeaturePlacementStatus::Unsupported
                            ? BoundaryConformingV4Status::Unsupported
                            : BoundaryConformingV4Status::Failed;
                    out.reason=placement.status==
                            ProtectedFeaturePlacementStatus::Rejected
                        ? BoundaryConformingV4Reason::ProviderRejected
                        : placement.status==
                              ProtectedFeaturePlacementStatus::Unsupported
                            ? BoundaryConformingV4Reason::ProviderUnsupported
                            : BoundaryConformingV4Reason::ProviderFailed;
                    out.phase=BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                Vec3 authoritative{};
                try {
                    authoritative=plan.curve
                        ? plan.geometryEdge->pointAt(plan.curveParameter)
                        : plan.geometryFace->pointAt(plan.faceParameter);
                } catch(const std::bad_alloc&) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::AllocationFailed;
                    out.phase=BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                } catch(...) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::ProviderFailed;
                    out.phase=BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                if(!finite(placement.coordinate)||!finite(authoritative)||
                   !exactPoint(placement.coordinate,authoritative)) {
                    out.status=BoundaryConformingV4Status::Failed;
                    out.reason=BoundaryConformingV4Reason::
                        PlacementContractViolation;
                    out.phase=BoundaryConformingV4Phase::EdgePlacement;
                    result->sourceOutcome=std::move(source);
                    out.retainedResult=std::move(result);
                    return out;
                }
                prepared.push_back({plan,placement});
            }

            // All provider work is complete. Prepare a bounded overlay only;
            // the live split/ledger/cache state is not changed until the next
            // source attempt has rebuilt and audited this exact proposal.
#if defined(PAACS_TOPO_TESTING)
            const t2c_detail::SurfaceCacheAllocationZoneScope allocationZone(
                t2c_detail::SurfaceCacheAllocationZone::
                    ProposedStateIntegration);
#endif
            pendingSplits.clear(); pendingFacePoints.clear();
            pendingSemanticStates.clear(); pendingLedger.clear();
            for(const auto& row:prepared) {
                pendingSemanticStates.push_back(
                    {row.plan.identity,row.plan.witness});
                const std::uint64_t ordinal=
                    static_cast<std::uint64_t>(
                        result->subdivisionLedger.size()+
                        pendingLedger.size()+1);
                if(row.plan.curve)
                    pendingSplits.push_back({row.plan.entityTag,
                        row.plan.segment.baseSegmentOrdinal,
                        row.plan.curveParameter,row.placement.coordinate,
                        row.plan.identity});
                else
                    pendingFacePoints.push_back({
                        row.plan.facePersistentId,
                        row.plan.faceParameter,row.placement.coordinate,
                        row.plan.identity});
                ProtectedFeatureChildIdentity ledger;
                ledger.parentEntityPersistentId=
                    row.plan.entityStamp.persistentId;
                ledger.childIdentity=row.plan.identity;
                ledger.splitOrdinal=ordinal;
                ledger.edgeTag=row.plan.curve?row.plan.entityTag:-1;
                ledger.faceTag=row.plan.curve?-1:row.plan.entityTag;
                ledger.parentParameterLo=row.plan.curve
                    ? row.plan.segment.parameterLo:0.0;
                ledger.parentParameterHi=row.plan.curve
                    ? row.plan.segment.parameterHi:0.0;
                ledger.midpointParameter=row.plan.curveParameter;
                ledger.faceParameter=row.plan.faceParameter;
                ledger.coordinate=row.placement.coordinate;
                ledger.classification={row.plan.curve?1:2,
                                       row.plan.entityTag};
                ledger.providerStamp=row.placement.providerStamp;
                pendingLedger.push_back(ledger);
            }
            source.reset();
        }
        out.status=BoundaryConformingV4Status::Incomplete;
        out.reason=BoundaryConformingV4Reason::SubdivisionBudgetReached;
        out.phase=BoundaryConformingV4Phase::SourceAttempt;
        out.retainedResult=std::move(result);
        return out;
    } catch(const std::bad_alloc&) {
        out.status=BoundaryConformingV4Status::Failed;
        out.reason=BoundaryConformingV4Reason::AllocationFailed;
        out.phase=activePhase;
        if(result) {
            if(cacheEnvelopePrepared) {
                auto& report=result->surfaceCacheReport;
                report.processed.clear(); report.unattempted.clear();
                report.reportComplete=true;
                report.hasFailedFace=false;
                report.failedFacePersistentId=0;
                SurfaceFailure failure;
                failure.reason=SurfaceReason::AllocationFailed;
                failure.allocation=SurfaceAllocationPhase::SetMapping;
                for(const auto face:cacheCanonicalFaces)
                    report.processed.push_back({face,
                        SurfaceCacheFaceDisposition::Stopped,
                        SurfaceFaceStatus::Failed,failure});
            }
            out.retainedResult=std::move(result);
        }
        return out;
    } catch(...) {
        out.status=BoundaryConformingV4Status::Failed;
        out.reason=BoundaryConformingV4Reason::InternalError;
        out.phase=activePhase;
        if(result) {
            if(cacheEnvelopePrepared) {
                auto& report=result->surfaceCacheReport;
                report.processed.clear(); report.unattempted.clear();
                report.reportComplete=true;
                report.hasFailedFace=false;
                report.failedFacePersistentId=0;
                SurfaceFailure failure;
                failure.reason=SurfaceReason::InternalError;
                for(const auto face:cacheCanonicalFaces)
                    report.processed.push_back({face,
                        SurfaceCacheFaceDisposition::Stopped,
                        SurfaceFaceStatus::Failed,failure});
            }
            out.retainedResult=std::move(result);
        }
        return out;
    }
}

}  // namespace paacs::algo
