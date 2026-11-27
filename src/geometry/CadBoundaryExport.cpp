// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0
#include "paacs/geometry/CadBoundaryExport.hpp"

#include <Standard_Failure.hxx>
#include <Standard_OutOfMemory.hxx>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <utility>

namespace paacs::geo {
namespace {
using Clock = std::chrono::steady_clock;
using Pair = std::array<Index, 2>;
struct Stop { CadBoundaryEdgeReason reason; };
void need(bool condition, CadBoundaryEdgeReason reason) {
    if (!condition) throw Stop{reason};
}
bool finite(Vec3 p) { return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z); }
double distance(Vec3 a, Vec3 b) { return std::hypot(a.x-b.x,a.y-b.y,a.z-b.z); }
Pair pair(Index a, Index b) { return {std::min(a,b),std::max(a,b)}; }
double segmentDistance(Vec3 p, Vec3 a, Vec3 b) {
    const auto d=b-a,q=p-a;
    const auto scale=std::max({std::abs(d.x),std::abs(d.y),std::abs(d.z),
                               std::abs(q.x),std::abs(q.y),std::abs(q.z)});
    if (!std::isfinite(scale)) return (std::numeric_limits<double>::infinity)();
    if (scale==0) return 0;
    const auto u=d/scale,v=q/scale;
    const auto norm=u.dot(u),dot=u.dot(v);
    const auto t=norm==0||dot<=0?0.:dot>=norm?1.:dot/norm;
    const auto delta=v-u*t;
    return std::hypot(delta.x,delta.y,delta.z)*scale;
}
struct Context {
    const CadBoundaryEdgeOptions& options;
    CadBoundaryEdgeOutcome& result;
    Clock::time_point started;
    std::uint64_t bytes=0;
    void tick(std::uint64_t amount=1) {
        need(amount<=options.maximumWork-result.work,CadBoundaryEdgeReason::BudgetExceeded);
        result.work+=amount;
        if (options.cancellation&&options.cancellation->load(std::memory_order_relaxed))
            throw Stop{CadBoundaryEdgeReason::Cancelled};
        need(std::chrono::duration<double>(Clock::now()-started).count()<=options.maximumSeconds,
             CadBoundaryEdgeReason::BudgetExceeded);
    }
    void storage(std::uint64_t count,std::uint64_t width) {
        need(width&&count<=(options.maximumWorkingBytes-bytes)/width,CadBoundaryEdgeReason::BudgetExceeded);
        bytes+=count*width;
    }
};
} // namespace

CadBoundaryEdgeOutcome exportCadBoundaryEdges(const OcctModel& model,int solidTag,
    const mesh::Mesh& boundary,const CadBoundaryEdgeOptions& options,mesh::Mesh& destination) noexcept {
    CadBoundaryEdgeOutcome result;result.solidTag=solidTag;
    const auto started=Clock::now();Context context{options,result,started};
    try {
        need(&boundary!=&destination&&model.solid(solidTag)&&options.maximumNodes&&options.maximumTriangles&&
            options.maximumSegments&&options.maximumWork&&options.maximumWorkingBytes&&options.maximumOutputBytes&&
            std::isfinite(options.maximumSeconds)&&options.maximumSeconds>0&&options.maximumSeconds<=86400&&
            std::isfinite(options.maximumNodeError)&&options.maximumNodeError>=0&&
            std::isfinite(options.maximumChordError)&&options.maximumChordError>0,CadBoundaryEdgeReason::InvalidRequest);
        context.tick();
        need(boundary.nodeCount()>0&&boundary.cellCount()>0,CadBoundaryEdgeReason::InvalidBoundary);
        const auto nodes=std::uint64_t(boundary.nodeCount()),triangles=std::uint64_t(boundary.cellCount());
        need(nodes<=options.maximumNodes&&triangles<=options.maximumTriangles,CadBoundaryEdgeReason::BudgetExceeded);
        result.nodes=nodes;result.triangles=triangles;
        // Includes node copies/indexes, graph tree nodes/incident-face sets,
        // per-edge ordered metadata and temporary sorting/storage overhead.
        context.storage(nodes,160);context.storage(triangles,1024);
        need(options.maximumOutputBytes>=40&&nodes<=(options.maximumOutputBytes-40)/32,
             CadBoundaryEdgeReason::BudgetExceeded);
        const auto outputBase=40+32*nodes;

        std::set<int> bodyFaces;
        std::map<int,std::set<int>> faceVertices;
        std::map<int,CadBoundaryEdgeRecord> records;
        for (const auto& face:model.faceUses(solidTag)) {
            context.tick();bodyFaces.insert(face.faceTag);
            std::uint64_t wireIndex=0;
            for (const auto& wire:model.boundaryUses(face.faceTag)) {
                for (const auto& use:wire.edges) {
                    context.tick();context.storage(1,256);
                    auto [position,inserted]=records.try_emplace(use.edgeTag);
                    auto& record=position->second;
                    if (inserted) {
                        context.storage(1,512);
                        record.edgeTag=use.edgeTag;
                        const auto* edge=model.edge(use.edgeTag);
                        need(edge!=nullptr,CadBoundaryEdgeReason::InvalidBoundary);
                        record.startVertexTag=edge->startVertexTag();record.endVertexTag=edge->endVertexTag();
                        record.periodic=edge->periodic();
                        record.closed=record.startVertexTag==record.endVertexTag;
                        need(std::isfinite(edge->length())&&edge->length()>=0,CadBoundaryEdgeReason::GeometryError);
                        record.degenerated=edge->length()==0;
                        const auto faces=model.edgeFaces(use.edgeTag);
                        record.faceTags.assign(faces.begin(),faces.end());
                        std::set<int> bodies;
                        for (const int tag:faces) for (const auto& owner:model.faceOwners(tag)) {
                            context.tick();context.storage(1,32);bodies.insert(owner.solidTag);
                        }
                        record.bodyTags.assign(bodies.begin(),bodies.end());
                    }
                    record.uses.push_back({face.faceTag,wireIndex,use.useId,use.orientation,use.sense,use.seam,use.degenerated});
                    faceVertices[face.faceTag].insert(use.startVertexTag);
                    faceVertices[face.faceTag].insert(use.endVertexTag);
                }
                ++wireIndex;
            }
        }
        need(!bodyFaces.empty()&&!records.empty(),CadBoundaryEdgeReason::InvalidBoundary);
        std::map<std::pair<int,int>,std::vector<Index>> classified;
        for (Index i=0;i<boundary.nodeCount();++i) {
            context.tick();const auto cls=boundary.nodeClass(i);const auto p=boundary.node(i);
            need(finite(p)&&cls.tag>0&&cls.dim>=0&&cls.dim<=2,CadBoundaryEdgeReason::InvalidBoundary);
            if (cls.dim==0) {
                const auto* vertex=model.vertex(cls.tag);
                need(vertex!=nullptr,CadBoundaryEdgeReason::InvalidBoundary);
                const auto error=distance(p,vertex->position());
                need(std::isfinite(error)&&error<=options.maximumNodeError+model.tolerance(0,cls.tag),
                     CadBoundaryEdgeReason::GeometryError);
                result.maximumNodeError=std::max(result.maximumNodeError,error);
            } else if (cls.dim==1) need(records.contains(cls.tag),CadBoundaryEdgeReason::InvalidBoundary);
            else need(bodyFaces.contains(cls.tag),CadBoundaryEdgeReason::InvalidBoundary);
            classified[{cls.dim,cls.tag}].push_back(i);
        }
        std::map<Pair,std::set<int>> graph;
        std::set<int> representedFaces;
        for (Index i=0;i<boundary.cellCount();++i) {
            context.tick();const auto cell=boundary.cell(i);
            need(cell.type==mesh::CellType::Tri3&&cell.onEntity.dim==2&&bodyFaces.contains(cell.onEntity.tag),
                 CadBoundaryEdgeReason::InvalidBoundary);
            representedFaces.insert(cell.onEntity.tag);
            const auto* ids=boundary.cellNodes(i);
            for (unsigned j=0;j<3;++j) {
                need(ids[j]>=0&&ids[j]<boundary.nodeCount()&&ids[j]!=ids[(j+1)%3],CadBoundaryEdgeReason::InvalidBoundary);
                const auto cls=boundary.nodeClass(ids[j]);
                need(cls.dim==2?cls.tag==cell.onEntity.tag:cls.dim==0?
                    faceVertices.contains(cell.onEntity.tag)&&faceVertices.at(cell.onEntity.tag).contains(cls.tag):
                    std::binary_search(model.edgeFaces(cls.tag).begin(),model.edgeFaces(cls.tag).end(),cell.onEntity.tag),
                    CadBoundaryEdgeReason::IncidenceMismatch);
                graph[pair(ids[j],ids[(j+1)%3])].insert(cell.onEntity.tag);
            }
        }
        need(representedFaces==bodyFaces,CadBoundaryEdgeReason::InvalidBoundary);
        std::map<Pair,int> claimed;
        std::vector<bool> used(nodes);
        for (auto& [tag,record]:records) {
            result.failureEdgeTag=tag;context.tick();const auto* edge=model.edge(tag);
            const auto anchor=[&](int vertexTag) {
                const auto at=classified.find({0,vertexTag});
                need(at!=classified.end()&&!at->second.empty(),CadBoundaryEdgeReason::MissingAnchor);
                need(at->second.size()==1,CadBoundaryEdgeReason::AmbiguousAnchor);
                return at->second.front();
            };
            if (record.degenerated) {
                ++result.degeneratedEdges;
                const auto at=classified.find({0,record.startVertexTag});
                if (at!=classified.end()) {
                    need(at->second.size()==1,CadBoundaryEdgeReason::AmbiguousAnchor);
                    record.naturalNodeOrder=at->second;record.naturalParameters={0.};
                }
                need(!classified.contains({1,tag}),CadBoundaryEdgeReason::IncidenceMismatch);
                continue;
            }
            double first,last;edge->parameterBounds(first,last);
            need(std::isfinite(first)&&std::isfinite(last)&&last>first,CadBoundaryEdgeReason::ParameterMappingFailed);
            const auto start=anchor(record.startVertexTag),end=anchor(record.endVertexTag);
            const double tolerance=options.maximumNodeError+model.tolerance(1,tag)+
                std::max(model.tolerance(0,record.startVertexTag),model.tolerance(0,record.endVertexTag));
            const auto startError=distance(boundary.node(start),edge->pointAt(first));
            const auto endError=distance(boundary.node(end),edge->pointAt(last));context.tick(2);
            need(std::isfinite(startError)&&std::isfinite(endError)&&std::max(startError,endError)<=tolerance,
                 CadBoundaryEdgeReason::GeometryError);
            result.maximumNodeError=std::max({result.maximumNodeError,startError,endError});
            const auto edgeNodes=classified.find({1,tag});
            const auto orderedCount=2+(edgeNodes==classified.end()?0:edgeNodes->second.size());
            context.storage(orderedCount,80);
            std::vector<std::pair<double,Index>> ordered;ordered.reserve(orderedCount);
            ordered.push_back({first,start});
            if (edgeNodes!=classified.end()) for (const auto id:edgeNodes->second) {
                context.tick(4);const auto parameter=edge->closestParameter(boundary.node(id));
                need(std::isfinite(parameter)&&parameter>first&&parameter<last,CadBoundaryEdgeReason::ParameterMappingFailed);
                const auto error=distance(boundary.node(id),edge->pointAt(parameter));
                need(std::isfinite(error)&&error<=tolerance,CadBoundaryEdgeReason::GeometryError);
                result.maximumNodeError=std::max(result.maximumNodeError,error);ordered.push_back({parameter,id});
            }
            ordered.push_back({last,end});
            context.tick(ordered.size()*std::uint64_t(std::max(1.,std::ceil(std::log2(double(ordered.size()))))));
            std::sort(ordered.begin(),ordered.end());
            need(!record.closed||ordered.size()>=4,CadBoundaryEdgeReason::MissingAnchor);
            std::set<int> expectedFaces;
            for (const int face:record.faceTags) if (bodyFaces.contains(face)) expectedFaces.insert(face);
            need(!expectedFaces.empty(),CadBoundaryEdgeReason::IncidenceMismatch);
            for (std::size_t i=0;i<ordered.size();++i) {
                record.naturalParameters.push_back(ordered[i].first);
                record.naturalNodeOrder.push_back(ordered[i].second);
                if (i==0) continue;
                context.tick();const auto [a,ia]=ordered[i-1];const auto [b,ib]=ordered[i];
                need(b>a&&ia!=ib,CadBoundaryEdgeReason::ParameterMappingFailed);
                const auto key=pair(ia,ib);const auto found=graph.find(key);
                need(found!=graph.end(),CadBoundaryEdgeReason::MissingSurfaceSegment);
                need(found->second==expectedFaces,CadBoundaryEdgeReason::IncidenceMismatch);
                need(claimed.emplace(key,tag).second,CadBoundaryEdgeReason::AmbiguousAnchor);
                for (const auto fraction:std::array{.25,.5,.75}) {
                    context.tick();const auto sample=edge->pointAt(std::lerp(a,b,fraction));
                    const auto error=segmentDistance(sample,boundary.node(ia),boundary.node(ib));
                    need(finite(sample)&&std::isfinite(error)&&error<=options.maximumChordError+model.tolerance(1,tag),
                         CadBoundaryEdgeReason::GeometryError);
                    result.maximumSampledChordError=std::max(result.maximumSampledChordError,error);
                }
                need(result.segments<options.maximumSegments&&result.segments<(options.maximumOutputBytes-outputBase)/28,
                     CadBoundaryEdgeReason::BudgetExceeded);
                record.naturalLineOrder.push_back(Index(result.segments++));used[ia]=used[ib]=true;
            }
        }
        need(result.segments>0,CadBoundaryEdgeReason::InvalidBoundary);
        need(result.segments<=std::uint64_t((std::numeric_limits<Index>::max)())/2,
             CadBoundaryEdgeReason::BudgetExceeded);
        mesh::Mesh output;output.reserveAppend(boundary.nodeCount(),Index(result.segments),Index(2*result.segments));
        for (Index i=0;i<boundary.nodeCount();++i) {
            context.tick();output.addNode(boundary.node(i),boundary.nodeClass(i));
            result.unusedNodes+=!used[std::size_t(i)];
        }
        result.edges.reserve(records.size());
        for (auto& [tag,record]:records) {
            for (std::size_t i=1;i<record.naturalNodeOrder.size()&&!record.degenerated;++i) {
                context.tick();output.addCell(mesh::CellType::Line2,
                    {record.naturalNodeOrder[i-1],record.naturalNodeOrder[i]},{1,tag});
            }
            result.edges.push_back(std::move(record));
        }
        context.tick();destination.swap(output);result.complete=true;result.destinationUnchanged=false;
        result.failureEdgeTag=-1;result.reason=CadBoundaryEdgeReason::None;
    } catch (const Stop& stop) { result.reason=stop.reason; }
      catch (const Standard_OutOfMemory&) { result.reason=CadBoundaryEdgeReason::AllocationFailed; }
      catch (const std::bad_alloc&) { result.reason=CadBoundaryEdgeReason::AllocationFailed; }
      catch (const Standard_Failure&) { result.reason=CadBoundaryEdgeReason::KernelFailure; }
      catch (...) { result.reason=CadBoundaryEdgeReason::InternalError; }
    result.seconds=std::chrono::duration<double>(Clock::now()-started).count();return result;
}

const char* cadBoundaryEdgeReasonName(CadBoundaryEdgeReason reason) noexcept {
    switch (reason) {
    case CadBoundaryEdgeReason::None:return "none";
    case CadBoundaryEdgeReason::InvalidRequest:return "invalid_request";
    case CadBoundaryEdgeReason::InvalidBoundary:return "invalid_boundary";
    case CadBoundaryEdgeReason::MissingAnchor:return "missing_cad_anchor";
    case CadBoundaryEdgeReason::AmbiguousAnchor:return "ambiguous_cad_anchor";
    case CadBoundaryEdgeReason::ParameterMappingFailed:return "cad_parameter_mapping_failed";
    case CadBoundaryEdgeReason::MissingSurfaceSegment:return "missing_surface_segment";
    case CadBoundaryEdgeReason::IncidenceMismatch:return "cad_incidence_mismatch";
    case CadBoundaryEdgeReason::GeometryError:return "geometry_error";
    case CadBoundaryEdgeReason::BudgetExceeded:return "budget_exceeded";
    case CadBoundaryEdgeReason::Cancelled:return "cancelled";
    case CadBoundaryEdgeReason::AllocationFailed:return "allocation_failed";
    case CadBoundaryEdgeReason::KernelFailure:return "kernel_failure";
    case CadBoundaryEdgeReason::InternalError:return "internal_error";
    }
    return "unknown";
}
} // namespace paacs::geo
