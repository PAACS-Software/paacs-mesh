// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctContactCylinder.hpp"
#include "paacs/geometry/OcctBoundary.hpp"
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <GProp_GProps.hxx>
#include <Message_ProgressIndicator.hxx>
#include <Standard_OutOfMemory.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <NCollection_List.hxx>
#include <gp_Circ.hxx>
#include <gp_Trsf.hxx>
#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <set>
#include <streambuf>
#include <ostream>

namespace paacs::geo::detail {
namespace {
using Reason = OcctContactCylinderReason;
struct Stop { Reason reason; };
void need(bool value, Reason reason) { if (!value) throw Stop{reason}; }
bool finite(Vec3 p) { return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z); }
Vec3 vec(gp_Pnt p) { return {p.X(),p.Y(),p.Z()}; }
struct Context {
    const OcctContactCylinderOptions& options;
    OcctContactCylinderOutcome& result;
    void tick(std::uint64_t count=1) {
        need(!options.cancellation||!options.cancellation->load(std::memory_order_relaxed),Reason::Cancelled);
        need(count<=options.maximumWork-result.workUsed,Reason::BudgetExceeded);
        result.workUsed+=count;
    }
};
class Progress final : public Message_ProgressIndicator {
public: explicit Progress(const std::atomic<bool>* flag):flag_(flag){}
protected:
    bool UserBreak() override { return flag_&&flag_->load(std::memory_order_relaxed); }
    void Show(const Message_ProgressScope&,bool) override {}
private: const std::atomic<bool>* flag_;
};
class Buffer final : public std::streambuf {
public: explicit Buffer(std::uint64_t maximum):maximum_(maximum){} std::string value;
protected:
    std::streamsize xsputn(const char* data,std::streamsize n) override {
        need(n>=0&&static_cast<std::uint64_t>(n)<=maximum_-value.size(),Reason::BudgetExceeded);
        value.append(data,static_cast<std::size_t>(n));return n;
    }
    int_type overflow(int_type c) override {
        if(traits_type::eq_int_type(c,traits_type::eof()))return traits_type::not_eof(c);
        const char b=traits_type::to_char_type(c);xsputn(&b,1);return c;
    }
private: std::uint64_t maximum_;
};
std::string snapshot(const TopoDS_Shape& shape,std::uint64_t maximum) {
    Buffer buffer(maximum);std::ostream stream(&buffer);stream.exceptions(std::ios::badbit|std::ios::failbit);
    BRepTools::Write(shape,stream);return std::move(buffer.value);
}
TopoDS_Shape common(const TopoDS_Shape& a,const TopoDS_Shape& b,Context& context) {
    context.tick();BRepAlgoAPI_Common operation;NCollection_List<TopoDS_Shape> args,tools;
    args.Append(a);tools.Append(b);operation.SetArguments(args);operation.SetTools(tools);
    operation.SetNonDestructive(true);operation.SetRunParallel(false);
    occ::handle<Progress> progress=new Progress(context.options.cancellation);
    operation.Build(progress->Start());context.tick();
    need(operation.IsDone()&&!operation.HasErrors(),Reason::NativeProofFailed);return operation.Shape();
}
double length(const TopoDS_Shape& shape) { GProp_GProps p;BRepGProp::LinearProperties(shape,p);return p.Mass(); }
double volume(const TopoDS_Shape& shape) { GProp_GProps p;BRepGProp::VolumeProperties(shape,p);return p.Mass(); }
bool has(const TopoDS_Shape& shape,TopAbs_ShapeEnum kind) { TopExp_Explorer e(shape,kind);return e.More(); }
struct Cap { int face=-1,edge=-1;gp_Circ circle; };
std::vector<Index> rim(const OcctModel& model,const mesh::Mesh& accepted,int face,int edge,Context& context) {
    const auto* curve=model.edge(edge);
    const auto belongs=[&](Index n){need(n<accepted.nodeCount(),Reason::InvalidRim);const auto c=accepted.nodeClass(n);
        return (c.dim==1&&c.tag==edge)||(c.dim==0&&(c.tag==curve->startVertexTag()||c.tag==curve->endVertexTag()));};
    std::map<std::pair<Index,Index>,unsigned> counts;
    for(Index i=0;i<accepted.cellCount();++i){context.tick();const auto& c=accepted.cell(i);
        if(c.onEntity!=mesh::Classification{2,face})continue;need(c.type==mesh::CellType::Tri3,Reason::InvalidRim);
        const auto* ids=accepted.cellNodes(i);for(int j=0;j<3;++j){auto a=ids[j],b=ids[(j+1)%3];if(belongs(a)&&belongs(b)){
            need(a!=b,Reason::InvalidRim);if(a>b)std::swap(a,b);++counts[{a,b}];}}
    }
    std::map<Index,std::vector<Index>> graph;
    for(const auto& [edgeIds,count]:counts){context.tick();need(count<=2,Reason::InvalidRim);if(count==1){
        graph[edgeIds.first].push_back(edgeIds.second);graph[edgeIds.second].push_back(edgeIds.first);}}
    if(graph.empty())return {};
    need(graph.size()>=3,Reason::InvalidRim);for(const auto& [id,adj]:graph){context.tick();need(adj.size()==2,Reason::InvalidRim);}
    std::vector<Index> loop;loop.reserve(graph.size());const auto start=graph.begin()->first;
    Index current=start,previous=(std::numeric_limits<Index>::max)();
    do {context.tick();loop.push_back(current);const auto& neighbors=graph.at(current);const auto next=neighbors[0]==previous?neighbors[1]:neighbors[0];
        previous=current;current=next;need(loop.size()<=graph.size(),Reason::InvalidRim);}while(current!=start);
    need(loop.size()==graph.size(),Reason::InvalidRim);return loop;
}
void admit(std::uint64_t nodes,std::uint64_t cells,const OcctContactCylinderOptions& options) {
    constexpr std::uint64_t nb=sizeof(Vec3)+sizeof(mesh::Classification),cb=sizeof(mesh::Cell)+4*sizeof(Index);
    need(nodes<=options.maximumNodes&&cells<=options.maximumCells&&nodes<=options.maximumOutputBytes/nb,
         Reason::BudgetExceeded);
    need(cells<=(options.maximumOutputBytes-nodes*nb)/cb,Reason::BudgetExceeded);
    // Boundary, staged fan, audit/certificate and graph container payload.
    need(nodes<=options.maximumWorkingBytes/512&&cells<=(options.maximumWorkingBytes-nodes*512)/512,
         Reason::BudgetExceeded);
}
} // namespace

const char* occtContactCylinderReasonName(OcctContactCylinderReason r) noexcept {
    constexpr const char* names[]={"none","invalid_request","unsupported_cylinder","contact_circle_not_found",
        "ambiguous_contact","native_proof_failed","invalid_rim","geometry_error","volume_mismatch","star_failed",
        "invalid_mesh","budget_exceeded","cancelled","allocation_failed","source_changed","internal_error"};
    const auto i=static_cast<unsigned>(r);return i<std::size(names)?names[i]:"internal_error";
}
OcctContactCylinderOutcome buildOcctContactCylinder(const OcctModel& model,int targetTag,
    const mesh::Mesh& accepted,const OcctContactCylinderOptions& options,mesh::Mesh& destination) noexcept {
    OcctContactCylinderOutcome result;
    try {
        Context context{options,result};context.tick();
        need(&accepted!=&destination&&model.solid(targetTag)&&accepted.nodeCount()>0&&accepted.cellCount()>0&&
             std::isfinite(options.maximumChordError)&&options.maximumChordError>0&&
             std::isfinite(options.maximumRelativeVolumeError)&&options.maximumRelativeVolumeError>=0,
             Reason::InvalidRequest);
        // Includes bounded source snapshots and conservative discovery-container
        // logical payload. The already accepted caller mesh is excluded.
        need(options.maximumSnapshotBytes<=options.maximumWorkingBytes/4,Reason::BudgetExceeded);
        auto remainingWorking=options.maximumWorkingBytes-4*options.maximumSnapshotBytes;
        const auto acceptedCells=static_cast<std::uint64_t>(accepted.cellCount());
        const auto acceptedNodes=static_cast<std::uint64_t>(accepted.nodeCount());
        need(acceptedCells<=remainingWorking/128,Reason::BudgetExceeded);
        remainingWorking-=acceptedCells*128;
        need(acceptedNodes<=remainingWorking/64,Reason::BudgetExceeded);
        remainingWorking-=acceptedNodes*64;
        const auto before=snapshot(OcctModelAccess::shape(model),options.maximumSnapshotBytes);
        BRepBuilderAPI_Copy copy(OcctModelAccess::shape(model),true,false);context.tick();
        const auto copiedEdge=[&](int tag){return TopoDS::Edge(copy.ModifiedShape(OcctModelAccess::edge(model,tag)));};
        std::vector<Cap> caps;int side=-1;std::set<int> targetEdges;
        need(model.faceUses(targetTag).size()==3,Reason::UnsupportedCylinder);
        for(const auto& use:model.faceUses(targetTag)) {
            context.tick();const auto* descriptor=model.surfaceDescriptor(use.faceTag);need(descriptor,Reason::UnsupportedCylinder);
            for(const auto& wire:model.boundaryUses(use.faceTag))for(const auto& edge:wire.edges){context.tick();targetEdges.insert(edge.edgeTag);}
            if(descriptor->kind==SurfaceKind::Cylinder){need(side<0,Reason::UnsupportedCylinder);side=use.faceTag;continue;}
            need(descriptor->kind==SurfaceKind::Plane,Reason::UnsupportedCylinder);
            const auto wires=model.boundaryUses(use.faceTag);need(wires.size()==1&&wires[0].outer&&wires[0].edges.size()==1,
                                                                 Reason::UnsupportedCylinder);
            const auto edgeTag=wires[0].edges[0].edgeTag;const auto* edge=model.edge(edgeTag);
            need(edge&&edge->kind()==CurveKind::Circle&&edge->startVertexTag()==edge->endVertexTag(),Reason::UnsupportedCylinder);
            BRepAdaptor_Curve curve(copiedEdge(edgeTag));need(curve.GetType()==GeomAbs_Circle&&
                std::abs(curve.LastParameter()-curve.FirstParameter()-2*std::numbers::pi)<=1e-8,Reason::UnsupportedCylinder);
            caps.push_back({use.faceTag,edgeTag,curve.Circle()});
        }
        need(side>0&&caps.size()==2&&targetEdges.size()==3,Reason::UnsupportedCylinder);
        std::set<int> acceptedFaces;
        for(Index i=0;i<accepted.cellCount();++i){context.tick();const auto& c=accepted.cell(i);if(c.onEntity.dim!=2)continue;
            need(model.face(c.onEntity.tag)&&c.type==mesh::CellType::Tri3,Reason::InvalidRim);
            const auto owners=model.faceOwners(c.onEntity.tag);need(!owners.empty(),Reason::InvalidRim);
            for(auto owner:owners)need(owner.solidTag!=targetTag,Reason::InvalidRequest);
            acceptedFaces.insert(c.onEntity.tag);
        }
        // Search only actual original edges incident to accepted boundary faces.
        std::uint64_t edgeUses=0;
        for(int face:acceptedFaces)for(const auto& wire:model.boundaryUses(face))for(const auto& edge:wire.edges){
            (void)edge;context.tick();need(edgeUses<remainingWorking/128,Reason::BudgetExceeded);++edgeUses;}
        remainingWorking-=edgeUses*128;
        auto stagedOptions=options;stagedOptions.maximumWorkingBytes=remainingWorking;
        std::map<int,std::set<int>> edgeFaces;
        for(int face:acceptedFaces)for(const auto& wire:model.boundaryUses(face))for(auto edge:wire.edges){context.tick();
            if(model.edge(edge.edgeTag)->kind()==CurveKind::Circle)edgeFaces[edge.edgeTag].insert(face);}
        std::vector<Index> selected;int selectedCap=-1;
        for(unsigned capIndex=0;capIndex<caps.size();++capIndex)for(const auto& [edgeTag,faces]:edgeFaces) {
            context.tick();const auto* edge=model.edge(edgeTag);if(edge->startVertexTag()!=edge->endVertexTag())continue;
            const auto source=copiedEdge(edgeTag),target=copiedEdge(caps[capIndex].edge);
            const double aLength=length(source),bLength=length(target);context.tick(2);
            const double allowance=std::max({8*model.tolerance(1,edgeTag),8*model.tolerance(1,caps[capIndex].edge),1e-9});
            if(!std::isfinite(aLength)||!std::isfinite(bLength)||std::abs(aLength-bLength)>allowance)continue;
            const auto overlap=common(source,target,context);if(!has(overlap,TopAbs_EDGE))continue;
            const double overlapLength=length(overlap);context.tick();
            if(!std::isfinite(overlapLength)||std::abs(overlapLength-aLength)>allowance||std::abs(overlapLength-bLength)>allowance)continue;
            for(int face:faces){auto loop=rim(model,accepted,face,edgeTag,context);if(loop.empty())continue;
                if(!selected.empty()) {auto a=loop,b=selected;std::sort(a.begin(),a.end());std::sort(b.begin(),b.end());
                    need(selectedCap==static_cast<int>(capIndex)&&result.sourceEdgeTag==edgeTag&&a==b,Reason::AmbiguousContact);continue;}
                selected=std::move(loop);selectedCap=static_cast<int>(capIndex);result.sourceEdgeTag=edgeTag;result.sourceFaceTag=face;
            }
        }
        need(!selected.empty(),Reason::ContactCircleNotFound);result.nativeCircleProved=true;
        const auto& bottom=caps[selectedCap];const auto& top=caps[1-selectedCap];result.targetEdgeTag=bottom.edge;result.contactCapFaceTag=bottom.face;
        const Vec3 origin=vec(bottom.circle.Location()),topCenter=vec(top.circle.Location()),delta=topCenter-origin;
        const double height=delta.norm();need(std::isfinite(height)&&height>0,Reason::UnsupportedCylinder);const auto axis=delta/height;
        need(std::abs(bottom.circle.Radius()-top.circle.Radius())<=std::max(8*model.tolerance(1,bottom.edge),1e-9),Reason::UnsupportedCylinder);
        const auto canonical=BRepPrimAPI_MakeCylinder(gp_Ax2(bottom.circle.Location(),gp_Dir(axis.x,axis.y,axis.z)),bottom.circle.Radius(),height).Shape();context.tick();
        result.cadVolume=model.volume(targetTag);const double canonicalVolume=volume(canonical);context.tick();
        const double overlapVolume=volume(common(copy.ModifiedShape(OcctModelAccess::solid(model,targetTag)),canonical,context));context.tick();
        need(std::isfinite(result.cadVolume)&&result.cadVolume>0&&std::isfinite(canonicalVolume)&&canonicalVolume>0&&std::isfinite(overlapVolume)&&
             std::abs(overlapVolume-result.cadVolume)<=1e-8*result.cadVolume&&
             std::abs(overlapVolume-canonicalVolume)<=1e-8*canonicalVolume,Reason::NativeProofFailed);result.nativeCylinderProved=true;
        const auto sourceVertex=TopExp::FirstVertex(copiedEdge(result.sourceEdgeTag),true);
        const auto bottomVertex=TopExp::FirstVertex(copiedEdge(bottom.edge),true),topVertex=TopExp::FirstVertex(copiedEdge(top.edge),true);
        need(has(common(sourceVertex,bottomVertex,context),TopAbs_VERTEX),Reason::NativeProofFailed);
        gp_Trsf transform;transform.SetTranslation(gp_Vec(delta.x,delta.y,delta.z));context.tick();
        const auto translated=BRepBuilderAPI_Transform(bottomVertex,transform,true).Shape();
        need(has(common(translated,topVertex,context),TopAbs_VERTEX),Reason::NativeProofFailed);
        result.ringNodes=selected.size();need(result.ringNodes<=UINT64_MAX/8,Reason::BudgetExceeded);
        admit(2*result.ringNodes+3,8*result.ringNodes,stagedOptions);
        double winding=0;for(std::size_t i=0;i<selected.size();++i){context.tick();
            winding+=(accepted.node(selected[i])-origin).cross(accepted.node(selected[(i+1)%selected.size()])-origin).dot(axis);}
        need(std::isfinite(winding)&&winding!=0,Reason::InvalidRim);if(winding<0)std::reverse(selected.begin(),selected.end());
        mesh::Mesh boundary;std::vector<Index> lower,upper;lower.reserve(selected.size());upper.reserve(selected.size());unsigned cornerCount=0;
        const auto* sourceEdge=model.edge(result.sourceEdgeTag);const auto* lowerEdge=model.edge(bottom.edge);const auto* upperEdge=model.edge(top.edge);
        for(Index sourceId:selected){context.tick();const auto cls=accepted.nodeClass(sourceId);mesh::Classification targetClass{1,bottom.edge};
            if(cls.dim==0){need(cls.tag==sourceEdge->startVertexTag(),Reason::InvalidRim);targetClass={0,lowerEdge->startVertexTag()};++cornerCount;}
            lower.push_back(boundary.addNode(accepted.node(sourceId),targetClass));}
        need(cornerCount==1,Reason::InvalidRim);
        for(Index id:lower){context.tick();const auto cls=boundary.nodeClass(id);upper.push_back(boundary.addNode(boundary.node(id)+delta,
            cls.dim==0?mesh::Classification{0,upperEdge->startVertexTag()}:mesh::Classification{1,top.edge}));}
        const auto lc=boundary.addNode(origin,{2,bottom.face}),uc=boundary.addNode(topCenter,{2,top.face});
        for(std::size_t i=0;i<selected.size();++i){context.tick();const auto j=(i+1)%selected.size();const auto a=lower[i],b=lower[j],c=upper[i],d=upper[j];
            boundary.addCell(mesh::CellType::Tri3,{lc,b,a},{2,bottom.face});boundary.addCell(mesh::CellType::Tri3,{uc,c,d},{2,top.face});
            boundary.addCell(mesh::CellType::Tri3,{a,b,d},{2,side});boundary.addCell(mesh::CellType::Tri3,{a,d,c},{2,side});}
        for(Index i=0;i<boundary.nodeCount();++i){context.tick();const auto p=boundary.node(i);const auto cls=boundary.nodeClass(i);Vec3 q;
            need(finite(p),Reason::GeometryError);if(cls.dim==0)q=model.vertex(cls.tag)->position();else if(cls.dim==1)q=model.edge(cls.tag)->closestPoint(p);
            else need(closestPointOnOcctFace(model,cls.tag,p,q),Reason::GeometryError);
            const double error=(p-q).norm();result.maximumBoundaryNodeError=std::max(result.maximumBoundaryNodeError,error);
            need(std::isfinite(error)&&error<=std::max(8*model.tolerance(cls.dim,cls.tag),1e-9)&&model.solid(targetTag)->isInside(p),Reason::GeometryError);}
        for(Index i=0;i<boundary.cellCount();++i){const auto* n=boundary.cellNodes(i);const auto a=boundary.node(n[0]),b=boundary.node(n[1]),c=boundary.node(n[2]);
            for(auto sample:{(a+b+c)/3.,(a+b)/2.,(b+c)/2.,(c+a)/2.}){context.tick();Vec3 q;
                need(closestPointOnOcctFace(model,boundary.cell(i).onEntity.tag,sample,q),Reason::GeometryError);const double error=(sample-q).norm();
                result.maximumBoundarySampleError=std::max(result.maximumBoundarySampleError,error);need(std::isfinite(error)&&error<=options.maximumChordError,Reason::GeometryError);}}
        auto starOptions=options.star;starOptions.maximumWork=std::min(starOptions.maximumWork,options.maximumWork-result.workUsed);
        starOptions.maximumNodes=std::min(starOptions.maximumNodes,options.maximumNodes);starOptions.maximumCells=std::min(starOptions.maximumCells,options.maximumCells);
        starOptions.maximumOutputBytes=std::min(starOptions.maximumOutputBytes,options.maximumOutputBytes);
        starOptions.maximumWorkingBytes=std::min(starOptions.maximumWorkingBytes,remainingWorking);
        starOptions.cancellation=options.cancellation;mesh::Mesh staged;const auto star=buildOcctStarVolume(boundary,targetTag,starOptions,staged);
        result.starReason=star.reason;context.tick(star.workUsed);if(!star.published){
            if(star.reason==OcctStarReason::Cancelled)throw Stop{Reason::Cancelled};
            if(star.reason==OcctStarReason::WorkBudget||star.reason==OcctStarReason::OutputBudget)throw Stop{Reason::BudgetExceeded};
            if(star.reason==OcctStarReason::AllocationFailed)throw Stop{Reason::AllocationFailed};throw Stop{Reason::StarFailed};}
        for(Index i=0;i<staged.nodeCount();++i){context.tick();need(model.solid(targetTag)->isInside(staged.node(i)),Reason::GeometryError);}
        for(Index i=0;i<staged.cellCount();++i){context.tick();if(staged.cell(i).onEntity.dim!=3)continue;const auto* n=staged.cellNodes(i);
            need(model.solid(targetTag)->isInside((staged.node(n[0])+staged.node(n[1])+staged.node(n[2])+staged.node(n[3]))/4.),Reason::GeometryError);}
        auto limits=options.star.auditLimits;limits.cancellation=options.cancellation;
        result.audit=mesh::auditLinearVolumeBoundary(staged,limits);context.tick(staged.cellCount());
        if(!result.audit.valid){if(result.audit.cancelled)throw Stop{Reason::Cancelled};if(result.audit.budgetExceeded)throw Stop{Reason::BudgetExceeded};
            if(result.audit.allocationFailed)throw Stop{Reason::AllocationFailed};throw Stop{Reason::InvalidMesh};}
        result.quality=algo::analyzeMixedMeshQuality(staged);context.tick(staged.cellCount());
        need(!result.quality.allocationFailed,Reason::AllocationFailed);need(result.quality.valid,Reason::InvalidMesh);
        result.meshVolume=result.quality.totalVolume;result.relativeVolumeError=std::abs(result.meshVolume-result.cadVolume)/result.cadVolume;
        need(std::isfinite(result.relativeVolumeError)&&result.relativeVolumeError<=options.maximumRelativeVolumeError,Reason::VolumeMismatch);
        context.tick();need(snapshot(OcctModelAccess::shape(model),options.maximumSnapshotBytes)==before,Reason::SourceChanged);
        result.sourceImmutabilityVerified=true;context.tick();destination.swap(staged);result.published=true;result.destinationUnchanged=false;result.reason=Reason::None;
    } catch(const Stop& s){result.reason=s.reason;}
      catch(const Standard_OutOfMemory&){result.reason=Reason::AllocationFailed;}
      catch(const std::bad_alloc&){result.reason=Reason::AllocationFailed;}
      catch(const Standard_Failure&){result.reason=Reason::NativeProofFailed;}
      catch(...){result.reason=Reason::InternalError;}
    return result;
}
} // namespace paacs::geo::detail
