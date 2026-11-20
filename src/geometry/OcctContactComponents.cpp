// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctContactComponents.hpp"
#include "OcctContactAssembly.hpp"
#include "OcctModelInternal.hpp"
#include "OcctNativeValidation.hpp"
#include "OcctStageTrace.hpp"
#include "OcctStarMesher.hpp"
#include "paacs/geometry/OcctBoundary.hpp"
#include "paacs/bnd/TetFacetGeometry.hpp"
#include <BRepAlgoAPI_Common.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepCheck_Analyzer.hxx>
#include "OcctSolidClassifier.hpp"
#include <BRepClass_FaceClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <NCollection_IndexedMap.hxx>
#include <Standard_OutOfMemory.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS_Compound.hxx>
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <map>
#include <numeric>
#include <set>
#include <string_view>
#include <vector>

namespace paacs::geo::detail {
namespace {
using ShapeMap=NCollection_IndexedMap<TopoDS_Shape,TopTools_ShapeMapHasher>;
using Tri=std::array<Index,3>;
using Clock=std::chrono::steady_clock;
struct Stop {NativeCadMeshReason reason;const char* diagnostic;};
void need(bool b,const char* s,NativeCadMeshReason reason=NativeCadMeshReason::GeometryError){if(!b)throw Stop{reason,s};}
struct Budget {
    const NativeCadMeshOptions& options;
    std::uint64_t limit,work=0,maximumWorking=0;
    std::uint64_t remaining() const{return limit-work;}
    void tick(std::uint64_t n=1){
        need(!options.geometry.cancellation||!options.geometry.cancellation->load(std::memory_order_relaxed),"cancelled",NativeCadMeshReason::Cancelled);
        need(n<=remaining(),"component_work_budget",NativeCadMeshReason::BudgetExceeded);work+=n;
    }
    void collect(std::uint64_t n){need(n<=remaining(),"component_work_budget",NativeCadMeshReason::BudgetExceeded);work+=n;tick(0);}
    void working(std::uint64_t bytes){maximumWorking=std::max(maximumWorking,bytes);
        need(bytes<=options.maximumPrimitiveWorkingBytes,"component_working_budget",NativeCadMeshReason::BudgetExceeded);tick(0);}
};
std::uint64_t add(std::uint64_t a,std::uint64_t b){need(b<=UINT64_MAX-a,"component_payload_budget",NativeCadMeshReason::BudgetExceeded);return a+b;}
std::uint64_t mul(std::uint64_t a,std::uint64_t b){need(!b||a<=UINT64_MAX/b,"component_payload_budget",NativeCadMeshReason::BudgetExceeded);return a*b;}
std::uint64_t payload(const mesh::Mesh& m){return add(add(mul(m.nodeCount(),sizeof(Vec3)+sizeof(mesh::Classification)),mul(m.cellCount(),sizeof(mesh::Cell))),mul(m.connectivityCount(),sizeof(Index)));}
void output(std::uint64_t nodes,std::uint64_t cells,std::uint64_t connectivity,Budget& b){
    const auto bytes=add(add(mul(nodes,sizeof(Vec3)+sizeof(mesh::Classification)),mul(cells,sizeof(mesh::Cell))),mul(connectivity,sizeof(Index)));
    need(nodes<=b.options.geometry.maximumNodes&&cells<=b.options.geometry.maximumCells&&bytes<=b.options.geometry.maximumOutputBytes,
         "component_output_budget",NativeCadMeshReason::BudgetExceeded);b.tick(0);
}
void audit(const mesh::Mesh& value,Budget& b,mesh::LinearVolumeAudit& result){
    mesh::LinearVolumeAuditLimits limits;limits.maximumTemporaryBytes=b.options.geometry.maximumAuditTemporaryBytes;
    limits.maximumFaceUses=b.options.geometry.maximumAuditFaceUses;limits.cancellation=b.options.geometry.cancellation;
    b.tick(value.cellCount());result=mesh::auditLinearVolumeBoundary(value,limits);
    need(!result.cancelled,"cancelled",NativeCadMeshReason::Cancelled);
    need(!result.budgetExceeded,"component_audit_budget",NativeCadMeshReason::BudgetExceeded);
    need(!result.allocationFailed,"allocation_failed",NativeCadMeshReason::AllocationFailed);
    need(result.valid,"component_boundary_audit_failed");
}
algo::MixedMeshQualityReport quality(const mesh::Mesh& value,Budget& b,std::uint64_t retained){
    // The analyzer reserves one double for every cell before sorting.
    b.working(add(retained,mul(value.cellCount(),sizeof(double))));
    b.tick(value.cellCount());auto q=algo::analyzeMixedMeshQuality(value);
    need(!q.allocationFailed,"allocation_failed",NativeCadMeshReason::AllocationFailed);
    need(q.complete&&q.valid,"component_nonpositive_or_incomplete_quality");return q;
}
bool same(Vec3 a,Vec3 b){return std::bit_cast<std::uint64_t>(a.x)==std::bit_cast<std::uint64_t>(b.x)&&
    std::bit_cast<std::uint64_t>(a.y)==std::bit_cast<std::uint64_t>(b.y)&&std::bit_cast<std::uint64_t>(a.z)==std::bit_cast<std::uint64_t>(b.z);}
Tri sorted(Tri t){std::sort(t.begin(),t.end());return t;}
TopoDS_Shape skin(const TopoDS_Shape& body){TopoDS_Compound shape;BRep_Builder builder;builder.MakeCompound(shape);
    ShapeMap faces;TopExp::MapShapes(body,TopAbs_FACE,faces);for(int i=1;i<=faces.Extent();++i)builder.Add(shape,faces(i));return shape;}
TopoDS_Shape common(const TopoDS_Shape& a,const TopoDS_Shape& b,Budget& budget){
    budget.tick();BRepAlgoAPI_Common op;NCollection_List<TopoDS_Shape> first,second;first.Append(a);second.Append(b);
    op.SetArguments(first);op.SetTools(second);op.SetNonDestructive(true);op.SetRunParallel(false);op.Build();budget.tick(0);
    need(op.IsDone()&&!op.HasErrors(),"component_contact_query_failed",NativeCadMeshReason::UnresolvedBodyInterface);return op.Shape();
}
double volume(const TopoDS_Shape& shape,Budget& b){
    b.tick(2);GProp_GProps first,second;
    const auto e1=BRepGProp::VolumePropertiesGK(shape,first,1e-8,false,true,false,false,false);
    const auto e2=BRepGProp::VolumePropertiesGK(shape,second,1e-10,false,true,false,false,false);b.tick(0);
    need(std::isfinite(e1)&&e1>=0&&std::isfinite(e2)&&e2>=0&&std::isfinite(first.Mass())&&first.Mass()>0&&std::isfinite(second.Mass())&&second.Mass()>0,
         "component_native_integration_failed");
    need(std::abs(first.Mass()-second.Mass())<=second.Mass()*1e-6,"component_native_integration_unstable");return second.Mass();
}
mesh::Classification remap(mesh::Classification c,const std::array<std::vector<int>,4>& tags){
    need(c.dim>=0&&c.dim<=3&&c.tag>0&&static_cast<std::size_t>(c.tag)<tags[c.dim].size(),"component_native_class_invalid");
    const int original=tags[c.dim][c.tag];need(original>0,"component_native_class_history_missing");return {c.dim,original};
}
struct Triangle {std::array<Vec3,3> p;Vec3 lo,hi;};
struct TreeNode {Vec3 lo,hi;std::size_t begin,end,left=0,right=0;};
double coordinate(Vec3 p,unsigned a){return a==0?p.x:a==1?p.y:p.z;}
bool overlaps(Vec3 a,Vec3 b,Vec3 c,Vec3 d){return a.x<=d.x&&c.x<=b.x&&a.y<=d.y&&c.y<=b.y&&a.z<=d.z&&c.z<=b.z;}
struct Tree {
    std::vector<Triangle> triangles;std::vector<TreeNode> nodes;
    std::size_t build(std::size_t begin,std::size_t end,Budget& b){
        b.tick(end-begin);Vec3 lo{1e300,1e300,1e300},hi{-1e300,-1e300,-1e300};
        for(auto i=begin;i<end;++i){const auto& t=triangles[i];lo={std::min(lo.x,t.lo.x),std::min(lo.y,t.lo.y),std::min(lo.z,t.lo.z)};hi={std::max(hi.x,t.hi.x),std::max(hi.y,t.hi.y),std::max(hi.z,t.hi.z)};}
        auto id=nodes.size();nodes.push_back({lo,hi,begin,end});if(end-begin>8){const auto span=hi-lo;unsigned axis=span.y>span.x?1:0;if(span.z>coordinate(span,axis))axis=2;
            const auto mid=begin+(end-begin)/2;std::nth_element(triangles.begin()+begin,triangles.begin()+mid,triangles.begin()+end,[&](const auto& x,const auto& y){b.tick();return coordinate(x.lo+x.hi,axis)<coordinate(y.lo+y.hi,axis);});
            const auto left=build(begin,mid,b),right=build(mid,end,b);nodes[id].left=left;nodes[id].right=right;}return id;
    }
};
Tree tree(const mesh::Mesh& value,Budget& b){Tree t;for(Index i=0;i<value.cellCount();++i){b.tick();if(value.cell(i).type!=mesh::CellType::Tri3)continue;
    const auto* ids=value.cellNodes(i);Triangle p{{value.node(ids[0]),value.node(ids[1]),value.node(ids[2])}};
    p.lo={std::min({p.p[0].x,p.p[1].x,p.p[2].x}),std::min({p.p[0].y,p.p[1].y,p.p[2].y}),std::min({p.p[0].z,p.p[1].z,p.p[2].z})};
    p.hi={std::max({p.p[0].x,p.p[1].x,p.p[2].x}),std::max({p.p[0].y,p.p[1].y,p.p[2].y}),std::max({p.p[0].z,p.p[1].z,p.p[2].z})};t.triangles.push_back(p);}
    need(!t.triangles.empty(),"component_body_boundary_missing");t.build(0,t.triangles.size(),b);return t;}
void contacts(const Tree& a,std::size_t ai,const Tree& b,std::size_t bi,bool nativeContact,Budget& budget,std::uint64_t& matched){
    budget.tick();const auto& x=a.nodes[ai];const auto& y=b.nodes[bi];if(!overlaps(x.lo,x.hi,y.lo,y.hi))return;
    if(x.left){contacts(a,x.left,b,bi,nativeContact,budget,matched);contacts(a,x.right,b,bi,nativeContact,budget,matched);return;}
    if(y.left){contacts(a,ai,b,y.left,nativeContact,budget,matched);contacts(a,ai,b,y.right,nativeContact,budget,matched);return;}
    for(auto i=x.begin;i<x.end;++i)for(auto j=y.begin;j<y.end;++j){budget.tick();const auto& p=a.triangles[i];const auto& q=b.triangles[j];if(!overlaps(p.lo,p.hi,q.lo,q.hi))continue;
        std::array<int,3> permutation{-1,-1,-1};for(int u=0;u<3;++u)for(int v=0;v<3;++v)if(same(p.p[u],q.p[v]))permutation[u]=v;
        if(permutation[0]>=0&&permutation[1]>=0&&permutation[2]>=0){const int inversions=(permutation[0]>permutation[1])+(permutation[0]>permutation[2])+(permutation[1]>permutation[2]);
            need(nativeContact&&inversions%2==1,"component_unproved_or_same_sense_contact",NativeCadMeshReason::UnresolvedBodyInterface);++matched;continue;}
        need(!bnd::tet_facet_detail::trianglesOpenIntersect(p.p,q.p)&&!bnd::tet_facet_detail::coplanarTrianglesHaveAreaOverlap(p.p,q.p),
             "component_unmatched_interbody_interface",NativeCadMeshReason::UnresolvedBodyInterface);
    }
}
void append(mesh::Mesh& joined,const mesh::Mesh& value,Budget& b){
    output(add(joined.nodeCount(),value.nodeCount()),add(joined.cellCount(),value.cellCount()),add(joined.connectivityCount(),value.connectivityCount()),b);
    const Index offset=joined.nodeCount();for(Index i=0;i<value.nodeCount();++i){b.tick();joined.addNode(value.node(i),value.nodeClass(i));}
    for(Index i=0;i<value.cellCount();++i){b.tick();const auto& c=value.cell(i);std::array<Index,4> ids{};const int count=mesh::nodesPerCell(c.type);
        need(count<=4,"component_unsupported_cell");for(int j=0;j<count;++j)ids[j]=value.cellNodes(i)[j]+offset;
        joined.addCell(c.type,std::span<const Index>(ids.data(),count),c.onEntity);}
}
void nativePoints(const OcctModel& original,int body,std::span<const Vec3> points,std::uint64_t retained,Budget& b,OcctNativeValidationOutcome& result,double maximumOutsideDistance=0){
    need(b.remaining()>0,"component_work_budget",NativeCadMeshReason::BudgetExceeded);
    OcctNativeValidationOptions options;options.maximumThreads=b.options.maximumNativeValidationThreads;options.maximumPoints=points.size();options.maximumWork=b.remaining();
    options.maximumSnapshotBytes=std::min<std::uint64_t>(64ULL*1024*1024,b.options.maximumPrimitiveWorkingBytes/8);
    options.maximumOutsideDistance=maximumOutsideDistance;
    b.working(retained);options.maximumWorkingBytes=b.options.maximumPrimitiveWorkingBytes-retained;options.cancellation=b.options.geometry.cancellation;
    result=validateOcctNativeInteriorPoints(original,body,points,options);b.collect(result.workUsed);
    b.working(add(retained,result.maximumWorkingPayloadEstimate));
    if(result.reason==OcctNativeValidationReason::BudgetExceeded)throw Stop{NativeCadMeshReason::BudgetExceeded,"component_native_query_budget"};
    if(result.reason==OcctNativeValidationReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled,"cancelled"};
    if(result.reason==OcctNativeValidationReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed,"allocation_failed"};
    need(result.completed&&result.sourceModelUnmodified,"component_original_native_node_or_centroid_outside");
}
}

NativeCadMeshOutcome buildOcctContactComponents(const OcctModel& original,const NativeCadMeshOptions& options,mesh::Mesh& destination,
                                                const std::string* sourceBefore) noexcept {
    NativeCadMeshOutcome out;const auto started=Clock::now();Budget budget{options,std::min(options.maximumValidationWork,options.geometry.maximumMeshingWork)};
    const auto snapshotLimit=std::min<std::uint64_t>(64ULL*1024*1024,options.maximumPrimitiveWorkingBytes/8);
    std::string before;
    const char* diagnostic="component_not_started";bool sourceUnchanged=false,sourceAudited=false;
    std::uint64_t originalFaceCoordinateConversions=0;
    try {
        // MSVC's mesh dedup container may allocate during construction. Keep
        // staging constructors within the allocation-failure transaction.
        mesh::Mesh joined;std::vector<mesh::Mesh> bodyMeshes;
        budget.tick(0);need(options.separateMatchedNativeFaceContacts&&options.elements==algo::PrimitiveElements::Tet4&&options.useNativeKernelAdmission,
                           "component_invalid_request",NativeCadMeshReason::InvalidRequest);
        need(options.faceSizes.empty()&&options.geometry.bodyControls.empty()&&!options.geometry.spatialField&&options.minimumSize==0&&
             options.geometry.strategy==algo::GeometryStrategy::Auto&&options.geometry.maximumElementSizeRatio==0&&options.geometry.warningElementSizeRatio==0,
             "component_unsupported_control",NativeCadMeshReason::UnsupportedControl);
        need(std::isfinite(options.geometry.targetSize)&&options.geometry.targetSize>0&&options.maximumNativeValidationThreads>=1&&options.maximumNativeValidationThreads<=8&&
             std::isfinite(options.maximumRelativeVolumeError)&&options.maximumRelativeVolumeError>0&&options.maximumRelativeVolumeError<1&&
             std::isfinite(options.maximumChordError)&&options.maximumChordError>=0&&
             std::isfinite(options.geometry.minimumElementQuality)&&options.geometry.minimumElementQuality>=0&&options.geometry.minimumElementQuality<=1&&
             std::isfinite(options.geometry.warningElementQuality)&&options.geometry.warningElementQuality>=0&&options.geometry.warningElementQuality<=1,
             "component_invalid_request",NativeCadMeshReason::InvalidRequest);
        const auto count=original.solids().size();need(count>1,"component_assembly_required",NativeCadMeshReason::UnsupportedBody);
        need(snapshotLimit>0&&budget.limit>0,"component_working_budget",NativeCadMeshReason::BudgetExceeded);
        auto sourceReason=snapshotOcctContactSource(original,snapshotLimit,before);budget.tick();need(sourceReason==NativeCadMeshReason::None,"component_source_audit_failed",sourceReason);
        if(sourceBefore)need(before==*sourceBefore,"component_original_source_changed",NativeCadMeshReason::InternalError);
        budget.working(add(mul(before.size(),sourceBefore?3:2),mul(count+1,sizeof(mesh::Mesh)+sizeof(NativeCadBodyMeshReport)+1024)));
        bodyMeshes.resize(count+1);const auto& source=OcctModelAccess::shape(original);
        std::vector<unsigned> parent(count);std::iota(parent.begin(),parent.end(),0);
        const auto find=[&](unsigned i){while(parent[i]!=i)i=parent[i];return i;};std::set<std::pair<int,int>> facePairs;
        OcctStageTrace discovery("contact_component_discovery");
        for(unsigned i=0;i<count;++i)for(unsigned j=0;j<i;++j){
            budget.tick();BRepBuilderAPI_Copy a(OcctModelAccess::solid(original,i+1),true,false),b(OcctModelAccess::solid(original,j+1),true,false);
            need(a.IsDone()&&b.IsDone(),"component_copy_failed");const auto material=common(a.Shape(),b.Shape(),budget);
            need(!TopExp_Explorer(material,TopAbs_SOLID).More(),"component_original_material_overlap",NativeCadMeshReason::UnresolvedBodyInterface);
            const auto contact=common(skin(a.Shape()),skin(b.Shape()),budget);if(TopExp_Explorer(contact,TopAbs_FACE).More()){
                GProp_GProps properties;BRepGProp::SurfaceProperties(contact,properties);budget.tick();need(std::isfinite(properties.Mass())&&properties.Mass()>0,"component_invalid_contact_area");
                facePairs.emplace(j+1,i+1);parent[find(i)]=find(j);
            }else if(TopExp_Explorer(contact,TopAbs_EDGE).More()||TopExp_Explorer(contact,TopAbs_VERTEX).More()){
                need(options.separateLowerDimensionalContacts,"component_lower_dimensional_contact_policy_required",NativeCadMeshReason::UnresolvedBodyInterface);
            }
        }
        discovery.finish();
        std::map<unsigned,std::vector<int>> components;for(unsigned i=0;i<count;++i)components[find(i)].push_back(i+1);
        bool contactComponent=false,isolated=false;for(const auto& [_,bodies]:components){contactComponent|=bodies.size()>1;isolated|=bodies.size()==1;}
        need(contactComponent&&isolated,"component_arrangement_not_applicable",NativeCadMeshReason::UnsupportedBody);
        ShapeMap originalMaps[4];const TopAbs_ShapeEnum kinds[]{TopAbs_VERTEX,TopAbs_EDGE,TopAbs_FACE,TopAbs_SOLID};
        for(unsigned d=0;d<4;++d)TopExp::MapShapes(source,kinds[d],originalMaps[d]);
        std::vector<NativeCadBodyMeshReport> reports(count+1);
        std::uint64_t retainedMetadata=mul(count+1,sizeof(mesh::Mesh)+sizeof(NativeCadBodyMeshReport)+1024);
        for(unsigned d=0;d<4;++d)retainedMetadata=add(retainedMetadata,mul(originalMaps[d].Extent(),256));
        const auto retained=[&](std::uint64_t active=0){auto bytes=add(add(mul(before.size(),sourceBefore?3:2),retainedMetadata),active);
            for(const auto& value:bodyMeshes)bytes=add(bytes,payload(value));bytes=add(bytes,payload(joined));budget.working(bytes);return bytes;};retained();
        for(const auto& [_,bodies]:components){
            OcctStageTrace component("contact_component_preparation",bodies.front());
            budget.tick();BRepBuilderAPI_Copy copy(source,true,false);need(copy.IsDone(),"component_copy_failed");
            TopoDS_Compound compound;BRep_Builder builder;builder.MakeCompound(compound);for(int tag:bodies)builder.Add(compound,copy.ModifiedShape(OcctModelAccess::solid(original,tag)));
            auto local=OcctModelAccess::fromShape(compound);need(local.imported(),"component_native_import_failed");
            std::array<std::vector<int>,4> tags;ShapeMap localMaps[4];std::uint64_t localMetadata=0;
            for(unsigned d=0;d<4;++d){TopExp::MapShapes(compound,kinds[d],localMaps[d]);tags[d].resize(localMaps[d].Extent()+1,-1);
                localMetadata=add(localMetadata,mul(tags[d].size(),sizeof(int)+256));
                for(int t=1;t<=originalMaps[d].Extent();++t){budget.tick();const int id=localMaps[d].FindIndex(copy.ModifiedShape(originalMaps[d](t)));if(id)tags[d][id]=t;}
                for(std::size_t i=1;i<tags[d].size();++i)need(tags[d][i]>0,"component_native_history_incomplete");}
            auto request=options;request.reuseNativeQueries=false;request.separateMatchedNativeFaceContacts=bodies.size()>1;
            need(budget.remaining()>0,"component_work_budget",NativeCadMeshReason::BudgetExceeded);
            request.useNativeKernelAdmission=bodies.size()>1;request.maximumValidationWork=budget.remaining();request.geometry.maximumMeshingWork=request.maximumValidationWork;
            request.maximumPrimitiveWorkingBytes=options.maximumPrimitiveWorkingBytes-retained(localMetadata);mesh::Mesh value;NativeCadMeshOutcome child;
            component.next("contact_component_meshing");
            if(bodies.size()>1){OcctContactAssemblyStrategy strategy;strategy.cuts=OcctContactCutPolicy::PrincipalFeatureAxes;strategy.radialCuts=2;
                strategy.maximumAdaptiveCuts=16;strategy.maximumPieceWork=10000000;strategy.controlDeflectionOnAllSurfaces=true;
                child=buildOcctContactAssembly(*local.model,request,strategy,value);
            }else child=meshNativeCad(*local.model,request,value);
            component.finish(nativeCadMeshReasonName(child.reason));
            component.next("partition_assembly");
            budget.collect(child.validationWork);out.preparationSeconds+=child.preparationSeconds;out.generationSeconds+=child.generationSeconds;
            out.validationSeconds+=child.validationSeconds;out.auditSeconds+=child.auditSeconds;out.qualitySeconds+=child.qualitySeconds;
            out.nativeContactPieces=add(out.nativeContactPieces,child.nativeContactPieces);out.nativeContactPlanes=add(out.nativeContactPlanes,child.nativeContactPlanes);
            retained(add(localMetadata,child.nativeContactMaximumWorkingPayloadEstimate));
            need(child.published,"component_child_not_published",child.reason);retained(add(localMetadata,payload(value)));
            for(auto report:child.bodies){const int body=tags[3].at(report.solidTag);report.solidTag=body;
                if(report.startFaceTag>0)report.startFaceTag=tags[2].at(report.startFaceTag);if(report.endFaceTag>0)report.endFaceTag=tags[2].at(report.endFaceTag);
                if(report.contactSourceEdgeTag>0)report.contactSourceEdgeTag=tags[1].at(report.contactSourceEdgeTag);if(report.contactTargetEdgeTag>0)report.contactTargetEdgeTag=tags[1].at(report.contactTargetEdgeTag);
                reports[body]=std::move(report);}
            // A boundary triangle's body comes from its adjacent Tet4, including
            // native faces used by more than one original body.
            std::map<Tri,int> facetBody;std::vector<int> cellBody(value.cellCount(),-1);
            for(Index i=0;i<value.cellCount();++i){budget.tick();const auto& c=value.cell(i);if(c.type!=mesh::CellType::Tet4)continue;
                const auto cls=remap(c.onEntity,tags);need(cls.dim==3,"component_tet_body_class_invalid");cellBody[i]=cls.tag;const auto* p=value.cellNodes(i);
                for(auto face:std::array<Tri,4>{Tri{p[0],p[1],p[2]},Tri{p[0],p[1],p[3]},Tri{p[0],p[2],p[3]},Tri{p[1],p[2],p[3]}}){
                    const auto [it,inserted]=facetBody.emplace(sorted(face),cls.tag);need(inserted||it->second==cls.tag,"component_shared_body_node_interface");}}
            const auto extractionMetadata=add(localMetadata,add(mul(facetBody.size(),96),mul(value.cellCount(),sizeof(int))));retained(add(extractionMetadata,payload(value)));
            std::vector<std::set<Index>> used(count+1);
            for(Index i=0;i<value.cellCount();++i){budget.tick();const auto& c=value.cell(i);const auto* p=value.cellNodes(i);
                if(c.type==mesh::CellType::Tri3){const auto it=facetBody.find(sorted(Tri{p[0],p[1],p[2]}));need(it!=facetBody.end(),"component_orphan_boundary_triangle");cellBody[i]=it->second;
                    const auto cls=remap(c.onEntity,tags);need(cls.dim==2,"component_boundary_class_invalid");bool owner=false;for(auto use:original.faceOwners(cls.tag))owner|=use.solidTag==cellBody[i];need(owner,"component_boundary_body_owner_mismatch");}
                else need(c.type==mesh::CellType::Tet4,"component_unsupported_cell");
                for(int j=0;j<mesh::nodesPerCell(c.type);++j)used.at(cellBody[i]).insert(p[j]);}
            std::uint64_t extractionSets=0;for(const auto& ids:used)extractionSets=add(extractionSets,mul(ids.size(),128));retained(add(add(extractionMetadata,extractionSets),payload(value)));
            for(int body:bodies){std::map<Index,Index> ids;for(Index i:used[body]){budget.tick();ids[i]=bodyMeshes[body].addNode(value.node(i),remap(value.nodeClass(i),tags));}
                for(Index i=0;i<value.cellCount();++i)if(cellBody[i]==body){budget.tick();std::array<Index,4> ns{};const auto& c=value.cell(i);const int n=mesh::nodesPerCell(c.type);
                    for(int j=0;j<n;++j)ns[j]=ids.at(value.cellNodes(i)[j]);bodyMeshes[body].addCell(c.type,std::span<const Index>(ns.data(),n),remap(c.onEntity,tags));}
                retained(add(add(extractionMetadata,extractionSets),payload(value)));}
        }
        for(const auto* body:original.solids()){
            OcctStageTrace validating("cad_validation",body->tag());
            const int bodyTag=body->tag();auto& value=bodyMeshes[bodyTag];auto& report=reports[bodyTag];need(report.solidTag==bodyTag,"component_missing_original_body");
            if(report.surfaceChordErrorLimit<=0){Vec3 lo,hi;body->bbox(lo,hi);
                report.surfaceChordErrorLimit=options.maximumChordError>0?options.maximumChordError:
                    std::min(options.geometry.targetSize*.05,(hi-lo).norm()*.0025);}
            // An owned-copy surface evaluation can differ from the original
            // evaluation by a few nanometres at large world coordinates. An
            // OUT face node may use the exact original trimmed-face evaluation
            // only within its existing classification allowance. Interior hubs,
            // curves and vertices are never moved. Every mesh/native/exterior
            // and interbody check below runs AFTER this private proposal.
            paacs::geo::detail::OcctSolidClassifier originalClassifier(OcctModelAccess::solid(original,bodyTag));
            const double bodyTolerance=original.tolerance(3,bodyTag);
            for(Index i=0;i<value.nodeCount();++i){budget.tick();const auto cls=value.nodeClass(i);if(cls.dim!=2)continue;
                const auto p=value.node(i);originalClassifier.Perform(gp_Pnt(p.x,p.y,p.z),bodyTolerance);
                if(originalClassifier.State()!=TopAbs_OUT)continue;
                bool owner=false;for(const auto use:original.faceOwners(cls.tag))owner|=use.solidTag==bodyTag;
                need(owner,"component_original_face_owner_mismatch");Vec3 q;budget.tick();
                need(closestPointOnOcctFace(original,cls.tag,p,q),"component_original_face_coordinate_query_failed");
                const double faceTolerance=original.tolerance(2,cls.tag);
                need(std::isfinite(q.x)&&std::isfinite(q.y)&&std::isfinite(q.z)&&(p-q).norm()<=std::max(8*faceTolerance,1e-9),
                     "component_original_face_coordinate_drift");
                BRepClass_FaceClassifier faceClassifier(OcctModelAccess::face(original,cls.tag),gp_Pnt(q.x,q.y,q.z),faceTolerance);budget.tick();
                need(faceClassifier.State()==TopAbs_IN||faceClassifier.State()==TopAbs_ON,"component_original_face_coordinate_outside_trim");
                originalClassifier.Perform(gp_Pnt(q.x,q.y,q.z),bodyTolerance);budget.tick();
                need(originalClassifier.State()==TopAbs_IN||originalClassifier.State()==TopAbs_ON,"component_original_face_coordinate_outside_body");
                value.setNode(i,q);++originalFaceCoordinateConversions;
            }
            mesh::LinearVolumeAudit bodyAudit;audit(value,budget,bodyAudit);const auto q=quality(value,budget,retained());const double cad=volume(OcctModelAccess::solid(original,bodyTag),budget);
            need(std::abs(q.totalVolume-cad)<=cad*options.maximumRelativeVolumeError,"component_original_body_volume_error");
            std::vector<Vec3> points;points.reserve(value.nodeCount());for(Index i=0;i<value.nodeCount();++i){budget.tick();points.push_back(value.node(i));}
            OcctNativeValidationOutcome nodes;nativePoints(original,bodyTag,points,retained(mul(points.capacity(),sizeof(Vec3))),budget,nodes);points.clear();
            std::set<int> vertices,faces;mesh::Mesh exterior;std::map<Index,Index> exteriorIds;
            for(Index i=0;i<value.nodeCount();++i){budget.tick();const auto cls=value.nodeClass(i);const auto p=value.node(i);Vec3 projected{};
                if(cls.dim==3){need(cls.tag==bodyTag,"component_original_interior_class_invalid");continue;}
                if(cls.dim==0){vertices.insert(cls.tag);projected=original.vertex(cls.tag)->position();}
                else if(cls.dim==1)projected=original.edge(cls.tag)->closestPoint(p);
                else need(cls.dim==2&&closestPointOnOcctFace(original,cls.tag,p,projected),"component_original_node_projection_failed");
                const auto error=(p-projected).norm();report.maximumBoundaryNodeError=std::max(report.maximumBoundaryNodeError,error);
                need(error<=std::max(8*original.tolerance(cls.dim,cls.tag),1e-9),"component_original_node_geometry_error");}
            for(Index i=0;i<value.cellCount();++i){budget.tick();const auto& c=value.cell(i);const auto* p=value.cellNodes(i);
                if(c.type==mesh::CellType::Tet4){Vec3 center{};for(int j=0;j<4;++j)center=center+value.node(p[j]);points.push_back(center*.25);}
                else {need(c.type==mesh::CellType::Tri3,"component_unsupported_boundary_cell");faces.insert(c.onEntity.tag);
                    const auto a=value.node(p[0]),b=value.node(p[1]),d=value.node(p[2]);for(auto sample:std::array{(a+b+d)/3.,(a+b)*.5,(b+d)*.5,(d+a)*.5}){
                        budget.tick();Vec3 projected;need(closestPointOnOcctFace(original,c.onEntity.tag,sample,projected),"component_original_facet_projection_failed");
                        const auto error=(sample-projected).norm();report.maximumBoundarySampleError=std::max(report.maximumBoundarySampleError,error);
                        need(report.surfaceChordErrorLimit>0&&error<=report.surfaceChordErrorLimit+original.tolerance(2,c.onEntity.tag),"component_original_facet_chord_error");}
                    Tri ns{};for(int j=0;j<3;++j){const auto [it,inserted]=exteriorIds.emplace(p[j],exterior.nodeCount());if(inserted)exterior.addNode(value.node(p[j]),value.nodeClass(p[j]));ns[j]=it->second;}
                    exterior.addCell(mesh::CellType::Tri3,ns,c.onEntity);}
            }
            const auto validationRetained=retained(add(mul(points.capacity(),sizeof(Vec3)),add(payload(exterior),mul(exteriorIds.size(),128))));
            OcctNativeValidationOutcome centroids;nativePoints(original,bodyTag,points,validationRetained,budget,centroids);
            std::set<int> expectedVertices,expectedFaces;for(const auto use:original.faceUses(bodyTag)){expectedFaces.insert(use.faceTag);
                for(const auto& wire:original.boundaryUses(use.faceTag))for(const auto& edge:wire.edges){expectedVertices.insert(edge.startVertexTag);expectedVertices.insert(edge.endVertexTag);}}
            need(vertices==expectedVertices&&faces==expectedFaces,"component_original_boundary_entities_missing");
            need(budget.remaining()>0,"component_work_budget",NativeCadMeshReason::BudgetExceeded);
            OcctStarOptions surface;surface.maximumWork=budget.remaining();surface.maximumBoundaryNodes=options.geometry.maximumNodes;
            surface.maximumBoundaryTriangles=options.geometry.maximumCells/2;surface.maximumWorkingBytes=options.maximumPrimitiveWorkingBytes-validationRetained;surface.cancellation=options.geometry.cancellation;
            const auto certificate=certifyOcctBoundarySurface(exterior,surface);budget.collect(certificate.workUsed);
            if(certificate.reason==OcctStarReason::WorkBudget||certificate.reason==OcctStarReason::OutputBudget)throw Stop{NativeCadMeshReason::BudgetExceeded,"component_exterior_budget"};
            if(certificate.reason==OcctStarReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled,"cancelled"};
            if(certificate.reason==OcctStarReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed,"allocation_failed"};
            need(certificate.certified,"component_original_exterior_not_certified");
            report.cadVolume=cad;report.meshVolume=q.totalVolume;report.relativeVolumeError=std::abs(q.totalVolume-cad)/cad;report.reason=NativeCadMeshReason::None;
            report.nativeCentroidThreads=centroids.threadsUsed;report.nativeCentroidsChecked=centroids.nativeQueries;report.nativeCentroidReferenceQueries=centroids.originalReferenceQueries;
            report.nativeCentroidGeometryPairs=centroids.copiedGeometryValuePairs;report.nativeCentroidSourceGuard=nodes.sourceModelUnmodified&&centroids.sourceModelUnmodified;
            report.boundaryEmbeddingCertified=true;
            retained(add(mul(points.capacity(),sizeof(Vec3)),add(payload(exterior),mul(exteriorIds.size(),128))));
            append(joined,value,budget);out.bodies.push_back(std::move(report));
        }
        OcctStageTrace validating("contact_validation");
        std::uint64_t treeBytes=0;for(const auto& value:bodyMeshes){std::uint64_t triangles=0;for(Index i=0;i<value.cellCount();++i){budget.tick();if(value.cell(i).type==mesh::CellType::Tri3)++triangles;}
            treeBytes=add(treeBytes,mul(triangles,sizeof(Triangle)+2*sizeof(TreeNode)));}retained(treeBytes);
        std::vector<Tree> trees;for(unsigned i=1;i<=count;++i)trees.push_back(tree(bodyMeshes[i],budget));std::uint64_t matched=0;
        for(unsigned i=0;i<count;++i)for(unsigned j=0;j<i;++j)contacts(trees[i],0,trees[j],0,facePairs.contains({j+1,i+1}),budget,matched);
        need(matched>0,"component_native_contact_triangles_missing",NativeCadMeshReason::UnresolvedBodyInterface);
        audit(joined,budget,out.audit);out.quality=quality(joined,budget,retained(treeBytes));
        need(out.quality.minimumQuality>=options.geometry.minimumElementQuality,"component_quality_requirement",NativeCadMeshReason::QualityRequirementNotMet);
        out.qualityWarning=out.quality.minimumQuality<options.geometry.warningElementQuality;
        std::string after;const auto finalSourceReason=snapshotOcctContactSource(original,snapshotLimit,after);budget.tick();need(finalSourceReason==NativeCadMeshReason::None,"component_source_audit_failed",finalSourceReason);
        sourceAudited=true;sourceUnchanged=before==after;need(sourceUnchanged,"component_original_source_changed",NativeCadMeshReason::InternalError);budget.tick(0);
        out.separateContactPairs=facePairs.size();out.matchedNativeContactTriangles=matched;
        destination.swap(joined);out.published=true;out.destinationUnchanged=false;out.reason=NativeCadMeshReason::None;
        diagnostic=originalFaceCoordinateConversions?"native_contact_components_original_face_coordinates":"native_contact_components";
    }catch(const Stop& s){out.reason=s.reason;diagnostic=s.diagnostic;}
    catch(const Standard_OutOfMemory&){out.reason=NativeCadMeshReason::AllocationFailed;diagnostic="allocation_failed";}
    catch(const std::bad_alloc&){out.reason=NativeCadMeshReason::AllocationFailed;diagnostic="allocation_failed";}
    catch(...){out.reason=NativeCadMeshReason::InternalError;diagnostic="component_internal_error";}
    if(!out.published&&!before.empty()){std::string after;const auto result=snapshotOcctContactSource(original,snapshotLimit,after);sourceAudited=result==NativeCadMeshReason::None;sourceUnchanged=sourceAudited&&before==after;
        if(sourceAudited&&!sourceUnchanged){out.reason=NativeCadMeshReason::InternalError;diagnostic="component_original_source_changed";}}
    if(!out.published)for(auto& report:out.bodies)report.reason=out.reason;
    out.validationWork=budget.work;out.nativeContactMaximumWorkingPayloadEstimate=budget.maximumWorking;
    out.nativeContactSourceGuard=sourceUnchanged;out.nativeContactSourceAuditAvailable=sourceAudited;out.nativeContactDiagnostic=diagnostic;
    out.totalSeconds=std::chrono::duration<double>(Clock::now()-started).count();return out;
}

NativeCadMeshOutcome buildOcctIndependentParts(const OcctModel& original,const NativeCadMeshOptions& options,mesh::Mesh& destination) noexcept {
    NativeCadMeshOutcome out;out.assemblyPolicy=NativeCadAssemblyPolicy::IndependentParts;
    const auto started=Clock::now();Budget budget{options,options.maximumValidationWork};std::string before;
    const auto snapshotLimit=std::min<std::uint64_t>(64ULL*1024*1024,options.maximumPrimitiveWorkingBytes/8);
    const char* diagnostic="independent_parts_not_started";
    try {
        budget.tick();need(options.assemblyPolicy==NativeCadAssemblyPolicy::IndependentParts&&options.elements==algo::PrimitiveElements::Tet4,"independent_parts_tet4_required",NativeCadMeshReason::UnsupportedControl);
        need(!original.solids().empty()&&snapshotLimit>0,"independent_parts_invalid_request",NativeCadMeshReason::InvalidRequest);
        need(std::isfinite(options.geometry.targetSize)&&options.geometry.targetSize>0&&options.maximumNativeValidationThreads>=1&&options.maximumNativeValidationThreads<=8&&std::isfinite(options.maximumRelativeVolumeError)&&options.maximumRelativeVolumeError>0&&options.maximumRelativeVolumeError<1&&std::isfinite(options.maximumChordError)&&options.maximumChordError>=0,"independent_parts_invalid_request",NativeCadMeshReason::InvalidRequest);
        need(std::isfinite(options.minimumSize)&&options.minimumSize>=0&&options.minimumSize<=options.geometry.targetSize,"independent_parts_minimum_size",NativeCadMeshReason::InvalidRequest);
        need(options.minimumSize==0||!options.faceSizes.empty(),"independent_parts_minimum_size_requires_face_control",NativeCadMeshReason::UnsupportedControl);
        std::set<int> selectedFaces,selectedBodies;for(const auto control:options.faceSizes)need(original.face(control.faceTag)&&!original.faceOwners(control.faceTag).empty()&&std::isfinite(control.size)&&control.size>0&&std::isfinite(control.growth)&&control.growth>=1&&selectedFaces.insert(control.faceTag).second,"independent_parts_face_control",NativeCadMeshReason::InvalidRequest);for(const auto& control:options.geometry.bodyControls)need(original.solid(control.solidTag)&&selectedBodies.insert(control.solidTag).second,"independent_parts_body_control",NativeCadMeshReason::InvalidRequest);
        const auto sourceReason=snapshotOcctContactSource(original,snapshotLimit,before);budget.tick();need(sourceReason==NativeCadMeshReason::None,"independent_parts_source_audit",sourceReason);
        mesh::Mesh joined;ShapeMap originalMaps[4];const TopAbs_ShapeEnum kinds[]{TopAbs_VERTEX,TopAbs_EDGE,TopAbs_FACE,TopAbs_SOLID};
        std::uint64_t metadata=0;for(unsigned d=0;d<4;++d){TopExp::MapShapes(OcctModelAccess::shape(original),kinds[d],originalMaps[d]);metadata=add(metadata,mul(originalMaps[d].Extent()+1,256));}
        const auto retained=[&](std::uint64_t active=0){const auto bytes=add(add(mul(before.size(),3),metadata),add(payload(joined),active));budget.working(bytes);return bytes;};retained();
        // Interference is descriptive in this explicit mode. An inconclusive
        // Boolean is reported as unresolved, never as proof of disjointness.
        OcctStageTrace interferenceTrace("contact_interference");
        for(std::size_t i=0;i<original.solids().size();++i)for(std::size_t j=0;j<i;++j){
            budget.tick();const auto* a=original.solids()[j];const auto* b=original.solids()[i];NativeCadInterference interference;interference.firstSolidTag=a->tag();interference.secondSolidTag=b->tag();
            Vec3 al,ah,bl,bh;a->bbox(al,ah);b->bbox(bl,bh);
            if(!overlaps(al,ah,bl,bh))interference.state=NativeCadInterferenceState::NoMaterialOverlap;
            else try {
                BRepBuilderAPI_Copy ca(OcctModelAccess::solid(original,a->tag()),true,false),cb(OcctModelAccess::solid(original,b->tag()),true,false);need(ca.IsDone()&&cb.IsDone(),"independent_parts_copy_failed");
                const auto material=common(ca.Shape(),cb.Shape(),budget);bool solid=false,valid=BRepCheck_Analyzer(material).IsValid();double mass=0;
                for(TopExp_Explorer k(material,TopAbs_SOLID);k.More();k.Next()){solid=true;if(valid)mass+=volume(k.Current(),budget);}
                if(valid&&!solid)interference.state=NativeCadInterferenceState::NoMaterialOverlap;
                else if(valid&&std::isfinite(mass)&&mass>0){interference.state=NativeCadInterferenceState::MaterialOverlap;interference.nativeOverlapVolume=mass;}
            }catch(const Stop& s){if(s.reason==NativeCadMeshReason::Cancelled||s.reason==NativeCadMeshReason::BudgetExceeded||s.reason==NativeCadMeshReason::AllocationFailed)throw;}
            catch(const Standard_OutOfMemory&){throw;}catch(const Standard_Failure&){}
            out.interferences.push_back(interference);metadata=add(metadata,sizeof(NativeCadInterference));retained();
        }
        interferenceTrace.finish();
        out.preparationSeconds=std::chrono::duration<double>(Clock::now()-started).count();NativeCadMeshReason failed=NativeCadMeshReason::None;
        for(std::size_t ordinal=0;ordinal<original.solids().size();++ordinal){
            OcctStageTrace part("contact_part_preparation",original.solids()[ordinal]->tag());
            const auto prepareStart=Clock::now();budget.tick();const int body=original.solids()[ordinal]->tag();BRepBuilderAPI_Copy copy(OcctModelAccess::solid(original,body),true,false);need(copy.IsDone(),"independent_parts_copy_failed");
            auto local=OcctModelAccess::fromShape(copy.Shape());need(local.imported()&&local.model->solids().size()==1,"independent_parts_import_failed");
            ShapeMap localMaps[4];std::array<std::vector<int>,4> tags;std::uint64_t localMetadata=0;
            for(unsigned d=0;d<4;++d){TopExp::MapShapes(copy.Shape(),kinds[d],localMaps[d]);tags[d].resize(localMaps[d].Extent()+1,-1);localMetadata=add(localMetadata,mul(tags[d].size(),256));
                ShapeMap bodyMap;TopExp::MapShapes(OcctModelAccess::solid(original,body),kinds[d],bodyMap);
                for(int t=1;t<=bodyMap.Extent();++t){budget.tick();const int originalTag=originalMaps[d].FindIndex(bodyMap(t));const auto mapped=copy.ModifiedShape(bodyMap(t));need(originalTag>0&&!mapped.IsNull(),"independent_parts_history_incomplete");const int id=localMaps[d].FindIndex(mapped);need(id>0,"independent_parts_history_incomplete");tags[d][id]=originalTag;}
                for(std::size_t t=1;t<tags[d].size();++t)need(tags[d][t]>0,"independent_parts_history_incomplete");}
            const auto live=retained(localMetadata);need(live<options.maximumPrimitiveWorkingBytes&&budget.remaining()>0,"independent_parts_working_budget",NativeCadMeshReason::BudgetExceeded);
            auto request=options;request.assemblyPolicy=NativeCadAssemblyPolicy::ConformingDomain;request.separateMatchedNativeFaceContacts=false;request.separateMatchedBoxFaceContacts=false;request.separateLowerDimensionalContacts=false;
            request.maximumValidationWork=std::min(options.geometry.maximumMeshingWork,budget.remaining()/(original.solids().size()-ordinal));request.geometry.maximumMeshingWork=request.maximumValidationWork;
            need(request.maximumValidationWork>0,"independent_parts_work_budget",NativeCadMeshReason::BudgetExceeded);request.maximumPrimitiveWorkingBytes=options.maximumPrimitiveWorkingBytes-live;
            // Selection identities are local only during the child call.
            request.faceSizes.clear();for(const auto control:options.faceSizes){for(std::size_t t=1;t<tags[2].size();++t)if(tags[2][t]==control.faceTag)request.faceSizes.push_back({static_cast<int>(t),control.size,control.growth});}
            // The clamp belongs to selected-face sizing. Unselected parts use
            // their unchanged global size and have no local field to clamp.
            if(request.faceSizes.empty())request.minimumSize=0;
            request.geometry.bodyControls.clear();for(auto control:options.geometry.bodyControls)if(control.solidTag==body){control.solidTag=1;request.geometry.bodyControls.push_back(control);}
            out.preparationSeconds+=std::chrono::duration<double>(Clock::now()-prepareStart).count();
            part.next("contact_part_meshing");
            mesh::Mesh value;auto child=meshNativeCad(*local.model,request,value);
            part.finish(nativeCadMeshReasonName(child.reason));
            budget.collect(child.validationWork);out.preparationSeconds+=child.preparationSeconds;out.generationSeconds+=child.generationSeconds;out.validationSeconds+=child.validationSeconds;out.auditSeconds+=child.auditSeconds;out.qualitySeconds+=child.qualitySeconds;
            if(child.bodies.empty()){NativeCadBodyMeshReport report;report.solidTag=body;report.reason=child.reason;out.bodies.push_back(report);}
            for(auto& report:child.bodies){report.solidTag=body;if(report.startFaceTag>0)report.startFaceTag=tags[2].at(report.startFaceTag);if(report.endFaceTag>0)report.endFaceTag=tags[2].at(report.endFaceTag);if(report.contactSourceEdgeTag>0)report.contactSourceEdgeTag=tags[1].at(report.contactSourceEdgeTag);if(report.contactTargetEdgeTag>0)report.contactTargetEdgeTag=tags[1].at(report.contactTargetEdgeTag);out.bodies.push_back(std::move(report));}
            if(!child.published){
                need(child.destinationUnchanged,"independent_parts_child_transaction",NativeCadMeshReason::InternalError);
                if(child.reason==NativeCadMeshReason::Cancelled||child.reason==NativeCadMeshReason::AllocationFailed||child.reason==NativeCadMeshReason::InternalError||child.reason==NativeCadMeshReason::InvalidRequest||child.reason==NativeCadMeshReason::QualityRequirementNotMet||child.reason==NativeCadMeshReason::SizeRequirementNotMet)throw Stop{child.reason,"independent_parts_child_refused"};
                if(failed==NativeCadMeshReason::None)failed=child.reason;continue;
            }
            part.next("cad_validation");
            const auto validateStart=Clock::now();need(child.bodies.size()==1&&!out.bodies.empty(),"independent_parts_child_report",NativeCadMeshReason::InternalError);auto& report=out.bodies.back();retained(add(localMetadata,mul(payload(value),2)));
            {mesh::Mesh classified;for(Index i=0;i<value.nodeCount();++i){budget.tick();classified.addNode(value.node(i),remap(value.nodeClass(i),tags));}
            for(Index i=0;i<value.cellCount();++i){budget.tick();const auto& c=value.cell(i);classified.addCell(c.type,std::span<const Index>(value.cellNodes(i),mesh::nodesPerCell(c.type)),remap(c.onEntity,tags));}value.swap(classified);}
            std::vector<Vec3> centers,interior;for(Index i=0;i<value.nodeCount();++i){budget.tick();const auto cls=value.nodeClass(i);if(cls.dim==3){need(cls.tag==body,"independent_parts_interior_identity");interior.push_back(value.node(i));}}
            for(Index i=0;i<value.cellCount();++i){budget.tick();const auto& c=value.cell(i);need(c.type==mesh::CellType::Tri3||c.type==mesh::CellType::Tet4,"independent_parts_cell_type");if(c.type==mesh::CellType::Tet4){need(c.onEntity.tag==body,"independent_parts_cell_identity");Vec3 p{};for(unsigned j=0;j<4;++j)p=p+value.node(value.cellNodes(i)[j])*.25;centers.push_back(p);}}
            const auto validationLive=retained(add(localMetadata,add(payload(value),mul(centers.capacity()+interior.capacity(),sizeof(Vec3)))));OcctNativeValidationOutcome nodes,centroids;
            const bool checkedPlcBand=options.allowCentroidsWithinChordTolerance&&std::string_view(report.route)=="native_plc_cavity"&&
                report.boundaryEmbeddingCertified&&report.centroidOutsideDistanceLimit>0&&
                report.centroidOutsideDistanceLimit<=report.surfaceChordErrorLimit&&report.nativeChordSamplesChecked==4*centers.size();
            const double allowance=checkedPlcBand?report.centroidOutsideDistanceLimit:0;
            if(!interior.empty())nativePoints(original,body,interior,validationLive,budget,nodes);nativePoints(original,body,centers,validationLive,budget,centroids,allowance);
            report.nativeCentroidThreads=centroids.threadsUsed;report.nativeCentroidsChecked=centroids.nativeQueries;report.nativeCentroidReferenceQueries=centroids.originalReferenceQueries;report.nativeCentroidGeometryPairs=centroids.copiedGeometryValuePairs;report.nativeCentroidSourceGuard=centroids.sourceModelUnmodified&&(interior.empty()||nodes.sourceModelUnmodified);
            if(checkedPlcBand){
                report.centroidsWithinChordTolerance=centroids.pointsWithinDistanceBand;report.nativeCentroidDistanceQueries=centroids.nativeDistanceQueries;report.maximumOutsideCentroidDistance=centroids.maximumOutsideDistance;
                std::vector<Vec3> samples;need(centers.size()<=(options.maximumPrimitiveWorkingBytes-validationLive)/(4*sizeof(Vec3)),"independent_parts_chord_sample_storage",NativeCadMeshReason::BudgetExceeded);samples.reserve(centers.size()*4);
                for(Index i=0;i<value.cellCount();++i)if(value.cell(i).type==mesh::CellType::Tet4){const auto* p=value.cellNodes(i);Vec3 center{};for(unsigned j=0;j<4;++j)center=center+value.node(p[j])*.25;
                    for(unsigned j=0;j<4;++j){budget.tick();samples.push_back(center*.5+value.node(p[j])*.5);}}
                OcctNativeValidationOutcome additional;nativePoints(original,body,samples,add(validationLive,mul(samples.capacity(),sizeof(Vec3))),budget,additional,allowance);
                report.nativeChordSamplesChecked=additional.nativeQueries;report.chordSamplesWithinTolerance=additional.pointsWithinDistanceBand;report.maximumOutsideChordSampleDistance=additional.maximumOutsideDistance;report.nativeCentroidSourceGuard=report.nativeCentroidSourceGuard&&additional.sourceModelUnmodified;
            }
            append(joined,value,budget);++out.validatedParts;retained();out.validationSeconds+=std::chrono::duration<double>(Clock::now()-validateStart).count();
        }
        OcctStageTrace validating("contact_validation");
        need(failed==NativeCadMeshReason::None,"independent_parts_child_refused",failed);need(out.validatedParts==original.solids().size(),"independent_parts_missing_part");
        const auto auditStart=Clock::now();audit(joined,budget,out.audit);out.auditSeconds+=std::chrono::duration<double>(Clock::now()-auditStart).count();const auto qualityStart=Clock::now();out.quality=quality(joined,budget,retained());out.qualitySeconds+=std::chrono::duration<double>(Clock::now()-qualityStart).count();need(out.quality.minimumQuality>=options.geometry.minimumElementQuality,"independent_parts_quality_requirement",NativeCadMeshReason::QualityRequirementNotMet);out.qualityWarning=out.quality.minimumQuality<options.geometry.warningElementQuality;
        std::string after;const auto end=snapshotOcctContactSource(original,snapshotLimit,after);budget.tick();need(end==NativeCadMeshReason::None,"independent_parts_source_audit",end);out.nativeContactSourceAuditAvailable=true;out.nativeContactSourceGuard=before==after;need(out.nativeContactSourceGuard,"independent_parts_source_changed",NativeCadMeshReason::InternalError);budget.tick(0);
        destination.swap(joined);out.published=true;out.destinationUnchanged=false;out.reason=NativeCadMeshReason::None;diagnostic="independent_parts";
    }catch(const Stop& s){out.reason=s.reason;diagnostic=s.diagnostic;}
    catch(const Standard_OutOfMemory&){out.reason=NativeCadMeshReason::AllocationFailed;diagnostic="allocation_failed";}
    catch(const std::bad_alloc&){out.reason=NativeCadMeshReason::AllocationFailed;diagnostic="allocation_failed";}
    catch(...){out.reason=NativeCadMeshReason::InternalError;diagnostic="independent_parts_internal_error";}
    if(!out.nativeContactSourceAuditAvailable&&!before.empty()){std::string after;const auto reason=snapshotOcctContactSource(original,snapshotLimit,after);out.nativeContactSourceAuditAvailable=reason==NativeCadMeshReason::None;out.nativeContactSourceGuard=out.nativeContactSourceAuditAvailable&&before==after;if(out.nativeContactSourceAuditAvailable&&!out.nativeContactSourceGuard){out.reason=NativeCadMeshReason::InternalError;diagnostic="independent_parts_source_changed";}}
    out.nativeContactDiagnostic=diagnostic;out.validationWork=budget.work;out.nativeContactMaximumWorkingPayloadEstimate=budget.maximumWorking;out.totalSeconds=std::chrono::duration<double>(Clock::now()-started).count();return out;
}
}
