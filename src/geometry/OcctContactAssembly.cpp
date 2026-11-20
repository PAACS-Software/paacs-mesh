// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// One copied native assembly supplies shared, authoritative interfaces.
#include "OcctModelInternal.hpp"
#include "OcctBoundaryInternal.hpp"
#include "OcctPlcFrontRecovery.hpp"
#include "OcctStarMesher.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "OcctContactAssembly.hpp"
#include "OcctStageTrace.hpp"
#include "paacs/bnd/TetFacetGeometry.hpp"
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include "OcctSolidClassifier.hpp"
#include <BRepGProp.hxx>
#include <BOPAlgo_PaveFiller.hxx>
#include <BOPDS_DS.hxx>
#include <BRepTools.hxx>
#include <GProp_GProps.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Pln.hxx>
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <map>
#include <numbers>
#include <numeric>
#include <set>
#include <ostream>
#include <streambuf>
#include <string_view>
#include <Standard_OutOfMemory.hxx>
#include <vector>
namespace paacs::geo::detail {
namespace {
using namespace paacs;
using ShapeMap=NCollection_IndexedMap<TopoDS_Shape,TopTools_ShapeMapHasher>;
using Class=mesh::Classification;
using Tri=std::array<Index,3>;
using Clock=std::chrono::steady_clock;
struct Stop {const char* diagnostic;};
void need(bool b,const char* s){if(!b)throw Stop{s};}
void checkStarResource(OcctStarReason reason){
    if(reason==OcctStarReason::WorkBudget||reason==OcctStarReason::OutputBudget)throw Stop{"child_budget"};
    if(reason==OcctStarReason::Cancelled)throw Stop{"cancelled"};
    if(reason==OcctStarReason::AllocationFailed)throw Stop{"allocation_failed"};
}
class SnapshotBuffer final:public std::streambuf {
    std::uint64_t maximum_;
public:
    std::string value;
    explicit SnapshotBuffer(std::uint64_t maximum):maximum_(maximum){}
protected:
    std::streamsize xsputn(const char* p,std::streamsize n) override {
        need(n>=0&&static_cast<std::uint64_t>(n)<=maximum_-value.size(),"snapshot_budget");
        value.append(p,static_cast<std::size_t>(n));return n;
    }
    int_type overflow(int_type c) override {
        if(traits_type::eq_int_type(c,traits_type::eof()))return traits_type::not_eof(c);
        const char byte=traits_type::to_char_type(c);xsputn(&byte,1);return c;
    }
};
std::string snapshot(const TopoDS_Shape& s,std::uint64_t maximum){
    SnapshotBuffer buffer(maximum);std::ostream stream(&buffer);
    stream.exceptions(std::ios::badbit|std::ios::failbit);BRepTools::Write(s,stream);return std::move(buffer.value);
}
double accurateMass(const TopoDS_Shape& s){
    GProp_GProps coarse,fine;
    const double ce=BRepGProp::VolumePropertiesGK(s,coarse,1e-8,false,true,false,false,false);
    const double fe=BRepGProp::VolumePropertiesGK(s,fine,1e-10,false,true,false,false,false);
    need(std::isfinite(ce)&&ce>=0&&std::isfinite(fe)&&fe>=0&&std::isfinite(coarse.Mass())&&coarse.Mass()>0&&std::isfinite(fine.Mass())&&fine.Mass()>0,"native_integration_failed");
    need(std::abs(coarse.Mass()-fine.Mass())<=fine.Mass()*1e-6,"native_integration_not_stable");
    return fine.Mass();
}
double mass(const TopoDS_Shape& s){GProp_GProps p;const auto error=BRepGProp::VolumeProperties(s,p,1e-10,false,false);need(std::isfinite(error)&&error>=0&&p.Mass()>0&&std::isfinite(p.Mass()),"native_integration_failed");return p.Mass();}
bool same(Vec3 a,Vec3 b){return std::bit_cast<std::uint64_t>(a.x)==std::bit_cast<std::uint64_t>(b.x)&&std::bit_cast<std::uint64_t>(a.y)==std::bit_cast<std::uint64_t>(b.y)&&std::bit_cast<std::uint64_t>(a.z)==std::bit_cast<std::uint64_t>(b.z);}
Tri ordered(Tri t){std::sort(t.begin(),t.end());return t;}
bool sameWinding(Tri a,Tri b){for(int i=0;i<3;++i)if(a[0]==b[i]&&a[1]==b[(i+1)%3]&&a[2]==b[(i+2)%3])return true;return false;}
struct Receipt{const NativeCadMeshOptions* options=nullptr;std::uint64_t limit=0;std::uint64_t pieces=0,planes=0,work=0,nodes=0,cells=0,interfaces=0,contactTriangles=0,matchedTriangles=0;
    int failedPiece=-1;bool published=false,sourceUnchanged=false,sourceAuditAvailable=false;const char* reason="not_started";
    std::uint64_t maximumWorkingPayloadEstimate=0;
    double meshVolume=0,cadVolume=0,minQuality=0,volumeChange=0,maxNodeError=0,maxSampleError=0;
    double partitionSeconds=0,surfaceSeconds=0,starSeconds=0,assemblySeconds=0,validationSeconds=0;};
struct Timer {double& value;Clock::time_point start=Clock::now();bool active=true;
    void finish(){if(active){value+=std::chrono::duration<double>(Clock::now()-start).count();active=false;}}
    ~Timer(){finish();}};
std::uint64_t remaining(const Receipt& r){return r.limit-r.work;}
void charge(Receipt& r,std::uint64_t n=1){
    need(!r.options->geometry.cancellation||!r.options->geometry.cancellation->load(std::memory_order_relaxed),"cancelled");
    need(n<=remaining(r),"work_budget");r.work+=n;
}
void collect(Receipt& r,std::uint64_t n){
    // The child already spent these visits, including on cancellation.
    need(n<=remaining(r),"work_budget");r.work+=n;charge(r,0);
}
std::uint64_t add(std::uint64_t a,std::uint64_t b){need(b<=UINT64_MAX-a,"payload_budget");return a+b;}
std::uint64_t mul(std::uint64_t a,std::uint64_t b){need(b==0||a<=UINT64_MAX/b,"payload_budget");return a*b;}
std::uint64_t payload(const mesh::Mesh& m){return add(add(mul(m.nodeCount(),sizeof(Vec3)+sizeof(Class)),mul(m.cellCount(),sizeof(mesh::Cell))),mul(m.connectivityCount(),sizeof(Index)));}
void output(const mesh::Mesh& m,const Receipt& r){
    need(m.nodeCount()<=r.options->geometry.maximumNodes&&m.cellCount()<=r.options->geometry.maximumCells&&payload(m)<=r.options->geometry.maximumOutputBytes,"output_budget");
}
void prospectiveOutput(std::uint64_t n,std::uint64_t c,std::uint64_t cn,const Receipt& r){
    const auto bytes=add(add(mul(n,sizeof(Vec3)+sizeof(Class)),mul(c,sizeof(mesh::Cell))),mul(cn,sizeof(Index)));
    need(n<=r.options->geometry.maximumNodes&&c<=r.options->geometry.maximumCells&&bytes<=r.options->geometry.maximumOutputBytes,"output_budget");
}
struct FaceUse{int piece,body,face;Tri canonical;std::array<Vec3,3> p;};
struct RetryCut {Vec3 low,high;};
struct RetrySurface {};
struct BoundaryTriangle{std::array<Vec3,3> p;Vec3 lo,hi;};
struct TreeNode{Vec3 lo,hi;std::size_t begin,end,left=0,right=0;};
bool overlaps(Vec3 a,Vec3 b,Vec3 c,Vec3 d){return a.x<=d.x&&c.x<=b.x&&a.y<=d.y&&c.y<=b.y&&a.z<=d.z&&c.z<=b.z;}
double coordinate(Vec3 p,int a){return a==0?p.x:a==1?p.y:p.z;}
struct Tree{
    std::vector<BoundaryTriangle> triangles;std::vector<TreeNode> nodes;
    std::size_t build(std::size_t begin,std::size_t end,Receipt& r){
        charge(r,end-begin);
        Vec3 lo{1e300,1e300,1e300},hi{-1e300,-1e300,-1e300};
        for(std::size_t i=begin;i<end;++i){const auto& t=triangles[i];lo={std::min(lo.x,t.lo.x),std::min(lo.y,t.lo.y),std::min(lo.z,t.lo.z)};hi={std::max(hi.x,t.hi.x),std::max(hi.y,t.hi.y),std::max(hi.z,t.hi.z)};}
        std::size_t id=nodes.size();nodes.push_back({lo,hi,begin,end});
        if(end-begin>8){auto span=hi-lo;int axis=span.y>span.x?1:0;if(span.z>coordinate(span,axis))axis=2;
            auto mid=begin+(end-begin)/2;std::nth_element(triangles.begin()+begin,triangles.begin()+mid,triangles.begin()+end,[axis](const auto& a,const auto& b){return coordinate(a.lo+a.hi,axis)<coordinate(b.lo+b.hi,axis);});
            auto left=build(begin,mid,r),right=build(mid,end,r);nodes[id].left=left;nodes[id].right=right;}
        return id;
    }
};
Tree tree(const mesh::Mesh& mesh,Receipt& receipt){Tree t;
    for(Index i=0;i<mesh.cellCount();++i){const auto& c=mesh.cell(i);if(c.type!=mesh::CellType::Tri3)continue;const auto* ids=mesh.cellNodes(i);
        BoundaryTriangle p{{mesh.node(ids[0]),mesh.node(ids[1]),mesh.node(ids[2])}};
        p.lo={std::min({p.p[0].x,p.p[1].x,p.p[2].x}),std::min({p.p[0].y,p.p[1].y,p.p[2].y}),std::min({p.p[0].z,p.p[1].z,p.p[2].z})};
        p.hi={std::max({p.p[0].x,p.p[1].x,p.p[2].x}),std::max({p.p[0].y,p.p[1].y,p.p[2].y}),std::max({p.p[0].z,p.p[1].z,p.p[2].z})};t.triangles.push_back(p);}
    need(!t.triangles.empty(),"missing_body_boundary");t.build(0,t.triangles.size(),receipt);return t;}
void auditContacts(const Tree& a,std::size_t ai,const Tree& b,std::size_t bi,Receipt& r){
    charge(r);const auto& x=a.nodes[ai];const auto& y=b.nodes[bi];if(!overlaps(x.lo,x.hi,y.lo,y.hi))return;
    if(x.left){auditContacts(a,x.left,b,bi,r);auditContacts(a,x.right,b,bi,r);return;}
    if(y.left){auditContacts(a,ai,b,y.left,r);auditContacts(a,ai,b,y.right,r);return;}
    for(auto i=x.begin;i<x.end;++i)for(auto j=y.begin;j<y.end;++j){charge(r);
        const auto& p=a.triangles[i];const auto& q=b.triangles[j];if(!overlaps(p.lo,p.hi,q.lo,q.hi))continue;
        std::array<int,3> permutation{-1,-1,-1};for(int u=0;u<3;++u)for(int v=0;v<3;++v)if(same(p.p[u],q.p[v]))permutation[u]=v;
        if(permutation[0]>=0&&permutation[1]>=0&&permutation[2]>=0){int inv=(permutation[0]>permutation[1])+(permutation[0]>permutation[2])+(permutation[1]>permutation[2]);need(inv%2==1,"same_sense_contact_triangle");continue;}
        need(!bnd::tet_facet_detail::trianglesOpenIntersect(p.p,q.p)&&!bnd::tet_facet_detail::coplanarTrianglesHaveAreaOverlap(p.p,q.p),"unmatched_body_interface");
    }
}

NativeCadMeshReason failureReason(std::string_view s){
    if(s=="cancelled")return NativeCadMeshReason::Cancelled;
    if(s.find("budget")!=s.npos)return NativeCadMeshReason::BudgetExceeded;
    if(s=="allocation_failed"||s=="allocation-failed")return NativeCadMeshReason::AllocationFailed;
    if(s=="quality_requirement")return NativeCadMeshReason::QualityRequirementNotMet;
    if(s=="unsupported_control")return NativeCadMeshReason::UnsupportedControl;
    if(s.find("source_changed")!=s.npos||s.find("source-changed")!=s.npos||s=="internal-error")return NativeCadMeshReason::InternalError;
    if(s=="original_material_overlap"||s.find("contact")!=s.npos||s=="overlap_query_failed")return NativeCadMeshReason::UnresolvedBodyInterface;
    if(s.find("original_")!=s.npos||s.find("volume")!=s.npos)return NativeCadMeshReason::GeometryError;
    if(s.find("request")!=s.npos||s=="assembly_required")return NativeCadMeshReason::InvalidRequest;
    return NativeCadMeshReason::MeshingFailed;
}
} // namespace
NativeCadMeshReason snapshotOcctContactSource(const OcctModel& model,std::uint64_t maximum,std::string& destination) noexcept {
    try{auto staged=snapshot(OcctModelAccess::shape(model),maximum);destination.swap(staged);return NativeCadMeshReason::None;}
    catch(const Stop& s){return failureReason(s.diagnostic);}
    catch(const Standard_OutOfMemory&){return NativeCadMeshReason::AllocationFailed;}
    catch(const std::bad_alloc&){return NativeCadMeshReason::AllocationFailed;}
    catch(...){return NativeCadMeshReason::InternalError;}
}
NativeCadMeshOutcome buildOcctContactAssembly(const OcctModel& original,const NativeCadMeshOptions& o,mesh::Mesh& destination,
                                             const std::string* sourceBefore) noexcept {
    return buildOcctContactAssembly(original,o,OcctContactAssemblyStrategy{},destination,sourceBefore);
}
NativeCadMeshOutcome buildOcctContactAssembly(const OcctModel& original,const NativeCadMeshOptions& o,
                                             const OcctContactAssemblyStrategy& strategy,mesh::Mesh& destination,
                                             const std::string* sourceBefore) noexcept {
    NativeCadMeshOutcome out;Receipt receipt;receipt.options=&o;
    receipt.limit=std::min(o.maximumValidationWork,o.geometry.maximumMeshingWork);
    const auto begin=Clock::now();std::string before;
    const std::uint64_t snapshotLimit=std::min<std::uint64_t>(64ULL*1024*1024,o.maximumPrimitiveWorkingBytes/8);
    mesh::LinearVolumeAuditLimits finalLimits;finalLimits.maximumTemporaryBytes=o.geometry.maximumAuditTemporaryBytes;
    finalLimits.maximumFaceUses=o.geometry.maximumAuditFaceUses;finalLimits.cancellation=o.geometry.cancellation;
    try {
        charge(receipt,0);
        need(o.separateMatchedNativeFaceContacts&&o.elements==algo::PrimitiveElements::Tet4,"invalid_request");
        need(o.geometry.targetSize>0&&std::isfinite(o.geometry.targetSize)&&o.maximumChordError>=0&&std::isfinite(o.maximumChordError)&&
             o.geometry.threads>=0&&o.geometry.threads<=1024&&
             o.maximumNativeValidationThreads>=1&&o.maximumNativeValidationThreads<=8&&
             o.geometry.relativeRecognitionTolerance>0&&o.geometry.relativeRecognitionTolerance<1e-3&&std::isfinite(o.geometry.relativeRecognitionTolerance)&&
             o.minimumSize>=0&&o.minimumSize<=o.geometry.targetSize&&std::isfinite(o.minimumSize)&&
             o.maximumRelativeVolumeError>0&&o.maximumRelativeVolumeError<1&&std::isfinite(o.maximumRelativeVolumeError)&&
             o.geometry.minimumElementQuality>=0&&o.geometry.minimumElementQuality<=1&&std::isfinite(o.geometry.minimumElementQuality)&&
             o.geometry.warningElementQuality>=0&&o.geometry.warningElementQuality<=1&&std::isfinite(o.geometry.warningElementQuality),"invalid_request");
        need(o.minimumSize==0&&o.geometry.bodyControls.empty()&&o.faceSizes.empty()&&!o.geometry.spatialField&&
             o.geometry.strategy==algo::GeometryStrategy::Auto&&o.geometry.maximumElementSizeRatio==0&&o.geometry.warningElementSizeRatio==0,"unsupported_control");
        need(receipt.limit>0&&o.geometry.maximumNodes>=5&&o.geometry.maximumCells>=8&&o.geometry.maximumOutputBytes>0&&snapshotLimit>0,"zero_admission_budget");
        need((strategy.cuts==OcctContactCutPolicy::AllFeatureAxes||strategy.cuts==OcctContactCutPolicy::PrincipalFeatureAxes)&&
             strategy.radialCuts>=2&&strategy.radialCuts<=16&&strategy.maximumAdaptiveCuts<=16&&
             strategy.maximumPieceWork>=500000&&strategy.maximumPieceWork<=10000000,"invalid_request");
        const unsigned density=strategy.radialCuts;
        const std::string_view mode=strategy.cuts==OcctContactCutPolicy::PrincipalFeatureAxes?"principal":"all";
        double chord=o.maximumChordError>0?o.maximumChordError:.05;
        for(const auto* body:original.solids()){Vec3 lo,hi;body->bbox(lo,hi);
            chord=std::min(chord,std::min(o.geometry.targetSize*.05,(hi-lo).norm()*.0025));}
        need(chord>0&&std::isfinite(chord),"invalid_request");double deflection=chord*.4;
        need(original.solids().size()>1,"assembly_required");
        const auto& source=geo::detail::OcctModelAccess::shape(original);
        need(!sourceBefore||sourceBefore->size()<=snapshotLimit,"snapshot_budget");
        before=sourceBefore?*sourceBefore:snapshot(source,snapshotLimit);
        if(sourceBefore){receipt.sourceAuditAvailable=true;receipt.sourceUnchanged=before==snapshot(source,snapshotLimit);need(receipt.sourceUnchanged,"original_source_changed");}
        auto tick=[&](std::uint64_t n=1){charge(receipt,n);};
        // Every source body remains material, and material overlap is forbidden.
        OcctStageTrace overlapTrace("contact_interference");
        for(std::size_t i=0;i<original.solids().size();++i)for(std::size_t j=0;j<i;++j){
            tick();BRepAlgoAPI_Common common;NCollection_List<TopoDS_Shape> a,b;
            BRepBuilderAPI_Copy ca(geo::detail::OcctModelAccess::solid(original,i+1),true,false),cb(geo::detail::OcctModelAccess::solid(original,j+1),true,false);
            a.Append(ca.Shape());b.Append(cb.Shape());common.SetArguments(a);common.SetTools(b);common.SetNonDestructive(true);common.Build();
            need(common.IsDone()&&!common.HasErrors(),"overlap_query_failed");need(!TopExp_Explorer(common.Shape(),TopAbs_SOLID).More(),"original_material_overlap");
        }
        overlapTrace.finish();
        std::vector<std::pair<Vec3,Vec3>> adaptivePlanes;
        unsigned adaptiveAttempts=0,surfaceRefinements=0,pieceRecoveryRefusals=0;
        std::uint64_t pieceRecoveryWork=0;
        const auto pieceRecoveryWorkLimit=std::min<std::uint64_t>(24'000'000,receipt.limit/4);
        for(;;){OcctStageTrace attemptTrace("contact_partition_attempt");try{
        OcctStageTrace partitionTrace("partition");
        Timer partitionTimer{receipt.partitionSeconds};
        receipt.cadVolume=0;receipt.volumeChange=0;receipt.interfaces=0;receipt.contactTriangles=0;receipt.matchedTriangles=0;
        BRepBuilderAPI_Copy copy(source,true,false);need(copy.IsDone(),"copy_failed");
        ShapeMap sourceMap[3],derivedMap[3],pieces;
        const TopAbs_ShapeEnum kinds[]{TopAbs_VERTEX,TopAbs_EDGE,TopAbs_FACE};
        for(int d=0;d<3;++d)TopExp::MapShapes(source,kinds[d],sourceMap[d]);
        std::vector<std::array<std::set<int>,3>> owned(original.solids().size()+1);
        NCollection_List<TopoDS_Shape> arguments,tools;
        Vec3 low{1e300,1e300,1e300},high{-1e300,-1e300,-1e300};
        for(const auto* body:original.solids()){
            tick();const auto native=geo::detail::OcctModelAccess::solid(original,body->tag());arguments.Append(copy.ModifiedShape(native));
            receipt.cadVolume+=original.volume(body->tag());Vec3 lo,hi;body->bbox(lo,hi);
            low={std::min(low.x,lo.x),std::min(low.y,lo.y),std::min(low.z,lo.z)};high={std::max(high.x,hi.x),std::max(high.y,hi.y),std::max(high.z,hi.z)};
            for(int d=0;d<3;++d){ShapeMap entities;TopExp::MapShapes(native,kinds[d],entities);
                for(int i=1;i<=entities.Extent();++i)owned[body->tag()][d].insert(sourceMap[d].FindIndex(entities(i)));}
        }
        const auto center=low+(high-low)*.5;const double extent=(high-low).norm()*2;
        std::vector<std::pair<Vec3,double>> equations;
        auto addPlane=[&](Vec3 origin,Vec3 normal){tick();normal=normal/normal.norm();const double offset=normal.dot(origin);
            for(const auto& [n,c]:equations){tick();const double dot=n.dot(normal);if(std::abs(dot)>1-1e-10&&std::abs(offset-(dot>0?c:-c))<=1e-8)return;}
            need(equations.size()<64,"plane_budget");auto p=center-normal*(normal.dot(center)-offset);
            BRepBuilderAPI_MakeFace face(gp_Pln(gp_Pnt(p.x,p.y,p.z),gp_Dir(normal.x,normal.y,normal.z)),-extent,extent,-extent,extent);need(face.IsDone(),"plane_failed");
            tools.Append(face.Face());equations.emplace_back(normal,offset);receipt.planes=equations.size();};
        // All body-derived cuts are applied once to the complete copied assembly.
        for(const auto* body:original.solids()){
            Vec3 lo,hi;body->bbox(lo,hi);const auto c=lo+(hi-lo)*.5;
            addPlane(c,{1,0,0});addPlane(c,{0,1,0});addPlane(c,{0,0,1});
            int principalTag=-1;double principalArea=0;Vec3 principalAxis{};
            if(mode=="principal")for(const auto& use:original.faceUses(body->tag())){
                tick();const auto* d=original.surfaceDescriptor(use.faceTag);
                if(!d||(d->kind!=geo::SurfaceKind::Cylinder&&d->kind!=geo::SurfaceKind::Torus))continue;
                GProp_GProps area;BRepGProp::SurfaceProperties(geo::detail::OcctModelAccess::face(original,use.faceTag),area);
                if(area.Mass()>principalArea){principalArea=area.Mass();principalTag=use.faceTag;principalAxis=d->axis;}
            }
            for(const auto& use:original.faceUses(body->tag())){
                tick();const auto* d=original.surfaceDescriptor(use.faceTag);need(d,"missing_descriptor");
                if(d->kind==geo::SurfaceKind::Cylinder||d->kind==geo::SurfaceKind::Cone||d->kind==geo::SurfaceKind::Torus){
                    if(mode=="cartesian")continue;
                    unsigned radial=density;
                    if(mode=="principal"&&use.faceTag!=principalTag){
                        if((d->kind!=geo::SurfaceKind::Cylinder&&d->kind!=geo::SurfaceKind::Cone)||use.sense!=-1)continue;
                        radial=2;
                    }
                    const auto x=d->xDirection,y=d->axis.cross(x);
                    for(unsigned i=0;i<radial;++i){double a=std::numbers::pi*i/radial;addPlane(d->origin,x*std::cos(a)+y*std::sin(a));}
                }else if(d->kind==geo::SurfaceKind::Plane){
                    if(mode=="cartesian"||(mode=="principal"&&(principalTag<0||std::abs(d->axis.dot(principalAxis))<1-1e-10)))continue;
                    double span=std::abs(d->axis.x)*(high.x-low.x)*.5+std::abs(d->axis.y)*(high.y-low.y)*.5+std::abs(d->axis.z)*(high.z-low.z)*.5;
                    if(std::abs(d->axis.dot(d->origin-center))<span-1e-7)addPlane(d->origin,d->axis);
                }
            }
        }
        for(const auto& [origin,normal]:adaptivePlanes)addPlane(origin,normal);
        BRepAlgoAPI_Splitter split;split.SetArguments(arguments);split.SetTools(tools);split.SetNonDestructive(true);
        split.SetRunParallel(false);split.SetToFillHistory(true);
        {OcctStageTrace phase("contact_split");split.Build();}
        need(split.IsDone()&&!split.HasErrors(),"split_failed");
        const auto result=split.Shape();TopExp::MapShapes(result,TopAbs_SOLID,pieces);receipt.pieces=pieces.Extent();
        need(receipt.pieces>0&&receipt.pieces<=o.maximumNativePartitionPieces,"piece_budget");
        for(int d=0;d<3;++d)TopExp::MapShapes(result,kinds[d],derivedMap[d]);
        std::vector<int> parent(pieces.Extent()+1,-1);
        for(const auto* body:original.solids()){
            const auto input=copy.ModifiedShape(geo::detail::OcctModelAccess::solid(original,body->tag()));
            auto assign=[&](const TopoDS_Shape& s){ShapeMap children;TopExp::MapShapes(s,TopAbs_SOLID,children);
                for(int i=1;i<=children.Extent();++i){tick();int tag=pieces.FindIndex(children(i));if(!tag)continue;
                    need(parent[tag]<0||parent[tag]==body->tag(),"ambiguous_body_parent");parent[tag]=body->tag();}};
            assign(input);for(const auto& s:split.Modified(input))assign(s);
        }
        std::array<std::vector<std::vector<Class>>,3> origins;
        for(int d=0;d<3;++d)origins[d].resize(derivedMap[d].Extent()+1);
        auto apply=[&](const TopoDS_Shape& s,Class cls,bool generated){
            for(int d=0;d<3;++d){if(d>cls.dim||(generated&&d==2))continue;ShapeMap entities;TopExp::MapShapes(s,kinds[d],entities);
                for(int i=1;i<=entities.Extent();++i){tick();const int t=derivedMap[d].FindIndex(entities(i));if(t)origins[d][t].push_back(cls);}}
        };
        for(int d=0;d<3;++d)for(int tag=1;tag<=sourceMap[d].Extent();++tag){
            const auto input=copy.ModifiedShape(sourceMap[d](tag));if(input.IsNull())continue;Class cls{d,tag};
            apply(input,cls,false);for(const auto& s:split.Modified(input))apply(s,cls,false);for(const auto& s:split.Generated(input))apply(s,cls,true);
        }
        need(split.DSFiller()!=nullptr,"missing_boolean_history");
        const auto& ds=split.DSFiller()->DS();
        for(int tag=1;tag<=sourceMap[0].Extent();++tag){tick();const auto input=copy.ModifiedShape(sourceMap[0](tag));if(input.IsNull())continue;
            const int index=ds.Index(input);if(index<0)continue;const int representative=ds.GetSameDomainIndex(index);
            need(representative>=0&&representative<ds.NbShapes(),"invalid_boolean_vertex_history");
            const auto& native=ds.Shape(representative);need(native.ShapeType()==TopAbs_VERTEX,"invalid_boolean_vertex_type");
            const int found=derivedMap[0].FindIndex(native);if(found)origins[0][found].push_back({0,tag});}
        auto originalClass=[&](Class c,int body){
            need(c.dim>=0&&c.dim<3&&c.tag>0&&c.tag<static_cast<int>(origins[c.dim].size()),"invalid_derived_class");
            Class selected{3,body};std::set<int> candidates;
            for(const auto p:origins[c.dim][c.tag])if(owned[body][p.dim].contains(p.tag)){
                if(p.dim<selected.dim){selected=p;candidates.clear();candidates.insert(p.tag);}
                else if(p.dim==selected.dim)candidates.insert(p.tag);}
            need(candidates.size()<=1,"ambiguous_original_entity");
            return selected;
        };
        std::vector<double> volumeByBody(original.solids().size()+1);
        std::vector<std::vector<int>> owners(derivedMap[2].Extent()+1);
        OcctStageTrace volumeTrace("contact_partition_volume");
        for(int t=1;t<=pieces.Extent();++t){tick();
            const char* refusal=parent[t]<=0?"unowned_partition_piece":
                !BRepCheck_Analyzer(pieces(t)).IsValid()?"invalid_partition_piece":nullptr;
            if(refusal){
                receipt.failedPiece=t;
                OcctStageTrace rejected("contact_piece_validation",t);
                rejected.finish(refusal);
                throw Stop{refusal};
            }
            double v=mass(pieces(t));need(v>0&&std::isfinite(v),"piece_volume_invalid");volumeByBody[parent[t]]+=v;
            for(TopExp_Explorer f(pieces(t),TopAbs_FACE);f.More();f.Next())owners[derivedMap[2].FindIndex(f.Current())].push_back(t);
        }
        for(const auto* body:original.solids()){
            const double sourceVolume=mass(geo::detail::OcctModelAccess::solid(original,body->tag()));
            double delta=std::abs(volumeByBody[body->tag()]-sourceVolume)/sourceVolume;
            if(delta>1e-6){
                double finePieces=0;for(int tag=1;tag<=pieces.Extent();++tag)if(parent[tag]==body->tag()){tick();finePieces+=accurateMass(pieces(tag));}
                const auto fineSource=accurateMass(geo::detail::OcctModelAccess::solid(original,body->tag()));
                delta=std::abs(finePieces-fineSource)/fineSource;
            }
            receipt.volumeChange=std::max(receipt.volumeChange,delta);need(delta<=1e-6,"partition_volume_mismatch");
        }
        volumeTrace.finish();
        partitionTimer.finish();
        partitionTrace.finish();
        auto derived=geo::detail::OcctModelAccess::fromShape(result);need(derived.imported(),"derived_import_failed");
        geo::OcctBoundaryOptions bo;bo.linearDeflection=deflection;bo.maximumChordError=chord;bo.maximumNodes=o.geometry.maximumNodes;
        bo.maximumTriangles=(o.geometry.maximumCells/2);bo.maximumOutputBytes=o.geometry.maximumOutputBytes;bo.maximumWork=remaining(receipt);bo.cancellation=o.geometry.cancellation;
        geo::OcctBoundaryOutcome prep;std::unique_ptr<geo::detail::OcctBoundarySession> session;
        {OcctStageTrace phase("surface");Timer timer{receipt.surfaceSeconds};session=geo::detail::OcctBoundarySession::prepare(*derived.model,bo,prep,
            geo::detail::OcctBoundaryAuditPolicy::Deferred,geo::detail::OcctBoundaryFacetPolicy::Immediate,
            strategy.controlDeflectionOnAllSurfaces);}
        collect(receipt,prep.work);need(session!=nullptr,geo::occtBoundaryReasonName(prep.reason));
        std::vector<mesh::Mesh> meshes(original.solids().size()+1);
        std::vector<std::map<geo::OcctBoundaryNodeKey,Index>> nodes(original.solids().size()+1);
        std::map<geo::OcctBoundaryNodeKey,Index> canonical;std::map<std::pair<int,Tri>,std::vector<FaceUse>> interfaces;
        const auto working=[&](std::uint64_t active=0){
            std::uint64_t bytes=add(mul(before.size(),sourceBefore?3:2),active);
            for(const auto& m:meshes)bytes=add(bytes,payload(m));
            for(const auto& n:nodes)bytes=add(bytes,mul(n.size(),128));
            bytes=add(bytes,add(mul(canonical.size(),128),mul(interfaces.size(),1024)));
            for(int d=0;d<3;++d)bytes=add(bytes,mul(sourceMap[d].Extent()+derivedMap[d].Extent(),256));
            bytes=add(bytes,mul(parent.size(),sizeof(int)));
            for(const auto& list:owners)bytes=add(bytes,add(sizeof(list),mul(list.size(),sizeof(int))));
            for(const auto& lists:origins)for(const auto& list:lists)bytes=add(bytes,add(sizeof(list),mul(list.size(),sizeof(Class))));
            for(const auto& dimensions:owned)for(const auto& ids:dimensions)bytes=add(bytes,mul(ids.size(),64));
            receipt.maximumWorkingPayloadEstimate=std::max(receipt.maximumWorkingPayloadEstimate,bytes);
            need(bytes<=o.maximumPrimitiveWorkingBytes,"working_payload_budget");
            return bytes;
        };working();
        std::uint64_t totalNodes=0,totalCells=0,totalConnectivity=0;
        for(const auto* piece:derived.model->solids()){
            const int tag=piece->tag(),body=parent[tag];receipt.failedPiece=tag;tick();
            mesh::Mesh boundary;OcctBoundaryOutcome surface;
            {OcctStageTrace phase("surface",body);Timer timer{receipt.surfaceSeconds};surface=session->boundary(tag,boundary,remaining(receipt));}
            collect(receipt,surface.work);
            if(!surface.published&&surface.reason==OcctBoundaryReason::GeometryError&&
               surface.maximumFacetSampleError>chord&&surfaceRefinements<4){
                const auto sourceAudit=session->auditSource(remaining(receipt));collect(receipt,sourceAudit.work);
                need(sourceAudit.reason==OcctBoundaryReason::None,occtBoundaryReasonName(sourceAudit.reason));
                need(sourceAudit.sourceModelUnmodified,"derived_source_changed");throw RetrySurface{};
            }
            need(surface.published,geo::occtBoundaryReasonName(surface.reason));
            geo::detail::OcctStarOptions so;so.maximumWork=std::min<std::uint64_t>(500000,remaining(receipt));
            so.maximumProposalSweeps=256;so.maximumWorkingBytes=o.maximumPrimitiveWorkingBytes;so.cancellation=o.geometry.cancellation;so.maximumBoundaryNodes=o.geometry.maximumNodes;so.maximumBoundaryTriangles=(o.geometry.maximumCells/2);
            so.maximumNodes=o.geometry.maximumNodes;so.maximumCells=o.geometry.maximumCells;so.maximumOutputBytes=o.geometry.maximumOutputBytes;so.auditLimits.maximumTemporaryBytes=o.geometry.maximumAuditTemporaryBytes;so.auditLimits.maximumFaceUses=o.geometry.maximumAuditFaceUses;
            const auto retained=working(payload(boundary));
            working(add(payload(boundary),mul(boundary.nodeCount()+boundary.cellCount(),1024)));
            so.maximumWorkingBytes=o.maximumPrimitiveWorkingBytes-retained;
            if(o.useNativeKernelAdmission){so.nativeAdmissionModel=&original;so.nativeAdmissionBodyTag=body;}
            so.auditLimits.cancellation=o.geometry.cancellation;
            mesh::Mesh fan;OcctStarOutcome star;
            {OcctStageTrace phase("volume_fill",body);Timer timer{receipt.starSeconds};star=geo::detail::buildOcctStarVolume(boundary,tag,so,fan);phase.finish(star.diagnostic);}
            collect(receipt,star.workUsed);
            // Large valid surfaces can exhaust the small per-piece allowance
            // during embedding certification, before any kernel is tested.
            // Retry the same authoritative boundary once; both attempts count
            // against the caller's remaining assembly budget.
            if(!star.published&&star.reason==OcctStarReason::WorkBudget&&
               std::string_view(star.diagnostic)=="work-budget"&&so.maximumWork<strategy.maximumPieceWork){
                so.maximumWork=std::min<std::uint64_t>(strategy.maximumPieceWork,remaining(receipt));
                {OcctStageTrace phase("volume_fill",body);Timer timer{receipt.starSeconds};star=geo::detail::buildOcctStarVolume(boundary,tag,so,fan);phase.finish(star.diagnostic);}
                collect(receipt,star.workUsed);
            }
            // A local search-work cap may be recovered by another cut. Other
            // storage/cancellation failures retain their typed child outcome.
            if(star.reason!=OcctStarReason::WorkBudget||std::string_view(star.diagnostic)!="work-budget")checkStarResource(star.reason);
            bool recoveredPiece=false;
            if(!star.published&&adaptiveAttempts<strategy.maximumAdaptiveCuts){
                bool certified=star.embedded;
                const bool hasGenus=star.reason==OcctStarReason::UnsupportedTopology&&
                    std::string_view(star.diagnostic)=="boundary-genus-not-zero";
                if(hasGenus){
                    so.maximumWork=std::min<std::uint64_t>(500000,remaining(receipt));
                    const auto certificate=geo::detail::certifyOcctBoundarySurface(boundary,so);collect(receipt,certificate.workUsed);certified=certificate.certified;
                    checkStarResource(certificate.reason);
                    if(certificate.reason==OcctStarReason::Cancelled)throw Stop{"cancelled"};
                    if(certificate.reason==OcctStarReason::AllocationFailed)throw Stop{"allocation_failed"};
                }
                const bool noKernel=star.reason==OcctStarReason::KernelNotFound;
                // A star kernel is a property of this fill method. An
                // embedded boundary can still admit a general Tet4 fill.
                // Try that immutable interface before cutting every body again.
                const bool requestRecovery=certified&&(noKernel||hasGenus)&&o.useAdaptiveTetRecovery&&
                    o.recoveryMode==NativeCadRecoveryMode::Advanced;
                const auto recoveryWorkRemaining=pieceRecoveryWorkLimit-pieceRecoveryWork;
                if(requestRecovery&&pieceRecoveryRefusals<4&&recoveryWorkRemaining&&remaining(receipt)){
                    auto limits=so;
                    limits.maximumWork=std::min<std::uint64_t>({2'000'000,remaining(receipt),recoveryWorkRemaining});
                    OcctPartitionRecoveryOptions recovery;
                    recovery.maximumPieceWork=limits.maximumWork;
                    recovery.maximumFrontSeconds=.75;
                    OcctPlcRecoveryOutcome candidate;
                    {
                        OcctStageTrace phase("contact_piece_recovery",tag);
                        Timer timer{receipt.starSeconds};
                        candidate=recoverOcctPlcFrontVolume(boundary,tag,*derived.model,original,body,limits,recovery,fan);
                        phase.finish(candidate.diagnostic);
                    }
                    collect(receipt,candidate.workUsed);
                    need(candidate.workUsed<=recoveryWorkRemaining,"piece_recovery_work_budget");
                    pieceRecoveryWork+=candidate.workUsed;
                    if(candidate.reason==OcctPlcRecoveryReason::Cancelled)throw Stop{"cancelled"};
                    if(candidate.reason==OcctPlcRecoveryReason::AllocationFailed)throw Stop{"allocation_failed"};
                    if(candidate.reason==OcctPlcRecoveryReason::InternalError)throw Stop{"internal-error"};
                    if(candidate.reason==OcctPlcRecoveryReason::InvalidRequest)throw Stop{"invalid_request"};
                    if(candidate.published){
                        need(fan.nodeCount()>=boundary.nodeCount(),"recovered_boundary_prefix_missing");
                        for(Index i=0;i<boundary.nodeCount();++i){
                            tick();
                            need(same(fan.node(i),boundary.node(i))&&fan.nodeClass(i)==boundary.nodeClass(i),
                                 "recovered_boundary_prefix_changed");
                        }
                        recoveredPiece=true;
                    }else ++pieceRecoveryRefusals;
                }else if(requestRecovery){
                    OcctStageTrace limit("contact_piece_recovery_limit",tag);
                    limit.finish(pieceRecoveryRefusals>=4?"piece_recovery_refusal_limit":"piece_recovery_work_limit");
                }
                if(!recoveredPiece&&certified&&(noKernel||star.reason==OcctStarReason::WorkBudget||hasGenus)){
                    const auto sourceAudit=session->auditSource(remaining(receipt));collect(receipt,sourceAudit.work);need(sourceAudit.reason==OcctBoundaryReason::None,occtBoundaryReasonName(sourceAudit.reason));need(sourceAudit.sourceModelUnmodified,"derived_source_changed");
                    Vec3 lo,hi;piece->bbox(lo,hi);throw RetryCut{lo,hi};
                }
            }
            if(!recoveredPiece){checkStarResource(star.reason);need(star.published,star.diagnostic);}
            OcctStageTrace assemblingTrace("partition_assembly",body);
            Timer assembling{receipt.assemblySeconds};
            auto& all=meshes[body];std::vector<Index> remap(fan.nodeCount()),global(boundary.nodeCount());
            for(Index i=0;i<boundary.nodeCount();++i){tick();const auto& key=surface.nodeKeys[i];const auto cls=originalClass(boundary.nodeClass(i),body);const auto p=boundary.node(i);
                const auto found=nodes[body].find(key);
                if(found==nodes[body].end()){
                    prospectiveOutput(totalNodes+1,totalCells,totalConnectivity,receipt);working(add(add(payload(boundary),payload(fan)),1024));
                    remap[i]=all.addNode(p,cls);nodes[body].emplace(key,remap[i]);++totalNodes;
                }else{need(same(all.node(found->second),p)&&all.nodeClass(found->second)==cls,"shared_node_mismatch");remap[i]=found->second;}
                auto [gi,newGlobal]=canonical.emplace(key,canonical.size());global[i]=gi->second;
            }
            for(Index i=boundary.nodeCount();i<fan.nodeCount();++i){tick();prospectiveOutput(totalNodes+1,totalCells,totalConnectivity,receipt);
                working(add(add(payload(boundary),payload(fan)),1024));remap[i]=all.addNode(fan.node(i),{3,body});++totalNodes;}
            for(Index i=0;i<fan.cellCount();++i){tick();const auto& cell=fan.cell(i);const auto* ids=fan.cellNodes(i);
                if(cell.type==mesh::CellType::Tet4){std::array<Index,4> ns{remap[ids[0]],remap[ids[1]],remap[ids[2]],remap[ids[3]]};
                    prospectiveOutput(totalNodes,totalCells+1,totalConnectivity+4,receipt);working(add(add(payload(boundary),payload(fan)),1024));
                    all.addCell(cell.type,ns,{3,body});++totalCells;totalConnectivity+=4;}
                else{need(cell.type==mesh::CellType::Tri3,"unexpected_cell_type");const int face=cell.onEntity.tag;Tri tri{global[ids[0]],global[ids[1]],global[ids[2]]};
                    const auto cls=originalClass({2,face},body);const auto& own=owners[face];
                    need(own.size()==1||own.size()==2,"face_ownership_invalid");
                    if(own.size()==1||parent[own[0]]!=parent[own[1]]){need(cls.dim==2,"boundary_face_parent_lost");
                        Tri ns{remap[ids[0]],remap[ids[1]],remap[ids[2]]};prospectiveOutput(totalNodes,totalCells+1,totalConnectivity+3,receipt);
                        working(add(add(payload(boundary),payload(fan)),1024));all.addCell(mesh::CellType::Tri3,ns,cls);++totalCells;totalConnectivity+=3;}
                    if(own.size()==2){working(add(add(payload(boundary),payload(fan)),2048));interfaces[{face,ordered(tri)}].push_back({tag,body,face,tri,{fan.node(ids[0]),fan.node(ids[1]),fan.node(ids[2])}});}
                }
                output(all,receipt);working(add(payload(boundary),payload(fan)));
            }
        }
        std::set<std::pair<int,int>> contactBodyPairs;
        for(const auto& [key,uses]:interfaces){tick();need(uses.size()==2&&uses[0].piece!=uses[1].piece&&!sameWinding(uses[0].canonical,uses[1].canonical),"interface_triangles_mismatch");
            ++receipt.interfaces;if(uses[0].body!=uses[1].body){++receipt.contactTriangles;++receipt.matchedTriangles;
                contactBodyPairs.emplace(std::min(uses[0].body,uses[1].body),std::max(uses[0].body,uses[1].body));}}
        need(receipt.contactTriangles>0,"no_face_contacts");
        mesh::Mesh joined;out.bodies.clear();
        for(const auto* body:original.solids()){
            OcctStageTrace validatingTrace("cad_validation",body->tag());
            Timer validating{receipt.validationSeconds};
            auto& mesh=meshes[body->tag()];need(mesh.nodeCount()>0&&mesh.cellCount()>0,"missing_body");
            mesh::Mesh exterior;std::map<Index,Index> exteriorIds;
            paacs::geo::detail::OcctSolidClassifier classifier(geo::detail::OcctModelAccess::solid(original,body->tag()));
            auto inside=[&](Vec3 p){tick();classifier.Perform(gp_Pnt(p.x,p.y,p.z),original.tolerance(3,body->tag()));return classifier.State()==TopAbs_IN||classifier.State()==TopAbs_ON;};
            std::set<int> retainedVertices,retainedFaces;
            for(Index i=0;i<mesh.nodeCount();++i){tick();const auto cls=mesh.nodeClass(i);const auto p=mesh.node(i);Vec3 q;
                if(cls.dim==0)retainedVertices.insert(cls.tag);
                if(cls.dim==3){need(inside(p),"interior_node_outside_original");continue;}
                if(cls.dim==0)q=original.vertex(cls.tag)->position();else if(cls.dim==1)q=original.edge(cls.tag)->closestPoint(p);
                else need(geo::closestPointOnOcctFace(original,cls.tag,p,q),"original_node_projection_failed");
                const double error=(p-q).norm();receipt.maxNodeError=std::max(receipt.maxNodeError,error);
                need(error<=std::max(8*original.tolerance(cls.dim,cls.tag),1e-9),"original_node_geometry_error");
            }
            for(Index i=0;i<mesh.cellCount();++i){tick();const auto& cell=mesh.cell(i);const auto* ids=mesh.cellNodes(i);
                if(cell.type==mesh::CellType::Tet4){Vec3 p{};for(int j=0;j<4;++j)p=p+mesh.node(ids[j]);need(inside(p/4.),"tet_centroid_outside_original");}
                else{const auto a=mesh.node(ids[0]),b=mesh.node(ids[1]),c=mesh.node(ids[2]);
                    retainedFaces.insert(cell.onEntity.tag);
                    for(const auto p:{(a+b+c)/3.,(a+b)/2.,(b+c)/2.,(c+a)/2.}){tick();Vec3 q;need(geo::closestPointOnOcctFace(original,cell.onEntity.tag,p,q),"original_facet_projection_failed");
                        const double error=(p-q).norm();receipt.maxSampleError=std::max(receipt.maxSampleError,error);need(error<=chord+original.tolerance(2,cell.onEntity.tag),"original_facet_error");}
                    working(add(add(payload(joined),payload(exterior)),add(mul(exteriorIds.size(),128),1024)));
                    Tri ns;for(int j=0;j<3;++j){auto [it,inserted]=exteriorIds.emplace(ids[j],exterior.nodeCount());if(inserted)exterior.addNode(mesh.node(ids[j]),mesh.nodeClass(ids[j]));ns[j]=it->second;}exterior.addCell(mesh::CellType::Tri3,ns,cell.onEntity);
                }
            }
            need(retainedVertices==owned[body->tag()][0]&&retainedFaces==owned[body->tag()][2],"original_boundary_entities_missing");
            geo::detail::OcctStarOptions audit; audit.maximumWork=remaining(receipt);audit.maximumBoundaryNodes=o.geometry.maximumNodes;audit.maximumBoundaryTriangles=(o.geometry.maximumCells/2);audit.maximumWorkingBytes=o.maximumPrimitiveWorkingBytes;audit.cancellation=o.geometry.cancellation;
            const auto retainedAudit=working(add(add(payload(joined),payload(exterior)),mul(exteriorIds.size(),128)));
            audit.maximumWorkingBytes=o.maximumPrimitiveWorkingBytes-retainedAudit;
            const auto certificate=geo::detail::certifyOcctBoundarySurface(exterior,audit);collect(receipt,certificate.workUsed);checkStarResource(certificate.reason);need(certificate.certified,certificate.diagnostic);
            mesh::LinearVolumeAuditLimits limits;limits.maximumTemporaryBytes=o.geometry.maximumAuditTemporaryBytes;limits.maximumFaceUses=o.geometry.maximumAuditFaceUses;limits.cancellation=o.geometry.cancellation;
            const auto valid=mesh::auditLinearVolumeBoundary(mesh,limits);need(!valid.cancelled,"cancelled");need(!valid.budgetExceeded,"audit_budget");need(!valid.allocationFailed,"allocation_failed");need(valid.valid,"body_boundary_audit_failed");
            const auto quality=algo::analyzeMixedMeshQuality(mesh);need(!quality.allocationFailed,"allocation_failed");need(quality.valid,"body_quality_invalid");
            need(std::abs(quality.totalVolume-original.volume(body->tag()))<=original.volume(body->tag())*o.maximumRelativeVolumeError,"body_mesh_volume_error");
            NativeCadBodyMeshReport report;report.solidTag=body->tag();report.route="native_contact_assembly";
            report.reason=NativeCadMeshReason::None;report.cadVolume=original.volume(body->tag());report.meshVolume=quality.totalVolume;
            report.relativeVolumeError=std::abs(report.meshVolume-report.cadVolume)/report.cadVolume;
            report.maximumBoundaryNodeError=receipt.maxNodeError;report.maximumBoundarySampleError=receipt.maxSampleError;
            report.surfaceDeflection=deflection;report.surfaceChordErrorLimit=chord;report.boundaryEmbeddingCertified=certificate.certified;
            report.surfaceAttempts=1+surfaceRefinements;
            report.nativeSurfaceReason="none";report.nativeVolumeReason="none";report.partitionReason="none";
            report.partitionPolicy="native_shared_assembly";report.partitionAttempts=1;report.partitionPieces=receipt.pieces;
            report.partitionInterfaces=receipt.interfaces;report.partitionVolumeChange=receipt.volumeChange;
            report.nativePartitionSeconds=receipt.partitionSeconds;report.nativeSurfaceSeconds=receipt.surfaceSeconds;
            report.nativeVolumeSeconds=receipt.starSeconds;report.nativeAssemblySeconds=receipt.assemblySeconds;
            for(Index i=0;i<mesh.cellCount();++i)if(mesh.cell(i).type==mesh::CellType::Tri3)++report.surfaceTriangles;
            out.bodies.push_back(std::move(report));
            prospectiveOutput(add(joined.nodeCount(),mesh.nodeCount()),add(joined.cellCount(),mesh.cellCount()),add(joined.connectivityCount(),mesh.connectivityCount()),receipt);
            working(add(add(payload(joined),payload(mesh)),add(payload(exterior),mul(exteriorIds.size(),128))));
            const Index offset=joined.nodeCount();for(Index i=0;i<mesh.nodeCount();++i){tick();joined.addNode(mesh.node(i),mesh.nodeClass(i));}
            for(Index i=0;i<mesh.cellCount();++i){tick();const auto& c=mesh.cell(i);std::array<Index,4> ns;const auto n=mesh::nodesPerCell(c.type);for(int j=0;j<n;++j)ns[j]=mesh.cellNodes(i)[j]+offset;joined.addCell(c.type,std::span<const Index>(ns.data(),n),c.onEntity);}
        }
        OcctStageTrace finalTrace("contact_validation");
        output(joined,receipt);working(payload(joined));
        std::uint64_t treePayload=0;
        for(const auto* body:original.solids()){std::uint64_t triangles=0;for(Index i=0;i<meshes[body->tag()].cellCount();++i){tick();if(meshes[body->tag()].cell(i).type==mesh::CellType::Tri3)++triangles;}
            treePayload=add(treePayload,mul(triangles,sizeof(BoundaryTriangle)+2*sizeof(TreeNode)));}
        working(add(payload(joined),treePayload));
        std::vector<Tree> boundaryTrees;for(const auto* body:original.solids())boundaryTrees.push_back(tree(meshes[body->tag()],receipt));
        for(std::size_t i=0;i<boundaryTrees.size();++i)for(std::size_t j=0;j<i;++j)auditContacts(boundaryTrees[i],0,boundaryTrees[j],0,receipt);
        const auto sourceAudit=session->auditSource(remaining(receipt));collect(receipt,sourceAudit.work);need(sourceAudit.reason==OcctBoundaryReason::None,occtBoundaryReasonName(sourceAudit.reason));need(sourceAudit.sourceModelUnmodified,"derived_source_changed");
        receipt.sourceAuditAvailable=true;receipt.sourceUnchanged=before==snapshot(source,snapshotLimit);need(receipt.sourceUnchanged,"original_source_changed");
        const auto quality=algo::analyzeMixedMeshQuality(joined);need(!quality.allocationFailed,"allocation_failed");need(quality.valid,"final_quality_invalid");receipt.meshVolume=quality.totalVolume;receipt.minQuality=quality.minimumQuality;
        receipt.nodes=joined.nodeCount();receipt.cells=joined.cellCount();need(receipt.nodes<=o.geometry.maximumNodes&&receipt.cells<=o.geometry.maximumCells,"final_count_budget");
        output(joined,receipt);
        out.audit=mesh::auditLinearVolumeBoundary(joined,finalLimits);
        need(!out.audit.cancelled,"cancelled");need(!out.audit.budgetExceeded,"audit_budget");need(!out.audit.allocationFailed,"allocation_failed");need(out.audit.valid,"final_boundary_audit_failed");
        out.quality=quality;need(quality.minimumQuality>=o.geometry.minimumElementQuality,"quality_requirement");
        out.qualityWarning=quality.minimumQuality<o.geometry.warningElementQuality;
        tick();need(before==snapshot(source,snapshotLimit),"original_source_changed");tick(0);
        out.separateContactPairs=contactBodyPairs.size();
        destination.swap(joined);out.published=true;out.destinationUnchanged=false;out.reason=NativeCadMeshReason::None;
        receipt.published=true;receipt.failedPiece=-1;receipt.reason="none";
        finalTrace.finish();
        attemptTrace.finish("none");
        break;
        }catch(const RetryCut& cut){
            attemptTrace.finish("adaptive_cut_retry");
            need(adaptiveAttempts<strategy.maximumAdaptiveCuts,"adaptive_cut_budget");++adaptiveAttempts;
            const auto center=cut.low+(cut.high-cut.low)*.5;
            adaptivePlanes.emplace_back(center,Vec3{1,0,0});
            adaptivePlanes.emplace_back(center,Vec3{0,1,0});
            adaptivePlanes.emplace_back(center,Vec3{0,0,1});
        }catch(const RetrySurface&){
            attemptTrace.finish("surface_refinement_retry");
            need(surfaceRefinements<4,"surface_refinement_budget");++surfaceRefinements;deflection*=.5;
        }}

    }catch(const Stop& s){receipt.reason=s.diagnostic;out.reason=failureReason(s.diagnostic);}
    catch(const Standard_OutOfMemory&){receipt.reason="allocation_failed";out.reason=NativeCadMeshReason::AllocationFailed;}
    catch(const std::bad_alloc&){receipt.reason="allocation_failed";out.reason=NativeCadMeshReason::AllocationFailed;}
    catch(...){receipt.reason="native_internal_error";out.reason=NativeCadMeshReason::InternalError;}
    if(!out.published&&!before.empty()){
        try{receipt.sourceUnchanged=before==snapshot(OcctModelAccess::shape(original),snapshotLimit);receipt.sourceAuditAvailable=true;}
        catch(...){receipt.sourceUnchanged=false;receipt.sourceAuditAvailable=false;}
        if(receipt.sourceAuditAvailable&&!receipt.sourceUnchanged){out.reason=NativeCadMeshReason::InternalError;receipt.reason="original_source_changed";}
    }
    if(!out.published)for(auto& report:out.bodies)report.reason=out.reason;
    out.nativeContactDiagnostic=receipt.reason;out.nativeContactSourceGuard=receipt.sourceUnchanged;
    out.nativeContactSourceAuditAvailable=receipt.sourceAuditAvailable;
    out.nativeContactMaximumWorkingPayloadEstimate=receipt.maximumWorkingPayloadEstimate;
    out.nativeContactPieces=receipt.pieces;out.nativeContactPlanes=receipt.planes;
    out.matchedNativeContactTriangles=receipt.matchedTriangles;
    out.validationWork=receipt.work;out.totalSeconds=std::chrono::duration<double>(Clock::now()-begin).count();
    out.preparationSeconds=receipt.partitionSeconds;out.generationSeconds=receipt.surfaceSeconds+receipt.starSeconds+receipt.assemblySeconds;
    out.validationSeconds=receipt.validationSeconds;
    return out;
}
} // namespace paacs::geo::detail
