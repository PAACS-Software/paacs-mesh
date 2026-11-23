// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/OcctMesher.hpp"
#include "OcctModelInternal.hpp"
#include "paacs/geometry/OcctBoundary.hpp"
#include "OcctStarMesher.hpp"
#include "OcctStageTrace.hpp"
#include "OcctPartitionMesher.hpp"
#include "OcctRecoveryMesher.hpp"
#include "OcctFrontMesher.hpp"
#include "OcctBoundaryInternal.hpp"
#include "OcctContactCylinder.hpp"
#include "OcctBoxContact.hpp"
#include "OcctQueryContext.hpp"
#include "OcctContactAssembly.hpp"
#include "OcctContactComponents.hpp"
#include "OcctNativeSweep.hpp"
#include "OcctPairedShell.hpp"
#include "paacs/bnd/TetFacetGeometry.hpp"
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <GProp_GProps.hxx>
#include <Standard_OutOfMemory.hxx>
#include <TopExp_Explorer.hxx>
#include <NCollection_List.hxx>
#include <TopoDS.hxx>
#include <gp_Ax2.hxx>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <map>
#include <string_view>
#include <limits>
#include <set>
#include <type_traits>

namespace paacs::geo {
namespace {
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point t) { return std::chrono::duration<double>(Clock::now()-t).count(); }
struct PhaseTimer {
    double& seconds;Clock::time_point started=Clock::now();bool active=true;
    explicit PhaseTimer(double& value):seconds(value){}
    void finish() noexcept {if(active){seconds+=elapsed(started);active=false;}}
    ~PhaseTimer(){finish();}
};
struct Stop { NativeCadMeshReason reason; };
void need(bool ok, NativeCadMeshReason reason) { if(!ok) throw Stop{reason}; }
bool finite(double x) { return std::isfinite(x); }
bool finite(Vec3 p) { return finite(p.x)&&finite(p.y)&&finite(p.z); }
double mass(const TopoDS_Shape& s) { GProp_GProps g; BRepGProp::VolumeProperties(s,g); return std::abs(g.Mass()); }
Vec3 vec(const gp_Pnt& p) { return {p.X(),p.Y(),p.Z()}; }
gp_Ax2 axes(const algo::PrimitiveShape& s) {
    return {gp_Pnt(s.origin.x,s.origin.y,s.origin.z), gp_Dir(s.axis.x,s.axis.y,s.axis.z),
            gp_Dir(s.xDirection.x,s.xDirection.y,s.xDirection.z)};
}
struct View final : Model {
    const OcctModel& source; const GSolid* body;
    std::vector<const GFace*> fs; std::vector<const GEdge*> es; std::vector<const GVertex*> vs;
    View(const OcctModel& m,const GSolid* s):source(m),body(s) {
        std::set<int> edges, vertices;
        for(int tag:s->boundaryFaceTags()) {
            fs.push_back(m.face(tag));
            for(const auto& loop:m.boundaryUses(tag)) for(const auto& use:loop.edges) {
                edges.insert(use.edgeTag); vertices.insert(use.startVertexTag); vertices.insert(use.endVertexTag);
            }
        }
        for(int tag:edges) es.push_back(m.edge(tag));
        for(int tag:vertices) if(tag>=0) vs.push_back(m.vertex(tag));
    }
    std::vector<const GVertex*> vertices() const override{return vs;}
    std::vector<const GEdge*> edges() const override{return es;}
    std::vector<const GFace*> faces() const override{return fs;}
    std::vector<const GSolid*> solids() const override{return {body};}
    const GVertex* vertex(int t) const override{return source.vertex(t);}
    const GEdge* edge(int t) const override{return source.edge(t);}
    const GFace* face(int t) const override{return source.face(t);}
    const GSolid* solid(int t) const override{return source.solid(t);}
};
struct Recognition {
    algo::PrimitiveShape shape;
    std::map<int,int> faces; // Primitive face identity -> imported face identity.
    bool analytic = false;
};
Recognition recognize(const OcctModel& m,const View& v,double relativeTolerance) {
    Recognition r;
    std::vector<const GFace*> curved, planes;
    for(auto f:v.fs) (f->kind()==SurfaceKind::Plane?planes:curved).push_back(f);
    if(curved.empty()) return r;
    need(curved.size()<=2,NativeCadMeshReason::UnsupportedBody);
    const auto* d=m.surfaceDescriptor(curved.front()->tag());
    need(d&&finite(d->origin)&&finite(d->axis)&&finite(d->xDirection),NativeCadMeshReason::UnsupportedBody);
    r.shape.origin=d->origin; r.shape.axis=d->axis; r.shape.xDirection=d->xDirection;
    r.shape.radius=d->radius; r.faces[1]=curved.front()->tag();
    if(d->kind==SurfaceKind::Torus) {
        need(curved.size()==1&&(planes.empty()||planes.size()==2),NativeCadMeshReason::UnsupportedBody);
        r.shape.kind=algo::PrimitiveKind::Torus; r.shape.majorRadius=d->majorRadius;
        if(!planes.empty()) {
            const double begin=d->parameterLower.u, angle=d->parameterUpper.u-begin;
            need(finite(begin)&&finite(angle)&&angle>relativeTolerance&&angle<2*std::numbers::pi-relativeTolerance,
                 NativeCadMeshReason::UnsupportedBody);
            const Vec3 y=d->axis.cross(d->xDirection);
            r.shape.xDirection=d->xDirection*std::cos(begin)+y*std::sin(begin);
            r.shape.kind=algo::PrimitiveKind::Elbow;r.shape.sweepAngle=angle;
            const Vec3 startTangent=d->axis.cross(r.shape.xDirection);
            const Vec3 endRadial=d->xDirection*std::cos(begin+angle)+y*std::sin(begin+angle);
            const Vec3 endTangent=d->axis.cross(endRadial);
            // Match both oriented cap planes by their complete major-circle
            // endpoint locations. The subsequent BRep common proves the trims.
            for(auto f:planes) {
                const auto* cap=m.surfaceDescriptor(f->tag());need(cap!=nullptr,NativeCadMeshReason::UnsupportedBody);
                const Vec3 startCenter=d->origin+r.shape.xDirection*d->majorRadius;
                const Vec3 endCenter=d->origin+endRadial*d->majorRadius;
                const double tolerance=std::max(m.tolerance(2,f->tag()),d->majorRadius*relativeTolerance);
                const bool start=std::abs(cap->axis.dot(startTangent))>1-relativeTolerance&&
                    std::abs((startCenter-cap->origin).dot(cap->axis))<=tolerance;
                const bool end=std::abs(cap->axis.dot(endTangent))>1-relativeTolerance&&
                    std::abs((endCenter-cap->origin).dot(cap->axis))<=tolerance;
                // At 180 degrees both planes coincide, so distinguish their
                // trimmed cap centers rather than making a plane-only guess.
                const Vec3 center=vec(BRepAdaptor_Surface(detail::OcctModelAccess::face(m,f->tag())).Value(
                    (cap->parameterLower.u+cap->parameterUpper.u)*.5,
                    (cap->parameterLower.v+cap->parameterUpper.v)*.5));
                const int role=start&&end?((center-startCenter).norm()<(center-endCenter).norm()?2:3):start?2:end?3:0;
                need(role!=0&&r.faces.emplace(role,f->tag()).second,NativeCadMeshReason::UnsupportedBody);
            }
        }
    } else {
        need(d->kind==SurfaceKind::Cylinder||d->kind==SurfaceKind::Cone,NativeCadMeshReason::UnsupportedBody);
        need(planes.size()==2||(planes.size()==1&&d->kind==SurfaceKind::Cone),NativeCadMeshReason::UnsupportedBody);
        double low=1e300, high=-1e300;
        for(auto vert:v.vs) { const double z=(vert->position()-d->origin).dot(d->axis); low=std::min(low,z); high=std::max(high,z); }
        need(finite(low)&&finite(high)&&high>low,NativeCadMeshReason::RecognitionFailed);
        r.shape.height=high-low; r.shape.origin=d->origin+d->axis*low;
        for(auto f:planes) {
            const auto* cap=m.surfaceDescriptor(f->tag());
            need(cap&&std::abs(cap->axis.dot(d->axis))>1-relativeTolerance,NativeCadMeshReason::UnsupportedBody);
            const double z=(cap->origin-d->origin).dot(d->axis);
            const int role=std::abs(z-low)<std::abs(z-high)?2:3;
            need(r.faces.emplace(role,f->tag()).second,NativeCadMeshReason::UnsupportedBody);
        }
        if(d->kind==SurfaceKind::Cylinder) {
            r.shape.kind=algo::PrimitiveKind::Cylinder;
            if(curved.size()==2) {
                const auto* inner=m.surfaceDescriptor(curved.back()->tag());
                need(inner&&inner->kind==SurfaceKind::Cylinder&&
                     std::abs(inner->axis.dot(d->axis))>1-relativeTolerance&&
                     (inner->origin-d->origin).cross(d->axis).norm()<std::max(m.tolerance(2,curved.front()->tag()),d->radius*relativeTolerance),
                     NativeCadMeshReason::UnsupportedBody);
                r.shape.kind=algo::PrimitiveKind::Tube;
                if(inner->radius>r.shape.radius) { r.shape.innerRadius=r.shape.radius; r.shape.radius=inner->radius;
                    r.faces[4]=r.faces[1]; r.faces[1]=curved.back()->tag(); }
                else { r.shape.innerRadius=inner->radius; r.faces[4]=curved.back()->tag(); }
            }
        } else {
            need(curved.size()==1,NativeCadMeshReason::UnsupportedBody);
            double a=d->radius+low*std::tan(d->semiAngle), b=d->radius+high*std::tan(d->semiAngle);
            if(a<b) { std::swap(a,b); r.shape.origin=d->origin+d->axis*high; r.shape.axis=d->axis*-1;
                std::swap(r.faces[2],r.faces[3]); }
            r.shape.radius=a; r.shape.endRadius=b;
            const double tol=std::max(m.tolerance(2,curved.front()->tag()),a*relativeTolerance);
            r.shape.kind=std::abs(b)<=tol?algo::PrimitiveKind::Cone:algo::PrimitiveKind::Frustum;
            if(r.shape.kind==algo::PrimitiveKind::Cone) r.faces.erase(3);
        }
    }
    TopoDS_Shape expected;
    if(r.shape.kind==algo::PrimitiveKind::Torus||r.shape.kind==algo::PrimitiveKind::Elbow)
        expected=r.shape.kind==algo::PrimitiveKind::Torus?
            BRepPrimAPI_MakeTorus(axes(r.shape),r.shape.majorRadius,r.shape.radius).Shape():
            BRepPrimAPI_MakeTorus(axes(r.shape),r.shape.majorRadius,r.shape.radius,r.shape.sweepAngle).Shape();
    else if(r.shape.kind==algo::PrimitiveKind::Cone||r.shape.kind==algo::PrimitiveKind::Frustum)
        expected=BRepPrimAPI_MakeCone(axes(r.shape),r.shape.radius,
            r.shape.kind==algo::PrimitiveKind::Cone?0:r.shape.endRadius,r.shape.height).Shape();
    else {
        expected=BRepPrimAPI_MakeCylinder(axes(r.shape),r.shape.radius,r.shape.height).Shape();
        if(r.shape.kind==algo::PrimitiveKind::Tube) {
            BRepAlgoAPI_Cut cut(expected,BRepPrimAPI_MakeCylinder(axes(r.shape),r.shape.innerRadius,r.shape.height).Shape());
            need(cut.IsDone(),NativeCadMeshReason::RecognitionFailed); expected=cut.Shape();
        }
    }
    const double volume=m.volume(v.body->tag()), expectedVolume=mass(expected);
    const double volumeTolerance=std::min(1e-8,100*relativeTolerance);
    need(volume>0&&std::abs(expectedVolume-volume)<=volume*volumeTolerance,NativeCadMeshReason::UnsupportedBody);
    // OCCT's convenience constructor builds in destructive mode. Recognition
    // must not add pcurves or change tolerances on the imported snapshot.
    BRepAlgoAPI_Common common;
    NCollection_List<TopoDS_Shape> arguments, tools;
    arguments.Append(detail::OcctModelAccess::solid(m,v.body->tag()));
    tools.Append(expected);
    common.SetArguments(arguments); common.SetTools(tools);
    common.SetNonDestructive(true); common.Build();
    need(common.IsDone()&&std::abs(mass(common.Shape())-volume)<=volume*volumeTolerance,NativeCadMeshReason::RecognitionFailed);
    r.analytic=true; return r;
}
struct Work {
    const NativeCadMeshOptions& o; std::uint64_t count=0;
    void tick(std::uint64_t n=1) { need(n<=o.maximumValidationWork-count,NativeCadMeshReason::BudgetExceeded);count+=n;
        if(o.geometry.cancellation&&o.geometry.cancellation->load(std::memory_order_relaxed)) throw Stop{NativeCadMeshReason::Cancelled}; }
};
struct BodyRange { int tag; Index begin, end; };
struct ContactTriangle { std::array<Vec3,3> p; Vec3 lo,hi; };
struct ContactBox { Vec3 lo,hi; std::size_t begin,end,left=0,right=0; };
bool boxesOverlap(Vec3 a,Vec3 b,Vec3 c,Vec3 d) {
    return a.x<=d.x&&c.x<=b.x&&a.y<=d.y&&c.y<=b.y&&a.z<=d.z&&c.z<=b.z;
}
double coord(Vec3 p,int axis){return axis==0?p.x:axis==1?p.y:p.z;}
struct ContactTree {
    std::vector<ContactTriangle> triangles; std::vector<ContactBox> boxes;
    std::size_t build(std::size_t begin,std::size_t end,Work& work) {
        work.tick(end-begin);ContactBox box{{1e300,1e300,1e300},{-1e300,-1e300,-1e300},begin,end};
        for(std::size_t i=begin;i<end;++i){const auto& t=triangles[i];
            box.lo={std::min(box.lo.x,t.lo.x),std::min(box.lo.y,t.lo.y),std::min(box.lo.z,t.lo.z)};
            box.hi={std::max(box.hi.x,t.hi.x),std::max(box.hi.y,t.hi.y),std::max(box.hi.z,t.hi.z)};}
        const auto id=boxes.size();boxes.push_back(box);
        if(end-begin>8){const auto span=box.hi-box.lo;int axis=span.y>span.x?1:0;if(span.z>coord(span,axis))axis=2;
            const auto mid=begin+(end-begin)/2;
            std::nth_element(triangles.begin()+begin,triangles.begin()+mid,triangles.begin()+end,
                [axis](const auto& a,const auto& b){return coord(a.lo+a.hi,axis)<coord(b.lo+b.hi,axis);});
            const auto left=build(begin,mid,work),right=build(mid,end,work);boxes[id].left=left;boxes[id].right=right;}
        return id;
    }
};
ContactTree contactTree(const mesh::Mesh& m,BodyRange range,Work& work) {
    std::uint64_t count=0;for(Index i=range.begin;i<range.end;++i){work.tick();const auto& c=m.cell(i);
        if(c.onEntity.dim==2){need(c.type==mesh::CellType::Tri3,NativeCadMeshReason::UnsupportedControl);++count;}}
    need(count>0&&count<=work.o.maximumPrimitiveWorkingBytes/(2*(sizeof(ContactTriangle)+2*sizeof(ContactBox))),NativeCadMeshReason::BudgetExceeded);
    ContactTree tree;tree.triangles.reserve(count);tree.boxes.reserve(2*count);
    for(Index i=range.begin;i<range.end;++i)if(m.cell(i).onEntity.dim==2){
        const auto* ids=m.cellNodes(i);ContactTriangle t{{m.node(ids[0]),m.node(ids[1]),m.node(ids[2])}};
        t.lo={std::min({t.p[0].x,t.p[1].x,t.p[2].x}),std::min({t.p[0].y,t.p[1].y,t.p[2].y}),std::min({t.p[0].z,t.p[1].z,t.p[2].z})};
        t.hi={std::max({t.p[0].x,t.p[1].x,t.p[2].x}),std::max({t.p[0].y,t.p[1].y,t.p[2].y}),std::max({t.p[0].z,t.p[1].z,t.p[2].z})};tree.triangles.push_back(t);}
    tree.build(0,count,work);return tree;
}
void auditSeparateContact(const ContactTree& a,std::size_t ai,const ContactTree& b,std::size_t bi,
                         Work& work,std::uint64_t& pairs,
                         std::vector<bool>* matchedA=nullptr,std::vector<bool>* matchedB=nullptr,
                         double* matchedArea=nullptr,std::uint64_t* matchedCount=nullptr) {
    work.tick();const auto& x=a.boxes[ai];const auto& y=b.boxes[bi];
    if(!boxesOverlap(x.lo,x.hi,y.lo,y.hi))return;
    if(x.left){auditSeparateContact(a,x.left,b,bi,work,pairs,matchedA,matchedB,matchedArea,matchedCount);auditSeparateContact(a,x.right,b,bi,work,pairs,matchedA,matchedB,matchedArea,matchedCount);return;}
    if(y.left){auditSeparateContact(a,ai,b,y.left,work,pairs,matchedA,matchedB,matchedArea,matchedCount);auditSeparateContact(a,ai,b,y.right,work,pairs,matchedA,matchedB,matchedArea,matchedCount);return;}
    for(std::size_t i=x.begin;i<x.end;++i)for(std::size_t j=y.begin;j<y.end;++j){work.tick();const auto& p=a.triangles[i];const auto& q=b.triangles[j];
        if(!boxesOverlap(p.lo,p.hi,q.lo,q.hi))continue;++pairs;
        if(matchedA){
            std::array<int,3> permutation{-1,-1,-1};
            for(int u=0;u<3;++u)for(int v=0;v<3;++v)
                if(p.p[u].x==q.p[v].x&&p.p[u].y==q.p[v].y&&p.p[u].z==q.p[v].z)permutation[u]=v;
            if(permutation[0]>=0&&permutation[1]>=0&&permutation[2]>=0){
                const int inversions=(permutation[0]>permutation[1])+(permutation[0]>permutation[2])+(permutation[1]>permutation[2]);
                need((inversions%2)==1&&!(*matchedA)[i]&&!(*matchedB)[j],NativeCadMeshReason::UnresolvedBodyInterface);
                (*matchedA)[i]=true;(*matchedB)[j]=true;
                const double area=(p.p[1]-p.p[0]).cross(p.p[2]-p.p[0]).norm()*.5;
                need(area>0&&finite(area),NativeCadMeshReason::GeometryError);
                *matchedArea+=area;++*matchedCount;continue;
            }
        }
        need(!bnd::tet_facet_detail::trianglesOpenIntersect(p.p,q.p)&&
             !bnd::tet_facet_detail::coplanarTrianglesHaveAreaOverlap(p.p,q.p),NativeCadMeshReason::UnresolvedBodyInterface);
    }
}
detail::OcctBoxContactOptions boxContactOptions(const NativeCadMeshOptions& o,const Work& work,
    detail::OcctBoxCertificateSession* certificates) {
    need(o.maximumPrimitiveWorkingBytes>=8,NativeCadMeshReason::BudgetExceeded);
    detail::OcctBoxContactOptions options;
    options.maximumWork=o.maximumValidationWork-work.count;
    options.maximumWorkingBytes=o.maximumPrimitiveWorkingBytes;
    options.maximumSnapshotBytes=std::min<std::uint64_t>(64ULL*1024*1024,o.maximumPrimitiveWorkingBytes/8);
    options.cancellation=o.geometry.cancellation;
    options.certificateSession=certificates;
    return options;
}
void checkBoxContactFailure(detail::OcctBoxContactReason reason) {
    if(reason==detail::OcctBoxContactReason::BudgetExceeded)throw Stop{NativeCadMeshReason::BudgetExceeded};
    if(reason==detail::OcctBoxContactReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled};
    if(reason==detail::OcctBoxContactReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed};
    if(reason==detail::OcctBoxContactReason::SourceChanged||reason==detail::OcctBoxContactReason::InternalError)throw Stop{NativeCadMeshReason::InternalError};
}
std::uint64_t payload(std::uint64_t nodes,std::uint64_t cells,std::uint64_t connectivity) {
    std::uint64_t total=0;
    for(const auto [count,width]:{std::pair{nodes,std::uint64_t(sizeof(Vec3)+sizeof(mesh::Classification))},
                                 std::pair{cells,std::uint64_t(sizeof(mesh::Cell))},std::pair{connectivity,std::uint64_t(sizeof(Index))}}) {
        need(count<=(UINT64_MAX-total)/width,NativeCadMeshReason::BudgetExceeded);total+=count*width;
    }
    return total;
}
int faceRole(const Recognition& r,int tag) {
    auto f=r.faces.find(tag); need(f!=r.faces.end()&&f->second>0,NativeCadMeshReason::RecognitionFailed); return f->second;
}
mesh::Classification classify(const OcctModel& m,const View& v,const Recognition& r,
                               Vec3 p,mesh::Classification c, Work& work,double& error,bool nativeSurface,
                               bool tryFaceWitness) {
    if(c.dim==3) return {3,v.body->tag()};
    if(!r.analytic||nativeSurface) {
        work.tick(); Vec3 q;
        const auto tolerance=std::max(8*m.tolerance(c.dim,c.tag),1e-9);
        if(c.dim==0) {const auto* vertex=m.vertex(c.tag);need(vertex!=nullptr,NativeCadMeshReason::GeometryError);q=vertex->position();}
        else if(c.dim==1) {const auto* edge=m.edge(c.tag);need(edge!=nullptr,NativeCadMeshReason::GeometryError);q=edge->closestPoint(p);}
        else {
            Vec2 uv;
            const auto* face=m.face(c.tag);
            need(c.dim==2&&face,NativeCadMeshReason::GeometryError);
            const bool witness=tryFaceWitness&&face->projectPoint(p,uv,q)&&(p-q).norm()<=tolerance;
            need(witness||(nativeSurface?closestPointOnOcctFace(m,c.tag,p,q):face->projectPoint(p,uv,q)),
                 NativeCadMeshReason::GeometryError);
        }
        const double e=(p-q).norm();error=std::max(error,e);
        need(finite(e)&&e<=tolerance,NativeCadMeshReason::GeometryError);return c;
    }
    // Resolve supported primitive boundary features against this body's actual
    // CAD incidence. This changes classifications only; no nodes are merged.
    if(c.dim<=2) {
        for(auto vert:v.vs) { work.tick();
            const double e=(vert->position()-p).norm();
            if(e<=std::max(m.tolerance(0,vert->tag()),1e-10)) {error=std::max(error,e);return {0,vert->tag()};}
        }
        for(auto edge:v.es) { work.tick();
            const double e=(edge->closestPoint(p)-p).norm();
            if(e<=std::max(m.tolerance(1,edge->tag()),1e-10)) {error=std::max(error,e);return {1,edge->tag()};}
        }
        if(c.dim<=1)need(false,NativeCadMeshReason::GeometryError);
    }
    const int tag=faceRole(r,c.tag); const auto* f=m.face(tag);
    Vec2 uv; Vec3 q; work.tick(); need(f->projectPoint(p,uv,q),NativeCadMeshReason::GeometryError);
    const double e=(p-q).norm(); error=std::max(error,e);
    need(e<=std::max(m.tolerance(2,tag)*8,1e-9),NativeCadMeshReason::GeometryError);
    return {2,tag};
}
void append(const OcctModel& m,const View& view,const Recognition& r,const mesh::Mesh& body,
            mesh::Mesh& all,NativeCadBodyMeshReport& report,Work& work,bool nativeSurface=false,
            bool originalFacetCentersCertified=false) {
    const auto& o=work.o.geometry;
    const auto n=static_cast<std::uint64_t>(all.nodeCount())+body.nodeCount();
    const auto c=static_cast<std::uint64_t>(all.cellCount())+body.cellCount();
    const auto cn=static_cast<std::uint64_t>(all.connectivityCount())+body.connectivityCount();
    need(n<=o.maximumNodes&&c<=o.maximumCells&&payload(n,c,cn)<=o.maximumOutputBytes,
         NativeCadMeshReason::BudgetExceeded);
    all.reserveAppend(body.nodeCount(),body.cellCount(),body.connectivityCount());
    const Index offset=all.nodeCount();
    const std::string_view route(report.route);
    const bool tryFaceWitness=route=="native_frontal_surface"||route=="native_frontal_quality"||
        route=="native_curvature_surface"||route=="native_spatial_sizing";
    const auto chord=std::max(report.surfaceChordErrorLimit*1.01,1e-7);
    detail::OcctStageTrace phase("assembly_boundary_node_checks",view.body->tag());
    for(Index i=0;i<body.nodeCount();++i) all.addNode(body.node(i),
        classify(m,view,r,body.node(i),body.nodeClass(i),work,report.maximumBoundaryNodeError,nativeSurface,tryFaceWitness));
    phase.next("assembly_cell_copy_and_facet_checks");
    std::array<Index,10> ids{};
    for(Index i=0;i<body.cellCount();++i) {
        work.tick(); const auto& cell=body.cell(i); const auto* nodes=body.cellNodes(i);
        auto cls=cell.onEntity;
        if(r.analytic&&!nativeSurface) cls=cls.dim==3?mesh::Classification{3,view.body->tag()}:mesh::Classification{2,faceRole(r,cls.tag)};
        const int width=mesh::nodesPerCell(cell.type);
        for(int j=0;j<width;++j) ids[j]=nodes[j]+offset;
        all.addCell(cell.type,std::span<const Index>(ids.data(),width),cls);
        if(cls.dim==2&&!originalFacetCentersCertified) {
            Vec3 center{}; for(int j=0;j<width;++j)center=center+body.node(nodes[j]); center=center/static_cast<double>(width);
            Vec2 uv;Vec3 projected;work.tick();
            // A point on the trimmed face is a direct upper bound for the
            // distance. Preserve the full shape-distance query when it misses.
            const bool witness=tryFaceWitness&&m.face(cls.tag)->projectPoint(center,uv,projected)&&
                (center-projected).norm()<=chord;
            need(witness||(nativeSurface?closestPointOnOcctFace(m,cls.tag,center,projected):
                m.face(cls.tag)->projectPoint(center,uv,projected)),NativeCadMeshReason::GeometryError);
            report.maximumBoundarySampleError=std::max(report.maximumBoundarySampleError,(center-projected).norm());
        }
    }
}
bool nativeStarFallback(const OcctModel& model,const View& view,const NativeCadMeshOptions& options,
                        mesh::Mesh& body,NativeCadBodyMeshReport& report,Work& work) {
    const auto& o=options.geometry;Vec3 lo,hi;view.body->bbox(lo,hi);const double diagonal=(hi-lo).norm();
    need(o.maximumMeshingWork>0&&o.maximumNodes>=5&&o.maximumCells>=8&&payload(5,8,28)<=o.maximumOutputBytes,NativeCadMeshReason::BudgetExceeded);
    const double chord=options.maximumChordError>0?options.maximumChordError:
        std::min(o.targetSize*.05,diagonal*.0025);
    need(chord>0&&finite(chord),NativeCadMeshReason::InvalidRequest);
    report.route="native_surface_star";report.surfaceChordErrorLimit=chord;
    double deflection=chord*.4;
    for(unsigned attempt=0;attempt<3;++attempt) {
        work.tick();++report.surfaceAttempts;report.surfaceDeflection=deflection;
        OcctBoundaryOptions bo;bo.linearDeflection=deflection;bo.maximumChordError=chord;
        bo.maximumNodes=o.maximumNodes;bo.maximumTriangles=o.maximumCells/2;
        bo.maximumOutputBytes=o.maximumOutputBytes;
        bo.maximumWork=std::min(o.maximumMeshingWork,options.maximumValidationWork-work.count);
        bo.cancellation=o.cancellation;mesh::Mesh boundary;
        const auto surfaceStarted=Clock::now();
        detail::OcctStageTrace phase("surface",view.body->tag());
        OcctBoundaryOutcome surface;
        std::unique_ptr<detail::OcctBoundarySession> candidateSession;
        if(options.deferDiscardedNativeFacetChecks) {
            OcctBoundaryOutcome preparation;
            candidateSession=detail::OcctBoundarySession::prepare(model,bo,preparation,
                detail::OcctBoundaryAuditPolicy::PerBoundary,detail::OcctBoundaryFacetPolicy::DeferredCandidate);
            if(candidateSession){surface=candidateSession->boundary(view.body->tag(),boundary);surface.work+=preparation.work;}
            else surface=std::move(preparation);
        } else surface=buildOcctBoundary(model,view.body->tag(),bo,boundary);
        phase.finish(occtBoundaryReasonName(surface.reason));
        report.nativeSurfaceSeconds+=elapsed(surfaceStarted);
        report.nativeSurfaceReason=occtBoundaryReasonName(surface.reason);report.nativeSurfaceWork+=surface.work;
        report.surfaceTriangles=surface.triangles;work.tick(surface.work);
        report.maximumBoundarySampleError=surface.maximumFacetSampleError;
        report.maximumBoundaryNodeError=surface.maximumNodeProjectionError;
        if(!surface.published) {
            if(surface.reason==OcctBoundaryReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled};
            if(surface.reason==OcctBoundaryReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed};
            if(surface.reason==OcctBoundaryReason::BudgetExceeded)throw Stop{NativeCadMeshReason::BudgetExceeded};
            if(surface.reason==OcctBoundaryReason::SourceChanged)throw Stop{NativeCadMeshReason::InternalError};
            if(surface.reason==OcctBoundaryReason::GeometryError&&attempt<2){deflection*=.5;continue;}
            return false;
        }
        if(surface.relativeVolumeError>options.maximumRelativeVolumeError){
            if(attempt<2){deflection*=.5;continue;}
            report.nativeVolumeReason="native-boundary-volume-error";return false;
        }
        detail::OcctStarOptions so;
        const auto remaining=std::min(o.maximumMeshingWork,options.maximumValidationWork-work.count);
        need(remaining>0,NativeCadMeshReason::BudgetExceeded);
        // A unsuccessful kernel search must leave work for another strategy.
        so.maximumWork=options.useNativePartitionFallback?
            std::min(remaining,std::max<std::uint64_t>(20000,o.maximumMeshingWork/5)):remaining;
        so.maximumBoundaryNodes=o.maximumNodes;so.maximumBoundaryTriangles=bo.maximumTriangles;
        so.maximumProposalSweeps=16;so.maximumNodes=o.maximumNodes;so.maximumCells=o.maximumCells;
        so.maximumOutputBytes=o.maximumOutputBytes;so.maximumWorkingBytes=options.maximumPrimitiveWorkingBytes;
        so.auditLimits.maximumTemporaryBytes=o.maximumAuditTemporaryBytes;so.auditLimits.maximumFaceUses=o.maximumAuditFaceUses;
        so.cancellation=o.cancellation;
        const auto volumeStarted=Clock::now();
        phase.next("volume_fill");
        const auto volume=detail::buildOcctStarVolume(boundary,view.body->tag(),so,body);
        phase.finish(volume.diagnostic);
        report.nativeVolumeSeconds+=elapsed(volumeStarted);
        report.nativeVolumeReason=volume.diagnostic;report.nativeVolumeWork+=volume.workUsed;
        report.exactSurfacePairs=volume.exactBoundaryPairs;report.kernelTrials=volume.kernelTrials;work.tick(volume.workUsed);
        if(volume.reason==detail::OcctStarReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled};
        if(volume.reason==detail::OcctStarReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed};
        if(volume.reason==detail::OcctStarReason::WorkBudget&&options.useNativePartitionFallback&&
           std::string_view(volume.diagnostic)=="work-budget"&&work.count<options.maximumValidationWork)return false;
        if(volume.reason==detail::OcctStarReason::WorkBudget||volume.reason==detail::OcctStarReason::OutputBudget)throw Stop{NativeCadMeshReason::BudgetExceeded};
        if(!volume.published)return false;
        phase.next("cad_validation");
        work.tick();bool inside=view.body->isInside(volume.kernel);
        for(Index i=0;inside&&i<body.cellCount();++i){work.tick();const auto& c=body.cell(i);
            if(c.type!=mesh::CellType::Tet4)continue;
            const auto ns=body.cellNodes(i);Vec3 center{};
            for(unsigned j=0;j<4;++j)center=center+body.node(ns[j]);
            inside=view.body->isInside(center/4.0);
        }
        if(!inside){report.nativeVolumeReason="native_tet_centroid_outside";
            if(attempt<2){deflection*=.5;continue;}return false;}
        if(options.deferDiscardedNativeFacetChecks) {
            bo.maximumWork=std::min(o.maximumMeshingWork,options.maximumValidationWork-work.count);
            const auto snapshotLimit=std::min<std::uint64_t>(64ULL*1024*1024,options.maximumPrimitiveWorkingBytes/8);
            need(bo.maximumWork&&snapshotLimit,NativeCadMeshReason::BudgetExceeded);
            const auto auditStarted=Clock::now();
            const auto audit=detail::auditOcctBoundaryFacets(model,view.body->tag(),boundary,bo,
                snapshotLimit);
            report.nativeSurfaceSeconds+=elapsed(auditStarted);
            report.nativeSurfaceWork+=audit.work;work.tick(audit.work);
            report.nativeSurfaceReason=occtBoundaryReasonName(audit.reason);
            report.maximumBoundarySampleError=audit.maximumFacetSampleError;
            if(!audit.completed) {
                if(audit.reason==OcctBoundaryReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled};
                if(audit.reason==OcctBoundaryReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed};
                if(audit.reason==OcctBoundaryReason::BudgetExceeded)throw Stop{NativeCadMeshReason::BudgetExceeded};
                if(audit.reason==OcctBoundaryReason::SourceChanged||audit.reason==OcctBoundaryReason::InternalError)
                    throw Stop{NativeCadMeshReason::InternalError};
                if(audit.reason==OcctBoundaryReason::GeometryError&&attempt<2){deflection*=.5;continue;}
                return false;
            }
            need(audit.sourceModelUnmodified&&audit.samples==4*boundary.cellCount(),NativeCadMeshReason::InternalError);
            // Compare with the original snapshot taken before copying, across
            // all volume feasibility/classification and deferred facet calls.
            const auto finalSourceStarted=Clock::now();
            const auto sourceAudit=candidateSession->auditSource(options.maximumValidationWork-work.count);
            report.nativeSurfaceSeconds+=elapsed(finalSourceStarted);
            report.nativeSurfaceWork+=sourceAudit.work;work.tick(sourceAudit.work);
            if(sourceAudit.reason==OcctBoundaryReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled};
            if(sourceAudit.reason==OcctBoundaryReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed};
            if(sourceAudit.reason==OcctBoundaryReason::BudgetExceeded)throw Stop{NativeCadMeshReason::BudgetExceeded};
            need(sourceAudit.reason==OcctBoundaryReason::None&&sourceAudit.sourceModelUnmodified,NativeCadMeshReason::InternalError);
        }
        report.route="native_surface_star";report.meshVolume=volume.boundaryVolume;report.boundaryEmbeddingCertified=volume.embedded;return true;
    }
    return false;
}
bool nativeContactCylinder(const OcctModel& model,const View& view,const mesh::Mesh& accepted,
                           const NativeCadMeshOptions& options,mesh::Mesh& body,
                           NativeCadBodyMeshReport& report,Work& work) {
    using namespace detail;const auto& o=options.geometry;
    Vec3 lo,hi;view.body->bbox(lo,hi);
    const double chord=options.maximumChordError>0?options.maximumChordError:
        std::min(o.targetSize*.05,(hi-lo).norm()*.0025);
    OcctContactCylinderOptions co;co.maximumChordError=chord;
    co.maximumRelativeVolumeError=options.maximumRelativeVolumeError;
    co.maximumWork=std::min(o.maximumMeshingWork,options.maximumValidationWork-work.count);
    need(co.maximumWork>0,NativeCadMeshReason::BudgetExceeded);
    co.maximumNodes=o.maximumNodes;co.maximumCells=o.maximumCells;
    co.maximumOutputBytes=o.maximumOutputBytes;co.maximumWorkingBytes=options.maximumPrimitiveWorkingBytes;
    co.maximumSnapshotBytes=std::min<std::uint64_t>(64ULL*1024*1024,options.maximumPrimitiveWorkingBytes/8);
    co.cancellation=o.cancellation;
    co.star.maximumWork=co.maximumWork;co.star.maximumProposalSweeps=16;
    co.star.maximumBoundaryNodes=o.maximumNodes;co.star.maximumBoundaryTriangles=o.maximumCells/2;
    co.star.maximumNodes=o.maximumNodes;co.star.maximumCells=o.maximumCells;
    co.star.maximumOutputBytes=o.maximumOutputBytes;co.star.maximumWorkingBytes=options.maximumPrimitiveWorkingBytes;
    co.star.auditLimits.maximumTemporaryBytes=o.maximumAuditTemporaryBytes;
    co.star.auditLimits.maximumFaceUses=o.maximumAuditFaceUses;co.star.auditLimits.cancellation=o.cancellation;
    const auto result=buildOcctContactCylinder(model,view.body->tag(),accepted,co,body);
    work.tick(result.workUsed);report.contactCylinderReason=occtContactCylinderReasonName(result.reason);
    report.contactSourceEdgeTag=result.sourceEdgeTag;report.contactTargetEdgeTag=result.targetEdgeTag;
    report.contactRingNodes=result.ringNodes;report.contactCylinderWork=result.workUsed;
    if(result.reason==OcctContactCylinderReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled};
    if(result.reason==OcctContactCylinderReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed};
    if(result.reason==OcctContactCylinderReason::BudgetExceeded)throw Stop{NativeCadMeshReason::BudgetExceeded};
    if(result.reason==OcctContactCylinderReason::SourceChanged||result.reason==OcctContactCylinderReason::InternalError)
        throw Stop{NativeCadMeshReason::InternalError};
    if(!result.published)return false;
    report.route="native_contact_cylinder";report.meshVolume=result.meshVolume;
    report.maximumBoundaryNodeError=result.maximumBoundaryNodeError;
    report.maximumBoundarySampleError=result.maximumBoundarySampleError;
    report.surfaceChordErrorLimit=chord;report.boundaryEmbeddingCertified=true;
    for(Index i=0;i<body.cellCount();++i){work.tick();if(body.cell(i).type==mesh::CellType::Tri3)++report.surfaceTriangles;}
    return true;
}
bool hasSplineRichBoundary(const OcctModel& model,const View& view) {
    unsigned splineFaces=0;
    for(const auto* face:view.fs){
        const auto kind=BRepAdaptor_Surface(detail::OcctModelAccess::face(model,face->tag()),false).GetType();
        if((kind==GeomAbs_BSplineSurface||kind==GeomAbs_BezierSurface)&&++splineFaces>=4)return true;
    }
    return false;
}
enum class NativePartitionPass {Full,Midplanes,ChordTrial};
bool nativePartitionFallback(const OcctModel& model,const View& view,const NativeCadMeshOptions& options,
                             mesh::Mesh& body,NativeCadBodyMeshReport& report,Work& work,
                             NativePartitionPass pass=NativePartitionPass::Full,bool preferFixedGrid=false) {
    using namespace detail;const auto& o=options.geometry;
    need(options.maximumNativePartitionPieces>0&&o.maximumMeshingWork>0,NativeCadMeshReason::BudgetExceeded);
    enum class Priority {Normal,ChordGrid,CoarseGrid,FixedGrid,FixedGridRepair,Concave};
    struct Attempt {OcctPartitionPolicy policy;unsigned density;const char* name;
        unsigned kernelSweeps=16;std::uint64_t pieceWork=200000;unsigned refinements=3;
        double deflectionFraction=1;bool containmentOnly=false;bool recoverPlc=false;
        bool featureRecovery=false,liveRecovery=false;Priority priority=Priority::Normal;};
    const std::array attempts{
        Attempt{OcctPartitionPolicy::CartesianMidplanes,2,"cartesian_midplanes"},
        Attempt{OcctPartitionPolicy::CartesianGrid,5,"cartesian_fifths_chord_trial",
                16,200000,1,1,false,false,false,false,Priority::ChordGrid},
        Attempt{OcctPartitionPolicy::CartesianGrid,3,"cartesian_thirds_local_recovery",
                16,1'000'000,1,1,false,false,false,true},
        Attempt{OcctPartitionPolicy::CartesianGrid,3,"cartesian_thirds_priority",
                16,1'000'000,1,1,false,false,false,true,Priority::CoarseGrid},
        Attempt{OcctPartitionPolicy::CartesianGrid,5,"cartesian_fifths_priority",
                16,200000,1,1,false,false,false,false,Priority::FixedGrid},
        Attempt{OcctPartitionPolicy::CartesianGrid,5,"cartesian_fifths_local_recovery",
                16,1'000'000,1,1,false,false,false,true,Priority::FixedGridRepair},
        Attempt{OcctPartitionPolicy::CoaxialAngular,8,"coaxial_angular"},
        Attempt{OcctPartitionPolicy::NativeFeatureAxes,4,"native_feature_axes_4"},
        Attempt{OcctPartitionPolicy::NativeFeatureAxes,4,"native_feature_axes_4_extended",256,500000},
        Attempt{OcctPartitionPolicy::NativeCircularFeatureAxes,4,"native_circular_feature_axes_4",256,500000,3,.007,true},
        Attempt{OcctPartitionPolicy::NativeCircularFeatureAxes,2,"native_circular_feature_axes_2",256,500000,3,.007,true},
        Attempt{OcctPartitionPolicy::NativeFeatureAxes,8,"native_feature_axes_8"},
        Attempt{OcctPartitionPolicy::NativeConcaveFeatureAxes,8,"native_concave_feature_axes_8_priority",
                256,500000,3,.2,false,false,false,false,Priority::Concave},
        Attempt{OcctPartitionPolicy::NativeCircularFeatureAxes,4,"native_circular_feature_axes_4",256,500000,3,.007},
        Attempt{OcctPartitionPolicy::NativeCircularFeatureAxes,2,"native_circular_feature_axes_2",256,500000,3,.007},
        Attempt{OcctPartitionPolicy::NativeFeatureAxes,2,"native_feature_axes_2"},
        Attempt{OcctPartitionPolicy::CartesianGrid,3,"cartesian_thirds"},
        Attempt{OcctPartitionPolicy::CartesianGrid,5,"cartesian_fifths"},
        Attempt{OcctPartitionPolicy::NativeConcaveFeatureAxes,8,"native_concave_feature_axes_8",256,500000,3,.2},
        Attempt{OcctPartitionPolicy::NativeConcaveFeatureAxes,16,"native_concave_feature_axes_16",256,500000,3,.2},
        Attempt{OcctPartitionPolicy::NativePrincipalFeatureAxes,8,"native_principal_feature_axes_8_fine",256,5000000,3,1./800},
        Attempt{OcctPartitionPolicy::NativeConcaveFeatureAxes,16,"native_concave_feature_axes_16_plc_recovery",256,1000000,3,.2,false,true},
        // Keep established successful policies ahead of these construction
        // attempts. Adaptive surface recovery can be expensive even when a
        // later established policy already produces the accepted mesh.
        Attempt{OcctPartitionPolicy::NativePrincipalFeatureAxes,4,"native_principal_feature_axes_4_surface_recovery",256,1000000,1,.4,false,false,true},
        Attempt{OcctPartitionPolicy::NativeConcaveFeatureAxes,8,"native_concave_feature_axes_8_surface_recovery",256,1000000,1,.4,false,false,true}};
    bool admissionRefused=false,extendedKernelRequired=false,circularContainmentRequired=false;
    unsigned circularPoliciesTried=0;double featureDeflection=0;
    std::vector<OcctPartitionPlane> retainedRecoveryPlanes;
    double retainedRecoveryDeflection=0;
    const auto reported=[&](std::string_view policy,std::string_view diagnostic){
        return std::ranges::any_of(report.partitionAttemptReports,[&](const auto& prior){
            return prior.policy==policy&&prior.diagnostic==diagnostic;
        });
    };
    bool triedPriorityConcave=false;
    const bool initialCartesianOnly=pass==NativePartitionPass::Midplanes;
    const bool chordTrial=pass==NativePartitionPass::ChordTrial;
    for(const auto& attempt:attempts){const auto policy=attempt.policy;
      if(initialCartesianOnly&&policy!=OcctPartitionPolicy::CartesianMidplanes)break;
      if(chordTrial!=(attempt.priority==Priority::ChordGrid))continue;
      if((attempt.priority==Priority::CoarseGrid||attempt.priority==Priority::FixedGrid)&&!preferFixedGrid)continue;
      if(attempt.priority==Priority::FixedGridRepair&&
         (!preferFixedGrid||!(reported("cartesian_fifths_priority","native-tet-centroid-outside")||
            reported("cartesian_fifths_priority","no-strict-star-kernel-found")||
            reported("cartesian_fifths_priority","boundary-genus-not-zero"))))continue;
      if(attempt.priority==Priority::Concave){
          if(!options.useAdaptiveTetRecovery||o.minimumElementQuality>0||
             !reported("cartesian_midplanes","work-budget")||
             !reported("native_feature_axes_4","partition-budget")||
             !reported("native_feature_axes_8","partition-budget"))continue;
          triedPriorityConcave=true;
      }else if(triedPriorityConcave&&policy==OcctPartitionPolicy::NativeConcaveFeatureAxes&&
               attempt.density==8&&!attempt.featureRecovery)continue;
      if(attempt.liveRecovery&&attempt.priority==Priority::Normal&&
         (!options.useLocalPartitionRecovery||o.minimumElementQuality>0||
         report.partitionAttemptReports.empty()||
         std::string_view(report.partitionAttemptReports.back().diagnostic)!="boundary-genus-not-zero"))continue;
      if(!initialCartesianOnly&&policy==OcctPartitionPolicy::CartesianMidplanes&&
         std::any_of(report.partitionAttemptReports.begin(),report.partitionAttemptReports.end(),
                     [](const auto& r){return std::string_view(r.policy)=="cartesian_midplanes";}))continue;
      // Fine native principal supports are an appended explicit-admission
      // route. Every previously accepted policy keeps its order and result.
      if(policy==OcctPartitionPolicy::NativePrincipalFeatureAxes&&!options.useNativeKernelAdmission)continue;
      if(attempt.featureRecovery&&!options.useNativeKernelAdmission)continue;
      if(attempt.recoverPlc&&!options.useNativeKernelAdmission)continue;
      const bool extended=policy==OcctPartitionPolicy::NativeFeatureAxes&&attempt.kernelSweeps>16;
      if(extended&&!extendedKernelRequired)continue;
      if(policy==OcctPartitionPolicy::NativeCircularFeatureAxes){
          if(attempt.containmentOnly&&!circularContainmentRequired)continue;
          const unsigned bit=attempt.density==4?1:2;
          if(circularPoliciesTried&bit)continue;
          circularPoliciesTried|=bit;
      }
      // Coarse facets can cut across curved material even when the native
      // partition is valid. Retry only failed geometric/surface gates, keeping
      // every earlier attempt charged to the same cumulative work admission.
      double deflection=std::min(report.surfaceDeflection,report.surfaceChordErrorLimit*attempt.deflectionFraction);
      if(chordTrial)deflection=report.surfaceChordErrorLimit;
      if(attempt.featureRecovery)deflection=report.surfaceChordErrorLimit*attempt.deflectionFraction;
      if(extended&&featureDeflection>0)deflection=featureDeflection;
      if(attempt.recoverPlc&&!retainedRecoveryPlanes.empty())deflection=retainedRecoveryDeflection;
      for(unsigned refinement=0;refinement<attempt.refinements;++refinement){
        work.tick();++report.partitionAttempts;
        OcctPartitionMeshOptions po;po.partition.policy=policy;
        po.maximumNativeValidationThreads=options.maximumNativeValidationThreads;
        po.partition.cartesianDivisions=attempt.density;po.partition.featureAxisPlanes=attempt.density;
        po.partition.maximumPieces=options.maximumNativePartitionPieces;
        po.partition.maximumHistoryWork=std::min(o.maximumMeshingWork,options.maximumValidationWork-work.count);
        po.partition.maximumSnapshotAuditBytes=options.maximumPrimitiveWorkingBytes;
        po.partition.retrySpanAwareVolumeReference=attempt.priority==Priority::CoarseGrid;
        if(policy==OcctPartitionPolicy::NativeConcaveFeatureAxes){
            po.useNativeKernelAdmission=options.useNativeKernelAdmission;
            po.maximumAdaptiveCuts=options.useNativeKernelAdmission&&attempt.density==16?16:4;
            // Initial feature cuts and bounded recovery cuts share this cap.
            // The native admission route can need more than 48 total planes.
            if(options.useNativeKernelAdmission)po.partition.maximumPlanes=64;
            if(attempt.recoverPlc&&!retainedRecoveryPlanes.empty()){
                po.partition.extraPlanes=retainedRecoveryPlanes;
                po.maximumAdaptiveCuts=0;
            }
        }
        if(policy==OcctPartitionPolicy::NativePrincipalFeatureAxes){
            po.useNativeKernelAdmission=true;po.maximumAdaptiveCuts=0;
            po.partition.maximumPlanes=64;
        }
        if(attempt.featureRecovery){
            po.useNativeKernelAdmission=true;po.repairSurfaceAliases=true;
            po.partition.meshConstructionCandidate=true;
            po.partition.includeCartesianMidplanes=false;
            // A disposable native split is a construction aid. Its small
            // integration/trim drift does not change the final original-CAD
            // node, facet, material, embedding or mesh-volume acceptance.
            po.partition.maximumRelativeVolumeChange=std::min(1e-4,options.maximumRelativeVolumeError*.01);
            po.maximumAdaptiveCuts=policy==OcctPartitionPolicy::NativeConcaveFeatureAxes?16:0;
        }
        po.boundary.linearDeflection=deflection;
        po.boundary.maximumChordError=report.surfaceChordErrorLimit;
        po.boundary.maximumNodes=o.maximumNodes;po.boundary.maximumTriangles=o.maximumCells/2;
        po.boundary.maximumOutputBytes=o.maximumOutputBytes;
        po.maximumWork=std::min(o.maximumMeshingWork,options.maximumValidationWork-work.count);
        if(attempt.liveRecovery){
            po.maximumWork=std::min<std::uint64_t>(po.maximumWork,90'000'000);
            po.useNativeKernelAdmission=true;
        }
        need(po.maximumWork>0,NativeCadMeshReason::BudgetExceeded);
        po.boundary.maximumWork=po.maximumWork;po.star.maximumWork=po.maximumWork;
        const auto pieceWork=options.useNativeKernelAdmission&&policy==OcctPartitionPolicy::NativeConcaveFeatureAxes
            ?std::max<std::uint64_t>(attempt.pieceWork,1000000):attempt.pieceWork;
        po.maximumPieceStarWork=std::min(pieceWork,po.maximumWork);
        po.star.maximumBoundaryNodes=o.maximumNodes;po.star.maximumBoundaryTriangles=o.maximumCells/2;
        po.star.maximumProposalSweeps=attempt.kernelSweeps;
        po.star.auditLimits.maximumTemporaryBytes=o.maximumAuditTemporaryBytes;
        po.star.auditLimits.maximumFaceUses=o.maximumAuditFaceUses;
        po.maximumNodes=o.maximumNodes;po.maximumCells=o.maximumCells;
        po.maximumOutputBytes=o.maximumOutputBytes;po.maximumWorkingBytes=options.maximumPrimitiveWorkingBytes;
        po.maximumRelativeVolumeError=options.maximumRelativeVolumeError;po.cancellation=o.cancellation;
        if(options.improveNativeTetQuality)po.tetQualityTargetSize=o.targetSize;
        OcctPartitionRecoveryOptions recovery;
        recovery.enabled=attempt.recoverPlc||attempt.liveRecovery;
        recovery.useLiveConstraints=attempt.liveRecovery;
        recovery.useAdaptiveRefill=attempt.liveRecovery&&options.useAdaptiveTetRecovery;
        if(attempt.liveRecovery)recovery.maximumPieceWork=12'000'000;
        recovery.maximumSteinerCandidates=options.maximumEdgeRecoverySteinerCandidates;
        recovery.maximumPieceWork=std::min(recovery.maximumPieceWork,po.maximumWork);
        OcctPartitionRecoveryReceipt recovered;
        OcctStageTrace trace(attempt.name,view.body->tag());
        const auto result=recovery.enabled
            ?buildOcctPartitionVolume(model,view.body->tag(),po,recovery,body,&recovered)
            :buildOcctPartitionVolume(model,view.body->tag(),po,body);
        trace.finish(result.diagnostic);
        if(policy==OcctPartitionPolicy::NativeFeatureAxes&&attempt.density==4&&attempt.kernelSweeps==16){
            extendedKernelRequired=(result.reason==OcctPartitionMeshReason::StarFailed&&result.starReason==OcctStarReason::KernelNotFound)||
                (result.reason==OcctPartitionMeshReason::BudgetExceeded&&std::string_view(result.diagnostic)=="work-budget");
            featureDeflection=deflection;
            circularContainmentRequired=result.reason==OcctPartitionMeshReason::GeometryError&&
                (std::string_view(result.diagnostic)=="native-tet-centroid-outside"||
                 std::string_view(result.diagnostic)=="interior-node-outside-original-solid");
        }
        report.recoverySteinerNodes+=recovered.steinerNodes;
        work.tick(result.workUsed);report.partitionWork+=result.partitionWork;
        // The earlier failed native attempt already paid for and checked these
        // original-body cuts. Replay that complete arrangement once with the
        // appended PLC recovery, retaining the earlier work charge and leaving
        // every established successful route unchanged.
        if(!attempt.recoverPlc&&policy==OcctPartitionPolicy::NativeConcaveFeatureAxes&&
           attempt.density==16&&options.useNativeKernelAdmission&&!result.published&&
           result.reason==OcctPartitionMeshReason::StarFailed&&
           result.starReason==OcctStarReason::KernelNotFound&&
           result.adaptiveCutReports.size()==po.maximumAdaptiveCuts&&po.maximumAdaptiveCuts>0){
            work.tick(result.adaptiveCutReports.size());
            retainedRecoveryPlanes.clear();
            retainedRecoveryPlanes.reserve(result.adaptiveCutReports.size());
            for(const auto& cut:result.adaptiveCutReports)retainedRecoveryPlanes.push_back(cut.plane);
            retainedRecoveryDeflection=deflection;
        }
        report.nativeCandidateQueries+=result.nativeCandidateQueries;
        report.surfaceDiagonalRepairs+=result.surfaceDiagonalRepairs;
        report.nativeSurfaceWork+=result.surfaceWork;report.nativeVolumeWork+=result.starWork;
        report.nativeSurfaceSeconds+=result.surfaceSeconds;report.nativeVolumeSeconds+=result.starSeconds;
        report.nativePartitionSeconds+=result.partitionSeconds;report.nativeAssemblySeconds+=result.assemblySeconds;
        report.partitionReason=result.diagnostic;
        report.partitionPolicy=attempt.name;
        report.partitionAttemptReports.push_back({report.partitionPolicy,result.diagnostic,result.failedPieceTag,
            result.pieces,result.interfaces,result.workUsed,result.totalSeconds,deflection});
        report.partitionPieces=result.pieces;report.partitionInterfaces=result.interfaces;
        report.partitionVolumeChange=result.partitionVolumeChange;
        report.partitionVolumeAdmission=po.partition.maximumRelativeVolumeChange;
        if(result.reason==OcctPartitionMeshReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled};
        if(result.reason==OcctPartitionMeshReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed};
        if(result.reason==OcctPartitionMeshReason::BudgetExceeded){
            // A local construction/search admission may fit a later strategy.
            // Its completed work stays charged to the shared cumulative cap.
            const auto why=std::string_view(result.diagnostic);
            if((attempt.liveRecovery||attempt.priority!=Priority::Normal)&&work.count<options.maximumValidationWork)break;
            if(work.count<options.maximumValidationWork&&(why=="work-budget"||why=="partition-budget")){
                admissionRefused=true;break;}
            throw Stop{NativeCadMeshReason::BudgetExceeded};
        }
        if(result.reason==OcctPartitionMeshReason::InternalError)throw Stop{NativeCadMeshReason::InternalError};
        if(!result.published){
            const bool refineSurface=result.reason==OcctPartitionMeshReason::GeometryError||
                result.reason==OcctPartitionMeshReason::VolumeMismatch||
                result.reason==OcctPartitionMeshReason::SurfaceFailed||
                (result.reason==OcctPartitionMeshReason::StarFailed&&result.starReason==OcctStarReason::BoundaryContact);
            if(refineSurface&&refinement+1<attempt.refinements){deflection*=.5;continue;}
            break;
        }
        report.route="native_partition_star";report.nativeSurfaceReason="none";report.nativeVolumeReason="none";
        if(recovered.piecesRecovered>0)report.route="native_partition_core";
        if(attempt.liveRecovery)report.route="native_partition_local_recovery";
        if(policy==OcctPartitionPolicy::NativePrincipalFeatureAxes||attempt.featureRecovery){
            report.cadVolume=result.cadVolume;report.surfaceDeflection=deflection;
        }
        if(recovered.piecesRecovered>0||attempt.priority==Priority::CoarseGrid)report.cadVolume=result.cadVolume;
        report.surfaceTriangles=result.surfaceTriangles;report.exactSurfacePairs=result.exactSurfacePairs;
        report.kernelTrials+=result.kernelTrials;report.meshVolume=result.meshVolume;
        report.nativeCentroidThreads=result.nativeCentroidValidation.threadsUsed;
        report.nativeCentroidsChecked=result.nativeCentroidValidation.nativeQueries;
        report.nativeCentroidReferenceQueries=result.nativeCentroidValidation.originalReferenceQueries;
        report.nativeCentroidGeometryPairs=result.nativeCentroidValidation.copiedGeometryValuePairs;
        report.nativeCentroidSourceGuard=result.nativeCentroidValidation.sourceModelUnmodified;
        report.maximumBoundaryNodeError=result.maximumBoundaryNodeError;
        report.maximumBoundarySampleError=result.maximumBoundarySampleError;
        report.surfaceDeflection=deflection;
        report.boundaryEmbeddingCertified=result.exteriorCertificate.certified;return true;
      }
    }
    if(admissionRefused&&!initialCartesianOnly)throw Stop{NativeCadMeshReason::BudgetExceeded};
    return false;
}
bool nativeLongitudinalFallback(const OcctModel& model,const View& view,const NativeCadMeshOptions& options,
    mesh::Mesh& body,NativeCadBodyMeshReport& report,Work& work) {
    using namespace detail;const auto& geometry=options.geometry;
    if(!options.useLongitudinalRecovery||!options.useNativeKernelAdmission||!options.useNativePartitionFallback)return false;
    const bool cartesian=std::none_of(view.fs.begin(),view.fs.end(),[](const GFace* face){return face->kind()==SurfaceKind::Cylinder;});
    Vec3 low,high;view.body->bbox(low,high);const auto span=high-low;
    unsigned axis=span.y>span.x?1:0;if(span.z>(axis==0?span.x:span.y))axis=2;
    const Vec3 normal=axis==0?Vec3{1,0,0}:axis==1?Vec3{0,1,0}:Vec3{0,0,1};
    struct Attempt {unsigned divisions;const char* nativeName;const char* cartesianName;};
    for(const auto attempt:std::array{Attempt{8,"longitudinal_8_front","cartesian_thirds_longitudinal_8_front"},
                                    Attempt{16,"longitudinal_16_front","cartesian_thirds_longitudinal_16_front"},
                                    Attempt{32,"longitudinal_32_front","cartesian_thirds_longitudinal_32_front"}}){
        const auto* name=cartesian?attempt.cartesianName:attempt.nativeName;
        work.tick();++report.partitionAttempts;
        OcctPartitionMeshOptions request;auto& partition=request.partition;
        // Noncylindrical solids have no native radial axis. Use a bounded
        // Cartesian subdivision rather than refusing that entire face family.
        partition.policy=cartesian?OcctPartitionPolicy::CartesianGrid:OcctPartitionPolicy::NativeConcaveFeatureAxes;
        partition.cartesianDivisions=3;partition.featureAxisPlanes=4;
        partition.includeCartesianMidplanes=false;partition.includePlanarSupportPlanes=true;
        partition.maximumPieces=std::min<std::uint64_t>(512,options.maximumNativePartitionPieces);
        partition.maximumPlanes=48;partition.maximumHistoryWork=10000000;
        // Conservation disagreements are checked by the existing independent
        // native trimming witness. No source tolerance/geometry is enlarged.
        partition.maximumRelativeVolumeChange=1e-6;partition.allowSampledTrimTolerance=!cartesian;
        partition.maximumPhysicalTrimDrift=report.surfaceChordErrorLimit;
        for(unsigned i=1;i<attempt.divisions;++i)partition.extraPlanes.push_back({low+span*(double(i)/attempt.divisions),normal});
        request.maximumWork=std::min(geometry.maximumMeshingWork,options.maximumValidationWork-work.count);
        need(request.maximumWork>0,NativeCadMeshReason::BudgetExceeded);
        request.maximumPieceStarWork=std::min<std::uint64_t>(5000000,request.maximumWork);
        request.maximumNodes=geometry.maximumNodes;request.maximumCells=geometry.maximumCells;
        request.maximumOutputBytes=geometry.maximumOutputBytes;request.maximumWorkingBytes=options.maximumPrimitiveWorkingBytes;
        request.maximumRelativeVolumeError=options.maximumRelativeVolumeError;
        request.maximumNativeValidationThreads=options.maximumNativeValidationThreads;request.cancellation=geometry.cancellation;
        request.useNativeKernelAdmission=true;request.repairSurfaceAliases=true;
        request.boundary.maximumChordError=report.surfaceChordErrorLimit;
        request.boundary.linearDeflection=report.surfaceChordErrorLimit*.4;
        request.boundary.maximumNodes=std::min<std::uint64_t>(50000,geometry.maximumNodes);
        request.boundary.maximumTriangles=std::min<std::uint64_t>(100000,geometry.maximumCells/2);
        request.boundary.maximumWork=request.maximumWork;request.boundary.maximumOutputBytes=geometry.maximumOutputBytes;
        request.star.maximumWork=request.maximumWork;request.star.maximumWorkingBytes=options.maximumPrimitiveWorkingBytes;
        request.star.maximumBoundaryNodes=request.boundary.maximumNodes;request.star.maximumBoundaryTriangles=request.boundary.maximumTriangles;
        request.star.maximumNodes=geometry.maximumNodes;request.star.maximumCells=geometry.maximumCells;
        request.star.maximumOutputBytes=geometry.maximumOutputBytes;request.star.maximumProposalSweeps=16;
        request.star.auditLimits.maximumTemporaryBytes=geometry.maximumAuditTemporaryBytes;
        request.star.auditLimits.maximumFaceUses=geometry.maximumAuditFaceUses;
        OcctPartitionRecoveryOptions recovery;recovery.enabled=true;recovery.useAdvancingFront=true;
        recovery.maximumFrontSeconds=cartesian?15.:3.;
        recovery.onlyFinalAdaptiveAttempt=false;recovery.maximumPieceWork=std::min<std::uint64_t>(3000000,request.maximumWork);
        OcctPartitionRecoveryReceipt recovered;
        const auto result=buildOcctPartitionVolume(model,view.body->tag(),request,recovery,body,&recovered);
        work.tick(result.workUsed);report.partitionWork+=result.partitionWork;
        report.nativeSurfaceWork+=result.surfaceWork;report.nativeVolumeWork+=result.starWork;
        report.nativeSurfaceSeconds+=result.surfaceSeconds;report.nativeVolumeSeconds+=result.starSeconds;
        report.nativePartitionSeconds+=result.partitionSeconds;report.nativeAssemblySeconds+=result.assemblySeconds;
        report.nativeCandidateQueries+=result.nativeCandidateQueries;report.surfaceDiagonalRepairs+=result.surfaceDiagonalRepairs;
        report.partitionReason=result.diagnostic;report.partitionPolicy=name;
        report.partitionPieces=result.pieces;report.partitionInterfaces=result.interfaces;
        report.partitionVolumeChange=result.partitionVolumeChange;report.partitionVolumeAdmission=partition.maximumRelativeVolumeChange;
        report.sampledPartitionTrimWitness=result.trimTolerance.admitted;
        report.partitionTrimSamples=result.trimTolerance.samples;
        report.maximumPartitionTrimDrift=result.trimTolerance.maximumPhysicalDrift;
        report.maximumPartitionTrimToleranceRatio=result.trimTolerance.maximumToleranceRatio;
        report.partitionAttemptReports.push_back({name,result.diagnostic,result.failedPieceTag,result.pieces,result.interfaces,result.workUsed,result.totalSeconds,request.boundary.linearDeflection});
        if(result.reason==OcctPartitionMeshReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled};
        if(result.reason==OcctPartitionMeshReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed};
        if(result.reason==OcctPartitionMeshReason::InternalError)throw Stop{NativeCadMeshReason::InternalError};
        if(!result.published)continue;
        report.route=cartesian?"native_cartesian_front":"native_longitudinal_front";report.nativeSurfaceReason="none";report.nativeVolumeReason="none";
        report.cadVolume=result.cadVolume;report.meshVolume=result.meshVolume;report.surfaceTriangles=result.surfaceTriangles;
        report.exactSurfacePairs=result.exactSurfacePairs;report.kernelTrials+=result.kernelTrials;
        report.maximumBoundaryNodeError=result.maximumBoundaryNodeError;report.maximumBoundarySampleError=result.maximumBoundarySampleError;
        report.surfaceDeflection=request.boundary.linearDeflection;report.boundaryEmbeddingCertified=result.exteriorCertificate.certified;
        report.nativeCentroidThreads=result.nativeCentroidValidation.threadsUsed;
        report.nativeCentroidsChecked=result.nativeCentroidValidation.nativeQueries;
        report.nativeCentroidReferenceQueries=result.nativeCentroidValidation.originalReferenceQueries;
        report.nativeCentroidGeometryPairs=result.nativeCentroidValidation.copiedGeometryValuePairs;
        report.nativeCentroidSourceGuard=result.nativeCentroidValidation.sourceModelUnmodified;
        return true;
    }
    return false;
}
// Cost counters concern both attempted routes. Geometry, quality and source
// certificates remain those of the route that is actually published.
bool accumulateInitialFrontTrial(NativeCadBodyMeshReport& result,const NativeCadBodyMeshReport& trial) noexcept {
    const std::array workCounters{
        &NativeCadBodyMeshReport::nativeSurfaceWork,&NativeCadBodyMeshReport::nativeVolumeWork,
        &NativeCadBodyMeshReport::partitionWork,&NativeCadBodyMeshReport::boxContactWork,
        &NativeCadBodyMeshReport::contactCylinderWork,&NativeCadBodyMeshReport::surfaceAttempts,
        &NativeCadBodyMeshReport::rebuiltSurfaceCharts,&NativeCadBodyMeshReport::exactSurfacePairs,
        &NativeCadBodyMeshReport::kernelTrials,&NativeCadBodyMeshReport::partitionAttempts,
        &NativeCadBodyMeshReport::partitionTrimSamples,&NativeCadBodyMeshReport::surfaceDiagonalRepairs,
        &NativeCadBodyMeshReport::nativeCentroidsChecked,&NativeCadBodyMeshReport::nativeCentroidReferenceQueries,
        &NativeCadBodyMeshReport::nativeCentroidExactReuses,&NativeCadBodyMeshReport::nativeCentroidGeometryPairs,
        &NativeCadBodyMeshReport::centroidsWithinChordTolerance,&NativeCadBodyMeshReport::nativeCentroidDistanceQueries,
        &NativeCadBodyMeshReport::nativeChordSamplesChecked,&NativeCadBodyMeshReport::chordSamplesWithinTolerance,
        &NativeCadBodyMeshReport::nativeCandidateQueries,&NativeCadBodyMeshReport::nativeCandidateRayAdmissions,
        &NativeCadBodyMeshReport::nativeFinalRayAdmissions};
    for(const auto member:workCounters){
        if(trial.*member>UINT64_MAX-result.*member)return false;
        result.*member+=trial.*member;
    }
    result.nativeSurfaceSeconds+=trial.nativeSurfaceSeconds;
    result.nativeVolumeSeconds+=trial.nativeVolumeSeconds;
    result.nativePartitionSeconds+=trial.nativePartitionSeconds;
    result.nativeAssemblySeconds+=trial.nativeAssemblySeconds;
    return true;
}
bool accumulateInitialFrontTrial(NativeCadMeshOutcome& result,const NativeCadMeshOutcome& trial) noexcept {
    const std::array queryCounters{
        &NativeCadMeshOutcome::nativeProjectionQueries,&NativeCadMeshOutcome::nativeProjectionInitializations,
        &NativeCadMeshOutcome::nativeProjectionReuses,&NativeCadMeshOutcome::nativeTrimmedProjectionQueries,
        &NativeCadMeshOutcome::nativeTrimmedProjectionMemoHits};
    for(const auto member:queryCounters){
        if(trial.*member>UINT64_MAX-result.*member)return false;
        result.*member+=trial.*member;
    }
    result.preparationSeconds+=trial.preparationSeconds;
    result.generationSeconds+=trial.generationSeconds;
    result.validationSeconds+=trial.validationSeconds;
    result.auditSeconds+=trial.auditSeconds;
    result.qualitySeconds+=trial.qualitySeconds;
    for(auto& body:result.bodies){
        const auto found=std::find_if(trial.bodies.begin(),trial.bodies.end(),[&](const auto& previous){
            return previous.solidTag==body.solidTag;
        });
        if(found!=trial.bodies.end()&&!accumulateInitialFrontTrial(body,*found))return false;
    }
    return true;
}
}
NativeCadMeshOutcome meshNativeCad(const OcctModel& model,const NativeCadMeshOptions& o,mesh::Mesh& destination) noexcept {
    // A positive allowance too small for even a snapshot header is resource
    // exhaustion. Do not turn it into an invalid derived zero-sized snapshot.
    // Preserve historical Bounded request handling and observed cancellation.
    if(o.construction.resourcePolicy==CadResourcePolicy::Normal&&
       o.maximumPrimitiveWorkingBytes>0&&o.maximumPrimitiveWorkingBytes<8){
        NativeCadMeshOutcome limited;
        limited.reason=o.geometry.cancellation&&o.geometry.cancellation->load(std::memory_order_relaxed)?
            NativeCadMeshReason::Cancelled:NativeCadMeshReason::BudgetExceeded;
        return limited;
    }
    if(!o.useAdvancingInteriorFront.has_value()){
        const auto automaticStarted=Clock::now();
        NativeCadMeshOutcome trial;
        NativeCadInitialFrontFallbackReport fallbackReport;
        // No allocation may occur while returning a result after publication.
        static_assert(std::is_nothrow_move_constructible_v<NativeCadMeshOutcome>);
        const auto failedAutomatic=[&](NativeCadMeshReason reason){
            trial.published=false;trial.destinationUnchanged=true;trial.reason=reason;
            trial.initialFrontFallback=fallbackReport;
            for(auto& body:trial.bodies){body.reason=reason;body.initialFrontFallback=fallbackReport;}
            trial.totalSeconds=elapsed(automaticStarted);
            return std::move(trial);
        };
        try {
            const auto solidCount=model.solids().size();
            const bool enabled=o.advancingInteriorFrontEnabled(solidCount);
            if(enabled&&o.maximumValidationWork!=o.geometry.maximumMeshingWork){
                // AUTO construction, source checks and final validation share
                // both original ceilings, including an independent wrapper's
                // own work. Preserve AUTO for the bounded child dispatch; equal
                // caps make this normalization re-entry happen at most once.
                auto bounded=o;
                const auto sharedWork=std::min(o.maximumValidationWork,o.geometry.maximumMeshingWork);
                bounded.maximumValidationWork=sharedWork;
                bounded.geometry.maximumMeshingWork=sharedWork;
                auto result=meshNativeCad(model,bounded,destination);
                result.totalSeconds=elapsed(automaticStarted);
                return result;
            }
            // A single independent part retains its wrapper/source audits. Its
            // conforming child resolves AUTO once; multi-part requests resolve
            // false here so their children retain the established dispatch.
            // The wrapper changes the child policy to ConformingDomain, so
            // deferral cannot recurse back into the independent-part wrapper.
            const bool deferToSinglePart=enabled&&solidCount==1&&
                o.assemblyPolicy==NativeCadAssemblyPolicy::IndependentParts;
            if(!deferToSinglePart){
                auto resolved=o;resolved.useAdvancingInteriorFront=enabled;
                const auto trialStarted=Clock::now();
                trial=meshNativeCad(model,resolved,destination);
                const auto trialSeconds=elapsed(trialStarted);
                const bool sourceRefusal=trial.reason==NativeCadMeshReason::GeometryError||
                    trial.reason==NativeCadMeshReason::UnsupportedBody||
                    trial.reason==NativeCadMeshReason::RecognitionFailed||
                    trial.reason==NativeCadMeshReason::MeshingFailed;
                const bool beforeInitialStage=!trial.bodies.empty()&&
                    std::all_of(trial.bodies.begin(),trial.bodies.end(),[](const auto& body){
                        return body.initialFront.diagnostic&&
                            std::string_view(body.initialFront.diagnostic)=="not_attempted";
                    });
                if(enabled&&!trial.published&&trial.destinationUnchanged&&sourceRefusal&&beforeInitialStage&&
                   trial.validationWork<o.maximumValidationWork&&trial.validationWork<o.geometry.maximumMeshingWork){
                    if(o.geometry.cancellation&&o.geometry.cancellation->load(std::memory_order_relaxed))
                        return failedAutomatic(NativeCadMeshReason::Cancelled);
                    fallbackReport={true,trial.reason,trial.validationWork,trialSeconds,
                        trial.bodies.front().nativeSurfaceReason?trial.bodies.front().nativeSurfaceReason:"not_attempted"};
                    auto legacy=o;legacy.useAdvancingInteriorFront=false;
                    // Work::tick enforces the validation admission while some
                    // legacy constructors separately enforce the meshing cap.
                    // Give execution their shared remaining allowance, rather
                    // than discover a meshing-cap overrun after the retry.
                    const auto retryWork=std::min(o.maximumValidationWork-trial.validationWork,
                        o.geometry.maximumMeshingWork-trial.validationWork);
                    legacy.maximumValidationWork=retryWork;
                    legacy.geometry.maximumMeshingWork=retryWork;
                    // The retry publishes into a private destination. Cumulative
                    // budgets, counters and cancellation are checked before the
                    // caller's destination can be changed.
                    mesh::Mesh candidate;
                    auto retry=meshNativeCad(model,legacy,candidate);
                    retry.destinationUnchanged=true;
                    retry.initialFrontFallback=fallbackReport;
                    for(auto& body:retry.bodies)body.initialFrontFallback=fallbackReport;
                    const bool withinBudget=retry.validationWork<=legacy.maximumValidationWork&&
                        retry.validationWork<=legacy.geometry.maximumMeshingWork;
                    if(!withinBudget){retry.published=false;retry.reason=NativeCadMeshReason::BudgetExceeded;}
                    if(trial.validationWork>UINT64_MAX-retry.validationWork){
                        retry.published=false;retry.reason=NativeCadMeshReason::InternalError;
                    }else retry.validationWork+=trial.validationWork;
                    if(!accumulateInitialFrontTrial(retry,trial)){
                        retry.published=false;retry.reason=NativeCadMeshReason::InternalError;
                    }
                    if(o.geometry.cancellation&&o.geometry.cancellation->load(std::memory_order_relaxed)){
                        retry.published=false;retry.reason=NativeCadMeshReason::Cancelled;
                    }
                    if(retry.published){destination.swap(candidate);retry.destinationUnchanged=false;}
                    else for(auto& body:retry.bodies)body.reason=retry.reason;
                    retry.totalSeconds=elapsed(automaticStarted);
                    return retry;
                }
                trial.totalSeconds=elapsed(automaticStarted);
                return trial;
            }
        }catch(const Standard_OutOfMemory&){return failedAutomatic(NativeCadMeshReason::AllocationFailed);}
         catch(const std::bad_alloc&){return failedAutomatic(NativeCadMeshReason::AllocationFailed);}
         catch(...){return failedAutomatic(NativeCadMeshReason::InternalError);}
    }
    if(o.useAdaptiveTetRecovery&&o.reuseNativeQueries){
        const auto retained=detail::OcctQueryContext::maximumRetainedPayloadBytes(true);
        if(retained<=o.maximumPrimitiveWorkingBytes/8){
            try {
                auto remaining=o;
                remaining.maximumPrimitiveWorkingBytes-=retained;
                remaining.reuseNativeQueries=false;
                detail::OcctQueryContext queries(true,true);
                auto result=meshNativeCad(model,remaining,destination);
                const auto counters=queries.counters();
                result.nativeProjectionQueries=counters.surfaceProjectionQueries;
                result.nativeProjectionInitializations=counters.surfaceProjectionInitializations;
                result.nativeProjectionReuses=counters.surfaceProjectionReuses;
                result.nativeTrimmedProjectionQueries=counters.trimmedProjectionQueries;
                result.nativeTrimmedProjectionMemoHits=counters.trimmedProjectionMemoHits;
                return result;
            }catch(const std::bad_alloc&){NativeCadMeshOutcome result;result.reason=NativeCadMeshReason::AllocationFailed;return result;}
             catch(...){NativeCadMeshOutcome result;result.reason=NativeCadMeshReason::InternalError;return result;}
        }
    }
    detail::OcctStageObserverScope observe(o.stageObserver);
    detail::OcctStageTrace meshTrace("native_mesh");
    const bool spatialSizing=o.useSpatialSizing||!o.sizeBoxes.empty()||o.useAdvancingInteriorFront.value_or(false);
    if(o.maximumRecoveryRandomRetries>4||o.maximumEdgeRecoverySteinerCandidates>64||o.maximumRecoveryPatchStates>65536||o.curvaturePointsPerTurn>64||
       (o.construction.effort!=CadConstructionEffort::Legacy&&o.construction.effort!=CadConstructionEffort::Automatic&&o.construction.effort!=CadConstructionEffort::Thorough)||
       (o.construction.resourcePolicy!=CadResourcePolicy::Bounded&&o.construction.resourcePolicy!=CadResourcePolicy::Normal)||
       (o.recoveryMode!=NativeCadRecoveryMode::Standard&&o.recoveryMode!=NativeCadRecoveryMode::Advanced)||
       (o.assemblyPolicy!=NativeCadAssemblyPolicy::ConformingDomain&&o.assemblyPolicy!=NativeCadAssemblyPolicy::IndependentParts)){NativeCadMeshOutcome invalid;invalid.reason=NativeCadMeshReason::InvalidRequest;return invalid;}
    if(spatialSizing&&(!o.useNativeSurfaceFallback||
        o.assemblyPolicy!=NativeCadAssemblyPolicy::ConformingDomain||o.separateMatchedNativeFaceContacts)){
        NativeCadMeshOutcome unsupported;
        unsupported.reason=NativeCadMeshReason::UnsupportedControl;
        return unsupported;
    }
    if(o.assemblyPolicy==NativeCadAssemblyPolicy::IndependentParts){
        detail::OcctStageTrace attempt("independent_parts");
        auto result=detail::buildOcctIndependentParts(model,o,destination);
        attempt.finish(nativeCadMeshReasonName(result.reason));
        return result;
    }
    if(o.separateMatchedNativeFaceContacts) {
        const auto began=Clock::now();NativeCadMeshOutcome result;
        std::unique_ptr<detail::OcctQueryContext> queries;
        std::uint64_t spent=0;
        try {
            const auto limit=std::min(o.maximumValidationWork,o.geometry.maximumMeshingWork);
            need(limit>=4,NativeCadMeshReason::BudgetExceeded);
            need(!o.geometry.cancellation||!o.geometry.cancellation->load(std::memory_order_relaxed),NativeCadMeshReason::Cancelled);
            if(o.reuseNativeQueries){
                need(detail::OcctQueryContext::maximumRetainedPayloadBytes()<=o.maximumPrimitiveWorkingBytes,NativeCadMeshReason::BudgetExceeded);
                queries=std::make_unique<detail::OcctQueryContext>(true);
            }
            std::string before;
            const auto snapshotLimit=std::min<std::uint64_t>(64ULL*1024*1024,o.maximumPrimitiveWorkingBytes/8);
            auto sourceReason=detail::snapshotOcctContactSource(model,snapshotLimit,before);++spent;
            need(sourceReason==NativeCadMeshReason::None,sourceReason);
            auto legacy=o;legacy.separateMatchedNativeFaceContacts=false;legacy.reuseNativeQueries=false;
            legacy.maximumValidationWork=limit-spent;legacy.geometry.maximumMeshingWork=limit-spent;
            mesh::Mesh candidate;result=meshNativeCad(model,legacy,candidate);
            need(result.validationWork<=limit-spent,NativeCadMeshReason::BudgetExceeded);
            spent+=result.validationWork;
            std::string after;
            need(spent<limit,NativeCadMeshReason::BudgetExceeded);
            sourceReason=detail::snapshotOcctContactSource(model,snapshotLimit,after);++spent;
            need(sourceReason==NativeCadMeshReason::None,sourceReason);
            need(before==after,NativeCadMeshReason::InternalError);
            need(!o.geometry.cancellation||!o.geometry.cancellation->load(std::memory_order_relaxed),NativeCadMeshReason::Cancelled);
            if(result.published) {destination.swap(candidate);result.destinationUnchanged=false;}
            else if(result.reason==NativeCadMeshReason::UnresolvedBodyInterface) {
                auto remaining=o;remaining.reuseNativeQueries=false;
                remaining.maximumValidationWork=limit-spent;remaining.geometry.maximumMeshingWork=limit-spent;
                need(remaining.maximumValidationWork>0,NativeCadMeshReason::BudgetExceeded);
                const double preparation=result.preparationSeconds,generation=result.generationSeconds,validation=result.validationSeconds;
                // Preserve an already published legacy assembly. The explicit
                // native-kernel route can mesh proper face-contact components
                // together while retaining each isolated body's existing path.
                if(o.useNativeKernelAdmission) {
                    detail::OcctStageTrace attempt("contact_components");
                    result=detail::buildOcctContactComponents(model,remaining,candidate,&before);
                    attempt.finish(result.nativeContactDiagnostic);
                    need(result.validationWork<=limit-spent,NativeCadMeshReason::BudgetExceeded);
                    spent+=result.validationWork;
                    const bool notApplicable=result.reason==NativeCadMeshReason::UnsupportedBody&&
                        std::string_view(result.nativeContactDiagnostic)=="component_arrangement_not_applicable";
                    // A bounded advanced-strategy refusal must not suppress an
                    // established whole-assembly route. Cancellation, allocation,
                    // source/internal failures, material overlap and requested
                    // quality/size refusals remain terminal. The whole route still
                    // validates every original node, Tet centroid and interface.
                    const bool childRefused=!result.published&&
                        std::string_view(result.nativeContactDiagnostic)=="component_child_not_published"&&
                        (result.reason==NativeCadMeshReason::BudgetExceeded||
                         result.reason==NativeCadMeshReason::UnsupportedBody||
                         result.reason==NativeCadMeshReason::MeshingFailed||
                         result.reason==NativeCadMeshReason::GeometryError);
                    if(notApplicable||childRefused) {
                        const double componentPreparation=result.preparationSeconds;
                        const double componentGeneration=result.generationSeconds;
                        const double componentValidation=result.validationSeconds;
                        const double componentAudit=result.auditSeconds,componentQuality=result.qualitySeconds;
                        if(childRefused) {
                            need(result.destinationUnchanged,NativeCadMeshReason::InternalError);
                            need(spent<limit,NativeCadMeshReason::BudgetExceeded);
                            sourceReason=detail::snapshotOcctContactSource(model,snapshotLimit,after);++spent;
                            need(sourceReason==NativeCadMeshReason::None,sourceReason);
                            need(before==after,NativeCadMeshReason::InternalError);
                            need(!o.geometry.cancellation||!o.geometry.cancellation->load(std::memory_order_relaxed),NativeCadMeshReason::Cancelled);
                            remaining.useNativeKernelAdmission=false;
                        }
                        remaining.maximumValidationWork=limit-spent;
                        remaining.geometry.maximumMeshingWork=limit-spent;
                        need(remaining.maximumValidationWork>0,NativeCadMeshReason::BudgetExceeded);
                        detail::OcctStageTrace whole("contact_assembly");
                        result=detail::buildOcctContactAssembly(model,remaining,candidate,&before);
                        whole.finish(result.nativeContactDiagnostic);
                        need(result.validationWork<=limit-spent,NativeCadMeshReason::BudgetExceeded);
                        spent+=result.validationWork;
                        result.preparationSeconds+=componentPreparation;
                        result.generationSeconds+=componentGeneration;
                        result.validationSeconds+=componentValidation;
                        result.auditSeconds+=componentAudit;result.qualitySeconds+=componentQuality;
                    }
                } else {
                    detail::OcctStageTrace whole("contact_assembly");
                    result=detail::buildOcctContactAssembly(model,remaining,candidate,&before);
                    whole.finish(result.nativeContactDiagnostic);
                    need(result.validationWork<=limit-spent,NativeCadMeshReason::BudgetExceeded);
                    spent+=result.validationWork;
                }
                result.preparationSeconds+=preparation;result.generationSeconds+=generation;result.validationSeconds+=validation;
                if(result.published) {
                    need(spent<limit,NativeCadMeshReason::BudgetExceeded);
                    sourceReason=detail::snapshotOcctContactSource(model,snapshotLimit,after);++spent;
                    need(sourceReason==NativeCadMeshReason::None,sourceReason);
                    need(before==after,NativeCadMeshReason::InternalError);
                    need(!o.geometry.cancellation||!o.geometry.cancellation->load(std::memory_order_relaxed),NativeCadMeshReason::Cancelled);
                    destination.swap(candidate);
                }
            }
        }catch(const Stop& s){result.published=false;result.destinationUnchanged=true;result.reason=s.reason;}
         catch(const Standard_OutOfMemory&){result.published=false;result.destinationUnchanged=true;result.reason=NativeCadMeshReason::AllocationFailed;}
         catch(const std::bad_alloc&){result.published=false;result.destinationUnchanged=true;result.reason=NativeCadMeshReason::AllocationFailed;}
         catch(...){result.published=false;result.destinationUnchanged=true;result.reason=NativeCadMeshReason::InternalError;}
        if(queries){const auto c=queries->counters();result.nativeProjectionQueries=c.surfaceProjectionQueries;
            result.nativeProjectionInitializations=c.surfaceProjectionInitializations;result.nativeProjectionReuses=c.surfaceProjectionReuses;
            result.nativeTrimmedProjectionQueries=c.trimmedProjectionQueries;
            result.nativeTrimmedProjectionMemoHits=c.trimmedProjectionMemoHits;}
        result.validationWork=spent;result.totalSeconds=elapsed(began);return result;
    }
    NativeCadMeshOutcome out; const auto started=Clock::now();Work work{o};
    std::unique_ptr<detail::OcctQueryContext> queryContext;
    try {
        if(o.reuseNativeQueries){
            need(detail::OcctQueryContext::maximumRetainedPayloadBytes()<=o.maximumPrimitiveWorkingBytes,NativeCadMeshReason::BudgetExceeded);
            queryContext=std::make_unique<detail::OcctQueryContext>(true);
        }
        need(o.geometry.targetSize>0&&finite(o.geometry.targetSize)&&o.geometry.threads>=0&&o.geometry.threads<=1024&&
             finite(o.geometry.relativeRecognitionTolerance)&&o.geometry.relativeRecognitionTolerance>0&&o.geometry.relativeRecognitionTolerance<1e-3&&
             o.maximumRelativeVolumeError>0&&o.maximumRelativeVolumeError<1&&finite(o.maximumRelativeVolumeError)&&
             o.maximumChordError>=0&&finite(o.maximumChordError)&&o.minimumSize>=0&&o.minimumSize<=o.geometry.targetSize&&finite(o.minimumSize)&&
             o.maximumNativeValidationThreads>=1&&o.maximumNativeValidationThreads<=8&&
             finite(o.geometry.minimumElementQuality)&&o.geometry.minimumElementQuality>=0&&o.geometry.minimumElementQuality<=1&&
             finite(o.geometry.warningElementQuality)&&o.geometry.warningElementQuality>=0&&o.geometry.warningElementQuality<=1&&
             (o.elements==algo::PrimitiveElements::Tet4||o.elements==algo::PrimitiveElements::HexDominant),NativeCadMeshReason::InvalidRequest);
        auto solids=model.solids(); need(!solids.empty(),NativeCadMeshReason::InvalidRequest);
        need(!spatialSizing||solids.size()==1,NativeCadMeshReason::UnsupportedControl);
        work.tick(); mesh::Mesh local;
        detail::OcctBoxCertificateSession boxCertificates;
        auto* cachedBoxCertificates=o.useAdaptiveTetRecovery?&boxCertificates:nullptr;
        PhaseTimer preparationTimer(out.preparationSeconds);
        detail::OcctStageTrace preparing("assembly_preparation");
        std::set<int> controlledBodies;
        for(const auto& control:o.geometry.bodyControls)need(model.solid(control.solidTag)&&controlledBodies.insert(control.solidTag).second,NativeCadMeshReason::InvalidRequest);
        need(o.minimumSize==0||!o.faceSizes.empty()||spatialSizing,NativeCadMeshReason::UnsupportedControl);
        std::set<int> selected;
        for(const auto& control:o.faceSizes) {
            need(model.face(control.faceTag)&&control.size>0&&finite(control.size)&&control.growth>=1&&
                 finite(control.growth)&&selected.insert(control.faceTag).second,NativeCadMeshReason::InvalidRequest);
            if(spatialSizing) {
                const auto owners=model.faceOwners(control.faceTag);
                need(control.size<=o.geometry.targetSize&&owners.size()==1&&
                     owners.front().solidTag==solids.front()->tag(),NativeCadMeshReason::InvalidRequest);
            }
        }
        std::vector<std::pair<int,int>> separateContacts;
        std::map<std::pair<int,int>,double> boxFaceContacts;
        std::set<int> boxContactBodies;
        std::vector<BodyRange> bodyRanges;
        for(std::size_t i=0;i<solids.size();++i) for(std::size_t j=0;j<i;++j) {
            work.tick(); Vec3 a,b,c,d;solids[i]->bbox(a,b);solids[j]->bbox(c,d);
            if(a.x>d.x||c.x>b.x||a.y>d.y||c.y>b.y||a.z>d.z||c.z>b.z)continue;
            BRepExtrema_DistShapeShape distance(detail::OcctModelAccess::solid(model,solids[i]->tag()),detail::OcctModelAccess::solid(model,solids[j]->tag()));
            need(distance.IsDone(),NativeCadMeshReason::UnresolvedBodyInterface);
            if(distance.Value()>1e-7)continue;
            need(o.separateLowerDimensionalContacts||o.separateMatchedBoxFaceContacts,NativeCadMeshReason::UnresolvedBodyInterface);
            // A contact policy is explicit: preserve distinct nodes/boundaries.
            // OCCT's non-destructive common must contain no material or faces.
            BRepAlgoAPI_Common common;NCollection_List<TopoDS_Shape> arguments,tools;
            arguments.Append(detail::OcctModelAccess::solid(model,solids[i]->tag()));
            tools.Append(detail::OcctModelAccess::solid(model,solids[j]->tag()));
            common.SetArguments(arguments);common.SetTools(tools);common.SetNonDestructive(true);common.Build();work.tick();
            need(common.IsDone(),NativeCadMeshReason::UnresolvedBodyInterface);
            need(!TopExp_Explorer(common.Shape(),TopAbs_SOLID).More(),NativeCadMeshReason::UnresolvedBodyInterface);
            // Solid Common alone discards even coincident face contacts.
            // Intersect compounds of the boundary faces to distinguish them.
            BRep_Builder builder;TopoDS_Compound af,bf;builder.MakeCompound(af);builder.MakeCompound(bf);
            for(TopExp_Explorer f(arguments.First(),TopAbs_FACE);f.More();f.Next()){work.tick();builder.Add(af,f.Current());}
            for(TopExp_Explorer f(tools.First(),TopAbs_FACE);f.More();f.Next()){work.tick();builder.Add(bf,f.Current());}
            BRepAlgoAPI_Common skinCommon;NCollection_List<TopoDS_Shape> skinsA,skinsB;skinsA.Append(af);skinsB.Append(bf);
            skinCommon.SetArguments(skinsA);skinCommon.SetTools(skinsB);skinCommon.SetNonDestructive(true);skinCommon.Build();work.tick();
            need(skinCommon.IsDone(),NativeCadMeshReason::UnresolvedBodyInterface);
            const std::pair<int,int> pair{solids[j]->tag(),solids[i]->tag()};
            if(TopExp_Explorer(common.Shape(),TopAbs_FACE).More()||TopExp_Explorer(skinCommon.Shape(),TopAbs_FACE).More()){
                need(o.separateMatchedBoxFaceContacts,NativeCadMeshReason::UnresolvedBodyInterface);
                const auto admission=detail::admitOcctBoxFaceContact(model,pair.first,pair.second,
                    boxContactOptions(o,work,cachedBoxCertificates));
                work.tick(admission.workUsed);checkBoxContactFailure(admission.reason);
                need(admission.admitted&&admission.sourceImmutabilityVerified,NativeCadMeshReason::UnresolvedBodyInterface);
                boxFaceContacts.emplace(pair,admission.nativeSkinArea);
                boxContactBodies.insert(pair.first);boxContactBodies.insert(pair.second);
            }else need(o.separateLowerDimensionalContacts,NativeCadMeshReason::UnresolvedBodyInterface);
            separateContacts.push_back(pair);
        }
        preparationTimer.finish();preparing.finish();
        for(auto solid:solids) {
            work.tick(); out.bodies.emplace_back(); auto& report=out.bodies.back(); report.solidTag=solid->tag();report.cadVolume=model.volume(solid->tag());
            PhaseTimer recognizing(out.preparationSeconds);View view(model,solid);
            need(view.fs.size()+view.es.size()+view.vs.size()<=o.geometry.maximumAnalysisWork,NativeCadMeshReason::BudgetExceeded);
            Recognition r;NativeCadMeshReason recognitionReason=NativeCadMeshReason::None;
            if(!spatialSizing) try{r=recognize(model,view,o.geometry.relativeRecognitionTolerance);}catch(const Stop& stopped){
                if(!o.useNativeSurfaceFallback||(stopped.reason!=NativeCadMeshReason::UnsupportedBody&&stopped.reason!=NativeCadMeshReason::RecognitionFailed))throw;
                recognitionReason=stopped.reason;
            }
            mesh::Mesh body;bool nativeSurface=false,originalFacetCentersCertified=false;
            recognizing.finish();
            std::vector<algo::GeometryBodyControl> controls;
            for(const auto& control:o.geometry.bodyControls)if(control.solidTag==solid->tag())controls.push_back(control);
            std::vector<CadFaceSizeControl> faceControls;
            for(const auto& control:o.faceSizes){const auto owners=model.faceOwners(control.faceTag);
                if(std::any_of(owners.begin(),owners.end(),[&](auto owner){return owner.solidTag==solid->tag();}))faceControls.push_back(control);}
            const bool earlierContact=std::any_of(separateContacts.begin(),separateContacts.end(),[&](const auto& contact){
                const int other=contact.first==solid->tag()?contact.second:contact.second==solid->tag()?contact.first:-1;
                return other>0&&std::any_of(bodyRanges.begin(),bodyRanges.end(),[other](const auto& range){return range.tag==other;});});
            if(earlierContact&&r.analytic&&r.shape.kind==algo::PrimitiveKind::Cylinder&&
               o.useNativeSurfaceFallback&&o.elements==algo::PrimitiveElements::Tet4&&controls.empty()&&faceControls.empty()&&
               !o.geometry.spatialField&&o.geometry.strategy==algo::GeometryStrategy::Auto&&
               o.geometry.maximumElementSizeRatio==0&&o.geometry.warningElementSizeRatio==0){
                PhaseTimer generating(out.generationSeconds);
                nativeSurface=nativeContactCylinder(model,view,local,o,body,report,work);
                if(r.faces.contains(2))report.startFaceTag=r.faces.at(2);if(r.faces.contains(3))report.endFaceTag=r.faces.at(3);
            }
            if(nativeSurface) {
                // This body already has original native classifications and a
                // certified boundary; final global contact/audit gates follow.
            } else if(!r.analytic) {
                need(spatialSizing||faceControls.empty(),NativeCadMeshReason::UnsupportedControl);
                auto geometry=o.geometry;geometry.bodyControls=controls;
                if(geometry.bodyControls.empty()) {
                    algo::GeometryBodyControl control;control.solidTag=solid->tag();
                    control.elements=o.elements==algo::PrimitiveElements::Tet4?algo::GeometryElements::Tet4:algo::GeometryElements::Native;
                    geometry.bodyControls.push_back(control);
                }
                report.route="planar_core";const auto generating=Clock::now();
                if(!spatialSizing&&recognitionReason==NativeCadMeshReason::None)report.planar=algo::meshGeometryModel(view,geometry,body);
                out.generationSeconds+=elapsed(generating);
                if(!report.planar.published) {
                    if(report.planar.reason==algo::GeometryReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled};
                    if(report.planar.reason==algo::GeometryReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed};
                    if(report.planar.reason==algo::GeometryReason::QualityRequirementNotMet)throw Stop{NativeCadMeshReason::QualityRequirementNotMet};
                    if(report.planar.reason==algo::GeometryReason::SizeRequirementNotMet)throw Stop{NativeCadMeshReason::SizeRequirementNotMet};
                    if(report.planar.reason==algo::GeometryReason::AnalysisBudget||report.planar.reason==algo::GeometryReason::MeshingBudget||report.planar.reason==algo::GeometryReason::OutputBudget)throw Stop{NativeCadMeshReason::BudgetExceeded};
                }
                if(!report.planar.published&&o.useNativeSurfaceFallback){
                    need(controls.empty()&&!o.geometry.spatialField&&o.geometry.strategy==algo::GeometryStrategy::Auto&&
                         o.elements==algo::PrimitiveElements::Tet4&&o.geometry.maximumElementSizeRatio==0&&o.geometry.warningElementSizeRatio==0,
                         NativeCadMeshReason::UnsupportedControl);
                    using detail::OcctFrontSurface;
                    const bool liveRecovery=o.useLiveMaterialRecovery||o.useLocalDelaunayRefinement;
                    const auto recoverFront=[&](std::uint64_t maximumWork,std::uint64_t earlyConstraintWork,
                                                OcctFrontSurface surfaceMode=OcctFrontSurface::Standard,
                                                bool liveMaterialTrial=false){
                        const bool spaceSurface=surfaceMode!=OcctFrontSurface::Standard;
                        const auto recovered=detail::buildOcctFrontVolume(model,solid->tag(),o,
                            maximumWork,body,earlyConstraintWork,surfaceMode,liveMaterialTrial);
                        work.tick(recovered.workUsed);report.nativeSurfaceWork+=recovered.surfaceWork;
                        report.nativeVolumeWork+=recovered.workUsed-recovered.surfaceWork;
                        report.nativeSurfaceSeconds+=recovered.surfaceSeconds;report.nativeVolumeSeconds+=recovered.volumeSeconds;
                        report.nativeCandidateQueries+=recovered.nativeQueries;
                        report.nativeCandidateRayAdmissions+=recovered.nativeRayAdmissions;
                        report.liveMaterialInitialRejects+=recovered.liveMaterialInitialRejects;
                        report.liveMaterialRemainingRejects+=recovered.liveMaterialRemainingRejects;
                        report.liveMaterialPatches+=recovered.liveMaterialPatches;
                        report.localDelaunayUsedGlobalFallback|=recovered.localDelaunayUsedGlobalFallback;
                        if(recovered.reason==NativeCadMeshReason::Cancelled||recovered.reason==NativeCadMeshReason::AllocationFailed||recovered.reason==NativeCadMeshReason::InternalError)throw Stop{recovered.reason};
                        if(earlyConstraintWork&&!recovered.published)return false;
                        if((spaceSurface||o.recoveryMode==NativeCadRecoveryMode::Advanced)&&recovered.surfaceDiagnostic){
                            report.nativeSurfaceReason=recovered.surfaceDiagnostic;
                            report.surfaceTriangles=recovered.surfaceTriangles;
                            report.surfaceChordErrorLimit=recovered.chord;
                            report.maximumBoundaryNodeError=recovered.maximumNodeError;
                            report.maximumBoundarySampleError=recovered.maximumFacetError;
                            report.rebuiltSurfaceCharts=recovered.rebuiltSurfaceCharts;
                        }
                        report.nativeVolumeReason=recovered.diagnostic;
                        report.recoveryRemovedInteriorNodes+=recovered.removedInteriorNodes;report.recoveryRemovedIncidentCells+=recovered.removedIncidentCells;report.recoveryRemovedCentroidCells+=recovered.removedCentroidCells;report.recoveryRandomRetries+=recovered.randomRetries;report.recoveryRandomCandidateProposals+=recovered.randomCandidateProposals;
                        report.recoveryMergedComponents+=recovered.mergedComponents;report.recoverySourceBlockedCandidates+=recovered.sourceBlockedCandidates;report.recoveryOtherComponentBlockedCandidates+=recovered.otherComponentBlockedCandidates;report.recoverySourceFacetNodes+=recovered.sourceFacetNodes;
                        report.recoveryRegions=recovered.regionRecovery;
                        report.constraintRecovery=recovered.constraintRecovery;
                        report.initialFront=recovered.initialFront;
                        report.surfaceSpacing=std::move(recovered.surfaceSpacing);
                        report.interiorFrontDiagnostic=recovered.interiorFrontDiagnostic;
                        report.interiorFrontNodes=recovered.interiorFrontNodes;
                        report.interiorFrontCandidates=recovered.interiorFrontCandidates;
                        report.interiorFrontWork=recovered.interiorFrontWork;
                        report.interiorFrontQueries=recovered.interiorFrontQueries;
                        report.interiorFrontRemainingOversized=recovered.interiorFrontRemainingOversized;
                        report.interiorFrontActivations=recovered.interiorFrontActivations;
                        report.interiorVoidActivations=recovered.interiorVoidActivations;
                        report.interiorFrontAdvances=recovered.interiorFrontAdvances;
                        report.interiorFrontSeconds=recovered.interiorFrontSeconds;
                        if(recovered.published){report.route=spatialSizing?"native_spatial_sizing":!spaceSurface?
                                (recovered.liveMaterialPatches?"native_plc_live_material":"native_plc_cavity"):
                                surfaceMode==OcctFrontSurface::Curvature?"native_curvature_surface":
                                surfaceMode==OcctFrontSurface::RegularQuality?"native_frontal_quality":
                                surfaceMode==OcctFrontSurface::Regular?"native_frontal_surface":"native_spaced_surface";
                            report.localDelaunaySupportNodes=recovered.localDelaunaySupportNodes;
                            report.localDelaunayRefinementNodes=recovered.localDelaunayRefinementNodes;
                            report.localDelaunayQualityEdits=recovered.localDelaunayQualityEdits;
                            originalFacetCentersCertified=recovered.originalFacetCentersCertified;
                            report.nativeSurfaceReason="none";
                            report.meshVolume=recovered.meshVolume;report.surfaceChordErrorLimit=recovered.chord;
                            report.surfaceTriangles=recovered.surfaceTriangles;
                            report.surfaceDeflection=spaceSurface?0.:recovered.chord*.4;
                            report.maximumBoundaryNodeError=recovered.maximumNodeError;report.maximumBoundarySampleError=recovered.maximumFacetError;
                            report.nativeCentroidThreads=recovered.nativeValidation.threadsUsed;report.nativeCentroidsChecked=recovered.nativeValidation.nativeQueries+recovered.nativeValidation.strictAdmissionReuses+recovered.nativeValidation.trimmedRayAdmissions;
                            report.nativeFinalRayAdmissions=recovered.nativeValidation.trimmedRayAdmissions;
                            report.nativeCentroidExactReuses=recovered.nativeValidation.strictAdmissionReuses;
                            report.nativeCentroidReferenceQueries=recovered.nativeValidation.originalReferenceQueries;report.nativeCentroidGeometryPairs=recovered.nativeValidation.copiedGeometryValuePairs;report.nativeCentroidSourceGuard=recovered.nativeValidation.sourceModelUnmodified;
                            report.centroidOutsideDistanceLimit=o.allowCentroidsWithinChordTolerance?recovered.chord:0;
                            report.centroidsWithinChordTolerance=recovered.nativeValidation.pointsWithinDistanceBand;
                            report.nativeCentroidDistanceQueries=recovered.nativeValidation.nativeDistanceQueries;
                            report.maximumOutsideCentroidDistance=recovered.nativeValidation.maximumOutsideDistance;
                            report.nativeChordSamplesChecked=recovered.chordSampleValidation.nativeQueries;
                            report.chordSamplesWithinTolerance=recovered.chordSampleValidation.pointsWithinDistanceBand;
                            report.maximumOutsideChordSampleDistance=recovered.chordSampleValidation.maximumOutsideDistance;
                            report.boundaryEmbeddingCertified=true;
                        }else if(spatialSizing||recovered.reason==NativeCadMeshReason::BudgetExceeded)throw Stop{recovered.reason};
                        return recovered.published;
                    };
                    PhaseTimer fallback(out.generationSeconds);
                    if(spatialSizing) {
                        const auto remaining=std::min(o.geometry.maximumMeshingWork,o.maximumValidationWork-work.count);
                        nativeSurface=recoverFront(remaining,0,OcctFrontSurface::SpatialSizing,
                            liveRecovery&&!o.allowCentroidsWithinChordTolerance);
                    } else {
                    nativeSurface=nativeStarFallback(model,view,o,body,report,work);
                    const bool improveRegularQuality=o.regularMeshQualityEnabled();
                    const bool regularFront=o.useRegularSurfaceFront||improveRegularQuality;
                    if(!nativeSurface&&(o.improveNativeSurfaceSpacing||regularFront)&&solids.size()==1&&model.faceUses(solid->tag()).size()<=32){
                        const auto trySurface=[&](OcctFrontSurface mode,const char* name){
                            const auto remaining=std::min(o.geometry.maximumMeshingWork,o.maximumValidationWork-work.count);
                            const auto trialWork=std::min<std::uint64_t>(10'000'000,remaining/5);
                            const auto constraintWork=std::min<std::uint64_t>(4'000'000,trialWork/2);
                            if(!constraintWork)return false;
                            detail::OcctStageTrace trial(name,solid->tag());
                            const bool accepted=recoverFront(trialWork,constraintWork,mode);
                            trial.finish(accepted?"accepted":"fallback");
                            return accepted;
                        };
                        if(improveRegularQuality)
                            nativeSurface=trySurface(OcctFrontSurface::RegularQuality,"regular_quality_trial");
                        if(!nativeSurface&&regularFront)
                            nativeSurface=trySurface(OcctFrontSurface::Regular,"regular_surface_trial");
                        if(!nativeSurface)
                            nativeSurface=trySurface(OcctFrontSurface::Spaced,"spaced_surface_trial");
                    }
                    const bool splineTrial=!nativeSurface&&o.useAdaptiveTetRecovery&&solids.size()==1&&
                        o.useNativePartitionFallback&&o.geometry.minimumElementQuality==0&&
                        std::string_view(report.nativeVolumeReason)=="no-strict-star-kernel-found"&&hasSplineRichBoundary(model,view);
                    if(!nativeSurface&&o.curvaturePointsPerTurn&&solids.size()==1&&
                       model.faceUses(solid->tag()).size()<=128&&o.recoveryMode==NativeCadRecoveryMode::Advanced){
                        const auto remaining=std::min(o.geometry.maximumMeshingWork,o.maximumValidationWork-work.count);
                        const auto trialWork=std::min<std::uint64_t>(80'000'000,remaining/3);
                        const auto constraintWork=std::min<std::uint64_t>(45'000'000,trialWork*2/3);
                        if(constraintWork){
                            detail::OcctStageTrace trial("curvature_surface_trial",solid->tag());
                            const bool live=liveRecovery&&!o.allowCentroidsWithinChordTolerance;
                            nativeSurface=recoverFront(trialWork,constraintWork,OcctFrontSurface::Curvature,live);
                            trial.finish(nativeSurface?"accepted":"fallback");
                        }
                    }
                    if(!nativeSurface&&liveRecovery&&solids.size()==1&&
                       o.recoveryMode==NativeCadRecoveryMode::Advanced&&!o.allowCentroidsWithinChordTolerance){
                        const auto remaining=std::min(o.geometry.maximumMeshingWork,o.maximumValidationWork-work.count);
                        const auto trialWork=std::min<std::uint64_t>(80'000'000,remaining/3);
                        const auto constraintWork=std::min<std::uint64_t>(45'000'000,trialWork*2/3);
                        if(constraintWork){
                            detail::OcctStageTrace trial("live_material_trial",solid->tag());
                            nativeSurface=recoverFront(trialWork,constraintWork,OcctFrontSurface::Standard,true);
                            trial.finish(nativeSurface?"accepted":"fallback");
                        }
                    }
                    if(!nativeSurface&&splineTrial)
                        nativeSurface=nativePartitionFallback(model,view,o,body,report,work,NativePartitionPass::ChordTrial);
                    // Keep the established successful coarse partition before
                    // new sweep recognition, so earlier meshes retain their
                    // coordinates/connectivity. Failed work remains charged.
                    if(!nativeSurface&&o.useNativePartitionFallback)
                        nativeSurface=nativePartitionFallback(model,view,o,body,report,work,NativePartitionPass::Midplanes);
                    if(!nativeSurface&&(o.useNativePartitionFallback||o.useNativeKernelAdmission)){
                        detail::OcctNativeSweepOptions sweep;
                        sweep.maximumWork=std::min(o.geometry.maximumMeshingWork,o.maximumValidationWork-work.count);
                        sweep.maximumLayerThickness=o.geometry.targetSize;
                        sweep.maximumNodes=o.geometry.maximumNodes;sweep.maximumCells=o.geometry.maximumCells;
                        sweep.maximumOutputBytes=o.geometry.maximumOutputBytes;sweep.maximumWorkingBytes=o.maximumPrimitiveWorkingBytes;
                        sweep.maximumSnapshotBytes=std::min<std::uint64_t>(64ULL*1024*1024,o.maximumPrimitiveWorkingBytes/8);
                        sweep.maximumRelativeVolumeError=o.maximumRelativeVolumeError;
                        sweep.minimumQuality=o.geometry.minimumElementQuality;sweep.warningQuality=o.geometry.warningElementQuality;
                        sweep.cancellation=o.geometry.cancellation;sweep.auditLimits.maximumTemporaryBytes=o.geometry.maximumAuditTemporaryBytes;
                        sweep.auditLimits.maximumFaceUses=o.geometry.maximumAuditFaceUses;sweep.auditLimits.cancellation=o.geometry.cancellation;
                        sweep.boundary.linearDeflection=report.surfaceDeflection;
                        sweep.boundary.maximumChordError=report.surfaceChordErrorLimit;
                        sweep.boundary.maximumNodes=o.geometry.maximumNodes;sweep.boundary.maximumTriangles=o.geometry.maximumCells/2;
                        sweep.boundary.maximumOutputBytes=o.geometry.maximumOutputBytes;sweep.boundary.maximumWork=sweep.maximumWork;
                        sweep.boundary.cancellation=o.geometry.cancellation;
                        detail::OcctStageTrace sweepTrace("cap_sweep",solid->tag());
                        const auto swept=detail::buildOcctNativeSweepVolume(model,solid->tag(),sweep,body);
                        sweepTrace.finish(swept.diagnostic);
                        work.tick(swept.workUsed);
                        if(swept.reason==detail::OcctNativeSweepReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled};
                        if(swept.reason==detail::OcctNativeSweepReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed};
                        if(swept.reason==detail::OcctNativeSweepReason::BudgetExceeded)throw Stop{NativeCadMeshReason::BudgetExceeded};
                        if(swept.reason==detail::OcctNativeSweepReason::QualityRequirement)throw Stop{NativeCadMeshReason::QualityRequirementNotMet};
                        if(swept.reason==detail::OcctNativeSweepReason::SourceMutation)throw Stop{NativeCadMeshReason::InternalError};
                        if(swept.published){nativeSurface=true;report.route="native_cap_sweep";
                            report.nativeVolumeReason=swept.diagnostic;report.cadVolume=swept.cadVolume;report.meshVolume=swept.meshVolume;
                            report.startFaceTag=swept.baseCapTag;report.endFaceTag=swept.targetCapTag;
                            report.maximumBoundaryNodeError=swept.maximumBoundaryNodeError;report.maximumBoundarySampleError=swept.maximumBoundarySampleError;
                            report.nativeSurfaceWork+=swept.surfaceWork;report.nativeVolumeWork+=swept.workUsed-swept.surfaceWork;
                            report.nativeSurfaceSeconds+=swept.surfaceSeconds;report.nativeVolumeSeconds+=swept.totalSeconds-swept.surfaceSeconds;
                            report.boundaryEmbeddingCertified=swept.exteriorCertificate.certified;
                        }
                        if(!nativeSurface){
                            detail::OcctPairedShellOptions paired;
                            static_cast<detail::OcctNativeSweepOptions&>(paired)=sweep;
                            paired.maximumWork=std::min(o.geometry.maximumMeshingWork,o.maximumValidationWork-work.count);
                            paired.boundary.maximumWork=paired.maximumWork;
                            paired.boundary.linearDeflection=paired.boundary.maximumChordError*.4;
                            paired.maximumNativeValidationThreads=o.maximumNativeValidationThreads;
                            detail::OcctStageTrace shellTrace("paired_shell",solid->tag());
                            const auto shell=detail::buildOcctPairedShellVolume(model,solid->tag(),paired,body);
                            shellTrace.finish(shell.diagnostic);
                            work.tick(shell.workUsed);
                            if(shell.reason==detail::OcctPairedShellReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled};
                            if(shell.reason==detail::OcctPairedShellReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed};
                            if(shell.reason==detail::OcctPairedShellReason::BudgetExceeded)throw Stop{NativeCadMeshReason::BudgetExceeded};
                            if(shell.reason==detail::OcctPairedShellReason::QualityRequirement)throw Stop{NativeCadMeshReason::QualityRequirementNotMet};
                            if(shell.reason==detail::OcctPairedShellReason::SourceMutation)throw Stop{NativeCadMeshReason::InternalError};
                            if(shell.published){nativeSurface=true;report.route="native_paired_shell";
                                report.nativeVolumeReason=shell.diagnostic;report.cadVolume=shell.cadVolume;report.meshVolume=shell.meshVolume;
                                report.surfaceDeflection=paired.boundary.linearDeflection;
                                report.surfaceTriangles=shell.cells-shell.quality.volumeCells;
                                report.startFaceTag=shell.baseCapTag;report.endFaceTag=shell.targetCapTag;
                                report.maximumBoundaryNodeError=shell.maximumBoundaryNodeError;report.maximumBoundarySampleError=shell.maximumBoundarySampleError;
                                report.nativeSurfaceWork+=shell.surfaceWork;report.nativeVolumeWork+=shell.workUsed-shell.surfaceWork;
                                report.nativeSurfaceSeconds+=shell.surfaceSeconds;report.nativeVolumeSeconds+=shell.totalSeconds-shell.surfaceSeconds;
                                report.boundaryEmbeddingCertified=shell.exteriorCertificate.certified;
                            }
                        }
                        if(!nativeSurface&&o.useLongitudinalRecovery)
                            nativeSurface=nativeLongitudinalFallback(model,view,o,body,report,work);
                        // When a single body's coarse partition still lacks a star kernel,
                        // try direct constraint recovery before the costly partition series.
                        // Charge its work even if the existing strategies need to continue.
                        const bool coarseKernelRefusal=!report.partitionAttemptReports.empty()&&
                            std::string_view(report.partitionAttemptReports.back().diagnostic)=="no-strict-star-kernel-found";
                        const bool preferFixedGrid=!nativeSurface&&o.useAdaptiveTetRecovery&&solids.size()==1&&
                            o.useNativePartitionFallback&&o.geometry.minimumElementQuality==0&&coarseKernelRefusal&&
                            hasSplineRichBoundary(model,view);
                        if(!nativeSurface&&!preferFixedGrid&&o.useNativePartitionFallback&&solids.size()==1&&
                           o.recoveryMode==NativeCadRecoveryMode::Advanced&&
                           o.geometry.minimumElementQuality==0&&coarseKernelRefusal){
                            const auto remaining=std::min(o.geometry.maximumMeshingWork,o.maximumValidationWork-work.count);
                            const auto trialWork=std::min<std::uint64_t>(70'000'000,remaining/4);
                            const auto constraintWork=std::min<std::uint64_t>(45'000'000,trialWork*2/3);
                            if(constraintWork)nativeSurface=recoverFront(trialWork,constraintWork);
                        }
                        if(!nativeSurface&&o.useNativePartitionFallback){
                            try{nativeSurface=nativePartitionFallback(model,view,o,body,report,work,NativePartitionPass::Full,preferFixedGrid);}
                            catch(const Stop& failure){
                                if(failure.reason!=NativeCadMeshReason::BudgetExceeded||work.count>=o.maximumValidationWork)throw;
                            }
                        }
                    }
                    if(!nativeSurface&&(o.useNativePartitionFallback||o.useNativeKernelAdmission||
                                       o.recoveryMode==NativeCadRecoveryMode::Advanced)){
                        nativeSurface=recoverFront(std::min(o.geometry.maximumMeshingWork,o.maximumValidationWork-work.count),0);
                    }
                    }
                }
                need(report.planar.published||nativeSurface,recognitionReason==NativeCadMeshReason::None?NativeCadMeshReason::MeshingFailed:recognitionReason);
                if(!nativeSurface){
                    if(report.planar.qualityAssessed)report.meshVolume=report.planar.quality.totalVolume;
                    else{const auto quality=algo::analyzeMixedMeshQuality(body);need(!quality.allocationFailed,NativeCadMeshReason::AllocationFailed);
                        need(quality.valid,NativeCadMeshReason::MeshingFailed);report.meshVolume=quality.totalVolume;}
                }
            } else {
                need(controls.empty()&&!o.geometry.spatialField&&o.geometry.strategy==algo::GeometryStrategy::Auto&&
                     o.geometry.maximumElementSizeRatio==0&&o.geometry.warningElementSizeRatio==0,
                     NativeCadMeshReason::UnsupportedControl);
                report.route=algo::primitiveKindName(r.shape.kind);
                if(r.faces.contains(2))report.startFaceTag=r.faces.at(2);if(r.faces.contains(3))report.endFaceTag=r.faces.at(3);
                algo::PrimitiveMeshingOptions po;po.targetSize=o.geometry.targetSize;po.elements=o.elements;
                po.maximumChordError=o.maximumChordError;po.maximumRelativeVolumeError=o.maximumRelativeVolumeError;
                po.minimumQuality=o.geometry.minimumElementQuality;po.threads=o.geometry.threads;
                po.warningQuality=o.geometry.warningElementQuality;
                po.maximumNodes=std::min(o.geometry.maximumNodes,o.geometry.maximumMeshingWork);
                po.maximumCells=std::min(o.geometry.maximumCells,o.geometry.maximumMeshingWork);
                po.maximumOutputBytes=o.geometry.maximumOutputBytes;
                po.maximumWorkingBytes=o.maximumPrimitiveWorkingBytes;
                po.maximumProfileSlabWork=o.geometry.maximumMeshingWork;
                po.profileLimits.maximumBytes=o.geometry.maximumProfileWorkingBytes;
                po.auditLimits.maximumTemporaryBytes=o.geometry.maximumAuditTemporaryBytes;po.auditLimits.maximumFaceUses=o.geometry.maximumAuditFaceUses;
                po.cancellation=o.geometry.cancellation;po.auditLimits.cancellation=o.geometry.cancellation;
                sizing::CappedSizeOutcome profile;
                if(!faceControls.empty()) {
                    PhaseTimer planning(out.preparationSeconds);
                    need(r.shape.kind==algo::PrimitiveKind::Cylinder||r.shape.kind==algo::PrimitiveKind::Tube,NativeCadMeshReason::UnsupportedControl);
                    sizing::CappedSizeSpecification spec;spec.length=r.shape.height;spec.globalSize=po.targetSize;spec.minimumSize=o.minimumSize;spec.cancellation=po.cancellation;
                    spec.maximumBytes=o.geometry.maximumProfileWorkingBytes;spec.maximumWork=o.geometry.maximumMeshingWork;
                    for(const auto& control:faceControls) {
                        need(control.faceTag==report.startFaceTag||control.faceTag==report.endFaceTag,NativeCadMeshReason::UnsupportedControl);
                        const double z=control.faceTag==report.startFaceTag?0:r.shape.height;
                        spec.controls.push_back({z,z,control.size,control.growth});
                    }
                    if(!spec.controls.empty()) {profile=sizing::compileCappedSizeProfile(spec);
                        if(profile.reason==sizing::ProfileReason::BudgetExceeded)throw Stop{NativeCadMeshReason::BudgetExceeded};
                        if(profile.reason==sizing::ProfileReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled};
                        if(profile.reason==sizing::ProfileReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed};
                        need(profile.compiled,NativeCadMeshReason::UnsupportedControl);
                        po.sizeProfile=&profile.profile;po.profileStrategy=r.shape.kind==algo::PrimitiveKind::Cylinder?algo::PrimitiveProfileStrategy::Volume:algo::PrimitiveProfileStrategy::Sweep;}
                }
                const bool optimize=o.optimizeProfileInterior&&po.profileStrategy==algo::PrimitiveProfileStrategy::Volume;
                if(optimize)po.minimumQuality=0; // Strict acceptance follows bounded improvement.
                const auto generating=Clock::now();report.primitive=algo::meshPrimitive(r.shape,po,body);
                out.generationSeconds+=elapsed(generating);
                if(!report.primitive.published) {
                    if(report.primitive.reason==algo::PrimitiveReason::Cancelled)throw Stop{NativeCadMeshReason::Cancelled};
                    if(report.primitive.reason==algo::PrimitiveReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed};
                    if(report.primitive.reason==algo::PrimitiveReason::QualityRequirement)throw Stop{NativeCadMeshReason::QualityRequirementNotMet};
                    if(report.primitive.reason==algo::PrimitiveReason::GeometryTolerance)throw Stop{NativeCadMeshReason::GeometryError};
                    if(report.primitive.reason==algo::PrimitiveReason::OutputBudget||report.primitive.reason==algo::PrimitiveReason::WorkingBudget)throw Stop{NativeCadMeshReason::BudgetExceeded};
                    throw Stop{NativeCadMeshReason::MeshingFailed};
                }
                report.meshVolume=report.primitive.meshVolume;
                if(optimize) {
                    PhaseTimer improving(out.qualitySeconds);report.optimizationAttempted=true;
                    algo::MixedMeshOptimizationRequest request;
                    request.maximumIterations=4;
                    request.maximumNodeTrials=std::min(o.maximumQualityNodeTrials,o.maximumValidationWork-work.count);
                    report.optimization=algo::optimizeMixedMeshQuality(body,request);
                    work.tick(report.optimization.nodeTrials);
                    if(report.optimization.reason==algo::MixedMeshOptimizationReason::AllocationFailed)throw Stop{NativeCadMeshReason::AllocationFailed};
                    if(report.optimization.reason==algo::MixedMeshOptimizationReason::BudgetExceeded)throw Stop{NativeCadMeshReason::BudgetExceeded};
                    need(report.optimization.completed,NativeCadMeshReason::MeshingFailed);
                    need(report.optimization.after.minimumQuality>=o.geometry.minimumElementQuality,NativeCadMeshReason::QualityRequirementNotMet);
                    report.meshVolume=report.optimization.after.totalVolume;
                    for(Index i=0;i<body.nodeCount();++i)if(body.nodeClass(i).dim==3){work.tick();need(solid->isInside(body.node(i)),NativeCadMeshReason::GeometryError);}
                }
            }
            if(boxContactBodies.contains(solid->tag())){
                PhaseTimer evaluating(out.validationSeconds);
                const auto evaluated=detail::evaluateOcctBoxContactCorners(model,solid->tag(),
                    boxContactOptions(o,work,cachedBoxCertificates),body);
                work.tick(evaluated.workUsed);report.boxContactWork=evaluated.workUsed;
                checkBoxContactFailure(evaluated.reason);
                need(evaluated.published&&evaluated.sourceImmutabilityVerified,NativeCadMeshReason::UnresolvedBodyInterface);
                report.boxContactVerticesEvaluated=evaluated.evaluatedVertices;
                report.boxContactMaximumEndpointCorrection=evaluated.maximumEndpointCorrection;
                report.meshVolume=evaluated.meshVolume;
                report.maximumBoundaryNodeError=std::max(report.maximumBoundaryNodeError,evaluated.maximumBoundaryNodeError);
                report.maximumBoundarySampleError=std::max(report.maximumBoundarySampleError,evaluated.maximumBoundarySampleError);
            }
            const auto nodes=static_cast<std::uint64_t>(body.nodeCount()), cells=static_cast<std::uint64_t>(body.cellCount());
            need(nodes<=o.geometry.maximumMeshingWork&&cells<=o.geometry.maximumMeshingWork-nodes,NativeCadMeshReason::BudgetExceeded);
            report.relativeVolumeError=std::abs(report.meshVolume-report.cadVolume)/report.cadVolume;
            need(report.relativeVolumeError<=o.maximumRelativeVolumeError,NativeCadMeshReason::GeometryError);
            PhaseTimer validating(out.validationSeconds);detail::OcctStageTrace validation("cad_validation",solid->tag());const auto bodyBegin=local.cellCount();append(model,view,r,body,local,report,work,nativeSurface,originalFacetCentersCertified);
            bodyRanges.push_back({solid->tag(),bodyBegin,local.cellCount()});
            const double chord=nativeSurface?report.surfaceChordErrorLimit:
                r.analytic?report.primitive.chordErrorBound:0;
            need(report.maximumBoundarySampleError<=std::max(chord*1.01,1e-7),NativeCadMeshReason::GeometryError);
            validating.finish();report.reason=NativeCadMeshReason::None;
        }
        {PhaseTimer contacts(out.validationSeconds);detail::OcctStageTrace contactValidation("contact_validation");
            for(const auto [a,b]:separateContacts){
                const auto ar=std::find_if(bodyRanges.begin(),bodyRanges.end(),[a](auto r){return r.tag==a;});
                const auto br=std::find_if(bodyRanges.begin(),bodyRanges.end(),[b](auto r){return r.tag==b;});
                need(ar!=bodyRanges.end()&&br!=bodyRanges.end(),NativeCadMeshReason::InternalError);
                const auto at=contactTree(local,*ar,work),bt=contactTree(local,*br,work);
                const auto contact=boxFaceContacts.find({a,b});
                if(contact==boxFaceContacts.end())auditSeparateContact(at,0,bt,0,work,out.contactTrianglePairs);
                else {
                    need(at.triangles.size()+bt.triangles.size()<=o.maximumPrimitiveWorkingBytes/4,NativeCadMeshReason::BudgetExceeded);
                    std::vector<bool> matchedA(at.triangles.size()),matchedB(bt.triangles.size());
                    double area=0;std::uint64_t matched=0;
                    auditSeparateContact(at,0,bt,0,work,out.contactTrianglePairs,&matchedA,&matchedB,&area,&matched);
                    need(matched>0&&finite(area)&&std::abs(area-contact->second)<=contact->second*1e-10,NativeCadMeshReason::UnresolvedBodyInterface);
                    ++out.separateBoxFaceContactPairs;out.matchedBoxFaceTriangles+=matched;
                }
                ++out.separateContactPairs;
            }
        }
        mesh::LinearVolumeAuditLimits limits;limits.maximumTemporaryBytes=o.geometry.maximumAuditTemporaryBytes;limits.maximumFaceUses=o.geometry.maximumAuditFaceUses;limits.cancellation=o.geometry.cancellation;
        const auto auditing=Clock::now();out.audit=mesh::auditLinearVolumeBoundary(local,limits);out.auditSeconds=elapsed(auditing);
        need(!out.audit.cancelled,NativeCadMeshReason::Cancelled);need(!out.audit.budgetExceeded,NativeCadMeshReason::BudgetExceeded);need(!out.audit.allocationFailed,NativeCadMeshReason::AllocationFailed);
        need(out.audit.valid,NativeCadMeshReason::MeshingFailed);
        const auto assessing=Clock::now();out.quality=algo::analyzeMixedMeshQuality(local);out.qualitySeconds+=elapsed(assessing);
        need(!out.quality.allocationFailed,NativeCadMeshReason::AllocationFailed);need(out.quality.valid,NativeCadMeshReason::MeshingFailed);
        need(out.quality.minimumQuality>=o.geometry.minimumElementQuality,NativeCadMeshReason::QualityRequirementNotMet);
        out.qualityWarning=out.quality.minimumQuality<o.geometry.warningElementQuality;
        work.tick(); out.validationWork=work.count;destination.swap(local);out.published=true;out.destinationUnchanged=false;out.reason=NativeCadMeshReason::None;
    } catch(const Stop& s){out.reason=s.reason;}catch(const Standard_OutOfMemory&){out.reason=NativeCadMeshReason::AllocationFailed;}catch(const std::bad_alloc&){out.reason=NativeCadMeshReason::AllocationFailed;}catch(...){out.reason=NativeCadMeshReason::InternalError;}
    if(!out.published&&!out.bodies.empty())out.bodies.back().reason=out.reason;
    if(queryContext){const auto queries=queryContext->counters();
        out.nativeProjectionQueries=queries.surfaceProjectionQueries;
        out.nativeProjectionInitializations=queries.surfaceProjectionInitializations;
        out.nativeProjectionReuses=queries.surfaceProjectionReuses;
        out.nativeTrimmedProjectionQueries=queries.trimmedProjectionQueries;
        out.nativeTrimmedProjectionMemoHits=queries.trimmedProjectionMemoHits;}
    out.validationWork=work.count;out.totalSeconds=elapsed(started);return out;
}
const char* nativeCadMeshReasonName(NativeCadMeshReason r) noexcept {
    constexpr const char* names[]={"none","invalid_request","unsupported_body","unresolved_body_interface","unsupported_control","recognition_failed","meshing_failed","geometry_error","budget_exceeded","cancelled","allocation_failed","internal_error","quality_requirement_not_met","size_requirement_not_met"};
    const auto i=static_cast<std::size_t>(r);return i<std::size(names)?names[i]:"unknown";
}
} // namespace paacs::geo
