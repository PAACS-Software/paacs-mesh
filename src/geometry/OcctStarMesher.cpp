// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Optional native-adapter volume fallback. CAD authority stays with the caller.
#include "paacs/bnd/TetFacetGeometry.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "OcctStarMesher.hpp"
#include "OcctModelInternal.hpp"
#include "OcctSolidClassifier.hpp"
#include <Standard_Failure.hxx>
#include <Standard_OutOfMemory.hxx>
#include <algorithm>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <cmath>
#include <limits>
namespace paacs::geo::detail {
namespace {
struct Stop { OcctStarReason reason; const char* diagnostic; };
void need(bool yes,const char* diagnostic, OcctStarReason reason=OcctStarReason::InvalidBoundary) {if(!yes)throw Stop{reason,diagnostic};}
struct Point {Vec3 p;mesh::Classification cls;};
struct Triangle {std::array<std::uint32_t,3> ids;int tag;Vec3 normal,lo,hi;};

bool finite(Vec3 p) {return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);}
using Edge=std::array<std::uint32_t,2>;
struct Use {std::size_t triangle;int sense;};
double coordinate(Vec3 p,int axis) {return axis==0?p.x:axis==1?p.y:p.z;}
struct Budget {
    const OcctStarOptions& options;
    std::uint64_t& work;
    std::uint64_t& pairs;
    void tick(std::uint64_t count=1) {
        need(!options.cancellation||!options.cancellation->load(std::memory_order_relaxed),"cancelled",OcctStarReason::Cancelled);
        need(count<=options.maximumWork-work,"work-budget",OcctStarReason::WorkBudget);work+=count;
    }
};
struct PreparedBoundary {
    std::vector<Point> points;
    std::vector<bnd::TetFacetVertex> vertices;
    std::vector<Triangle> faces;
    Vec3 lo{},hi{},mean{},centroid{};
    double scale=0,volume=0;
    std::uint64_t shells=0;
    std::int64_t euler=0;
};
struct Box {Vec3 lo{},hi{};};
bool overlaps(Box a,Box b) {
    return !(a.lo.x>b.hi.x||b.lo.x>a.hi.x||a.lo.y>b.hi.y||b.lo.y>a.hi.y||a.lo.z>b.hi.z||b.lo.z>a.hi.z);
}
Box united(Box a,Box b) {
    return {{std::min(a.lo.x,b.lo.x),std::min(a.lo.y,b.lo.y),std::min(a.lo.z,b.lo.z)},
            {std::max(a.hi.x,b.hi.x),std::max(a.hi.y,b.hi.y),std::max(a.hi.z,b.hi.z)}};
}
struct BvhNode {Box box;std::size_t begin=0,end=0,left=0,right=0;bool leaf=true;};
void certifyEmbedding(const std::vector<Point>& points,const std::vector<bnd::TetFacetVertex>& vertices,
                      const std::vector<Triangle>& faces,Budget& budget) {
    constexpr std::array<std::array<int,2>,3> triEdges{{{0,1},{1,2},{2,0}}};
    constexpr std::array<std::array<int,3>,1> triFaces{{{0,1,2}}};
    const bnd::TetFacetComplexView view{vertices,{}, {}, {}, {}, {3,0}};
    std::vector<std::size_t> order(faces.size());std::iota(order.begin(),order.end(),0);
    std::vector<BvhNode> tree;tree.reserve(faces.size()*2);
    const auto build=[&](auto&& self,std::size_t begin,std::size_t end)->std::size_t {
        budget.tick();const auto id=tree.size();tree.push_back({});
        Box box{faces[order[begin]].lo,faces[order[begin]].hi};
        for(auto i=begin+1;i<end;++i){budget.tick();box=united(box,{faces[order[i]].lo,faces[order[i]].hi});}
        tree[id]={box,begin,end,0,0,true};
        if(end-begin>8) {
            const auto extent=box.hi-box.lo;int axis=0;if(extent.y>extent.x)axis=1;if(extent.z>coordinate(extent,axis))axis=2;
            const auto middle=begin+(end-begin)/2;
            std::nth_element(order.begin()+begin,order.begin()+middle,order.begin()+end,[&](auto a,auto b){
                const auto ca=coordinate(faces[a].lo*.5+faces[a].hi*.5,axis);
                const auto cb=coordinate(faces[b].lo*.5+faces[b].hi*.5,axis);
                return std::pair{ca,a}<std::pair{cb,b};
            });
            const auto left=self(self,begin,middle),right=self(self,middle,end);
            tree[id].left=left;tree[id].right=right;tree[id].leaf=false;
        }
        return id;
    };
    build(build,0,order.size());
    const auto pair=[&](std::size_t ai,std::size_t bi) {
        budget.tick();const auto& a=faces[ai];const auto& b=faces[bi];
        if(!overlaps({a.lo,a.hi},{b.lo,b.hi}))return;
        ++budget.pairs;
        const auto contact=bnd::tet_facet_detail::auditIdentityAwareClosedContact(view,a.ids,triEdges,triFaces,b.ids,triEdges,triFaces);
        const std::array<Vec3,3> ap{points[a.ids[0]].p,points[a.ids[1]].p,points[a.ids[2]].p},
                                 bp{points[b.ids[0]].p,points[b.ids[1]].p,points[b.ids[2]].p};
        need(!contact.illegal()&&!bnd::tet_facet_detail::trianglesOpenIntersect(ap,bp),
             "boundary-geometric-contact",OcctStarReason::BoundaryContact);
    };
    const auto visit=[&](auto&& self,std::size_t ai,std::size_t bi)->void {
        budget.tick();const auto a=tree[ai],b=tree[bi];
        if(!overlaps(a.box,b.box))return;
        if(ai==bi) {
            if(a.leaf)for(auto i=a.begin;i<a.end;++i)for(auto j=i+1;j<a.end;++j)pair(order[i],order[j]);
            else {self(self,a.left,a.left);self(self,a.left,a.right);self(self,a.right,a.right);}
        }else if(a.leaf&&b.leaf) {
            for(auto i=a.begin;i<a.end;++i)for(auto j=b.begin;j<b.end;++j)pair(order[i],order[j]);
        }else if(!a.leaf&&(b.leaf||a.end-a.begin>=b.end-b.begin)) {
            self(self,a.left,bi);self(self,a.right,bi);
        }else {self(self,ai,b.left);self(self,ai,b.right);}
    };
    visit(visit,0,0);
}
PreparedBoundary prepareBoundary(const mesh::Mesh& boundary,const OcctStarOptions& options,
                                 Budget& budget,OcctBoundaryCertificate& certificate,
                                 bool requireStarTopology=false) {
        const std::uint64_t nodes=static_cast<std::uint64_t>(boundary.nodeCount());
        const std::uint64_t triangles=static_cast<std::uint64_t>(boundary.cellCount());
        certificate.nodes=nodes;certificate.triangles=triangles;
        const auto tick=[&](std::uint64_t count=1){budget.tick(count);};
        need(nodes>=4&&triangles>=4&&nodes<=UINT32_MAX,"input-count");
        need(nodes<=options.maximumBoundaryNodes&&triangles<=options.maximumBoundaryTriangles,
             "boundary-count-budget",OcctStarReason::WorkBudget);
        need(nodes<=UINT64_MAX-triangles&&(nodes+triangles)<=options.maximumWorkingBytes/1024,
             "working-payload-admission",OcctStarReason::WorkBudget);
        tick(0);
        std::vector<Point> points(nodes);std::vector<bnd::TetFacetVertex> vertices;
        vertices.reserve(nodes);Vec3 lo{1e300,1e300,1e300},hi{-1e300,-1e300,-1e300},mean{};
        for(std::size_t i=0;i<nodes;++i) {
            auto& p=points[i];p.p=boundary.node(static_cast<Index>(i));p.cls=boundary.nodeClass(static_cast<Index>(i));need(finite(p.p)&&
                p.cls.dim>=0&&p.cls.dim<=2&&p.cls.tag>=0,"point-payload");tick();
            vertices.push_back({p.p,i});mean=mean+p.p;
            lo={std::min(lo.x,p.p.x),std::min(lo.y,p.p.y),std::min(lo.z,p.p.z)};
            hi={std::max(hi.x,p.p.x),std::max(hi.y,p.p.y),std::max(hi.z,p.p.z)};
        }
        const double scale=(hi-lo).norm();need(scale>0&&std::isfinite(scale),"coordinate-range");
        mean=mean/static_cast<double>(nodes);
        std::vector<Triangle> faces(triangles);std::vector<std::vector<std::size_t>> dual(triangles);
        std::map<Edge,std::vector<Use>> edges;std::set<std::array<std::uint32_t,3>> faceKeys;
        std::vector<std::vector<Edge>> links(nodes);std::vector<bool> used(nodes);
        long double volume=0;Vec3 centroid{};
        for(std::size_t i=0;i<triangles;++i) {
            auto& f=faces[i];const auto& cell=boundary.cell(static_cast<Index>(i));
            need(cell.type==mesh::CellType::Tri3&&cell.onEntity.dim==2&&cell.onEntity.tag>=0,"triangle-payload");
            const auto ids=boundary.cellNodes(static_cast<Index>(i));
            for(unsigned j=0;j<3;++j){need(ids[j]>=0&&static_cast<std::uint64_t>(ids[j])<nodes,"triangle-identity");f.ids[j]=static_cast<std::uint32_t>(ids[j]);}f.tag=cell.onEntity.tag;
            auto key=f.ids;std::sort(key.begin(),key.end());
            need(key[0]!=key[1]&&key[1]!=key[2]&&key[2]<nodes&&faceKeys.insert(key).second,"triangle-identity");tick();
            for(std::size_t e=0;e<3;++e) {
                tick();
                const auto a=f.ids[e],b=f.ids[(e+1)%3];const Edge k{std::min(a,b),std::max(a,b)};
                edges[k].push_back({i,a<b?1:-1});used[a]=true;
                links[a].push_back({f.ids[(e+1)%3],f.ids[(e+2)%3]});
            }
            const auto a=points[f.ids[0]].p,b=points[f.ids[1]].p,c=points[f.ids[2]].p;
            need(bnd::tet_facet_detail::projectionForTriangle(a,b,c)!=bnd::tet_facet_detail::Projection::Invalid,"degenerate-triangle");
            f.normal=(b-a).cross(c-a);const double norm=f.normal.norm();need(norm>0&&std::isfinite(norm),"normal-range");f.normal=f.normal/norm;
            f.lo={std::min({a.x,b.x,c.x}),std::min({a.y,b.y,c.y}),std::min({a.z,b.z,c.z})};
            f.hi={std::max({a.x,b.x,c.x}),std::max({a.y,b.y,c.y}),std::max({a.z,b.z,c.z})};
            const long double v=static_cast<long double>((a-lo).dot((b-lo).cross(c-lo)))/6;
            volume+=v;centroid=centroid+(lo+a+b+c)*static_cast<double>(v*.25L);
        }
        need(std::find(used.begin(),used.end(),false)==used.end(),"orphan-boundary-point");
        for(const auto& [edge,owners]:edges) {
            tick();need(owners.size()==2&&owners[0].sense!=owners[1].sense,"open-or-nonorientable-boundary",OcctStarReason::UnsupportedTopology);
            dual[owners[0].triangle].push_back(owners[1].triangle);dual[owners[1].triangle].push_back(owners[0].triangle);
        }
        const std::int64_t euler=static_cast<std::int64_t>(nodes)-static_cast<std::int64_t>(edges.size())+static_cast<std::int64_t>(triangles);certificate.eulerCharacteristic=euler;
        std::vector<bool> seen(triangles);std::uint64_t shells=0;
        for(std::size_t first=0;first<triangles;++first)if(!seen[first]) {
            ++shells;std::queue<std::size_t> pending;pending.push(first);seen[first]=true;
            while(!pending.empty()) {auto at=pending.front();pending.pop();tick();for(auto next:dual[at])if(!seen[next]){seen[next]=true;pending.push(next);}}
        }
        certificate.connectedShells=shells;
        for(const auto& link:links) {
            std::map<std::uint32_t,std::vector<std::uint32_t>> graph;
            for(const auto& e:link){tick();graph[e[0]].push_back(e[1]);graph[e[1]].push_back(e[0]);}
            for(auto& [id,next]:graph){std::sort(next.begin(),next.end());need(next.size()==2&&next[0]!=next[1],"pinched-boundary-vertex",OcctStarReason::UnsupportedTopology);}
            std::set<std::uint32_t> visited;std::queue<std::uint32_t> q;q.push(graph.begin()->first);visited.insert(q.front());
            while(!q.empty()){auto at=q.front();q.pop();tick();for(auto next:graph[at])if(visited.insert(next).second)q.push(next);}
            need(visited.size()==graph.size(),"pinched-boundary-vertex",OcctStarReason::UnsupportedTopology);
        }
        certificate.closedOrientedVertexManifold=true;
        const double plcVolume=static_cast<double>(volume);need(plcVolume>0&&std::isfinite(plcVolume),"nonpositive-boundary-volume");
        centroid=centroid/plcVolume;certificate.signedVolume=plcVolume;
        if(requireStarTopology) {
            need(shells==1,"disconnected-boundary",OcctStarReason::UnsupportedTopology);
            need(euler==2,"boundary-genus-not-zero",OcctStarReason::UnsupportedTopology);
        }
        certifyEmbedding(points,vertices,faces,budget);
        certificate.embedded=true;
        return {std::move(points),std::move(vertices),std::move(faces),lo,hi,mean,centroid,scale,plcVolume,shells,euler};
}

struct ProposalPlane { Vec3 normal; double offset; std::uint64_t order; };
// A bijective integer mixer fixes the constraint order without mutable state.
std::uint64_t proposalOrder(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}
// Numerical halfspace feasibility proposes a point only. The original exact
// orient3d tests still decide whether that point is a strict fan kernel.
// A violated 3D plane reduces the preceding constraints to its 2D plane;
// a violated 2D line reduces preceding constraints to a bounded 1D interval.
bool halfspaceProposal(const PreparedBoundary& prepared, double margin,
                       Vec3& candidate, Budget& budget, bool retryParallelCancellation=false) {
    budget.tick(0);
    need(prepared.faces.size() <= UINT64_MAX - 6, "proposal-plane-count-overflow", OcctStarReason::WorkBudget);
    const auto planeCount = static_cast<std::uint64_t>(prepared.faces.size()) + 6;
    // Faces were admitted against maximumBoundaryTriangles by prepareBoundary.
    // The six box constraints and this array have a separate payload admission,
    // in addition to the existing conservative boundary/container estimate.
    const auto boundaryBytes = (static_cast<std::uint64_t>(prepared.points.size()) +
                                static_cast<std::uint64_t>(prepared.faces.size())) * 1024;
    need(boundaryBytes <= budget.options.maximumWorkingBytes &&
         planeCount <= (budget.options.maximumWorkingBytes - boundaryBytes) / sizeof(ProposalPlane),
         "proposal-working-byte-budget", OcctStarReason::WorkBudget);
    std::vector<ProposalPlane> planes;
    need(planeCount <= planes.max_size(), "proposal-plane-count-budget", OcctStarReason::WorkBudget);
    planes.reserve(static_cast<std::size_t>(planeCount));
    const auto origin = prepared.lo * .5 + prepared.hi * .5;
    const auto scale = prepared.scale;
    const auto low = (prepared.lo - origin) / scale, high = (prepared.hi - origin) / scale;
    planes.push_back({{1, 0, 0}, high.x, 0}); planes.push_back({{-1, 0, 0}, -low.x, 1});
    planes.push_back({{0, 1, 0}, high.y, 2}); planes.push_back({{0, -1, 0}, -low.y, 3});
    planes.push_back({{0, 0, 1}, high.z, 4}); planes.push_back({{0, 0, -1}, -low.z, 5});
    for (std::size_t i = 0; i < prepared.faces.size(); ++i) {
        budget.tick(); const auto& face = prepared.faces[i];
        planes.push_back({face.normal, face.normal.dot((prepared.points[face.ids[0]].p - origin) / scale) - margin,
                          proposalOrder(static_cast<std::uint64_t>(i))});
    }
    std::sort(planes.begin() + 6, planes.end(), [&](const auto& a, const auto& b) {
        budget.tick(); return a.order < b.order;
    });
    Vec3 point{};
    for (std::size_t i = 0; i < planes.size(); ++i) {
        budget.tick(); const auto& plane = planes[i];
        const auto excess = plane.normal.dot(point) - plane.offset;
        if (!(excess > 0)) continue;
        const Vec3 axis = std::abs(plane.normal.x) <= std::abs(plane.normal.y) && std::abs(plane.normal.x) <= std::abs(plane.normal.z) ? Vec3{1, 0, 0} :
                          std::abs(plane.normal.y) <= std::abs(plane.normal.z) ? Vec3{0, 1, 0} : Vec3{0, 0, 1};
        auto u = plane.normal.cross(axis);
        const auto norm = u.norm(); if (!(norm > 0 && std::isfinite(norm))) return false;
        u = u / norm;
        const auto v = plane.normal.cross(u), base = point - plane.normal * excess;
        double s = 0, t = 0;
        for (std::size_t j = 0; j < i; ++j) {
            budget.tick(); const auto& previous = planes[j];
            const auto a = previous.normal.dot(u), b = previous.normal.dot(v), c = previous.offset - previous.normal.dot(base);
            if (!(a * s + b * t > c)) continue;
            const auto length = std::hypot(a, b); if (!(length > 0 && std::isfinite(length))) return false;
            const auto nx = a / length, ny = b / length, height = c / length;
            const auto sx = nx * height, sy = ny * height, dx = -ny, dy = nx;
            double lower = -std::numeric_limits<double>::infinity(), upper = std::numeric_limits<double>::infinity();
            for (std::size_t k = 0; k < j; ++k) {
                budget.tick(); const auto& prior = planes[k];
                const auto ra = prior.normal.dot(u), rb = prior.normal.dot(v), rc = prior.offset - prior.normal.dot(base);
                const auto coefficient = ra * dx + rb * dy, rhs = rc - ra * sx - rb * sy;
                // Proposal-only retry for cancellation in nearly parallel
                // projected constraints. Frozen exact acceptance follows.
                if(retryParallelCancellation&&std::abs(coefficient)<=32*std::numeric_limits<double>::epsilon()){
                    const auto uncertainty=64*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(rc),std::abs(ra*sx),std::abs(rb*sy)});
                    if(rhs < -uncertainty)return false;
                    continue;
                }
                if (coefficient == 0) { if (rhs < 0) return false; continue; }
                const auto bound = rhs / coefficient; if (!std::isfinite(bound)) return false;
                if (coefficient > 0) upper = std::min(upper, bound); else lower = std::max(lower, bound);
                if (lower > upper) return false;
            }
            const auto parameter = std::clamp(0., lower, upper); if (!std::isfinite(parameter)) return false;
            s = sx + dx * parameter; t = sy + dy * parameter;
            if (!std::isfinite(s) || !std::isfinite(t)) return false;
        }
        point = base + u * s + v * t; if (!finite(point)) return false;
    }
    candidate = origin + point * scale;
    return finite(candidate);
}
} // namespace
OcctBoundaryCertificate certifyOcctBoundarySurface(const mesh::Mesh& boundary,
                                                   const OcctStarOptions& options) noexcept {
    OcctBoundaryCertificate outcome;Budget budget{options,outcome.workUsed,outcome.exactBoundaryPairs};
    try {
        prepareBoundary(boundary,options,budget,outcome);
        outcome.certified=true;outcome.reason=OcctStarReason::None;outcome.diagnostic="none";
    }catch(const Stop& stop){outcome.reason=stop.reason;outcome.diagnostic=stop.diagnostic;}
     catch(const std::bad_alloc&){outcome.reason=OcctStarReason::AllocationFailed;outcome.diagnostic="allocation-failed";}
     catch(...){outcome.reason=OcctStarReason::InternalError;outcome.diagnostic="internal-error";}
    return outcome;
}
OcctStarOutcome buildOcctStarVolume(const mesh::Mesh& boundary,int body,
                                   const OcctStarOptions& options,mesh::Mesh& destination) noexcept {
    OcctStarOutcome outcome;
    const std::uint64_t nodes=static_cast<std::uint64_t>(boundary.nodeCount());
    const std::uint64_t triangles=static_cast<std::uint64_t>(boundary.cellCount());
    auto& work=outcome.workUsed;auto& exactPairs=outcome.exactBoundaryPairs;
    auto& kernelTrials=outcome.kernelTrials;auto& plcVolume=outcome.boundaryVolume;
    auto& kernel=outcome.kernel;auto& kernelFeasible=outcome.kernelFeasible;
    auto& embedded=outcome.embedded;auto& minQuality=outcome.minimumQuality;
    Budget budget{options,work,exactPairs};OcctBoundaryCertificate certificate;
    const auto tick=[&](std::uint64_t count=1){budget.tick(count);};
    try {
        need(body>=0&&nodes>=4&&triangles>=4&&nodes<=UINT32_MAX,"input-count");
        need((options.nativeAdmissionModel!=nullptr)==(options.nativeAdmissionBodyTag>=0),
             "invalid-native-kernel-admission");
        need(nodes<=options.maximumBoundaryNodes&&triangles<=options.maximumBoundaryTriangles,
             "boundary-count-budget",OcctStarReason::WorkBudget);
        need(nodes+1<=options.maximumNodes&&triangles<=options.maximumCells/2,
             "output-count-budget",OcctStarReason::OutputBudget);
        const std::uint64_t nodeBytes=sizeof(Vec3)+sizeof(mesh::Classification);
        const std::uint64_t cellBytes=2*sizeof(mesh::Cell)+7*sizeof(Index);
        need(nodes+1<=UINT64_MAX/nodeBytes&&triangles<=UINT64_MAX/cellBytes,
             "output-size-overflow",OcctStarReason::OutputBudget);
        const auto nb=(nodes+1)*nodeBytes,cb=triangles*cellBytes;
        need(nb<=options.maximumOutputBytes&&cb<=options.maximumOutputBytes-nb,
             "output-byte-budget",OcctStarReason::OutputBudget);
        auto prepared=prepareBoundary(boundary,options,budget,certificate,true);
        embedded=certificate.embedded;plcVolume=prepared.volume;
        const auto& points=prepared.points;const auto& faces=prepared.faces;
        const auto lo=prepared.lo,hi=prepared.hi,mean=prepared.mean,centroid=prepared.centroid;
        const auto scale=prepared.scale;
        paacs::geo::detail::OcctSolidClassifier nativeClassifier;
        double nativeTolerance=0;
        if(options.nativeAdmissionModel){
            need(options.nativeAdmissionModel->solid(options.nativeAdmissionBodyTag)!=nullptr,
                 "invalid-native-admission-body");
            nativeTolerance=options.nativeAdmissionModel->tolerance(3,options.nativeAdmissionBodyTag);
            need(std::isfinite(nativeTolerance)&&nativeTolerance>=0,"invalid-native-admission-tolerance");
            tick();nativeClassifier.Load(OcctModelAccess::solid(*options.nativeAdmissionModel,
                                                               options.nativeAdmissionBodyTag));
        }
        const auto accepts=[&](Vec3 candidate) {
            ++kernelTrials;if(!finite(candidate))return false;
            for(const auto& f:faces) {tick();if(predicates::orient3d(points[f.ids[0]].p,points[f.ids[1]].p,points[f.ids[2]].p,candidate)!=predicates::Sign::Negative)return false;}
            if(options.nativeAdmissionModel){
                const auto admits=[&](Vec3 p){tick();++outcome.nativeCandidateQueries;
                    need(finite(p),"nonfinite-native-proposal",OcctStarReason::NativeQueryFailed);
                    nativeClassifier.Perform(gp_Pnt(p.x,p.y,p.z),nativeTolerance);
                    const auto state=nativeClassifier.State();
                    need(state!=TopAbs_UNKNOWN,"native-kernel-query-failed",OcctStarReason::NativeQueryFailed);
                    return state==TopAbs_IN||state==TopAbs_ON;};
                if(!admits(candidate))return false;
                for(const auto& f:faces){
                    const auto p=(candidate+points[f.ids[0]].p+points[f.ids[1]].p+points[f.ids[2]].p)*.25;
                    if(!admits(p))return false;
                }
            }
            kernel=candidate;return true;
        };
        std::vector<Vec3> starts{centroid,lo*.5+hi*.5,mean};
        for(auto seed:starts) {
            if(accepts(seed)){kernelFeasible=true;break;}
            Vec3 x=seed;
            for(unsigned sweep=0;sweep<options.maximumProposalSweeps;++sweep) {
                const double margin=scale*1e-8*std::pow(.5,static_cast<double>(sweep/32));
                for(const auto& f:faces) {tick();const double d=(x-points[f.ids[0]].p).dot(f.normal)+margin;if(d>0)x=x-f.normal*d;}
                if(accepts(x)){kernelFeasible=true;break;}
            }
            if(kernelFeasible)break;
        }
        // Append proposals after the existing bounded search so previously
        // accepted kernels, input identities and published fans stay identical.
        if (!kernelFeasible) for (const double margin : std::array{1e-10, 1e-8, 1e-12}) {
            Vec3 proposal;
            if (halfspaceProposal(prepared, margin, proposal, budget) && accepts(proposal)) {
                kernelFeasible = true; break;
            }
        }
        if (!kernelFeasible) for (const double margin : std::array{1e-10, 1e-8, 1e-12}) {
            Vec3 proposal;
            if (halfspaceProposal(prepared, margin, proposal, budget, true) && accepts(proposal)) {
                kernelFeasible = true; break;
            }
        }
        // Native admission may reject a mathematically valid PLC fan. Try a
        // finite set of deeper numerical kernels and interior lattice points;
        // every proposal still passes every exact side and native fan test.
        if(!kernelFeasible&&options.nativeAdmissionModel){
            for(const double margin:std::array{.05,.02,.01,.005,.001,.0001}){
                Vec3 p;if(halfspaceProposal(prepared,margin,p,budget,true)&&accepts(p)){
                    kernelFeasible=true;break;
                }
            }
            if(!kernelFeasible)for(const double x:std::array{.25,.5,.75}){
                for(const double y:std::array{.25,.5,.75}){
                    for(const double z:std::array{.25,.5,.75}){
                        if(accepts({lo.x+(hi.x-lo.x)*x,lo.y+(hi.y-lo.y)*y,lo.z+(hi.z-lo.z)*z})){
                            kernelFeasible=true;break;
                        }
                    }if(kernelFeasible)break;
                }if(kernelFeasible)break;
            }
        }
        need(kernelFeasible,"no-strict-star-kernel-found",OcctStarReason::KernelNotFound);
        mesh::Mesh result;
        result.reserveAppend(static_cast<std::size_t>(nodes+1),static_cast<std::size_t>(2*triangles),static_cast<std::size_t>(7*triangles));
        for(const auto& p:points){tick();result.addNode(p.p,p.cls);}
        const Index center=result.addNode(kernel,{3,body});
        for(const auto& f:faces) {tick();const std::array<Index,4> ids{center,f.ids[0],f.ids[1],f.ids[2]};
            need(predicates::orient3d(kernel,points[f.ids[0]].p,points[f.ids[1]].p,points[f.ids[2]].p)==predicates::Sign::Positive,"nonpositive-fan-tet",OcctStarReason::InvalidFan);
            result.addCell(mesh::CellType::Tet4,ids,{3,body});}
        for(const auto& f:faces){tick();const std::array<Index,3> ids{f.ids[0],f.ids[1],f.ids[2]};result.addCell(mesh::CellType::Tri3,ids,{2,f.tag});}
        auto auditLimits=options.auditLimits;auditLimits.cancellation=options.cancellation;
        const auto audited=mesh::auditLinearVolumeBoundary(result,auditLimits);
        need(!audited.cancelled,"cancelled",OcctStarReason::Cancelled);
        need(!audited.budgetExceeded,"audit-budget",OcctStarReason::WorkBudget);
        need(!audited.allocationFailed,"allocation-failed",OcctStarReason::AllocationFailed);
        need(audited.valid,"fan-boundary-audit",OcctStarReason::InvalidFan);
        tick(triangles);
        const auto quality=algo::analyzeMixedMeshQuality(result);
        need(!quality.allocationFailed,"fan-quality-allocation-failed",OcctStarReason::AllocationFailed);
        need(quality.valid,"fan-quality-audit",OcctStarReason::InvalidFan);minQuality=quality.minimumQuality;
        need(std::abs(quality.totalVolume-plcVolume)<=plcVolume*1e-10,"fan-plc-volume-mismatch",OcctStarReason::InvalidFan);
        tick(0);
        destination.swap(result);
        outcome.published=true;outcome.destinationUnchanged=false;
        outcome.reason=OcctStarReason::None;outcome.diagnostic="none";
    }catch(const Stop& stop){outcome.reason=stop.reason;outcome.diagnostic=stop.diagnostic;}
     catch(const Standard_OutOfMemory&){outcome.reason=OcctStarReason::AllocationFailed;outcome.diagnostic="native-allocation-failed";}
     catch(const std::bad_alloc&){outcome.reason=OcctStarReason::AllocationFailed;outcome.diagnostic="allocation-failed";}
     catch(const Standard_Failure&){outcome.reason=OcctStarReason::NativeQueryFailed;outcome.diagnostic="native-kernel-query-failed";}
     catch(...){outcome.reason=OcctStarReason::InternalError;outcome.diagnostic="internal-error";}
    outcome.embedded=certificate.embedded;
    return outcome;
}
} // namespace paacs::geo::detail
