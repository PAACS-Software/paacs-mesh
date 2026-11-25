// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/OcctBoundary.hpp"
#include "OcctModelInternal.hpp"
#include "OcctBoundaryInternal.hpp"
#include "OcctChartTrim.hpp"
#include "OcctQueryContext.hpp"
#include "OcctSurfaceAliasRepair.hpp"
#include "OcctUvRefinement.hpp"
#include "OcctSurfaceOppositePairs.hpp"
#include "OcctWireChart.hpp"
#include "OcctStageTrace.hpp"
#include "OcctSurfaceSpacing.hpp"
#include "OcctSurfaceRelaxation.hpp"
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <Message_ProgressIndicator.hxx>
#include <Poly_PolygonOnTriangulation.hxx>
#include <Poly_Triangulation.hxx>
#include <Standard_OutOfMemory.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <ostream>
#include <set>
#include <sstream>
#include <streambuf>
#include <string_view>

namespace paacs::geo {
namespace {
using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;
bool finite(Vec3 p) { return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z); }
Vec3 vec(const gp_Pnt& p) { return {p.X(),p.Y(),p.Z()}; }
bool closestSolution(const BRepExtrema_DistShapeShape& distance, const Vec3& sample, Vec3& closest) {
    if (!distance.IsDone() || distance.NbSolution() < 1) return false;
    const double value = distance.Value();
    if (!std::isfinite(value) || value < 0) return false;
    // OCCT may report several tolerance-equivalent supports: its first edge
    // solution can differ from Value while a later face solution realizes it.
    // Preserve the physical consistency threshold and inspect every solution.
    Vec3 best{};
    double bestDistance = (std::numeric_limits<double>::infinity)();
    for (int i = 1; i <= distance.NbSolution(); ++i) {
        const auto candidate = vec(distance.PointOnShape2(i));
        const double measured = (sample-candidate).norm();
        if (!finite(candidate) || !std::isfinite(measured) ||
            std::abs(measured-value) > 1e-10*std::max({1.,measured,value})) continue;
        if (measured < bestDistance) { best = candidate; bestDistance = measured; }
    }
    if (!std::isfinite(bestDistance)) return false;
    closest = best;
    return true;
}
bool closestOnFace(const TopoDS_Face& face, const Vec3& sample, Vec3& closest) {
    const BRepExtrema_DistShapeShape distance(
        BRepBuilderAPI_MakeVertex(gp_Pnt(sample.x,sample.y,sample.z)).Vertex(),face);
    return closestSolution(distance,sample,closest);
}
struct Failure { OcctBoundaryReason reason; };
struct Context {
    const OcctBoundaryOptions& options;
    OcctBoundaryOutcome& result;
    void require(bool value, OcctBoundaryReason reason) const { if (!value) throw Failure{reason}; }
    void tick(std::uint64_t amount = 1) {
        if (options.cancellation && options.cancellation->load(std::memory_order_relaxed))
            throw Failure{OcctBoundaryReason::Cancelled};
        if (amount > options.maximumWork-result.work) throw Failure{OcctBoundaryReason::BudgetExceeded};
        result.work += amount;
    }
    void admit(std::uint64_t nodes, std::uint64_t triangles) const {
        constexpr std::uint64_t nodeBytes = sizeof(Vec3)+sizeof(mesh::Classification)+sizeof(OcctBoundaryNodeKey);
        constexpr std::uint64_t triangleBytes = sizeof(mesh::Cell)+3*sizeof(Index);
        require(nodes <= options.maximumNodes && triangles <= options.maximumTriangles &&
                nodes < static_cast<std::uint64_t>((std::numeric_limits<int>::max)()) &&
                triangles <= static_cast<std::uint64_t>((std::numeric_limits<int>::max)()),
                OcctBoundaryReason::BudgetExceeded);
        require(nodes <= options.maximumOutputBytes/nodeBytes, OcctBoundaryReason::BudgetExceeded);
        const auto remaining = options.maximumOutputBytes-nodes*nodeBytes;
        require(triangles <= remaining/triangleBytes, OcctBoundaryReason::BudgetExceeded);
    }
};
class CancellationProgress final : public Message_ProgressIndicator {
public:
    explicit CancellationProgress(const std::atomic<bool>* cancellation) : cancellation_(cancellation) {}
protected:
    bool UserBreak() override { return cancellation_ && cancellation_->load(std::memory_order_relaxed); }
    void Show(const Message_ProgressScope&, bool) override {}
private:
    const std::atomic<bool>* cancellation_;
};
struct Sample { Vec3 point; mesh::Classification classification; OcctBoundaryNodeKey key; };
struct EdgeSamples { std::vector<int> ids; std::vector<double> parameters; std::uint64_t uses = 0; };
struct FaceSamples {
    int tag = 0;
    TopoDS_Face face;
    TopLoc_Location location;
    occ::handle<Poly_Triangulation> triangulation;
    std::vector<int> ids;
    std::vector<bool> wireNodes;
    std::vector<std::pair<int,Vec2>> interior;
    std::size_t triangleFirst=0,triangleLast=0;
    bool rebuild=false;
};
std::string snapshot(const TopoDS_Shape& shape) {
    std::ostringstream stream; BRepTools::Write(shape,stream);
    if (!stream.good()) throw Failure{OcctBoundaryReason::KernelFailure};
    return stream.str();
}
std::string facetSnapshot(const TopoDS_Shape&,std::uint64_t maximum);

mesh::Mesh extract(const OcctModel& model, int bodyTag, BRepBuilderAPI_Copy* copy,
                   const ShapeMap& copiedEdges, const std::vector<int>& copiedEdgeTags, Context& context,
                   detail::OcctBoundaryFacetPolicy facetPolicy,
                   const detail::OcctBoundaryAliasOptions& aliasOptions,
                   detail::OcctBoundaryAliasReceipt& aliasReceipt,bool regularFront=false,bool improveQuality=false) {
    auto& result = context.result;
    const bool repairCharts=aliasOptions.maximumRebuiltCharts!=0;
    const bool curveMetric=repairCharts&&aliasOptions.curveMetricSizeCap>0&&model.solids().size()==1;
    const bool spaceSurface=copy==nullptr;
    const bool normalResources=context.options.construction.resourcePolicy==CadResourcePolicy::Normal;
    if(normalResources&&repairCharts)
        context.require(model.faceUses(bodyTag).size()<=aliasOptions.maximumWorkingBytes/1024,
                        OcctBoundaryReason::BudgetExceeded);
    std::vector<Sample> points;
    std::vector<std::array<int,3>> triangles;
    std::vector<int> faceTags;
    std::map<int,int> vertexIds;
    std::map<int,EdgeSamples> edges;
    const auto appendPoint = [&](Vec3 p, mesh::Classification classification, std::uint64_t ordinal) {
        context.tick(); context.require(finite(p),OcctBoundaryReason::GeometryError);
        context.admit(points.size()+1,triangles.size());
        const int id = static_cast<int>(points.size());
        points.push_back({p,classification,{classification.dim,classification.tag,ordinal}});
        result.nodes = points.size(); return id;
    };
    const auto vertex = [&](int tag) {
        context.tick();
        if (auto found = vertexIds.find(tag); found != vertexIds.end()) return found->second;
        const auto* source = model.vertex(tag);
        context.require(source != nullptr,OcctBoundaryReason::UnsupportedTopology);
        const int id = appendPoint(source->position(),{0,tag},0); vertexIds[tag] = id; return id;
    };
    std::vector<FaceSamples> faces;
    std::uint64_t localNodeCount = 0, rawTriangleCount = 0;
    for (const auto& use : model.faceUses(bodyTag)) {
        context.tick();
        context.require(use.orientation == CadOrientation::Forward || use.orientation == CadOrientation::Reversed,
                        OcctBoundaryReason::UnsupportedTopology);
        FaceSamples face; face.tag = use.faceTag;
        const auto original=detail::OcctModelAccess::face(model,face.tag);
        face.face = copy?TopoDS::Face(copy->ModifiedShape(original)):original;
        // ModifiedShape can normalize orientation; restore the actual owning
        // solid's occurrence independently of the shared GFace normal.
        face.face.Orientation(use.orientation == CadOrientation::Forward ? TopAbs_FORWARD : TopAbs_REVERSED);
        if(spaceSurface){
            face.rebuild=true;
            ++result.faces;
            for(const auto& wire:model.boundaryUses(face.tag))for(const auto& edgeUse:wire.edges){
                context.tick();
                ++result.edgeUses;
                result.seamUses+=edgeUse.seam;
                result.degeneratedUses+=edgeUse.degenerated;
                auto [chain,inserted]=edges.try_emplace(edgeUse.edgeTag);
                if(!inserted)continue;
                const auto* edge=model.edge(edgeUse.edgeTag);
                const auto start=vertex(edge->startVertexTag()),end=vertex(edge->endVertexTag());
                context.require(!edgeUse.degenerated||start==end,OcctBoundaryReason::UnsupportedTopology);
                chain->second.ids={start,end};
                chain->second.parameters={edgeUse.firstParameter,edgeUse.lastParameter};
            }
            faces.push_back(std::move(face));
            continue;
        }
        face.triangulation = BRep_Tool::Triangulation(face.face,face.location);
        context.require(!face.triangulation.IsNull() && face.triangulation->HasUVNodes(),
                        OcctBoundaryReason::MissingTriangulation);
        const auto n = static_cast<std::uint64_t>(face.triangulation->NbNodes());
        const auto t = static_cast<std::uint64_t>(face.triangulation->NbTriangles());
        context.require(n > 0 && t > 0 && n <= context.options.maximumWork-localNodeCount &&
                        t <= context.options.maximumTriangles-rawTriangleCount,
                        OcctBoundaryReason::BudgetExceeded);
        localNodeCount += n; rawTriangleCount += t;
        face.ids.assign(static_cast<std::size_t>(n)+1,-1); ++result.faces;
        if(curveMetric)face.wireNodes.resize(static_cast<std::size_t>(n)+1);
        // On a forward chart, the edge-use orientation selects OCCT's first
        // or second polygon branch at a periodic seam.
        const auto forward = TopoDS::Face(face.face.Oriented(TopAbs_FORWARD));
        for (TopExp_Explorer explorer(forward,TopAbs_EDGE); explorer.More(); explorer.Next()) {
            context.tick();
            const auto used = TopoDS::Edge(explorer.Current());
            const int copiedIndex = copiedEdges.FindIndex(used);
            context.require(copiedIndex > 0,OcctBoundaryReason::UnsupportedTopology);
            const int edgeTag = copiedEdgeTags[static_cast<std::size_t>(copiedIndex)-1];
            const auto* edge = model.edge(edgeTag);
            const bool pole = BRep_Tool::Degenerated(used);
            const auto polygon = BRep_Tool::PolygonOnTriangulation(used,face.triangulation,face.location);
            context.require(!polygon.IsNull() && polygon->NbNodes() >= 2,OcctBoundaryReason::MissingTriangulation);
            ++result.edgeUses;
            if (BRep_Tool::IsClosed(used,face.triangulation,face.location)) ++result.seamUses;
            if (pole) ++result.degeneratedUses;
            const double tolerance = 8*std::max({model.tolerance(1,edgeTag),model.tolerance(2,face.tag),1e-10});
            auto& canonical = edges[edgeTag];
            const auto bindSample=[&](int local,int id,double error,OcctBoundaryReason reason){
                if(curveMetric)face.wireNodes[local]=true;
                const bool fits=std::isfinite(error)&&error<=tolerance;
                const bool distinct=face.ids[local]>=0&&face.ids[local]!=id;
                if(repairCharts&&(!fits||distinct)){
                    face.rebuild=true;
                    return;
                }
                if(!pole)result.maximumStitchError=std::max(result.maximumStitchError,error);
                result.maximumNodeProjectionError=std::max(result.maximumNodeProjectionError,error);
                context.require(fits,reason);
                context.require(!distinct,OcctBoundaryReason::NonconformingEdge);
                face.ids[local]=id;
            };
            if (pole) {
                context.require(edge->startVertexTag() == edge->endVertexTag(),OcctBoundaryReason::UnsupportedTopology);
                const int id = vertex(edge->startVertexTag());
                if(repairCharts&&canonical.uses==0)canonical.ids={id,id};
                for (int j=1; j<=polygon->NbNodes(); ++j) {
                    context.tick(); const int local = polygon->Node(j);
                    context.require(local >= 1 && local <= face.triangulation->NbNodes(),OcctBoundaryReason::InvalidSurface);
                    const auto raw = vec(face.triangulation->Node(local).Transformed(face.location.Transformation()));
                    bindSample(local,id,(raw-points[id].point).norm(),OcctBoundaryReason::GeometryError);
                }
                ++canonical.uses; continue;
            }
            context.require(polygon->HasParameters(),OcctBoundaryReason::UnsupportedTopology);
            double first,last; edge->parameterBounds(first,last);
            context.require(std::isfinite(first) && std::isfinite(last) && last > first,OcctBoundaryReason::GeometryError);
            const double parameterTolerance = 1e-10*std::max(1.,last-first);
            const bool reverse = polygon->Parameter(1) > polygon->Parameter(polygon->NbNodes());
            context.require(std::abs(polygon->Parameter(reverse?polygon->NbNodes():1)-first) <= parameterTolerance &&
                            std::abs(polygon->Parameter(reverse?1:polygon->NbNodes())-last) <= parameterTolerance,
                            OcctBoundaryReason::NonconformingEdge);
            if (canonical.uses > 0)
                context.require(canonical.ids.size() == static_cast<std::size_t>(polygon->NbNodes()),
                                OcctBoundaryReason::NonconformingEdge);
            for (int j=0; j<polygon->NbNodes(); ++j) {
                context.tick(); const int k = reverse?polygon->NbNodes()-j:j+1, local = polygon->Node(k);
                context.require(local >= 1 && local <= face.triangulation->NbNodes(),OcctBoundaryReason::InvalidSurface);
                const double parameter = polygon->Parameter(k);
                context.require(std::isfinite(parameter) && parameter >= first-parameterTolerance && parameter <= last+parameterTolerance,
                                OcctBoundaryReason::GeometryError);
                if (j > 0) context.require(parameter > (canonical.uses?canonical.parameters[j-1]:canonical.parameters.back()),
                                           OcctBoundaryReason::NonconformingEdge);
                const auto raw = vec(face.triangulation->Node(local).Transformed(face.location.Transformation()));
                int id;
                if (canonical.uses == 0) {
                    if (j == 0) id = vertex(edge->startVertexTag());
                    else if (j+1 == polygon->NbNodes()) id = vertex(edge->endVertexTag());
                    else id = appendPoint(edge->pointAt(parameter),{1,edgeTag},static_cast<std::uint64_t>(j));
                    canonical.ids.push_back(id); canonical.parameters.push_back(parameter);
                } else {
                    context.require(std::abs(parameter-canonical.parameters[j]) <= parameterTolerance,
                                    OcctBoundaryReason::NonconformingEdge);
                    id = canonical.ids[j];
                }
                bindSample(local,id,(raw-points[id].point).norm(),OcctBoundaryReason::NonconformingEdge);
            }
            ++canonical.uses;
        }
        faces.push_back(std::move(face));
    }
    Vec3 lower,upper; model.solid(bodyTag)->bbox(lower,upper);
    context.require(finite(lower) && finite(upper),OcctBoundaryReason::GeometryError);
    const Vec3 origin = lower+(upper-lower)*.5;
    context.require(finite(origin),OcctBoundaryReason::GeometryError);
    double sum = 0, compensation = 0;
    result.minimumTriangleArea = (std::numeric_limits<double>::infinity)();
    int midpointFace=-1;
    std::map<std::array<int,2>,std::pair<Vec3,double>> midpointWitnesses;
    const auto accumulate=[&](Vec3 p,Vec3 q,Vec3 r){
        const double area=(q-p).cross(r-p).norm()*.5;
        context.require(std::isfinite(area)&&area>0,OcctBoundaryReason::InvalidSurface);
        result.minimumTriangleArea=std::min(result.minimumTriangleArea,area);
        const double volume=(p-origin).dot((q-origin).cross(r-origin))/6.;
        context.require(std::isfinite(volume),OcctBoundaryReason::GeometryError);
        const double next=sum+volume;
        compensation+=std::abs(sum)>=std::abs(volume)?(sum-next)+volume:(volume-next)+sum;
        sum=next;
    };
    const auto facetFits=[&](const FaceSamples& face,const std::array<int,3>& nodes,
                             const std::array<Vec2,3>* hint){
        if(facetPolicy!=detail::OcctBoundaryFacetPolicy::Immediate)return true;
        const auto* source=model.face(face.tag);
        const auto chord=context.options.maximumChordError;
        const auto tolerance=model.tolerance(2,face.tag);
        const auto p=points[nodes[0]].point,q=points[nodes[1]].point,r=points[nodes[2]].point;
        const std::array samples{(p+q+r)/3.,(p+q)/2.,(q+r)/2.,(r+p)/2.};
        std::array<Vec2,4> parameters{};
        if(hint){
            if(midpointFace!=face.tag){
                midpointWitnesses.clear();
                midpointFace=face.tag;
            }
            const auto [a,b,c]=*hint;
            parameters={{{(a.u+b.u+c.u)/3.,(a.v+b.v+c.v)/3.},
                         {(a.u+b.u)/2.,(a.v+b.v)/2.},
                         {(b.u+c.u)/2.,(b.v+c.v)/2.},
                         {(c.u+a.u)/2.,(c.v+a.v)/2.}}};
        }
        for(std::size_t i=0;i<samples.size();++i){
            context.tick();
            const auto sample=samples[i];
            Vec2 uv;
            Vec3 closest;
            double error=(std::numeric_limits<double>::infinity)();
            std::array<int,2> midpoint{};
            if(hint&&i!=0){
                midpoint={nodes[i-1],nodes[i%3]};
                if(midpoint[1]<midpoint[0])std::swap(midpoint[0],midpoint[1]);
                if(auto found=midpointWitnesses.find(midpoint);found!=midpointWitnesses.end()){
                    closest=found->second.first;
                    error=found->second.second;
                    ++aliasReceipt.reusedMidpoints;
                }
            }
            // A point on the trimmed CAD face is a direct distance witness.
            // Fall back to the regular physical query when the hint misses.
            if(hint&&!std::isfinite(error)){
                const auto witness=[&](Vec2 parameter){
                    const auto point=source->pointAt(parameter);
                    const auto distance=(sample-point).norm();
                    if(!(distance<=chord)||!source->containsParam(parameter))return false;
                    closest=point;
                    error=distance;
                    return true;
                };
                if(!witness(parameters[i])&&i!=0){
                    // The chord in UV can pass just outside a curved trim.
                    // A nearby point toward the triangle center can witness
                    // the distance without running a global CAD projection.
                    const auto edge=parameters[i],center=parameters[0];
                    for(double inset:{1./64,1./8,1./2}){
                        context.tick();
                        if(witness({edge.u+(center.u-edge.u)*inset,
                                    edge.v+(center.v-edge.v)*inset}))break;
                    }
                }
            }
            if(!std::isfinite(error)||error>chord){
                const auto began=hint?std::chrono::steady_clock::now():std::chrono::steady_clock::time_point{};
                if(source->projectPoint(sample,uv,closest))error=(sample-closest).norm();
                if(!std::isfinite(error)||error>chord){
                    ++result.closestTrimmedFaceSamples;
                    context.require(closestOnFace(face.face,sample,closest),OcctBoundaryReason::GeometryError);
                    error=(sample-closest).norm();
                }
                if(hint)aliasReceipt.projectionSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
            }
            context.require(std::isfinite(error),OcctBoundaryReason::GeometryError);
            if(error>result.maximumFacetSampleError){
                result.maximumFacetSampleError=error;
                result.maximumErrorFace=face.tag;
                result.maximumErrorSample=sample;
                result.maximumErrorClosest=closest;
            }
            if(error>chord+tolerance)return false;
            if(hint&&i!=0&&error<=chord&&midpointWitnesses.size()<16384)
                midpointWitnesses.try_emplace(midpoint,closest,error);
        }
        return true;
    };
    for (auto& face : faces) {
        context.tick(); const auto* source = model.face(face.tag);
        if(spaceSurface)continue;
        const auto interiorStart=points.size();
        for (int i=1; i<=face.triangulation->NbNodes()&&(!face.rebuild||curveMetric); ++i) {
            context.tick(); if (face.ids[i] >= 0||(curveMetric&&face.wireNodes[i])) continue;
            const auto uv = face.triangulation->UVNode(i);
            const Vec2 parameter{uv.X(),uv.Y()};
            if(!source->containsParam(parameter)){
                context.require(repairCharts,OcctBoundaryReason::GeometryError);
                face.rebuild=true;
                if(curveMetric)continue;
                break;
            }
            const auto p = source->pointAt(parameter);
            const auto raw = vec(face.triangulation->Node(i).Transformed(face.location.Transformation()));
            const double error = (p-raw).norm();
            const bool fits=finite(p)&&std::isfinite(error)&&error<=8*std::max(model.tolerance(2,face.tag),1e-10);
            if(repairCharts&&!fits){face.rebuild=true;if(curveMetric)continue;break;}
            result.maximumNodeProjectionError = std::max(result.maximumNodeProjectionError,error);
            context.require(fits,OcctBoundaryReason::GeometryError);
            face.ids[i] = appendPoint(p,{2,face.tag},static_cast<std::uint64_t>(i));
            if(repairCharts)face.interior.emplace_back(face.ids[i],parameter);
        }
        if(face.rebuild){
            if(!curveMetric){
                points.resize(interiorStart);
                result.nodes=points.size();
                face.interior.clear();
            }
            continue;
        }
        // Keep native chart connectivity. Repair only a face whose actual
        // facet samples fail; shared CAD-wire segments remain protected.
        std::vector<std::array<int,3>> refinedTriangles;
        const auto faceStart=triangles.size();const auto savedSum=sum,savedCompensation=compensation;
        face.triangleFirst=faceStart;
        const auto savedMinimumArea=result.minimumTriangleArea,savedFacetError=result.maximumFacetSampleError;
        const auto savedErrorFace=result.maximumErrorFace;const auto savedSample=result.maximumErrorSample,savedClosest=result.maximumErrorClosest;
        const auto emitFace=[&](){
        for (int t=1; t<=(refinedTriangles.empty()?face.triangulation->NbTriangles():int(refinedTriangles.size())); ++t) {
            context.tick(); int a=1,b=1,c=1; if(refinedTriangles.empty())face.triangulation->Triangle(t).Get(a,b,c);
            context.require(a >= 1 && b >= 1 && c >= 1 && a <= face.triangulation->NbNodes() &&
                            b <= face.triangulation->NbNodes() && c <= face.triangulation->NbNodes(),
                            OcctBoundaryReason::InvalidSurface);
            std::array<int,3> triangle=refinedTriangles.empty()?std::array<int,3>{face.ids[a],face.ids[b],face.ids[c]}:refinedTriangles[t-1];
            if (triangle[0] == triangle[1] || triangle[1] == triangle[2] || triangle[2] == triangle[0]) {
                bool pole = false;
                for (int i=0; i<3; ++i) for (int j=i+1; j<3; ++j)
                    if (triangle[i] == triangle[j]) pole |= points[triangle[i]].classification.dim == 0;
                context.require(pole,OcctBoundaryReason::InvalidSurface);
                ++result.collapsedPoleTriangles; continue;
            }
            std::array<Vec2,3> hint;
            const bool useHint=repairCharts&&refinedTriangles.empty();
            if(useHint){
                const std::array local{a,b,c};
                for(unsigned i=0;i<3;++i){
                    const auto uv=face.triangulation->UVNode(local[i]);
                    hint[i]={uv.X(),uv.Y()};
                }
            }
            if(face.face.Orientation()==TopAbs_REVERSED){
                std::swap(triangle[1],triangle[2]);
                if(useHint)std::swap(hint[1],hint[2]);
            }
            const auto p = points[triangle[0]].point, q = points[triangle[1]].point, r = points[triangle[2]].point;
            accumulate(p,q,r);
            if(!facetFits(face,triangle,useHint?&hint:nullptr))return false;
            context.admit(points.size(),triangles.size()+1);
            triangles.push_back(triangle); faceTags.push_back(face.tag); result.triangles = triangles.size();
        }
        return true;};
        if(!emitFace()){
            triangles.resize(faceStart);faceTags.resize(faceStart);result.triangles=triangles.size();
            sum=savedSum;compensation=savedCompensation;result.minimumTriangleArea=savedMinimumArea;
            std::vector<Vec2> uv;std::vector<int> global;std::vector<algo::Tri> nativeTriangles;
            for(int i=1;i<=face.triangulation->NbNodes();++i){context.tick();const auto p=face.triangulation->UVNode(i);uv.push_back({p.X(),p.Y()});global.push_back(face.ids[i]);}
            std::map<std::array<int,2>,unsigned> uses;
            for(int i=1;i<=face.triangulation->NbTriangles();++i){int a,b,c;face.triangulation->Triangle(i).Get(a,b,c);nativeTriangles.push_back({{a-1,b-1,c-1}});for(auto edge:std::array{std::array{a-1,b-1},std::array{b-1,c-1},std::array{c-1,a-1}}){context.tick();std::sort(edge.begin(),edge.end());++uses[edge];}}
            std::set<std::array<int,2>> protectedUv;for(auto [edge,count]:uses)if(count==1)protectedUv.insert(edge);
            // Interior insertions cannot repair a chord on a shared CAD
            // wire. Refuse that case without spending its point allowance;
            // the caller can refine the common wire on every incident face.
            for(auto edge:protectedUv){
                context.tick();const auto sample=(points[global[edge[0]]].point+points[global[edge[1]]].point)*.5;
                Vec2 parameter;Vec3 closest;double error=(std::numeric_limits<double>::infinity)();
                if(source->projectPoint(sample,parameter,closest))error=(sample-closest).norm();
                if(!std::isfinite(error)||error>context.options.maximumChordError){++result.closestTrimmedFaceSamples;context.require(closestOnFace(face.face,sample,closest),OcctBoundaryReason::GeometryError);error=(sample-closest).norm();}
                context.require(std::isfinite(error)&&error<=context.options.maximumChordError+model.tolerance(2,face.tag),OcctBoundaryReason::GeometryError);
            }
            const auto originalUv=uv;const auto originalGlobal=global;const auto originalTriangles=nativeTriangles;
            const auto pointStart=points.size();const auto room=context.options.maximumNodes-pointStart;
            const unsigned maximumAdded=unsigned(std::min<std::uint64_t>(1024,room));
            bool accepted=false,bounded=false;
            // Closest-distance guidance avoids overrefining distorted charts.
            // If it stalls, restart this face with the stricter parametric
            // correspondence bound. Both attempts protect the same wires and
            // must pass the unchanged physical sample audit before acceptance.
            for(unsigned attempt=0;attempt<2&&!accepted;++attempt){
                triangles.resize(faceStart);faceTags.resize(faceStart);result.triangles=triangles.size();
                points.resize(pointStart);result.nodes=points.size();sum=savedSum;compensation=savedCompensation;result.minimumTriangleArea=savedMinimumArea;
                uv=originalUv;global=originalGlobal;nativeTriangles=originalTriangles;refinedTriangles.clear();
                const auto position=[&](int id){return points[id].point;};const auto evaluate=[&](Vec2 p){return source->pointAt(p);};
                const auto inside=[&](Vec2 p){return source->containsParam(p);};
                const auto append=[&](Vec2,Vec3 p){return appendPoint(p,{2,face.tag},std::uint64_t(face.triangulation->NbNodes())+uv.size()+1);};
                const auto charge=[&](std::uint64_t n){context.tick(n);};
                const auto deviation=[&](Vec3 p){context.tick();Vec3 closest;return closestOnFace(face.face,p,closest)?(p-closest).norm():std::numeric_limits<double>::infinity();};
                const auto refined=attempt==0?
                    detail::refineOcctUvTriangles(uv,global,nativeTriangles,protectedUv,context.options.maximumChordError*.9,maximumAdded,position,evaluate,inside,append,charge,deviation):
                    detail::refineOcctUvTriangles(uv,global,nativeTriangles,protectedUv,context.options.maximumChordError*.9,maximumAdded,position,evaluate,inside,append,charge);
                bounded=refined.added==maximumAdded||uv.size()>=4096;
                if(!refined.completed||!refined.added)continue;
                for(auto t:nativeTriangles)refinedTriangles.push_back({global[t.v[0]],global[t.v[1]],global[t.v[2]]});
                result.maximumFacetSampleError=savedFacetError;result.maximumErrorFace=savedErrorFace;result.maximumErrorSample=savedSample;result.maximumErrorClosest=savedClosest;
                accepted=emitFace();
            }
            context.require(accepted,bounded?OcctBoundaryReason::BudgetExceeded:OcctBoundaryReason::GeometryError);
        }
        face.triangleLast=triangles.size();
    }
    std::set<int> retiredWireNodes;
    if(repairCharts){
        using Edge=std::array<int,2>;
        const auto toggle=[](std::set<Edge>& boundary,int a,int b){
            if(a==b)return;
            const Edge edge{std::min(a,b),std::max(a,b)};
            if(!boundary.insert(edge).second)boundary.erase(edge);
        };
        for(auto& face:faces){
            context.tick();
            if(face.rebuild)continue;
            // Interior uses cancel in pairs. The remaining edges must match
            // the CAD wire, including all segments incident to a pole.
            std::set<Edge> discrepancy;
            for(auto i=face.triangleFirst;i<face.triangleLast;++i){
                context.tick(3);
                const auto t=triangles[i];
                for(unsigned j=0;j<3;++j)toggle(discrepancy,t[j],t[(j+1)%3]);
            }
            for(const auto& wire:model.boundaryUses(face.tag))for(const auto& use:wire.edges){
                context.tick();
                const auto& chain=edges.at(use.edgeTag).ids;
                for(std::size_t i=1;i<chain.size();++i){
                    context.tick();
                    toggle(discrepancy,chain[i-1],chain[i]);
                }
            }
            face.rebuild=!discrepancy.empty();
        }
        const auto position=[&](int id){return points[id].point;};
        const auto tick=[&](std::uint64_t n){context.tick(n);};
        const auto remaining=[&]{return context.options.maximumWork-result.work;};
        const auto noRefine=[](auto&,auto&,auto&,const auto&){return true;};
        const auto admitRepair=[&](std::uint64_t reserve){
            const auto topology=normalResources?std::uint64_t(faces.size())+edges.size():0;
            const auto count=std::uint64_t(points.size())+triangles.size()+topology+reserve;
            context.require(count<=aliasOptions.maximumWorkingBytes/1024,OcctBoundaryReason::BudgetExceeded);
            aliasReceipt.maximumWorkingPayloadEstimate=std::max(aliasReceipt.maximumWorkingPayloadEstimate,count*1024);
        };
        if(curveMetric){
            detail::OcctStageTrace curveTrace("surface_shared_edges",bodyTag);
            std::set<int> crossed;
            if(spaceSurface){
                for(const auto& face:faces)for(const auto& wire:model.boundaryUses(face.tag))
                    for(const auto& use:wire.edges){
                        context.tick();
                        if(!use.degenerated)crossed.insert(use.edgeTag);
                    }
            }else{
                for(const auto& face:faces)if(face.rebuild){
                    std::vector<detail::OcctChartTriangle> unused;
                    const auto complete=detail::rebuildOcctWireChart(model,face.tag,edges,face.interior,
                        position,tick,remaining,unused,&crossed,false,noRefine,false,4096,&aliasReceipt.witness);
                    context.require(complete,OcctBoundaryReason::NonconformingEdge);
                }
            }
            const std::size_t maximumRephasedEdges=normalResources?edges.size():
                spaceSurface?(aliasOptions.curvaturePointsPerTurn?256:96):32;
            context.require(crossed.size()<=maximumRephasedEdges,OcctBoundaryReason::BudgetExceeded);
            if(normalResources){
                admitRepair(0);
                context.require(crossed.size()<=remaining(),OcctBoundaryReason::BudgetExceeded);
            }
            if(!crossed.empty()){
                detail::OcctCurveUses uses;
                detail::OcctVertexEdges incident;
                for(const auto& [tag,chain]:edges){
                    context.tick();
                    const auto* edge=model.edge(tag);
                    incident[edge->startVertexTag()].push_back(tag);
                    if(edge->endVertexTag()!=edge->startVertexTag())incident[edge->endVertexTag()].push_back(tag);
                }
                for(auto& face:faces)for(const auto& wire:model.boundaryUses(face.tag))for(const auto& use:wire.edges){
                    context.tick();
                    uses[use.edgeTag].push_back({face.tag,use.useId});
                    if(crossed.contains(use.edgeTag))face.rebuild=true;
                }
                for(int tag:crossed){
                    auto& chain=edges.at(tag);
                    detail::OcctCurveSpacing spacing;
                    if(spaceSurface)spacing={aliasOptions.curvaturePointsPerTurn?aliasOptions.curvaturePointsPerTurn:16,
                                             context.options.maximumChordError*.8};
                    detail::OcctCurveSamplingReport sampling;
                    spacing.report=&sampling;
                    if(normalResources){
                        admitRepair(0);
                        constexpr std::uint64_t nodeBytes=sizeof(Vec3)+sizeof(mesh::Classification)+sizeof(OcctBoundaryNodeKey);
                        constexpr std::uint64_t triangleBytes=sizeof(mesh::Cell)+3*sizeof(Index);
                        const auto retained=std::uint64_t(points.size())+triangles.size()+faces.size()+edges.size();
                        const auto curveBytes=aliasOptions.maximumWorkingBytes-retained*1024;
                        const auto outputUsed=std::uint64_t(points.size())*nodeBytes+std::uint64_t(triangles.size())*triangleBytes;
                        context.require(outputUsed<=context.options.maximumOutputBytes,OcctBoundaryReason::BudgetExceeded);
                        spacing.resourcePolicy=CadResourcePolicy::Normal;
                        spacing.maximumSegments=std::min({context.options.maximumNodes-std::uint64_t(points.size())+1,
                            std::uint64_t((std::numeric_limits<int>::max)())-points.size(),
                            (context.options.maximumOutputBytes-outputUsed)/nodeBytes+1,remaining()});
                        spacing.maximumIntegrationIntervals=std::min<std::uint64_t>(curveBytes/4096,remaining());
                        spacing.maximumWorkingBytes=curveBytes;
                    }
                    const auto schedule=detail::sampleOcctCurveMetric(model,tag,aliasOptions.curveMetricSizeCap,
                        uses,incident,tick,spacing,aliasOptions.sizeField);
                    if(normalResources&&sampling.budgetLimited){
                        aliasReceipt.witness={"shared_edge_sampling",sampling.diagnostic,tag,-1,std::nullopt,false};
                        context.require(false,OcctBoundaryReason::BudgetExceeded);
                    }
                    if(schedule.size()<2){
                        const auto samplingReason=std::string_view(sampling.diagnostic);
                        const bool geometric=samplingReason=="curve_parameter_resolution"||
                            samplingReason=="curve_metric_resolution";
                        aliasReceipt.witness={"shared_edge_sampling",
                            samplingReason=="not_attempted"?"curve_metric_schedule_incomplete":sampling.diagnostic,
                            tag,-1,std::nullopt,geometric};
                    }
                    context.require(schedule.size()>=2,OcctBoundaryReason::NonconformingEdge);
                    admitRepair(schedule.size()+4096);
                    for(std::size_t j=1;j+1<chain.ids.size();++j){context.tick();retiredWireNodes.insert(chain.ids[j]);}
                    std::vector<int> ids{chain.ids.front()};
                    for(std::size_t j=1;j+1<schedule.size();++j)
                        ids.push_back(appendPoint(model.edge(tag)->pointAt(schedule[j]),{1,tag},j));
                    ids.push_back(chain.ids.back());
                    chain.ids.swap(ids);
                    chain.parameters=schedule;
                    ++aliasReceipt.rephasedEdges;
                }
            }
        }
        if(spaceSurface&&curveMetric){
            // A geometrically accurate coarse curve polygon can still cross
            // another trim chord across a narrow feature. Repair the shared
            // canonical interval before emitting any incident face triangles.
            // Existing samples stay fixed; unrelated curves are not resampled.
            detail::OcctStageTrace topologyTrace("surface_wire_topology",bodyTag);
            std::set<int> pending;
            for(const auto& face:faces)if(face.rebuild)pending.insert(face.tag);
            constexpr unsigned maximumRounds=std::numeric_limits<double>::digits;
            for(unsigned round=0;!pending.empty();++round){
                std::set<detail::OcctChartSplitInterval> intervals;
                for(int tag:pending){
                    std::uint64_t nodes=0;
                    for(const auto& wire:model.boundaryUses(tag))for(const auto& use:wire.edges){
                        context.tick();
                        const auto count=edges.at(use.edgeTag).ids.size()-1;
                        context.require(count<=std::uint64_t((std::numeric_limits<int>::max)())-nodes,
                                        OcctBoundaryReason::BudgetExceeded);
                        nodes+=count;
                    }
                    admitRepair(std::max<std::uint64_t>(4096,nodes*2));
                    std::vector<detail::OcctChartTriangle> unused;
                    const auto complete=detail::rebuildOcctWireChart(model,tag,edges,{},
                        position,tick,remaining,unused,nullptr,false,noRefine,true,
                        std::size_t(nodes),&aliasReceipt.witness,&intervals);
                    context.require(complete,OcctBoundaryReason::NonconformingEdge);
                }
                if(intervals.empty())break;
                if(round>=maximumRounds){
                    aliasReceipt.witness={"shared_chart_sampling","chart_crossing_round_limit",
                                          -1,-1,std::nullopt,true};
                    context.require(false,OcctBoundaryReason::NonconformingEdge);
                }
                const auto append=[&](int tag,double at,std::uint64_t ordinal){
                    return appendPoint(model.edge(tag)->pointAt(at),{1,tag},ordinal);
                };
                const auto complete=detail::bisectOcctChartIntervals(edges,intervals,append,tick,
                    [&](std::size_t reserve){admitRepair(std::uint64_t(reserve)*2);},&aliasReceipt.witness);
                context.require(complete,OcctBoundaryReason::NonconformingEdge);
                pending.clear();
                for(auto first=intervals.begin();first!=intervals.end();){
                    const int tag=first->edgeTag;
                    const auto& chain=edges.at(tag);
                    // Node IDs/coordinates/classifications stay unchanged.
                    // Keys retain their documented increasing edge ordinal.
                    for(std::size_t i=1;i+1<chain.ids.size();++i){
                        context.tick();points[chain.ids[i]].key.ordinal=i;
                    }
                    for(int face:model.edgeFaces(tag)){
                        context.tick();
                        for(const auto& owner:model.faceOwners(face))
                            if(owner.solidTag==bodyTag){pending.insert(face);break;}
                    }
                    do {context.tick();++first;}while(first!=intervals.end()&&first->edgeTag==tag);
                }
            }
        }
        std::vector<FaceSamples*> rebuild;
        for(auto& face:faces)if(face.rebuild)rebuild.push_back(&face);
        context.require(rebuild.size()<=aliasOptions.maximumRebuiltCharts,OcctBoundaryReason::BudgetExceeded);
        if(!rebuild.empty()){
            admitRepair(4096);
            std::size_t kept=0;
            for(const auto& face:faces)if(!face.rebuild){
                for(auto i=face.triangleFirst;i<face.triangleLast;++i){
                    context.tick();
                    triangles[kept]=triangles[i];
                    faceTags[kept++]=faceTags[i];
                }
            }
            triangles.resize(kept);
            faceTags.resize(kept);
            std::sort(rebuild.begin(),rebuild.end(),[](const auto* a,const auto* b){return a->tag<b->tag;});
            for(const auto* face:rebuild){
                if(aliasReceipt.rephasedEdges)admitRepair(4096);
                const auto owners=model.faceOwners(face->tag);
                context.require(owners.size()==1&&owners.front().solidTag==bodyTag,OcctBoundaryReason::UnsupportedTopology);
                std::vector<detail::OcctChartTriangle> replacements;
                detail::OcctStageTrace chartTrace("surface_face_chart",bodyTag);
                const auto& construction=context.options.construction;
                const bool progressiveConstruction=construction.resourcePolicy==CadResourcePolicy::Normal;
                const bool legacyConstruction=!progressiveConstruction&&construction.effort==CadConstructionEffort::Legacy&&
                    !construction.surfacePointsPerFace&&!construction.surfaceChordPointsPerFace;
                // An inserted chart point adds two live surface triangles.
                // This allowance only plans effort; each allocation, append,
                // output facet and work unit keeps its existing admission.
                const auto pointCapacity=[&](std::uint64_t localTriangles){
                    constexpr std::uint64_t nodeBytes=sizeof(Vec3)+sizeof(mesh::Classification)+sizeof(OcctBoundaryNodeKey);
                    constexpr std::uint64_t triangleBytes=sizeof(mesh::Cell)+3*sizeof(Index);
                    const auto nodes=std::uint64_t(points.size());
                    const auto facets=std::uint64_t(triangles.size())+localTriangles;
                    if(nodes>=context.options.maximumNodes||facets>=context.options.maximumTriangles||
                       nodes>=std::uint64_t((std::numeric_limits<int>::max)())||
                       nodes>context.options.maximumOutputBytes/nodeBytes)return std::uint64_t{0};
                    auto output=context.options.maximumOutputBytes-nodes*nodeBytes;
                    if(facets>output/triangleBytes)return std::uint64_t{0};
                    output-=facets*triangleBytes;
                    const auto topology=normalResources?std::uint64_t(faces.size())+edges.size():0;
                    const auto units=aliasOptions.maximumWorkingBytes/1024;
                    const auto memory=nodes+facets+topology<units?(units-nodes-facets-topology)/8:0;
                    return std::min({context.options.maximumNodes-nodes,
                        (context.options.maximumTriangles-facets)/2,output/(nodeBytes+2*triangleBytes),
                        memory,context.options.maximumWork-result.work,
                        std::uint64_t((std::numeric_limits<int>::max)())-nodes});
                };
                std::size_t chartNodeLimit=4096;
                if(!legacyConstruction){
                    std::uint64_t initial=face->interior.size();
                    context.require(initial<=std::uint64_t((std::numeric_limits<int>::max)()),
                                    OcctBoundaryReason::BudgetExceeded);
                    for(const auto& wire:model.boundaryUses(face->tag))for(const auto& use:wire.edges){
                        const auto samples=edges.at(use.edgeTag).ids.size()-1;
                        context.require(samples<=std::uint64_t((std::numeric_limits<int>::max)())-initial,
                                        OcctBoundaryReason::BudgetExceeded);
                        initial+=samples;
                    }
                    const auto room=pointCapacity(0);
                    const auto spacing=construction.surfacePointsPerFace?
                        std::min(construction.surfacePointsPerFace,room):
                        !progressiveConstruction&&construction.effort==CadConstructionEffort::Legacy?std::min<std::uint64_t>(256,room):room;
                    const auto chord=construction.surfaceChordPointsPerFace?
                        std::min(construction.surfaceChordPointsPerFace,room):
                        !progressiveConstruction&&construction.effort==CadConstructionEffort::Legacy?std::min<std::uint64_t>(regularFront?96:48,room):room;
                    const auto extra=std::min(room,spacing+chord);
                    const auto units=aliasOptions.maximumWorkingBytes/1024;
                    const auto retained=std::uint64_t(points.size())+triangles.size()+
                        (normalResources?std::uint64_t(faces.size())+edges.size():0);
                    const auto memoryLimit=retained<units?(units-retained)/4:0;
                    chartNodeLimit=std::size_t(std::min({std::max<std::uint64_t>(4096,initial+extra),
                        memoryLimit,std::uint64_t((std::numeric_limits<int>::max)())}));
                    context.require(chartNodeLimit>=initial,OcctBoundaryReason::BudgetExceeded);
                    admitRepair(chartNodeLimit*2);
                }
                std::optional<detail::OcctChartTrim> chartTrim;
                if(spaceSurface&&aliasOptions.curvaturePointsPerTurn)
                    chartTrim.emplace(*model.face(face->tag),face->face,model.tolerance(2,face->tag));
                const auto refine=[&](auto& uv,auto& ids,auto& cells,const auto& protectedEdges){
                    const auto* native=model.face(face->tag);
                    const auto evaluate=[&](Vec2 p){return native->pointAt(p);};
                    const auto inside=[&](Vec2 p){
                        return chartTrim?chartTrim->contains(p):native->containsParam(p);
                    };
                    const auto originalNodes=spaceSurface?0:face->triangulation->NbNodes();
                    const auto append=[&](Vec2,Vec3 p){
                        return appendPoint(p,{2,face->tag},std::uint64_t(originalNodes)+uv.size()+1);
                    };
                    const auto deviation=[&](Vec3 p){
                        context.tick();
                        Vec2 uv;Vec3 closest;
                        if(native->projectPoint(p,uv,closest)||closestOnFace(face->face,p,closest))
                            return (p-closest).norm();
                        return std::numeric_limits<double>::infinity();
                    };
                    if(spaceSurface){
                        const auto chord=context.options.maximumChordError;
                        aliasReceipt.surfaceSpacing.emplace_back();
                        auto& spacingReport=aliasReceipt.surfaceSpacing.back();
                        spacingReport.faceTag=face->tag;
                        spacingReport.chartNodeLimit=chartNodeLimit;
                        std::optional<std::size_t> chordStart;
                        try {
                        detail::OcctStageTrace chartStep("surface_front_insertion",bodyTag);
                        if(!detail::spaceOcctChart(model,face->tag,uv,ids,cells,protectedEdges,
                            aliasOptions.curveMetricSizeCap,chord,position,append,tick,
                            [&](std::uint64_t cells){admitRepair(std::max<std::uint64_t>(4096,cells));},
                            regularFront,aliasOptions.curvaturePointsPerTurn?aliasOptions.curvaturePointsPerTurn:16,
                            chartTrim?&*chartTrim:nullptr,aliasOptions.sizeField,construction,
                            legacyConstruction?(std::numeric_limits<std::uint64_t>::max)():pointCapacity(cells.size()),
                            chartNodeLimit,&spacingReport)){
                            spacingReport.termination="invalid_chart";
                            return false;
                        }
                        // Keep a short, trim-constrained chord pass for edges
                        // whose local Delaunay cavities cannot accept a point.
                        chartStep.next("surface_chord_refinement");
                        auto chordLimit=construction.surfaceChordPointsPerFace?construction.surfaceChordPointsPerFace:
                            std::uint64_t(regularFront?96:48);
                        if(progressiveConstruction&&!construction.surfaceChordPointsPerFace)
                            chordLimit=pointCapacity(cells.size());
                        if(construction.effort!=CadConstructionEffort::Legacy&&!construction.surfaceChordPointsPerFace)
                            chordLimit=std::max(chordLimit,spacingReport.pointLimit/2);
                        if(!legacyConstruction)chordLimit=std::min(chordLimit,pointCapacity(cells.size()));
                        chordLimit=std::min<std::uint64_t>(chordLimit,(std::numeric_limits<unsigned>::max)());
                        spacingReport.chordPointLimit=chordLimit;
                        chordStart=uv.size();
                        const auto refined=detail::refineOcctUvTriangles(uv,ids,cells,protectedEdges,
                            chord*.96,unsigned(chordLimit),position,evaluate,inside,append,tick,deviation,chartNodeLimit);
                        spacingReport.chordAddedPoints=refined.added;
                        spacingReport.pendingEntries+=refined.pendingTriangles;
                        spacingReport.queueExhausted=spacingReport.queueExhausted&&refined.completed;
                        if(!refined.completed){
                            const bool limited=refined.added==chordLimit||uv.size()>=chartNodeLimit;
                            if(!spacingReport.budgetLimited)
                                spacingReport.termination=uv.size()>=chartNodeLimit?"chord_chart_node_limit":
                                    refined.added==chordLimit?"chord_point_limit":"chord_stalled";
                            spacingReport.budgetLimited=spacingReport.budgetLimited||limited;
                        }
                        chartStep.next("surface_vertex_relaxation");
                        if(regularFront)
                            detail::relaxOcctSurface(model,face->tag,uv,ids,cells,protectedEdges,chord,
                                position,[&](int id,Vec3 p){points[id].point=p;},tick,false,chartTrim?&*chartTrim:nullptr);
                        if(improveQuality){
                            detail::improveOcctSurfaceDiagonals(*native,uv,ids,cells,protectedEdges,chord,position,tick);
                            detail::relaxOcctSurface(model,face->tag,uv,ids,cells,protectedEdges,chord,
                                position,[&](int id,Vec3 p){points[id].point=p;},tick,true,chartTrim?&*chartTrim:nullptr);
                        }
                        return true; // Facet witnesses below make the final decision.
                        } catch(const Failure& failure) {
                            if(chordStart)spacingReport.chordAddedPoints=uv.size()-*chordStart;
                            spacingReport.budgetLimited=spacingReport.budgetLimited||failure.reason==OcctBoundaryReason::BudgetExceeded;
                            spacingReport.termination=failure.reason==OcctBoundaryReason::BudgetExceeded?"global_budget":
                                failure.reason==OcctBoundaryReason::Cancelled?"cancelled":"construction_refused";
                            spacingReport.queueExhausted=false;
                            throw;
                        } catch(const std::bad_alloc&) {
                            if(chordStart)spacingReport.chordAddedPoints=uv.size()-*chordStart;
                            spacingReport.termination="allocation_failed";
                            spacingReport.queueExhausted=false;
                            throw;
                        } catch(const Standard_Failure&) {
                            if(chordStart)spacingReport.chordAddedPoints=uv.size()-*chordStart;
                            spacingReport.termination="kernel_failure";
                            spacingReport.queueExhausted=false;
                            throw;
                        }
                    }
                    const auto refined=detail::refineOcctUvTriangles(uv,ids,cells,protectedEdges,
                        context.options.maximumChordError*.9,256,position,evaluate,inside,append,tick,deviation);
                    return refined.completed;
                };
                const auto complete=detail::rebuildOcctWireChart(model,face->tag,edges,face->interior,
                    position,tick,remaining,replacements,nullptr,spaceSurface||aliasReceipt.rephasedEdges!=0,refine,
                    // Parameter correspondence can need bounded recovery at
                    // any spacing. Curvature sampling is an independent size
                    // request, not permission to fit an interior pcurve point.
                    spaceSurface,chartNodeLimit,&aliasReceipt.witness);
                context.require(complete,OcctBoundaryReason::NonconformingEdge);
                chartTrace.next("surface_facet_witnesses");
                // Only new interior chart vertices used the fast trim
                // path. Check their final parameters on the original trim
                // before accepting the chart.
                std::set<int> checkedInterior;
                for(auto [t,uv]:replacements){
                    if(chartTrim)for(unsigned i=0;i<3;++i){
                        if(points[t[i]].classification.dim!=2||!checkedInterior.insert(t[i]).second)continue;
                        context.tick();
                        context.require(model.face(face->tag)->containsParam(uv[i]),OcctBoundaryReason::GeometryError);
                    }
                    if(face->face.Orientation()==TopAbs_REVERSED){
                        std::swap(t[1],t[2]);
                        std::swap(uv[1],uv[2]);
                    }
                    context.require(facetFits(*face,t,&uv),
                                    OcctBoundaryReason::GeometryError);
                    context.admit(points.size(),triangles.size()+1);
                    triangles.push_back(t);
                    faceTags.push_back(face->tag);
                }
                ++aliasReceipt.rebuiltCharts;
            }
            result.triangles=triangles.size();
        }
    }
    if(aliasOptions.maximumFlips&&!aliasOptions.faceTags.empty()){
        context.require(aliasOptions.maximumFlips<=256,OcctBoundaryReason::InvalidRequest);
        const auto count=std::uint64_t(points.size())+triangles.size()+aliasOptions.faceTags.size();
        context.require(count<=aliasOptions.maximumWorkingBytes/1024,OcctBoundaryReason::BudgetExceeded);
        aliasReceipt.maximumWorkingPayloadEstimate=count*1024;
        std::set<detail::OcctAliasEdge> protectedEdges;
        for(const auto& [tag,edge]:edges){(void)tag;for(std::size_t j=1;j<edge.ids.size();++j){context.tick();const auto a=edge.ids[j-1],b=edge.ids[j];protectedEdges.insert({std::min(a,b),std::max(a,b)});}}
        const auto position=[&](int i){return points[i].point;};
        const auto admission=[&](int faceTag,auto first,auto second){
            const auto face=detail::OcctModelAccess::face(model,faceTag);
            for(auto f:std::array{first,second}){const auto a=position(f[0]),b=position(f[1]),c=position(f[2]);
                for(auto p:std::array{(a+b+c)/3.,(a+b)*.5,(b+c)*.5,(c+a)*.5}){context.tick();Vec3 q;
                    if(!closestOnFace(face,p,q)||(p-q).norm()>context.options.maximumChordError+model.tolerance(2,faceTag))return false;}}
            return true;
        };
        aliasReceipt.flips=detail::repairOcctSurfaceAliases(triangles,faceTags,protectedEdges,aliasOptions.faceTags,aliasOptions.maximumFlips,position,[&](std::uint64_t n){context.tick(n);},admission);
        if(aliasReceipt.flips){
            sum=0;compensation=0;result.minimumTriangleArea=(std::numeric_limits<double>::infinity)();
            for(auto f:triangles){context.tick();const auto p=position(f[0]),q=position(f[1]),r=position(f[2]);
                result.minimumTriangleArea=std::min(result.minimumTriangleArea,(q-p).cross(r-p).norm()*.5);
                const double volume=(p-origin).dot((q-origin).cross(r-origin))/6.,next=sum+volume;
                compensation+=std::abs(sum)>=std::abs(volume)?(sum-next)+volume:(volume-next)+sum;sum=next;}
        }
    }
    // A collapsed parametric pole can leave two oppositely oriented
    // copies of the same physical triangle. Neutralizing the exact pair
    // on one CAD face preserves its oriented boundary and signed volume.
    detail::OcctStageTrace topologyTrace("surface_topology",bodyTag);
    std::map<int,unsigned> faceOccurrences;for(const auto& face:faces){context.tick();++faceOccurrences[face.tag];}
    const auto cancelledPairs=detail::cancelOcctOppositeFacePairs(triangles,faceTags,[&](std::uint64_t n){context.tick(n);},[&](int face){context.tick();return faceOccurrences.at(face)==1;});
    if(cancelledPairs||aliasReceipt.rebuiltCharts){
        std::vector<bool> used(points.size());for(auto t:triangles)for(int id:t)used[id]=true;
        std::vector<int> remap(points.size(),-1);std::vector<Sample> compact;compact.reserve(points.size());
        for(std::size_t i=0;i<points.size();++i){
            context.tick();
            const bool retired=retiredWireNodes.contains(int(i));
            if(!used[i]){
                context.require(points[i].classification.dim==2||retired,OcctBoundaryReason::NonconformingEdge);
                continue;
            }
            context.require(!retired,OcctBoundaryReason::NonconformingEdge);
            remap[i]=int(compact.size());compact.push_back(points[i]);
        }
        for(auto& t:triangles)for(int& id:t)id=remap[id];points.swap(compact);result.nodes=points.size();result.triangles=triangles.size();
        sum=0;compensation=0;result.minimumTriangleArea=(std::numeric_limits<double>::infinity)();
        for(auto t:triangles){context.tick();const auto p=points[t[0]].point,q=points[t[1]].point,r=points[t[2]].point;
            result.minimumTriangleArea=std::min(result.minimumTriangleArea,(q-p).cross(r-p).norm()*.5);
            const double volume=(p-origin).dot((q-origin).cross(r-origin))/6.,next=sum+volume;
            compensation+=std::abs(sum)>=std::abs(volume)?(sum-next)+volume:(volume-next)+sum;sum=next;
        }
    }
    result.signedSurfaceVolume = sum+compensation;
    result.cadVolume = model.volume(bodyTag);
    context.require(std::isfinite(result.signedSurfaceVolume) && result.signedSurfaceVolume > 0 &&
                    std::isfinite(result.cadVolume) && result.cadVolume > 0,OcctBoundaryReason::InvalidSurface);
    result.relativeVolumeError = std::abs(result.signedSurfaceVolume-result.cadVolume)/result.cadVolume;
    std::map<std::array<int,2>,std::vector<std::pair<int,int>>> incidence;
    std::vector<std::vector<int>> adjacency(triangles.size());
    std::vector<std::map<int,std::vector<int>>> links(points.size());
    for (int t=0; t<static_cast<int>(triangles.size()); ++t) {
        const auto triangle = triangles[t];
        for (int j=0; j<3; ++j) {
            context.tick(); int a = triangle[j], b = triangle[(j+1)%3]; const int sign = a<b?1:-1;
            if (a>b) std::swap(a,b); auto& uses = incidence[{a,b}];
            context.require(uses.size() < 2,OcctBoundaryReason::Nonmanifold); uses.push_back({t,sign});
            const int linkVertex = triangle[j], x = triangle[(j+1)%3], y = triangle[(j+2)%3];
            links[linkVertex][x].push_back(y); links[linkVertex][y].push_back(x);
        }
    }
    for (const auto& [edge,uses] : incidence) {
        (void)edge; context.tick();
        context.require(uses.size() == 2,OcctBoundaryReason::OpenBoundary);
        context.require(uses[0].second+uses[1].second == 0,OcctBoundaryReason::InvalidSurface);
        adjacency[uses[0].first].push_back(uses[1].first); adjacency[uses[1].first].push_back(uses[0].first);
    }
    result.closed = true; result.oriented = true;
    std::vector<bool> seen(triangles.size());
    for (int t=0; t<static_cast<int>(triangles.size()); ++t) {
        context.tick(); if (seen[t]) continue; ++result.connectedShells;
        std::vector<int> pending{t}; seen[t] = true;
        for (std::size_t i=0; i<pending.size(); ++i) for (int next : adjacency[pending[i]]) {
            context.tick(); if (!seen[next]) { seen[next] = true; pending.push_back(next); }
        }
    }
    for (const auto& link : links) {
        context.tick(); context.require(link.size() >= 3,OcctBoundaryReason::Nonmanifold);
        for (const auto& [linkVertex,neighbors] : link) {
            (void)linkVertex; context.tick();
            context.require(neighbors.size() == 2 && neighbors[0] != neighbors[1],OcctBoundaryReason::Nonmanifold);
        }
        std::set<int> visited; std::vector<int> pending{link.begin()->first}; visited.insert(pending[0]);
        for (std::size_t i=0; i<pending.size(); ++i) for (int next : link.at(pending[i])) {
            context.tick(); if (visited.insert(next).second) pending.push_back(next);
        }
        context.require(visited.size() == link.size(),OcctBoundaryReason::Nonmanifold);
    }
    result.vertexManifold = true;
    mesh::Mesh boundary; boundary.reserveAppend(points.size(),triangles.size(),3*triangles.size());
    result.nodeKeys.reserve(points.size());
    for (const auto& point : points) {
        context.tick(); boundary.addNode(point.point,point.classification); result.nodeKeys.push_back(point.key);
    }
    for (std::size_t i=0; i<triangles.size(); ++i) {
        context.tick(); const auto t = triangles[i];
        const std::array<Index,3> nodes{t[0],t[1],t[2]};
        boundary.addCell(mesh::CellType::Tri3,nodes,{2,faceTags[i]});
    }
    return boundary;
}
} // namespace

bool closestPointOnOcctFace(const OcctModel& model, int tag, const Vec3& sample, Vec3& closest) {
    if (!finite(sample) || !model.face(tag)) return false;
    try {
        const auto face = detail::OcctModelAccess::face(model,tag);
        if (auto* queries = detail::OcctQueryContext::active()) {
            // Only immutable model faces are cached. Native construction may
            // retessellate private copies that also use closestOnFace().
            if (queries->findTrimmedFacePoint(face,sample,closest)) return true;
            const auto& distance = queries->trimmedFaceDistance(face,gp_Pnt(sample.x,sample.y,sample.z));
            if (!closestSolution(distance,sample,closest)) return false;
            queries->rememberTrimmedFacePoint(face,sample,closest);
            return true;
        }
        return closestOnFace(face,sample,closest);
    }
    catch (const Standard_OutOfMemory&) { throw std::bad_alloc{}; }
    catch (const Standard_Failure&) { return false; }
}

struct detail::OcctBoundarySession::Impl {
    const OcctModel* model = nullptr;
    OcctBoundaryOptions options;
    BRepBuilderAPI_Copy copy;
    ShapeMap copiedEdges;
    std::vector<int> copiedEdgeTags;
    std::string sourceSnapshot;
    OcctBoundaryAuditPolicy auditPolicy = OcctBoundaryAuditPolicy::PerBoundary;
    OcctBoundaryFacetPolicy facetPolicy = OcctBoundaryFacetPolicy::Immediate;
    std::map<OcctBoundaryNodeKey,Vec3> nodeCoordinates;
    std::uint64_t work = 0;
};
detail::OcctBoundarySession::OcctBoundarySession(std::unique_ptr<Impl> impl) : impl_(std::move(impl)) {}
detail::OcctBoundarySession::~OcctBoundarySession() = default;
std::uint64_t detail::OcctBoundarySession::work() const noexcept { return impl_->work; }

std::unique_ptr<detail::OcctBoundarySession> detail::OcctBoundarySession::prepare(
    const OcctModel& model,const OcctBoundaryOptions& options,OcctBoundaryOutcome& preparation,
    OcctBoundaryAuditPolicy auditPolicy) noexcept {
    return prepare(model,options,preparation,auditPolicy,OcctBoundaryFacetPolicy::Immediate);
}

std::unique_ptr<detail::OcctBoundarySession> detail::OcctBoundarySession::prepare(
    const OcctModel& model, const OcctBoundaryOptions& options, OcctBoundaryOutcome& preparation,
    OcctBoundaryAuditPolicy auditPolicy, OcctBoundaryFacetPolicy facetPolicy) noexcept {
    return prepare(model,options,preparation,auditPolicy,facetPolicy,false);
}

std::unique_ptr<detail::OcctBoundarySession> detail::OcctBoundarySession::prepare(
    const OcctModel& model, const OcctBoundaryOptions& options, OcctBoundaryOutcome& preparation,
    OcctBoundaryAuditPolicy auditPolicy, OcctBoundaryFacetPolicy facetPolicy,
    bool controlDeflectionOnAllSurfaces) noexcept {
    return prepare(model,options,preparation,auditPolicy,facetPolicy,
                   controlDeflectionOnAllSurfaces,nullptr);
}

std::unique_ptr<detail::OcctBoundarySession> detail::OcctBoundarySession::prepare(
    const OcctModel& model, const OcctBoundaryOptions& options, OcctBoundaryOutcome& preparation,
    OcctBoundaryAuditPolicy auditPolicy, OcctBoundaryFacetPolicy facetPolicy,
    bool controlDeflectionOnAllSurfaces, const OcctBoundaryConstructionOptions* construction) noexcept {
    preparation = OcctBoundaryOutcome{};
    try {
        Context context{options,preparation};
        context.require((auditPolicy == OcctBoundaryAuditPolicy::PerBoundary || auditPolicy == OcctBoundaryAuditPolicy::Deferred) &&
                        (facetPolicy == OcctBoundaryFacetPolicy::Immediate || facetPolicy == OcctBoundaryFacetPolicy::DeferredCandidate) &&
                        !model.solids().empty() && std::isfinite(options.linearDeflection) &&
                        options.linearDeflection >= Precision::Confusion() && std::isfinite(options.angularDeflection) &&
                        options.angularDeflection >= Precision::Angular() && options.angularDeflection < std::numbers::pi &&
                        std::isfinite(options.maximumChordError) && options.maximumChordError > 0 &&
                        options.maximumNodes && options.maximumTriangles && options.maximumOutputBytes && options.maximumWork,
                        OcctBoundaryReason::InvalidRequest);
        context.tick();
        auto impl = std::make_unique<Impl>(); impl->model = &model; impl->options = options;
        impl->auditPolicy = auditPolicy;
        impl->facetPolicy = facetPolicy;
        impl->sourceSnapshot = snapshot(OcctModelAccess::shape(model));
        impl->copy.Perform(OcctModelAccess::shape(model),true,false);
        context.require(impl->copy.IsDone(),OcctBoundaryReason::KernelFailure); context.tick();
        occ::handle<Message_ProgressIndicator> progress = new CancellationProgress(options.cancellation);
        IMeshTools_Parameters parameters; parameters.Deflection = options.linearDeflection;
        parameters.Angle = options.angularDeflection; parameters.Relative = false; parameters.InParallel = false;
        parameters.EnableControlSurfaceDeflectionAllSurfaces = controlDeflectionOnAllSurfaces;
        if (construction && !construction->refinedEdges.empty()) {
            context.require(construction->maximumRefinedEdges && construction->maximumIncidentFaces &&
                            std::isfinite(construction->localLinearDeflection) &&
                            construction->localLinearDeflection >= Precision::Confusion() &&
                            construction->localLinearDeflection <= options.linearDeflection,
                            OcctBoundaryReason::InvalidRequest);
            context.require(construction->refinedEdges.size() <= construction->maximumRefinedEdges,
                            OcctBoundaryReason::BudgetExceeded);
            std::set<int> seen;
            std::vector<TopoDS_Shape> selected;
            for (const int tag : construction->refinedEdges) {
                context.tick();
                context.require(model.edge(tag) && seen.insert(tag).second,OcctBoundaryReason::InvalidRequest);
                const auto edge = impl->copy.ModifiedShape(OcctModelAccess::edge(model,tag));
                context.require(!edge.IsNull(),OcctBoundaryReason::UnsupportedTopology);
                selected.push_back(edge);
            }
            TopoDS_Compound patch; BRep_Builder builder; builder.MakeCompound(patch);
            unsigned incidentFaces = 0;
            for (const auto* face : model.faces()) {
                context.tick();
                const auto copied = impl->copy.ModifiedShape(OcctModelAccess::face(model,face->tag()));
                context.require(!copied.IsNull(),OcctBoundaryReason::UnsupportedTopology);
                bool incident = false;
                for (TopExp_Explorer edges(copied,TopAbs_EDGE);edges.More();edges.Next())
                    for (const auto& edge : selected) {
                        context.tick();
                        incident |= edges.Current().IsSame(edge);
                    }
                if (!incident) continue;
                context.require(++incidentFaces <= construction->maximumIncidentFaces,
                                OcctBoundaryReason::BudgetExceeded);
                builder.Add(patch,copied); context.tick();
            }
            context.require(incidentFaces != 0,OcctBoundaryReason::UnsupportedTopology);
            auto localParameters = parameters;
            localParameters.Deflection = construction->localLinearDeflection;
            localParameters.AllowQualityDecrease = false;
            BRepMesh_IncrementalMesh local(patch,localParameters,progress->Start()); context.tick();
            context.require(local.IsDone(),OcctBoundaryReason::KernelFailure);
            // Retain the refined shared handles during whole-compound meshing.
            parameters.AllowQualityDecrease = false;
        }
        BRepMesh_IncrementalMesh mesher(impl->copy.Shape(),parameters,progress->Start()); context.tick();
        context.require(mesher.IsDone(),OcctBoundaryReason::KernelFailure);
        // The whole copied compound retains exact edge identity across bodies.
        // Build its tag lookup once instead of scanning all source edges for
        // every extracted piece. Face occurrences still choose their owner
        // orientation and seam branch locally.
        for (const auto* edge : model.edges()) {
            context.tick();
            const auto original = OcctModelAccess::edge(model,edge->tag());
            const auto copied = impl->copy.ModifiedShape(original);
            context.require(!copied.IsNull(),OcctBoundaryReason::UnsupportedTopology);
            const int index = impl->copiedEdges.Add(copied);
            context.require(index == static_cast<int>(impl->copiedEdgeTags.size())+1,
                            OcctBoundaryReason::UnsupportedTopology);
            impl->copiedEdgeTags.push_back(edge->tag());
        }
        context.require(impl->sourceSnapshot == snapshot(OcctModelAccess::shape(model)),OcctBoundaryReason::SourceChanged);
        preparation.sourceModelUnmodified = true; preparation.reason = OcctBoundaryReason::None;
        impl->work = preparation.work;
        return std::unique_ptr<OcctBoundarySession>(new OcctBoundarySession(std::move(impl)));
    } catch (const Failure& failure) { preparation.reason = failure.reason; }
      catch (const Standard_OutOfMemory&) { preparation.reason = OcctBoundaryReason::AllocationFailed; }
      catch (const std::bad_alloc&) { preparation.reason = OcctBoundaryReason::AllocationFailed; }
      catch (const Standard_Failure&) { preparation.reason = OcctBoundaryReason::KernelFailure; }
      catch (...) { preparation.reason = OcctBoundaryReason::InternalError; }
    return nullptr;
}
OcctBoundaryOutcome detail::OcctBoundarySession::boundary(int tag, mesh::Mesh& destination,
                                                        std::uint64_t maximumAdditionalWork) noexcept {
    OcctBoundaryAliasOptions disabled;disabled.maximumFlips=0;OcctBoundaryAliasReceipt ignored;
    return boundary(tag,destination,disabled,ignored,maximumAdditionalWork);
}
OcctBoundaryOutcome detail::OcctBoundarySession::boundary(int tag, mesh::Mesh& destination,
    const OcctBoundaryAliasOptions& aliases,OcctBoundaryAliasReceipt& aliasReceipt,
    std::uint64_t maximumAdditionalWork) noexcept {
    aliasReceipt={};
    OcctBoundaryOutcome result; result.solidTag = tag;
    try {
        auto options = impl_->options;
        options.maximumWork = std::min(options.maximumWork-impl_->work,maximumAdditionalWork);
        Context context{options,result};
        context.require(impl_->model->solid(tag) != nullptr,OcctBoundaryReason::InvalidRequest); context.tick();
        context.require(aliases.maximumFlips<=256,OcctBoundaryReason::InvalidRequest);
        context.require(options.construction.resourcePolicy==CadResourcePolicy::Normal||
                        aliases.maximumRebuiltCharts<=64,OcctBoundaryReason::InvalidRequest);
        context.require(std::isfinite(aliases.curveMetricSizeCap)&&aliases.curveMetricSizeCap>=0,
                        OcctBoundaryReason::InvalidRequest);
        for(int face:aliases.faceTags){context.tick();context.require(impl_->model->face(face)!=nullptr,OcctBoundaryReason::InvalidRequest);}
        auto boundary = extract(*impl_->model,tag,&impl_->copy,impl_->copiedEdges,impl_->copiedEdgeTags,context,impl_->facetPolicy,aliases,aliasReceipt); context.tick();
        if (impl_->auditPolicy == OcctBoundaryAuditPolicy::PerBoundary) {
            context.require(impl_->sourceSnapshot == snapshot(OcctModelAccess::shape(*impl_->model)),OcctBoundaryReason::SourceChanged);
            result.sourceModelUnmodified = true;
        }
        // Shared native entities are joined by session identity only. Even a
        // tolerance-close coordinate mismatch is refused before publication.
        for (Index i=0; i<boundary.nodeCount(); ++i) {
            context.tick(); const auto& key = result.nodeKeys[static_cast<std::size_t>(i)]; const auto p = boundary.node(i);
            auto [found,inserted] = impl_->nodeCoordinates.emplace(key,p);
            const auto q = found->second;
            context.require(inserted ||
                (std::bit_cast<std::uint64_t>(p.x) == std::bit_cast<std::uint64_t>(q.x) &&
                 std::bit_cast<std::uint64_t>(p.y) == std::bit_cast<std::uint64_t>(q.y) &&
                 std::bit_cast<std::uint64_t>(p.z) == std::bit_cast<std::uint64_t>(q.z)),
                 OcctBoundaryReason::NonconformingEdge);
        }
        context.tick(); destination.swap(boundary);
        result.published = true; result.destinationUnchanged = false; result.reason = OcctBoundaryReason::None;
    } catch (const Failure& failure) { result.reason = failure.reason; }
      catch (const Standard_OutOfMemory&) { result.reason = OcctBoundaryReason::AllocationFailed; }
      catch (const std::bad_alloc&) { result.reason = OcctBoundaryReason::AllocationFailed; }
      catch (const Standard_Failure&) { result.reason = OcctBoundaryReason::KernelFailure; }
      catch (...) { result.reason = OcctBoundaryReason::InternalError; }
    impl_->work += result.work;
    if (!result.published) result.nodeKeys.clear();
    return result;
}
OcctBoundaryOutcome detail::buildSpacedOcctBoundary(const OcctModel& model,int tag,
    const OcctBoundaryOptions& options,double targetSize,std::uint64_t maximumWorkingBytes,
    mesh::Mesh& destination,OcctBoundaryAliasReceipt& receipt,bool regularFront,bool improveQuality,
    unsigned curvaturePointsPerTurn,const OcctSpatialSizing* sizeField) noexcept {
    OcctBoundaryOutcome result;
    result.solidTag=tag;
    receipt={};
    try {
        Context context{options,result};
        context.tick();
        context.require(model.solid(tag)&&model.solids().size()==1&&
                        std::isfinite(targetSize)&&targetSize>0&&std::isfinite(options.maximumChordError)&&
                        options.maximumChordError>0,OcctBoundaryReason::InvalidRequest);
        const auto faceCount=model.faceUses(tag).size();
        const bool normal=options.construction.resourcePolicy==CadResourcePolicy::Normal;
        const auto chartAllowance=normal?std::min({options.maximumWork-result.work,maximumWorkingBytes/1024,
            std::uint64_t((std::numeric_limits<unsigned>::max)())}):std::uint64_t(curvaturePointsPerTurn?128:32);
        context.require(faceCount<=chartAllowance&&maximumWorkingBytes>=4096*1024,
                        OcctBoundaryReason::BudgetExceeded);
        const unsigned maximumCharts=normal?static_cast<unsigned>(faceCount):(curvaturePointsPerTurn?128:32);
        const auto& shape=OcctModelAccess::shape(model);
        const auto snapshotLimit=std::min<std::uint64_t>(64ULL*1024*1024,maximumWorkingBytes/8);
        const auto before=facetSnapshot(shape,snapshotLimit);
        OcctBoundaryAliasOptions spacing;
        spacing.maximumFlips=0;
        spacing.maximumRebuiltCharts=maximumCharts;
        spacing.maximumWorkingBytes=maximumWorkingBytes-2*snapshotLimit;
        spacing.curveMetricSizeCap=targetSize;
        spacing.curvaturePointsPerTurn=curvaturePointsPerTurn;
        spacing.sizeField=sizeField;
        auto boundary=extract(model,tag,nullptr,{},{},context,OcctBoundaryFacetPolicy::Immediate,spacing,receipt,
                              regularFront,regularFront&&improveQuality);
        context.tick();
        context.require(before==facetSnapshot(shape,snapshotLimit),OcctBoundaryReason::SourceChanged);
        result.sourceModelUnmodified=true;
        destination.swap(boundary);
        result.published=true;
        result.destinationUnchanged=false;
        result.reason=OcctBoundaryReason::None;
    } catch(const Failure& failure) { result.reason=failure.reason; }
      catch(const Standard_OutOfMemory&) { result.reason=OcctBoundaryReason::AllocationFailed; }
      catch(const std::bad_alloc&) { result.reason=OcctBoundaryReason::AllocationFailed; }
      catch(const Standard_Failure&) { result.reason=OcctBoundaryReason::KernelFailure; }
      catch(...) { result.reason=OcctBoundaryReason::InternalError; }
    if(!result.published)result.nodeKeys.clear();
    return result;
}
OcctBoundaryOutcome detail::OcctBoundarySession::auditSource(std::uint64_t maximumAdditionalWork) noexcept {
    OcctBoundaryOutcome result;
    try {
        auto options = impl_->options;
        options.maximumWork = std::min(options.maximumWork-impl_->work,maximumAdditionalWork);
        Context context{options,result}; context.tick();
        context.require(impl_->sourceSnapshot == snapshot(OcctModelAccess::shape(*impl_->model)),OcctBoundaryReason::SourceChanged);
        context.tick(); result.sourceModelUnmodified = true; result.reason = OcctBoundaryReason::None;
    } catch (const Failure& failure) { result.reason = failure.reason; }
      catch (const Standard_OutOfMemory&) { result.reason = OcctBoundaryReason::AllocationFailed; }
      catch (const std::bad_alloc&) { result.reason = OcctBoundaryReason::AllocationFailed; }
      catch (const Standard_Failure&) { result.reason = OcctBoundaryReason::KernelFailure; }
      catch (...) { result.reason = OcctBoundaryReason::InternalError; }
    impl_->work += result.work;
    return result;
}
OcctBoundaryOutcome buildOcctBoundary(const OcctModel& model, int tag, const OcctBoundaryOptions& options,
                                     mesh::Mesh& destination) noexcept {
    if (!model.solid(tag)) { OcctBoundaryOutcome result; result.solidTag = tag; result.reason = OcctBoundaryReason::InvalidRequest; return result; }
    OcctBoundaryOutcome preparation;
    auto session = detail::OcctBoundarySession::prepare(model,options,preparation);
    if (!session) { preparation.solidTag = tag; return preparation; }
    auto result = session->boundary(tag,destination); result.work += preparation.work; return result;
}
OcctBoundaryOutcome detail::buildOcctBoundaryCandidate(const OcctModel& model,int tag,
    const OcctBoundaryOptions& options,mesh::Mesh& destination) noexcept {
    if(!model.solid(tag)){OcctBoundaryOutcome result;result.solidTag=tag;result.reason=OcctBoundaryReason::InvalidRequest;return result;}
    OcctBoundaryOutcome preparation;
    auto session=OcctBoundarySession::prepare(model,options,preparation,
        OcctBoundaryAuditPolicy::PerBoundary,OcctBoundaryFacetPolicy::DeferredCandidate);
    if(!session){preparation.solidTag=tag;return preparation;}
    auto result=session->boundary(tag,destination);result.work+=preparation.work;return result;
}
namespace {
class BoundedFacetSnapshot final : public std::streambuf {
public:
    explicit BoundedFacetSnapshot(std::uint64_t maximum):maximum_(maximum){}
    std::string value;bool exceeded=false;
protected:
    std::streamsize xsputn(const char* data,std::streamsize size)override {
        if(size<0||static_cast<std::uint64_t>(size)>maximum_-value.size()){exceeded=true;return 0;}
        value.append(data,static_cast<std::size_t>(size));return size;
    }
    int_type overflow(int_type c)override {
        if(traits_type::eq_int_type(c,traits_type::eof()))return traits_type::not_eof(c);
        const char byte=traits_type::to_char_type(c);return xsputn(&byte,1)==1?c:traits_type::eof();
    }
private:std::uint64_t maximum_;
};
std::string facetSnapshot(const TopoDS_Shape& shape,std::uint64_t maximum) {
    BoundedFacetSnapshot buffer(maximum);std::ostream stream(&buffer);
    stream.exceptions(std::ios::badbit|std::ios::failbit);
    try{BRepTools::Write(shape,stream);}
    catch(const std::ios_base::failure&){throw Failure{buffer.exceeded?OcctBoundaryReason::BudgetExceeded:OcctBoundaryReason::KernelFailure};}
    if(buffer.exceeded)throw Failure{OcctBoundaryReason::BudgetExceeded};
    if(!stream.good())throw Failure{OcctBoundaryReason::KernelFailure};
    return std::move(buffer.value);
}
}
detail::OcctBoundaryFacetAuditOutcome detail::auditOcctBoundaryFacets(const OcctModel& model,
    int tag,const mesh::Mesh& boundary,const OcctBoundaryOptions& options,
    std::uint64_t maximumSnapshotBytes) noexcept {
    OcctBoundaryFacetAuditOutcome result;OcctBoundaryOutcome counters;
    try {
        Context context{options,counters};
        context.require(model.solid(tag)&&boundary.nodeCount()&&boundary.cellCount()&&
            maximumSnapshotBytes&&options.maximumWork&&std::isfinite(options.maximumChordError)&&
            options.maximumChordError>0,OcctBoundaryReason::InvalidRequest);
        context.admit(boundary.nodeCount(),boundary.cellCount());
        context.tick();const auto before=facetSnapshot(OcctModelAccess::shape(model),maximumSnapshotBytes);
        std::set<int> ownedFaces;
        for(const auto& use:model.faceUses(tag)){context.tick();ownedFaces.insert(use.faceTag);}
        for(Index i=0;i<boundary.cellCount();++i) {
            context.tick();const auto& cell=boundary.cell(i);
            context.require(cell.type==mesh::CellType::Tri3&&cell.onEntity.dim==2&&
                ownedFaces.contains(cell.onEntity.tag),OcctBoundaryReason::InvalidSurface);
            const auto* ids=boundary.cellNodes(i);const auto p=boundary.node(ids[0]),q=boundary.node(ids[1]),r=boundary.node(ids[2]);
            context.require(finite(p)&&finite(q)&&finite(r),OcctBoundaryReason::GeometryError);
            const auto* face=model.face(cell.onEntity.tag);
            const auto tolerance=model.tolerance(2,cell.onEntity.tag);
            context.require(face&&std::isfinite(tolerance)&&tolerance>=0,OcctBoundaryReason::InvalidSurface);
            for(const auto sample:std::array{(p+q+r)/3.,(p+q)/2.,(q+r)/2.,(r+p)/2.}) {
                context.tick();++result.samples;Vec2 uv;Vec3 closest;
                double error=(std::numeric_limits<double>::infinity)();
                if(face->projectPoint(sample,uv,closest))error=(sample-closest).norm();
                if(!std::isfinite(error)||error>options.maximumChordError) {
                    ++result.closestTrimmedFaceSamples;
                    context.require(closestOnFace(OcctModelAccess::face(model,cell.onEntity.tag),sample,closest),OcctBoundaryReason::GeometryError);
                    error=(sample-closest).norm();
                }
                context.require(std::isfinite(error),OcctBoundaryReason::GeometryError);
                if(error>result.maximumFacetSampleError){result.maximumFacetSampleError=error;
                    result.maximumErrorFace=cell.onEntity.tag;result.maximumErrorSample=sample;result.maximumErrorClosest=closest;}
                context.require(error<=options.maximumChordError+tolerance,OcctBoundaryReason::GeometryError);
            }
        }
        context.tick();context.require(before==facetSnapshot(OcctModelAccess::shape(model),maximumSnapshotBytes),OcctBoundaryReason::SourceChanged);
        result.sourceModelUnmodified=true;result.completed=true;result.reason=OcctBoundaryReason::None;
    } catch(const Failure& failure){result.reason=failure.reason;}
      catch(const Standard_OutOfMemory&){result.reason=OcctBoundaryReason::AllocationFailed;}
      catch(const std::bad_alloc&){result.reason=OcctBoundaryReason::AllocationFailed;}
      catch(const Standard_Failure&){result.reason=OcctBoundaryReason::KernelFailure;}
      catch(...){result.reason=OcctBoundaryReason::InternalError;}
    result.work=counters.work;return result;
}
const char* occtBoundaryReasonName(OcctBoundaryReason reason) noexcept {
    switch (reason) {
    case OcctBoundaryReason::None:return "none";
    case OcctBoundaryReason::InvalidRequest:return "invalid_request";
    case OcctBoundaryReason::UnsupportedTopology:return "unsupported_topology";
    case OcctBoundaryReason::MissingTriangulation:return "missing_triangulation";
    case OcctBoundaryReason::NonconformingEdge:return "nonconforming_edge";
    case OcctBoundaryReason::GeometryError:return "geometry_error";
    case OcctBoundaryReason::InvalidSurface:return "invalid_surface";
    case OcctBoundaryReason::OpenBoundary:return "open_boundary";
    case OcctBoundaryReason::Nonmanifold:return "nonmanifold";
    case OcctBoundaryReason::BudgetExceeded:return "budget_exceeded";
    case OcctBoundaryReason::Cancelled:return "cancelled";
    case OcctBoundaryReason::AllocationFailed:return "allocation_failed";
    case OcctBoundaryReason::KernelFailure:return "kernel_failure";
    case OcctBoundaryReason::SourceChanged:return "source_changed";
    case OcctBoundaryReason::InternalError:return "internal_error";
    }
    return "internal_error";
}
} // namespace paacs::geo
