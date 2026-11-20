// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "test_harness.hpp"

#include "paacs/algo/MultiRegionMesher.hpp"
#include "paacs/algo/Meshers.hpp"
#include "paacs/algo/ProtectedFeatureRefinement.hpp"
#if defined(PAACS_T5_SCALING_MAIN) || defined(PAACS_T2C_SCALING_MAIN)
#include "paacs/algo/ProductionSurfaceMesher.hpp"
#include "paacs/algo/TetCandidate.hpp"
#endif

#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/mesh/Quality.hpp"
#include "paacs/sizing/SizeField.hpp"
#include "../src/algo/ProductionSurfaceMesherDetail.hpp"
#include "../src/algo/T2CSurfaceCache.hpp"

#include <array>
#include <atomic>
#include <chrono>
#include <cfenv>
#if defined(PAACS_T5_SCALING_MAIN) || defined(PAACS_T2C_SCALING_MAIN)
#include <cstdlib>
#endif
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string_view>
#include <thread>
#include <tuple>

namespace {
using namespace paacs;

class Ownership final : public algo::VolumeBoundaryOwnershipProvider {
public:
    mutable algo::VolumeBoundaryProviderStamp current{15001,1};
    std::map<algo::FacePersistentId,std::vector<bnd::VolumeRegionUse>> rows;
    mutable std::size_t stampCalls=0,ownershipCalls=0,envelopeCalls=0;
    std::size_t mutateAtStampCall=0;
    algo::VolumeBoundaryProviderStamp stamp() const noexcept override {
        ++stampCalls;
        if(mutateAtStampCall!=0&&stampCalls==mutateAtStampCall)
            ++current.revision;
        return current;
    }
    algo::VolumeFaceOwnership ownership(
            algo::FacePersistentId face) const override {
        ++ownershipCalls;
        algo::VolumeFaceOwnership out;out.providerStamp=current;
        out.facePersistentId=face;
        const auto found=rows.find(face);
        if(found!=rows.end())out.owners=found->second;
        return out;
    }
    algo::ProviderRecoveryEnvelope recoveryEnvelope(
            std::span<const Vec3>) const override {
        ++envelopeCalls;
        algo::ProviderRecoveryEnvelope out;out.providerStamp=current;
        out.availability=
            algo::VolumeProviderAvailability::RecoveryEnvelopeUnavailable;
        return out;
    }
};

class CurvePlacement final : public algo::ProtectedFeaturePlacementProvider {
public:
    enum class FacetFault : std::uint8_t {
        None,WrongEntity,WrongParameter,OffEntity,Rejected,Throw,BadAlloc
    };
    mutable algo::ProtectedFeatureProviderStamp current{15002,1};
    mutable std::size_t stampCalls=0,placementCalls=0;
    mutable std::size_t facetPlacementCalls=0;
    mutable std::vector<algo::ProtectedFeaturePlacementRequest> requests;
    std::size_t mutateAtStampCall=0;
    algo::ProtectedFeaturePlacementStatus status =
        algo::ProtectedFeaturePlacementStatus::Placed;
    bool throwPlacement=false;
    bool allocatePlacement=false;
    bool wrongEntity=false;
    bool wrongParameter=false;
    bool offEntity=false;
    std::size_t facetFaultAt=0;
    FacetFault facetFault=FacetFault::None;

    algo::ProtectedFeatureProviderStamp stamp() const noexcept override {
        ++stampCalls;
        if(mutateAtStampCall!=0&&stampCalls==mutateAtStampCall)
            ++current.revision;
        return current;
    }
    algo::ProtectedFeaturePlacement placeOnEdge(
            const geo::GEdge& edge,
            const algo::ProtectedFeaturePlacementRequest& request) const override {
        ++placementCalls;
        requests.push_back(request);
        if(allocatePlacement)throw std::bad_alloc{};
        if(throwPlacement)throw std::runtime_error("placement failure");
        auto entity=request.entityStamp;
        auto identity=request.requestedIdentity;
        auto coordinate=edge.pointAt(request.curveParameter);
        if(wrongEntity)++entity.revision;
        if(wrongParameter)identity.parameterUBits^=UINT64_C(1);
        if(offEntity)coordinate.x+=1.0;
        return {current,entity,identity,status,coordinate};
    }
    algo::ProtectedFeaturePlacement placeOnFace(
            const geo::GFace& face,
            const algo::ProtectedFeaturePlacementRequest& request) const override {
        ++placementCalls;
        requests.push_back(request);
        const bool facet=request.kind==
            algo::ProtectedFeatureEntityKind::FacetInterior;
        if(facet)++facetPlacementCalls;
        const bool facetFaultActive=facet&&facetFaultAt!=0&&
            facetPlacementCalls==facetFaultAt;
        if(allocatePlacement||
           (facetFaultActive&&facetFault==FacetFault::BadAlloc))
            throw std::bad_alloc{};
        if(throwPlacement||
           (facetFaultActive&&facetFault==FacetFault::Throw))
            throw std::runtime_error("placement failure");
        auto entity=request.entityStamp;
        auto identity=request.requestedIdentity;
        auto coordinate=face.pointAt(request.faceParameter);
        if(wrongEntity||
           (facetFaultActive&&facetFault==FacetFault::WrongEntity))
            ++entity.revision;
        if(wrongParameter||
           (facetFaultActive&&facetFault==FacetFault::WrongParameter))
            identity.parameterUBits^=UINT64_C(1);
        if(offEntity||
           (facetFaultActive&&facetFault==FacetFault::OffEntity))
            coordinate.x+=1.0;
        const auto placementStatus=facetFaultActive&&
                facetFault==FacetFault::Rejected
            ? algo::ProtectedFeaturePlacementStatus::Rejected:status;
        return {current,entity,identity,placementStatus,coordinate};
    }
};

class InterfaceEdge final : public geo::GEdge {
public:
    InterfaceEdge(int tag,const geo::GVertex* a,const geo::GVertex* b)
        : GEdge(tag),a_(a),b_(b) {}
    Vec3 pointAt(double t) const override {
        return a_->position()+(b_->position()-a_->position())*t;
    }
    Vec3 tangentAt(double) const override {
        return (b_->position()-a_->position()).normalized();
    }
    double length() const override {
        return (b_->position()-a_->position()).norm();
    }
    int startVertexTag() const override{return a_->tag();}
    int endVertexTag() const override{return b_->tag();}
private:
    const geo::GVertex* a_=nullptr;const geo::GVertex* b_=nullptr;
};

class InterfaceFace final : public geo::GFace {
public:
    InterfaceFace(int tag,Vec3 origin,Vec3 u,Vec3 v,std::vector<int> loop)
        : GFace(tag),origin_(origin),u_(u),v_(v),loop_(std::move(loop)) {}
    Vec3 pointAt(const Vec2& p) const override{return origin_+u_*p.u+v_*p.v;}
    Vec2 paramAt(const Vec3& p) const override {
        const auto d=p-origin_;
        const double uu=u_.norm2(),uv=u_.dot(v_),vv=v_.norm2();
        const double determinant=uu*vv-uv*uv;
        return {(d.dot(u_)*vv-d.dot(v_)*uv)/determinant,
                (d.dot(v_)*uu-d.dot(u_)*uv)/determinant};
    }
    Vec3 normalAt(const Vec2&) const override{return u_.cross(v_).normalized();}
    void paramBounds(Vec2& lo,Vec2& hi) const override{lo={0,0};hi={1,1};}
    std::vector<int> boundaryEdgeTags() const override{return loop_;}
    double curvatureAt(const Vec2&) const override{return 0;}
    Vec2 metricScale() const override{return {u_.norm(),v_.norm()};}
private:
    Vec3 origin_{},u_{},v_{};std::vector<int> loop_;
};

class LayeredModel final : public geo::Model {
public:
    explicit LayeredModel(int regionCount=2,bool independentPairs=false,
                          double crossSectionHalfExtent=1.0,
                          double layerSpacing=1.0,
                          bool reversePresentation=false)
        :regionCount_(regionCount),reversePresentation_(reversePresentation) {
        if(independentPairs) {
            buildIndependentPairs(regionCount/2);
            return;
        }
        for(int x=0;x<=regionCount_;++x)
            for(int z=0;z<2;++z)for(int y=0;y<2;++y)
            {
                Vec3 point{
                    layerSpacing*static_cast<double>(x),
                    crossSectionHalfExtent*static_cast<double>(2*y-1),
                    crossSectionHalfExtent*static_cast<double>(2*z-1)};
                vertices_.push_back(std::make_unique<geo::GVertex>(
                    static_cast<int>(vertices_.size()+1),point));
            }
        addFace({v(0,0,0),v(0,0,1),v(0,1,1),v(0,1,0)},{{0,1}});
        addFace({v(regionCount_,0,0),v(regionCount_,1,0),
                 v(regionCount_,1,1),v(regionCount_,0,1)},
                {{regionCount_-1,1}});
        for(int interval=0;interval<regionCount_;++interval) {
            const int lo=interval,hi=interval+1;
            addFace({v(lo,0,0),v(hi,0,0),v(hi,0,1),v(lo,0,1)},
                    {{interval,1}});
            addFace({v(lo,1,0),v(lo,1,1),v(hi,1,1),v(hi,1,0)},
                    {{interval,1}});
            addFace({v(lo,0,0),v(lo,1,0),v(hi,1,0),v(hi,0,0)},
                    {{interval,1}});
            addFace({v(lo,0,1),v(hi,0,1),v(hi,1,1),v(lo,1,1)},
                    {{interval,1}});
        }
        for(int plane=1;plane<regionCount_;++plane)
            addFace({v(plane,0,0),v(plane,1,0),v(plane,1,1),v(plane,0,1)},
                    {{plane-1,1},{plane,-1}});
    }
    const std::vector<bnd::VolumeRegionUse>& owners(int tag) const {
        return owners_.at(static_cast<std::size_t>(tag-1));
    }
    std::vector<const geo::GVertex*> vertices() const override {
        std::vector<const geo::GVertex*> out;
        for(const auto& x:vertices_)out.push_back(x.get());return out;
    }
    std::vector<const geo::GEdge*> edges() const override {
        std::vector<const geo::GEdge*> out;
        for(const auto& x:edges_)out.push_back(x.get());return out;
    }
    std::vector<const geo::GFace*> faces() const override {
        std::vector<const geo::GFace*> out;
        for(const auto& x:faces_)out.push_back(x.get());
        if(reversePresentation_)std::reverse(out.begin(),out.end());
        return out;
    }
    std::vector<const geo::GSolid*> solids() const override{return {};}
    const geo::GVertex* vertex(int tag) const override {
        for(const auto& x:vertices_)if(x->tag()==tag)return x.get();return nullptr;
    }
    const geo::GEdge* edge(int tag) const override {
        for(const auto& x:edges_)if(x->tag()==tag)return x.get();return nullptr;
    }
    const geo::GFace* face(int tag) const override {
        for(const auto& x:faces_)if(x->tag()==tag)return x.get();return nullptr;
    }
    const geo::GSolid* solid(int) const override{return nullptr;}
private:
    std::size_t v(int x,int y,int z) const noexcept {
        return static_cast<std::size_t>(x*4+z*2+y);
    }
    void addFace(std::array<std::size_t,4> loop,
                 std::vector<bnd::VolumeRegionUse> owners) {
        std::vector<int> faceEdges;
        for(std::size_t i=0;i<4;++i) {
            auto key=std::array<std::size_t,2>{loop[i],loop[(i+1)%4]};
            if(key[1]<key[0])std::swap(key[0],key[1]);
            auto found=edgeTags_.find(key);
            if(found==edgeTags_.end()) {
                const int tag=static_cast<int>(edges_.size()+1);
                edges_.push_back(std::make_unique<InterfaceEdge>(
                    tag,vertices_[key[0]].get(),vertices_[key[1]].get()));
                found=edgeTags_.emplace(key,tag).first;
            }
            faceEdges.push_back(found->second);
        }
        const auto p0=vertices_[loop[0]]->position();
        const auto p1=vertices_[loop[1]]->position();
        const auto p3=vertices_[loop[3]]->position();
        faces_.push_back(std::make_unique<InterfaceFace>(
            static_cast<int>(faces_.size()+1),p0,p1-p0,p3-p0,
            std::move(faceEdges)));
        owners_.push_back(std::move(owners));
    }
    void buildIndependentPairs(int pairCount) {
        int side=1;
        while(side*side*side<pairCount)++side;
        for(int pair=0;pair<pairCount;++pair) {
            const int px=pair%side;
            const int py=(pair/side)%side;
            const int pz=pair/(side*side);
            const double offset=std::ldexp(static_cast<double>(pair),-20);
            const Vec3 center{static_cast<double>(4*px)+1.0+offset,
                              static_cast<double>(4*py)+2*offset,
                              static_cast<double>(4*pz)+3*offset};
            const auto perturb=[pair](Vec3 p) noexcept {
                const double a=static_cast<double>(pair%7-3)/128.0;
                const double b=static_cast<double>((pair/7)%7-3)/128.0;
                const double c=static_cast<double>((pair/49)%7-3)/128.0;
                return Vec3{p.x+a*p.y+b*p.z,p.y+c*p.z,p.z};
            };
            constexpr double half=0.99;
            for(int plane=0;plane<3;++plane)
                for(int z=0;z<2;++z)for(int y=0;y<2;++y)
                    vertices_.push_back(std::make_unique<geo::GVertex>(
                        static_cast<int>(vertices_.size()+1),
                        center+perturb(Vec3{
                            half*static_cast<double>(plane-1),
                            half*static_cast<double>(2*y-1),
                            half*static_cast<double>(2*z-1)})));
        }
        const auto pv=[](int pair,int plane,int y,int z) noexcept {
            return static_cast<std::size_t>(
                pair*12+plane*4+z*2+y);
        };
        for(int pair=0;pair<pairCount;++pair) {
            const int leftRegion=2*pair,rightRegion=leftRegion+1;
            addFace({pv(pair,0,0,0),pv(pair,0,0,1),
                     pv(pair,0,1,1),pv(pair,0,1,0)},
                    {{leftRegion,1}});
            addFace({pv(pair,2,0,0),pv(pair,2,1,0),
                     pv(pair,2,1,1),pv(pair,2,0,1)},
                    {{rightRegion,1}});
            for(int interval=0;interval<2;++interval) {
                const int region=leftRegion+interval;
                addFace({pv(pair,interval,0,0),
                         pv(pair,interval+1,0,0),
                         pv(pair,interval+1,0,1),
                         pv(pair,interval,0,1)},{{region,1}});
                addFace({pv(pair,interval,1,0),
                         pv(pair,interval,1,1),
                         pv(pair,interval+1,1,1),
                         pv(pair,interval+1,1,0)},{{region,1}});
                addFace({pv(pair,interval,0,0),
                         pv(pair,interval,1,0),
                         pv(pair,interval+1,1,0),
                         pv(pair,interval+1,0,0)},{{region,1}});
                addFace({pv(pair,interval,0,1),
                         pv(pair,interval+1,0,1),
                         pv(pair,interval+1,1,1),
                         pv(pair,interval,1,1)},{{region,1}});
            }
            addFace({pv(pair,1,0,0),pv(pair,1,1,0),
                     pv(pair,1,1,1),pv(pair,1,0,1)},
                    {{leftRegion,1},{rightRegion,-1}});
        }
    }
    std::vector<std::unique_ptr<geo::GVertex>> vertices_;
    std::vector<std::unique_ptr<InterfaceEdge>> edges_;
    std::vector<std::unique_ptr<InterfaceFace>> faces_;
    std::vector<std::vector<bnd::VolumeRegionUse>> owners_;
    std::map<std::array<std::size_t,2>,int> edgeTags_;
    int regionCount_=2;
    bool reversePresentation_=false;
};

class LayerQualificationField final : public sizing::SizeField {
public:
    LayerQualificationField(double firstInterface,double sizeScale) noexcept
        : firstInterface_(firstInterface),sizeScale_(sizeScale) {}
    double eval(const Vec3& point) const override {
        return sizeScale_*(
            std::abs(point.x-firstInterface_)<=firstInterface_*0.125
                ? 0.5 : 0.75);
    }
private:
    double firstInterface_=1.0;
    double sizeScale_=1.0;
};

void configure(const geo::Model& model,Ownership& provider) {
    Vec3 lo{},hi{};model.solids().front()->bbox(lo,hi);
    const Vec3 center=(lo+hi)*0.5;
    for(const auto* face:model.faces()) {
        Vec2 uvLo{},uvHi{};face->paramBounds(uvLo,uvHi);
        const Vec2 uv=(uvLo+uvHi)*0.5;
        const int sense=face->normalAt(uv).dot(
            face->pointAt(uv)-center)>0?1:-1;
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            {{0,sense}};
    }
}

void requireSurfaceCacheReportXor(
        const geo::Model& model,const algo::SurfaceCacheReport& report) {
    CHECK(report.reportComplete);
    std::vector<algo::FacePersistentId> expected;
    expected.reserve(model.faces().size());
    for(const auto* face:model.faces()) {
        CHECK(face!=nullptr);
        CHECK(face&&face->tag()>=0);
        if(face&&face->tag()>=0)
            expected.push_back(
                static_cast<algo::FacePersistentId>(face->tag())+1U);
    }
    std::sort(expected.begin(),expected.end());
    CHECK(std::adjacent_find(expected.begin(),expected.end())==expected.end());

    std::vector<algo::FacePersistentId> observed;
    observed.reserve(report.processed.size()+report.unattempted.size());
    std::map<algo::FacePersistentId,std::array<std::size_t,2>> multiplicity;
    algo::FacePersistentId previous=0;
    bool first=true;
    for(const auto& row:report.processed) {
        CHECK(first||previous<row.facePersistentId);
        first=false;
        previous=row.facePersistentId;
        observed.push_back(row.facePersistentId);
        ++multiplicity[row.facePersistentId][0];
    }
    previous=0;
    first=true;
    for(const auto face:report.unattempted) {
        CHECK(first||previous<face);
        first=false;
        previous=face;
        observed.push_back(face);
        ++multiplicity[face][1];
    }
    if(!report.processed.empty()&&!report.unattempted.empty())
        CHECK(report.processed.back().facePersistentId<
              report.unattempted.front());
    CHECK(observed.size()==expected.size());
    std::sort(observed.begin(),observed.end());
    CHECK(observed==expected);
    for(const auto face:expected) {
        const auto found=multiplicity.find(face);
        CHECK(found!=multiplicity.end());
        if(found==multiplicity.end())continue;
        CHECK(found->second[0]+found->second[1]==1);
        CHECK(found->second[0]==0||found->second[1]==0);
    }
    if(report.hasFailedFace) {
        CHECK(report.failedFacePersistentId!=0);
        const auto found=multiplicity.find(report.failedFacePersistentId);
        CHECK(found!=multiplicity.end());
        if(found!=multiplicity.end()) {
            CHECK(found->second[0]==1);
            CHECK(found->second[1]==0);
        }
        std::size_t failedRows=0;
        for(const auto& row:report.processed)
            if(row.facePersistentId==report.failedFacePersistentId) {
                ++failedRows;
                CHECK(row.surfaceStatus!=algo::SurfaceFaceStatus::Ready);
            }
        CHECK(failedRows==1);
    } else {
        CHECK(report.failedFacePersistentId==0);
    }
}

bool sameClassification(mesh::Classification a,
                        mesh::Classification b) noexcept {
    return a.dim==b.dim&&a.tag==b.tag;
}

bool sameVec2Exact(const Vec2& a,const Vec2& b) noexcept {
    return std::bit_cast<std::uint64_t>(a.u==0.0?0.0:a.u)==
               std::bit_cast<std::uint64_t>(b.u==0.0?0.0:b.u)&&
           std::bit_cast<std::uint64_t>(a.v==0.0?0.0:a.v)==
               std::bit_cast<std::uint64_t>(b.v==0.0?0.0:b.v);
}

bool sameVec3Exact(const Vec3& a,const Vec3& b) noexcept {
    return std::bit_cast<std::uint64_t>(a.x==0.0?0.0:a.x)==
               std::bit_cast<std::uint64_t>(b.x==0.0?0.0:b.x)&&
           std::bit_cast<std::uint64_t>(a.y==0.0?0.0:a.y)==
               std::bit_cast<std::uint64_t>(b.y==0.0?0.0:b.y)&&
           std::bit_cast<std::uint64_t>(a.z==0.0?0.0:a.z)==
               std::bit_cast<std::uint64_t>(b.z==0.0?0.0:b.z);
}

bool sameParam(const bnd::AuthoritativeParam& a,
               const bnd::AuthoritativeParam& b) noexcept {
    if(a.index()!=b.index())return false;
    if(const auto* stored=std::get_if<bnd::StoredParam>(&a)) {
        const auto* other=std::get_if<bnd::StoredParam>(&b);
        return other&&sameVec2Exact(stored->uv,other->uv);
    }
    const auto& x=std::get<bnd::AffineParam>(a);
    const auto& y=std::get<bnd::AffineParam>(b);
    return x.lo==y.lo&&x.hi==y.hi&&x.num==y.num&&x.level==y.level;
}

bool sameChain(const algo::SurfaceChain& a,const algo::SurfaceChain& b) {
    if(a.chainPersistentId!=b.chainPersistentId||a.role!=b.role||
       !sameClassification(a.classification,b.classification)||
       a.closed!=b.closed||a.segmentEntityIds!=b.segmentEntityIds||
       a.owningCurve!=b.owningCurve||a.uses.size()!=b.uses.size())return false;
    for(std::size_t i=0;i<a.uses.size();++i)
        if(a.uses[i].faceLocalVertexId!=b.uses[i].faceLocalVertexId||
           a.uses[i].sharedKey!=b.uses[i].sharedKey)return false;
    return true;
}

bool sameCachedReport(
        const algo::t2c_detail::CachedSteinerReport& a,
        const algo::t2c_detail::CachedSteinerReport& b) {
    if(a.unattempted!=b.unattempted||a.processed.size()!=b.processed.size())
        return false;
    for(std::size_t i=0;i<a.processed.size();++i) {
        const auto& x=a.processed[i];const auto& y=b.processed[i];
        if(x.entity!=y.entity||x.result.status!=y.result.status||
           x.result.reject!=y.result.reject||x.result.refuse!=y.result.refuse||
           x.result.fail!=y.result.fail||x.result.revision!=y.result.revision||
           x.result.recoveredChain!=y.result.recoveredChain||
           x.result.createdParams.size()!=y.result.createdParams.size())
            return false;
        for(std::size_t j=0;j<x.result.createdParams.size();++j)
            if(x.result.createdParams[j].first!=
                    y.result.createdParams[j].first||
               !sameParam(x.result.createdParams[j].second,
                          y.result.createdParams[j].second))return false;
    }
    return true;
}

bool sameAudit(const algo::SurfaceTopologyAudit& a,
               const algo::SurfaceTopologyAudit& b) noexcept {
    return a.workspaceReady==b.workspaceReady&&
        a.connectedGenusZero==b.connectedGenusZero&&
        a.boundaryCyclesExact==b.boundaryCyclesExact&&
        a.nodeLinksValid==b.nodeLinksValid&&
        a.constraintsPersistent==b.constraintsPersistent&&
        a.classificationsPreserved==b.classificationsPreserved&&
        a.sharedMappingsPreserved==b.sharedMappingsPreserved&&
        a.noFrameSurvivors==b.noFrameSurvivors&&
        a.noDuplicateTriangles==b.noDuplicateTriangles&&
        a.physicalOrientation==b.physicalOrientation&&
        a.boundaryCycleCount==b.boundaryCycleCount&&
        a.eulerCharacteristic==b.eulerCharacteristic;
}

bool sameConstructionMetrics(const algo::Delaunay2DScalableMetrics& a,
                             const algo::Delaunay2DScalableMetrics& b)
        noexcept {
    return a.inputPoints==b.inputPoints&&
        a.coordinateSortComparisons==b.coordinateSortComparisons&&
        a.insertionOrderComparisons==b.insertionOrderComparisons&&
        a.bootstrapPoints==b.bootstrapPoints&&
        a.hullOrientationTests==b.hullOrientationTests&&
        a.bootstrapLegalizationEdgeVisits==
            b.bootstrapLegalizationEdgeVisits&&
        a.bootstrapLegalizationOrientationTests==
            b.bootstrapLegalizationOrientationTests&&
        a.bootstrapLegalizationIncircleTests==
            b.bootstrapLegalizationIncircleTests&&
        a.bootstrapCocircularTieTests==b.bootstrapCocircularTieTests&&
        a.bootstrapLegalizationFlips==b.bootstrapLegalizationFlips&&
        a.bootstrapWorklistInsertions==b.bootstrapWorklistInsertions&&
        a.bootstrapWorklistUpdates==b.bootstrapWorklistUpdates&&
        a.bootstrapStaleEdgeVisits==b.bootstrapStaleEdgeVisits&&
        a.bootstrapWorklistHighWater==b.bootstrapWorklistHighWater&&
        a.pointsInserted==b.pointsInserted&&
        a.locationTriangleTests==b.locationTriangleTests&&
        a.walkSteps==b.walkSteps&&
        a.cavityPredicateTests==b.cavityPredicateTests&&
        a.cavityMembers==b.cavityMembers&&
        a.boundaryEdges==b.boundaryEdges&&
        a.trianglesCommitted==b.trianglesCommitted&&
        a.fullScanFallbacks==b.fullScanFallbacks&&
        a.fullScanTriangleTests==b.fullScanTriangleTests&&
        a.peakLiveTriangles==b.peakLiveTriangles&&
        a.materializedTriangleSlots==b.materializedTriangleSlots&&
        a.peakCavitySize==b.peakCavitySize&&
        a.peakBoundarySize==b.peakBoundarySize&&
        a.scratchHighWaterBytes==b.scratchHighWaterBytes&&
        a.accountedCapacityHighWaterBytes==
            b.accountedCapacityHighWaterBytes&&
        a.adjacencyAuditPassed==b.adjacencyAuditPassed&&
        a.runComplete==b.runComplete;
}

bool sameFaceMetrics(const algo::ProductionSurfaceFaceMetrics& a,
                     const algo::ProductionSurfaceFaceMetrics& b) noexcept {
    return a.facePersistentId==b.facePersistentId&&
        a.boundaryNodes==b.boundaryNodes&&
        a.interiorCandidates==b.interiorCandidates&&
        a.spacingQueries==b.spacingQueries&&
        a.spacingBucketsVisited==b.spacingBucketsVisited&&
        a.spacingDistanceTests==b.spacingDistanceTests&&
        a.spacingIndexHighWater==b.spacingIndexHighWater&&
        a.constraintPointTests==b.constraintPointTests&&
        a.protectedSegments==b.protectedSegments&&
        a.recoveredProtectedSegments==b.recoveredProtectedSegments&&
        sameConstructionMetrics(a.construction,b.construction);
}

bool sameCachedArtifact(
        const algo::t2c_detail::CachedSurfaceFaceArtifact& a,
        const algo::t2c_detail::CachedSurfaceFaceArtifact& b) {
    if(a.facePersistentId!=b.facePersistentId||
       a.invariantAuthorityBytes!=b.invariantAuthorityBytes||
       a.authorityBytes!=b.authorityBytes||
       a.sizingSamples.size()!=b.sizingSamples.size()||
       a.faceStamp!=b.faceStamp||a.chart!=b.chart||
       !sameClassification(a.faceClass,b.faceClass)||
       a.desiredNormalSign!=b.desiredNormalSign||
       a.physicalPolicy.minPhysicalDoubleArea!=
            b.physicalPolicy.minPhysicalDoubleArea||
       a.physicalPolicy.maxVertexMappingError!=
            b.physicalPolicy.maxVertexMappingError||
       a.physicalPolicy.maxMappingDeviation!=
            b.physicalPolicy.maxMappingDeviation||
       a.physicalPolicy.minNormalAlignment!=
            b.physicalPolicy.minNormalAlignment||
       a.physicalPolicy.maxArtifactPositionError!=
            b.physicalPolicy.maxArtifactPositionError||
       !(a.expectedGeometryStamp==b.expectedGeometryStamp)||
       a.chartCertificate.providerStamp!=b.chartCertificate.providerStamp||
       a.chartCertificate.chart!=b.chartCertificate.chart||
       a.chartCertificate.certificateStamp!=
            b.chartCertificate.certificateStamp||
       a.chartCertificate.chartPositiveToNormalSign!=
            b.chartCertificate.chartPositiveToNormalSign||
       a.chartCertificate.regular!=b.chartCertificate.regular||
       a.chartCertificate.periodic!=b.chartCertificate.periodic||
       a.chartCertificate.seam!=b.chartCertificate.seam||
       a.chartCertificate.certifiedRegionBoundary.size()!=
            b.chartCertificate.certifiedRegionBoundary.size()||
       a.localNodes.size()!=b.localNodes.size()||
       a.localParams.size()!=b.localParams.size()||
       a.localMappingParams.size()!=b.localMappingParams.size()||
       a.refinedIdentityByLocalNode!=b.refinedIdentityByLocalNode||
       a.localCells.size()!=b.localCells.size()||
       a.sharedKeysByLocalNode!=b.sharedKeysByLocalNode||
       a.boundaryLoops!=b.boundaryLoops||
       a.normalizedChains.size()!=b.normalizedChains.size()||
       a.segmentOrigins.size()!=b.segmentOrigins.size()||
       a.persistentEdges.size()!=b.persistentEdges.size()||
       !sameCachedReport(a.protectedRecovery,b.protectedRecovery)||
       !sameCachedReport(a.localRecovery,b.localRecovery)||
       !sameAudit(a.audit,b.audit)||!sameFaceMetrics(a.metrics,b.metrics)||
       a.complete!=b.complete)return false;
    for(std::size_t i=0;i<a.sizingSamples.size();++i) {
        const auto& x=a.sizingSamples[i];const auto& y=b.sizingSamples[i];
        if(x.face!=y.face||x.ordinal!=y.ordinal||
           x.globalOrdinal!=y.globalOrdinal||x.status!=y.status||
           !sameVec3Exact(x.point,y.point)||x.curvature!=y.curvature||
           x.result!=y.result)return false;
    }
    for(std::size_t i=0;i<a.chartCertificate.certifiedRegionBoundary.size();++i)
        if(!sameVec2Exact(a.chartCertificate.certifiedRegionBoundary[i],
                          b.chartCertificate.certifiedRegionBoundary[i]))
            return false;
    for(std::size_t i=0;i<a.localNodes.size();++i)
        if(!sameVec3Exact(a.localNodes[i].coord,b.localNodes[i].coord)||
           !sameClassification(a.localNodes[i].cls,b.localNodes[i].cls))
            return false;
    for(std::size_t i=0;i<a.localParams.size();++i)
        if(!sameParam(a.localParams[i],b.localParams[i]))return false;
    for(std::size_t i=0;i<a.localMappingParams.size();++i)
        if(!sameVec2Exact(a.localMappingParams[i],b.localMappingParams[i]))
            return false;
    for(std::size_t i=0;i<a.localCells.size();++i)
        if(a.localCells[i].nodes!=b.localCells[i].nodes||
           !sameClassification(a.localCells[i].classification,
                               b.localCells[i].classification))return false;
    for(std::size_t i=0;i<a.normalizedChains.size();++i)
        if(!sameChain(a.normalizedChains[i],b.normalizedChains[i]))return false;
    auto ao=a.segmentOrigins.begin();auto bo=b.segmentOrigins.begin();
    for(;ao!=a.segmentOrigins.end();++ao,++bo) {
        const auto& x=ao->second;const auto& y=bo->second;
        if(ao->first!=bo->first||x.bndEntityId!=y.bndEntityId||
           x.chainPersistentId!=y.chainPersistentId||
           x.originalSegmentOrdinal!=y.originalSegmentOrdinal||
           x.normalizedSegmentOrdinal!=y.normalizedSegmentOrdinal||
           x.endpoints[0].faceLocalVertexId!=
                y.endpoints[0].faceLocalVertexId||
           x.endpoints[0].sharedKey!=y.endpoints[0].sharedKey||
           x.endpoints[1].faceLocalVertexId!=
                y.endpoints[1].faceLocalVertexId||
           x.endpoints[1].sharedKey!=y.endpoints[1].sharedKey||
           !sameClassification(x.classification,y.classification))return false;
    }
    for(std::size_t i=0;i<a.persistentEdges.size();++i) {
        const auto& x=a.persistentEdges[i];const auto& y=b.persistentEdges[i];
        if(x.faceStamp!=y.faceStamp||x.bndEntityId!=y.bndEntityId||
           x.chainPersistentId!=y.chainPersistentId||
           x.originalSegmentOrdinal!=y.originalSegmentOrdinal||
           x.normalizedSegmentOrdinal!=y.normalizedSegmentOrdinal||
           x.role!=y.role||x.localEndpoints!=y.localEndpoints||
           x.endpointSharedKeys!=y.endpointSharedKeys||
           !sameClassification(x.classification,y.classification)||
           x.required!=y.required)return false;
    }
    // Capacity/accounting is intentionally excluded: immutable logical rows
    // must match while allocator reserve choices may differ.
    return true;
}

bool sameSurfaceFailure(const algo::SurfaceFailure& a,
                        const algo::SurfaceFailure& b) noexcept {
    return a.reason==b.reason&&a.candidate==b.candidate&&
        a.domain==b.domain&&a.allocation==b.allocation&&
        a.bndReject==b.bndReject&&a.bndRefuse==b.bndRefuse&&
        a.bndFail==b.bndFail;
}

bool sameSharedLedger(
        const std::map<algo::SharedCurveNodeKey,
                       algo::SharedCurveNodeBinding>& a,
        const std::map<algo::SharedCurveNodeKey,
                       algo::SharedCurveNodeBinding>& b) noexcept {
    if(a.size()!=b.size())return false;
    auto ai=a.begin();auto bi=b.begin();
    for(;ai!=a.end();++ai,++bi)
        if(ai->first!=bi->first||
           ai->second.existingGlobalIndex!=bi->second.existingGlobalIndex||
           !sameClassification(ai->second.classification,
                               bi->second.classification))return false;
    return true;
}

template <class Handle>
bool sameCanonicalHandle(const Handle& a,const Handle& b) noexcept {
    // Pool ids are attempt-local ownership tokens and are deliberately absent
    // from every deterministic/public identity domain.
    return a.slot==b.slot&&a.gen==b.gen;
}

bool sameSteinerReport(const bnd::SteinerSegmentsReport& a,
                       const bnd::SteinerSegmentsReport& b) {
    if(a.unattempted!=b.unattempted||a.processed.size()!=b.processed.size())
        return false;
    for(std::size_t i=0;i<a.processed.size();++i) {
        const auto& x=a.processed[i];const auto& y=b.processed[i];
        if(x.entity!=y.entity||x.result.status!=y.result.status||
           x.result.reject!=y.result.reject||
           x.result.refuse!=y.result.refuse||x.result.fail!=y.result.fail||
           x.result.revision!=y.result.revision||
           x.result.recoveredChain.size()!=y.result.recoveredChain.size()||
           x.result.createdParams.size()!=y.result.createdParams.size())
            return false;
        for(std::size_t j=0;j<x.result.recoveredChain.size();++j)
            if(!sameCanonicalHandle(x.result.recoveredChain[j],
                                    y.result.recoveredChain[j]))return false;
        for(std::size_t j=0;j<x.result.createdParams.size();++j)
            if(!sameCanonicalHandle(x.result.createdParams[j].first,
                                    y.result.createdParams[j].first)||
               !sameParam(x.result.createdParams[j].second,
                          y.result.createdParams[j].second))return false;
    }
    return true;
}

bool sameRetainedSurfaceFace(const algo::SurfaceFaceResult& a,
                             const algo::SurfaceFaceResult& b) {
    const auto samePolicy=[](const algo::SurfacePhysicalPolicy& x,
                             const algo::SurfacePhysicalPolicy& y) noexcept {
        return x.minPhysicalDoubleArea==y.minPhysicalDoubleArea&&
            x.maxVertexMappingError==y.maxVertexMappingError&&
            x.maxMappingDeviation==y.maxMappingDeviation&&
            x.minNormalAlignment==y.minNormalAlignment&&
            x.maxArtifactPositionError==y.maxArtifactPositionError;
    };
    if(a.faceStamp!=b.faceStamp||a.chart!=b.chart||
       !sameClassification(a.faceClass,b.faceClass)||
       a.desiredNormalSign!=b.desiredNormalSign||
       !samePolicy(a.physicalPolicy,b.physicalPolicy)||
       a.expectedGeometryStamp!=b.expectedGeometryStamp||
       a.chartCertificate.providerStamp!=b.chartCertificate.providerStamp||
       a.chartCertificate.chart!=b.chartCertificate.chart||
       a.chartCertificate.certificateStamp!=
            b.chartCertificate.certificateStamp||
       a.chartCertificate.chartPositiveToNormalSign!=
            b.chartCertificate.chartPositiveToNormalSign||
       a.chartCertificate.regular!=b.chartCertificate.regular||
       a.chartCertificate.periodic!=b.chartCertificate.periodic||
       a.chartCertificate.seam!=b.chartCertificate.seam||
       a.chartCertificate.certifiedRegionBoundary.size()!=
            b.chartCertificate.certifiedRegionBoundary.size()||
       a.localNodes.size()!=b.localNodes.size()||
       a.localParams.size()!=b.localParams.size()||
       a.localMappingParams.size()!=b.localMappingParams.size()||
       a.refinedIdentityByLocalNode!=b.refinedIdentityByLocalNode||
       a.localHandles.size()!=b.localHandles.size()||
       a.localCells.size()!=b.localCells.size()||
       a.existingGlobalIndex!=b.existingGlobalIndex||
       a.publishedGlobalByLocalNode!=b.publishedGlobalByLocalNode||
       a.publishedGlobalCellByLocalCell!=b.publishedGlobalCellByLocalCell||
       a.boundaryLoops!=b.boundaryLoops||
       a.normalizedChains.size()!=b.normalizedChains.size()||
       a.segmentOrigins.size()!=b.segmentOrigins.size()||
       a.persistentEdges.size()!=b.persistentEdges.size()||
       !sameSteinerReport(a.protectedRecovery,b.protectedRecovery)||
       !sameSteinerReport(a.localRecovery,b.localRecovery)||
       !sameAudit(a.audit,b.audit)||
       static_cast<bool>(a.localWorkspace)!=
            static_cast<bool>(b.localWorkspace))return false;
    for(std::size_t i=0;i<a.chartCertificate.certifiedRegionBoundary.size();
            ++i)
        if(!sameVec2Exact(a.chartCertificate.certifiedRegionBoundary[i],
                          b.chartCertificate.certifiedRegionBoundary[i]))
            return false;
    for(std::size_t i=0;i<a.localNodes.size();++i)
        if(!sameVec3Exact(a.localNodes[i].coord,b.localNodes[i].coord)||
           !sameClassification(a.localNodes[i].cls,b.localNodes[i].cls))
            return false;
    for(std::size_t i=0;i<a.localParams.size();++i)
        if(!sameParam(a.localParams[i],b.localParams[i]))return false;
    for(std::size_t i=0;i<a.localMappingParams.size();++i)
        if(!sameVec2Exact(a.localMappingParams[i],b.localMappingParams[i]))
            return false;
    for(std::size_t i=0;i<a.localHandles.size();++i)
        if(!sameCanonicalHandle(a.localHandles[i],b.localHandles[i]))
            return false;
    for(std::size_t i=0;i<a.localCells.size();++i)
        if(a.localCells[i].nodes!=b.localCells[i].nodes||
           !sameClassification(a.localCells[i].classification,
                               b.localCells[i].classification))return false;
    for(std::size_t i=0;i<a.normalizedChains.size();++i)
        if(!sameChain(a.normalizedChains[i],b.normalizedChains[i]))return false;
    auto ao=a.segmentOrigins.begin();auto bo=b.segmentOrigins.begin();
    for(;ao!=a.segmentOrigins.end();++ao,++bo) {
        const auto& x=ao->second;const auto& y=bo->second;
        if(ao->first!=bo->first||x.bndEntityId!=y.bndEntityId||
           x.chainPersistentId!=y.chainPersistentId||
           x.originalSegmentOrdinal!=y.originalSegmentOrdinal||
           x.normalizedSegmentOrdinal!=y.normalizedSegmentOrdinal||
           x.endpoints[0].faceLocalVertexId!=
                y.endpoints[0].faceLocalVertexId||
           x.endpoints[0].sharedKey!=y.endpoints[0].sharedKey||
           x.endpoints[1].faceLocalVertexId!=
                y.endpoints[1].faceLocalVertexId||
           x.endpoints[1].sharedKey!=y.endpoints[1].sharedKey||
           !sameClassification(x.classification,y.classification))return false;
    }
    for(std::size_t i=0;i<a.persistentEdges.size();++i) {
        const auto& x=a.persistentEdges[i];const auto& y=b.persistentEdges[i];
        if(x.faceStamp!=y.faceStamp||x.bndEntityId!=y.bndEntityId||
           x.chainPersistentId!=y.chainPersistentId||
           x.originalSegmentOrdinal!=y.originalSegmentOrdinal||
           x.normalizedSegmentOrdinal!=y.normalizedSegmentOrdinal||
           x.role!=y.role||
           !sameCanonicalHandle(x.localEndpoints[0],y.localEndpoints[0])||
           !sameCanonicalHandle(x.localEndpoints[1],y.localEndpoints[1])||
           x.endpointSharedKeys!=y.endpointSharedKeys||
           !sameClassification(x.classification,y.classification)||
           x.required!=y.required||
           x.publishedGlobalEndpoints!=y.publishedGlobalEndpoints)return false;
    }
    if(a.localWorkspace&&
       (a.localWorkspace->semanticTopologyHash()!=
            b.localWorkspace->semanticTopologyHash()||
        a.localWorkspace->strictWorkspaceHash()!=
            b.localWorkspace->strictWorkspaceHash()))return false;
    return true;
}

bool sameRetainedSurfaceSet(const algo::SurfaceSetOutcome& a,
                            const algo::SurfaceSetOutcome& b) {
    if(a.setStatus!=b.setStatus||!sameSurfaceFailure(a.setFailure,b.setFailure)||
       a.reportComplete!=b.reportComplete||
       a.unattempted!=b.unattempted||
       a.processed.size()!=b.processed.size()||
       static_cast<bool>(a.retainedArtifact)!=
            static_cast<bool>(b.retainedArtifact))return false;
    for(std::size_t i=0;i<a.processed.size();++i) {
        const auto& x=a.processed[i];const auto& y=b.processed[i];
        if(x.facePersistentId!=y.facePersistentId||x.status!=y.status||
           !sameSurfaceFailure(x.failure,y.failure)||
           static_cast<bool>(x.readyResult)!=
                static_cast<bool>(y.readyResult))return false;
        if(x.readyResult&&
           !sameRetainedSurfaceFace(*x.readyResult,*y.readyResult))return false;
    }
    if(!a.retainedArtifact)return true;
    const auto& x=*a.retainedArtifact;const auto& y=*b.retainedArtifact;
    if(x.setStamp!=y.setStamp||
       x.publishedGlobalByRefinedIdentity!=
            y.publishedGlobalByRefinedIdentity||
       !sameSharedLedger(x.sharedLedger,y.sharedLedger)||
       x.faces.size()!=y.faces.size())return false;
    for(std::size_t i=0;i<x.faces.size();++i) {
        if(static_cast<bool>(x.faces[i])!=static_cast<bool>(y.faces[i]))
            return false;
        if(x.faces[i]&&
           !sameRetainedSurfaceFace(*x.faces[i],*y.faces[i]))return false;
    }
    return true;
}

bool sameSurfaceCacheReport(const algo::SurfaceCacheReport& a,
                            const algo::SurfaceCacheReport& b) noexcept {
    if(a.reportComplete!=b.reportComplete||
       a.hasFailedFace!=b.hasFailedFace||
       a.failedFacePersistentId!=b.failedFacePersistentId||
       a.unattempted!=b.unattempted||a.processed.size()!=b.processed.size())
        return false;
    for(std::size_t i=0;i<a.processed.size();++i)
        if(a.processed[i].facePersistentId!=b.processed[i].facePersistentId||
           a.processed[i].disposition!=b.processed[i].disposition||
           a.processed[i].surfaceStatus!=b.processed[i].surfaceStatus||
           !sameSurfaceFailure(a.processed[i].originalSurfaceFailure,
                               b.processed[i].originalSurfaceFailure))
            return false;
    return true;
}

bool sameSubdivisionLedger(
        const std::vector<algo::ProtectedFeatureChildIdentity>& a,
        const std::vector<algo::ProtectedFeatureChildIdentity>& b) noexcept {
    if(a.size()!=b.size())return false;
    for(std::size_t i=0;i<a.size();++i) {
        const auto& x=a[i];const auto& y=b[i];
        if(x.parentEntityPersistentId!=y.parentEntityPersistentId||
           x.childIdentity!=y.childIdentity||x.splitOrdinal!=y.splitOrdinal||
           x.edgeTag!=y.edgeTag||x.faceTag!=y.faceTag||
           std::bit_cast<std::uint64_t>(x.parentParameterLo)!=
                std::bit_cast<std::uint64_t>(y.parentParameterLo)||
           std::bit_cast<std::uint64_t>(x.parentParameterHi)!=
                std::bit_cast<std::uint64_t>(y.parentParameterHi)||
           std::bit_cast<std::uint64_t>(x.midpointParameter)!=
                std::bit_cast<std::uint64_t>(y.midpointParameter)||
           !sameVec2Exact(x.faceParameter,y.faceParameter)||
           !sameVec3Exact(x.coordinate,y.coordinate)||
           !sameClassification(x.classification,y.classification)||
           x.providerStamp!=y.providerStamp)return false;
    }
    return true;
}

bool sameProviderRequests(
        const std::vector<algo::ProtectedFeaturePlacementRequest>& a,
        const std::vector<algo::ProtectedFeaturePlacementRequest>& b)
        noexcept {
    if(a.size()!=b.size())return false;
    for(std::size_t i=0;i<a.size();++i)
        if(a[i].kind!=b[i].kind||a[i].entityStamp!=b[i].entityStamp||
           a[i].requestedIdentity!=b[i].requestedIdentity||
           a[i].entityTag!=b[i].entityTag||
           std::bit_cast<std::uint64_t>(a[i].curveParameter)!=
                std::bit_cast<std::uint64_t>(b[i].curveParameter)||
           !sameVec2Exact(a[i].faceParameter,b[i].faceParameter))return false;
    return true;
}

bool sameV4Metrics(const algo::BoundaryConformingV4Metrics& a,
                   const algo::BoundaryConformingV4Metrics& b) noexcept {
    // Wall-clock phase samples are intentionally non-canonical. Every count,
    // work value, capacity/accounting value and completion bit is pinned.
    return std::tie(
        a.rounds,a.edgeSplits,a.facetSplits,a.providerStampCalls,
        a.providerPlacementCalls,a.rebuiltSurfaceNodes,
        a.rebuiltSurfaceSegments,a.candidateTetCells,a.semanticStates,
        a.frontierRows,a.uniquePlacements,a.deferredFrontierRows,
        a.nestedConstructionWorkUnits,a.nestedRecoveryWorkUnits,
        a.totalMeasuredWorkUnits,a.coordinatorWorkUnits,
        a.surfaceCacheFacesBuilt,a.surfaceCacheFacesReused,
        a.surfaceCacheFacesReplayed,a.surfaceCacheFacesRejected,
        a.surfaceCacheKeyRowsVisited,a.surfaceCacheDirtyPropagationVisits,
        a.surfaceCacheSizingReplayCalls,a.surfaceCachePreparationWorkUnits,
        a.surfaceCacheReplayWorkUnits,a.surfaceCacheSetAuditWorkUnits,
        a.surfaceCacheCommits,a.surfaceCacheCleanNodeRows,
        a.surfaceCacheDirtyNodeRows,a.surfaceCacheCleanCellRows,
        a.surfaceCacheDirtyCellRows,a.surfaceCacheCleanConstraintRows,
        a.surfaceCacheDirtyConstraintRows,a.surfaceCacheCleanMappingRows,
        a.surfaceCacheDirtyMappingRows,a.surfaceCacheCopiedBytes,
        a.surfaceCacheHighWaterSlots,a.surfaceCacheHighWaterRows,
        a.surfaceCacheAllFaceScansPerDirtyFace,a.facetConflictRounds,
        a.facetConflictSccMembers,a.facetConflictDependencyRows,
        a.facetConflictProviderCalls,a.facetConflictCoordinatorWorkUnits,
        a.surfaceCacheAccountedBytes,a.peakAccountedBytes,a.runComplete)==
           std::tie(
        b.rounds,b.edgeSplits,b.facetSplits,b.providerStampCalls,
        b.providerPlacementCalls,b.rebuiltSurfaceNodes,
        b.rebuiltSurfaceSegments,b.candidateTetCells,b.semanticStates,
        b.frontierRows,b.uniquePlacements,b.deferredFrontierRows,
        b.nestedConstructionWorkUnits,b.nestedRecoveryWorkUnits,
        b.totalMeasuredWorkUnits,b.coordinatorWorkUnits,
        b.surfaceCacheFacesBuilt,b.surfaceCacheFacesReused,
        b.surfaceCacheFacesReplayed,b.surfaceCacheFacesRejected,
        b.surfaceCacheKeyRowsVisited,b.surfaceCacheDirtyPropagationVisits,
        b.surfaceCacheSizingReplayCalls,b.surfaceCachePreparationWorkUnits,
        b.surfaceCacheReplayWorkUnits,b.surfaceCacheSetAuditWorkUnits,
        b.surfaceCacheCommits,b.surfaceCacheCleanNodeRows,
        b.surfaceCacheDirtyNodeRows,b.surfaceCacheCleanCellRows,
        b.surfaceCacheDirtyCellRows,b.surfaceCacheCleanConstraintRows,
        b.surfaceCacheDirtyConstraintRows,b.surfaceCacheCleanMappingRows,
        b.surfaceCacheDirtyMappingRows,b.surfaceCacheCopiedBytes,
        b.surfaceCacheHighWaterSlots,b.surfaceCacheHighWaterRows,
        b.surfaceCacheAllFaceScansPerDirtyFace,b.facetConflictRounds,
        b.facetConflictSccMembers,b.facetConflictDependencyRows,
        b.facetConflictProviderCalls,b.facetConflictCoordinatorWorkUnits,
        b.surfaceCacheAccountedBytes,b.peakAccountedBytes,b.runComplete)&&
        a.mixedPresentScheduling==b.mixedPresentScheduling&&
        a.latestSourceMixedPresentScheduling==
            b.latestSourceMixedPresentScheduling;
}

bool sameRecoveryMetrics(const algo::TetRecoveryMetrics& a,
                         const algo::TetRecoveryMetrics& b) noexcept {
    return std::tie(
        a.contextNodes,a.contextCells,a.contextEdges,a.contextFaces,
        a.contextAuditVisits,a.indexBuildComparisons,a.indexBytes,
        a.edgeAabbQueries,a.facetAabbQueries,a.indexNodesVisited,
        a.returnedTetCandidates,a.segmentTetTests,a.targetTetTests,
        a.candidateSurvivorTests,a.segmentWalkCells,a.linkCycleHighWater,
        a.edgeCavityCommits,a.predictedEdgeMatches,a.completeStarQueries,
        a.closureInsertions,a.apexPairs,a.facetCommits,a.indexRemovals,
        a.indexInsertions,a.indexRebalances,a.indexHighWater,
        a.routeAssessments,a.dependencyEdges,a.dependencyInvalidations,
        a.selectedRebuilds,a.batchPreflightAssessments,
        a.alreadyPresentBatchEntities,a.protectedImportFallbacks,
        a.protectedConflicts,a.conflictSccVisits,a.fullStageRescans,
        a.productionGlobalScanFallbacks,a.scratchHighWaterBytes,
        a.runComplete)==std::tie(
        b.contextNodes,b.contextCells,b.contextEdges,b.contextFaces,
        b.contextAuditVisits,b.indexBuildComparisons,b.indexBytes,
        b.edgeAabbQueries,b.facetAabbQueries,b.indexNodesVisited,
        b.returnedTetCandidates,b.segmentTetTests,b.targetTetTests,
        b.candidateSurvivorTests,b.segmentWalkCells,b.linkCycleHighWater,
        b.edgeCavityCommits,b.predictedEdgeMatches,b.completeStarQueries,
        b.closureInsertions,b.apexPairs,b.facetCommits,b.indexRemovals,
        b.indexInsertions,b.indexRebalances,b.indexHighWater,
        b.routeAssessments,b.dependencyEdges,b.dependencyInvalidations,
        b.selectedRebuilds,b.batchPreflightAssessments,
        b.alreadyPresentBatchEntities,b.protectedImportFallbacks,
        b.protectedConflicts,b.conflictSccVisits,b.fullStageRescans,
        b.productionGlobalScanFallbacks,b.scratchHighWaterBytes,
        b.runComplete);
}

bool sameTetEdgeOutcome(const bnd::TetEdgeRecoveryOutcome& a,
                        const bnd::TetEdgeRecoveryOutcome& b) noexcept {
    return a.status==b.status&&a.reason==b.reason&&a.phase==b.phase&&
        sameCanonicalHandle(a.edge,b.edge)&&
        a.geometryReason==b.geometryReason&&a.counters==b.counters&&
        a.capacityResiduePermitted==b.capacityResiduePermitted;
}

bool sameTetFacetOutcome(const bnd::TetFacetRecoveryOutcome& a,
                         const bnd::TetFacetRecoveryOutcome& b) noexcept {
    return a.status==b.status&&a.reason==b.reason&&a.phase==b.phase&&
        sameCanonicalHandle(a.face,b.face)&&
        a.geometryReason==b.geometryReason&&a.counters==b.counters&&
        a.capacityResiduePermitted==b.capacityResiduePermitted;
}

bool sameV3WorkKey(const algo::BoundaryConformingV3WorkKey& a,
                   const algo::BoundaryConformingV3WorkKey& b) noexcept {
    return a.kind==b.kind&&a.edgeKey==b.edgeKey&&a.facetKey==b.facetKey;
}

bool sameV3Failure(const algo::BoundaryConformingV3Failure& a,
                   const algo::BoundaryConformingV3Failure& b) noexcept {
    return a.reason==b.reason&&a.phase==b.phase&&
        a.providerCheckpoint==b.providerCheckpoint&&
        a.workKind==b.workKind&&a.edgeKey==b.edgeKey&&
        a.facetKey==b.facetKey&&a.scheduleOrdinal==b.scheduleOrdinal&&
        sameTetEdgeOutcome(a.edgeOutcome,b.edgeOutcome)&&
        sameTetFacetOutcome(a.facetOutcome,b.facetOutcome)&&
        a.exactWitnessDigest==b.exactWitnessDigest;
}

bool sameCompleteV3Report(const algo::BoundaryConformingV3Outcome& a,
                          const algo::BoundaryConformingV3Outcome& b) {
    if(a.reportComplete!=b.reportComplete||
       static_cast<bool>(a.report)!=static_cast<bool>(b.report))return false;
    if(!a.report)return true;
    if(a.report->reportComplete!=b.report->reportComplete||
       a.report->processed.size()!=b.report->processed.size()||
       a.report->unattempted.size()!=b.report->unattempted.size())return false;
    for(std::size_t i=0;i<a.report->processed.size();++i) {
        const auto& x=a.report->processed[i];const auto& y=b.report->processed[i];
        if(!sameV3WorkKey(x.key,y.key)||x.status!=y.status||
           x.scheduleOrdinal!=y.scheduleOrdinal||
           !sameV3Failure(x.failure,y.failure))return false;
    }
    for(std::size_t i=0;i<a.report->unattempted.size();++i)
        if(!sameV3WorkKey(a.report->unattempted[i],
                          b.report->unattempted[i]))return false;
    return true;
}

void capacitySnapshotContract() {
    TEST("t5 mesh capacity snapshot is exact stable and read-only");
    mesh::Mesh mesh;
    const auto empty=mesh.capacitySnapshot();
    CHECK(empty.nodeCoordinateCapacity==0);
    CHECK(empty.nodeClassificationCapacity==0);
    CHECK(empty.cellCapacity==0);
    CHECK(empty.connectivityCapacity==0);
    CHECK(empty.dedupEntryCount==0);
    CHECK(empty.dedupIndexCapacity==0);
    CHECK(empty.complete);
    mesh.reserveAppend(8,3,10);
    const auto reserved=mesh.capacitySnapshot();
    CHECK(reserved.nodeCoordinateCapacity>=8);
    CHECK(reserved.nodeClassificationCapacity>=8);
    CHECK(reserved.cellCapacity>=3);
    CHECK(reserved.connectivityCapacity>=10);
    const auto a=mesh.addNode({0,0,0},{0,1});
    const auto b=mesh.addNode({1,0,0},{0,2});
    const auto c=mesh.addNode({0,1,0},{0,3});
    mesh.addCell(mesh::CellType::Tri3,{a,b,c},{2,7});
    const auto populated=mesh.capacitySnapshot();
    CHECK(populated.nodeCoordinateCapacity==reserved.nodeCoordinateCapacity);
    CHECK(populated.nodeClassificationCapacity==
          reserved.nodeClassificationCapacity);
    CHECK(populated.cellCapacity==reserved.cellCapacity);
    CHECK(populated.connectivityCapacity==reserved.connectivityCapacity);
    CHECK(populated.dedupEntryCount==0);
    CHECK(populated.dedupIndexCapacity==0);
    const auto hash=mesh::meshHash(mesh);
    std::atomic<bool> good{true};
    std::array<std::thread,4> readers;
    for(auto& reader:readers)reader=std::thread([&]{
        for(int i=0;i<1000;++i) {
            const auto now=mesh.capacitySnapshot();
            if(now.nodeCoordinateCapacity!=populated.nodeCoordinateCapacity||
               now.connectivityCapacity!=populated.connectivityCapacity||
               !now.complete)good=false;
        }
    });
    for(auto& reader:readers)reader.join();
    CHECK(good.load());
    CHECK(mesh::meshHash(mesh)==hash);
}

void capacitySnapshotDedupAndOverflow() {
    TEST("t5 mesh capacity snapshot observes nested dedup and saturates");
    mesh::Mesh mesh;mesh.setDedupTolerance(1e-6);
    (void)mesh.addNodeUnique({0,0,0},{0,1});
    (void)mesh.addNodeUnique({1,0,0},{0,2});
    const auto snapshot=mesh.capacitySnapshot();
    CHECK(snapshot.complete);
    CHECK(snapshot.dedupEntryCount==2);
    CHECK(snapshot.dedupIndexCapacity>=2);
    std::size_t total=(std::numeric_limits<std::size_t>::max)()-3;
    CHECK(!mesh::meshCapacityCheckedAdd(total,4));
    CHECK(total==(std::numeric_limits<std::size_t>::max)());
}

void sha256KnownAnswer() {
    TEST("t5 local sha256 matches the standard abc vector");
    const std::array<std::uint8_t,3> input{{'a','b','c'}};
    const auto digest=algo::t5_detail::sha256(input);
    constexpr std::array<std::uint8_t,32> expected{{
        0xba,0x78,0x16,0xbf,0x8f,0x01,0xcf,0xea,
        0x41,0x41,0x40,0xde,0x5d,0xae,0x22,0x23,
        0xb0,0x03,0x61,0xa3,0x96,0x17,0x7a,0x9c,
        0xb4,0x10,0xff,0x61,0xf2,0x00,0x15,0xad}};
    CHECK(digest==expected);
}

std::int64_t integerSixVolume(const std::array<Vec3,4>& p) {
    const auto value=[](double x) {
        return static_cast<std::int64_t>(x);
    };
    const std::int64_t ax=value(p[1].x)-value(p[0].x);
    const std::int64_t ay=value(p[1].y)-value(p[0].y);
    const std::int64_t az=value(p[1].z)-value(p[0].z);
    const std::int64_t bx=value(p[2].x)-value(p[0].x);
    const std::int64_t by=value(p[2].y)-value(p[0].y);
    const std::int64_t bz=value(p[2].z)-value(p[0].z);
    const std::int64_t cx=value(p[3].x)-value(p[0].x);
    const std::int64_t cy=value(p[3].y)-value(p[0].y);
    const std::int64_t cz=value(p[3].z)-value(p[0].z);
    return ax*(by*cz-bz*cy)+ay*(bz*cx-bx*cz)+az*(bx*cy-by*cx);
}

void exactVolumeContract() {
    TEST("t5 exact volume pins public bound translation and hard cases");
    static_assert(algo::kT5SixVolumeExpansionTermsPerTet==192);
    const std::array<Vec3,4> unit{{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}};
    const auto unitProbe=algo::t5_detail::exactSixVolumeForTesting(unit);
    CHECK(unitProbe.complete);
    CHECK(unitProbe.termCount==unitProbe.workUnits);
    CHECK(unitProbe.termCount<=algo::kT5SixVolumeExpansionTermsPerTet);
    CHECK(unitProbe.volumeLower<=1.0/6.0);
    CHECK(unitProbe.volumeUpper>=1.0/6.0);

    constexpr double translation=1125899906842624.0; // 2^50, exact.
    const std::array<Vec3,4> translated{{
        {translation,translation,translation},
        {translation+3,translation+1,translation},
        {translation+1,translation+4,translation+1},
        {translation,translation+1,translation+5}}};
    const auto translatedProbe=
        algo::t5_detail::exactSixVolumeForTesting(translated);
    CHECK(translatedProbe.complete);
    const auto exact=integerSixVolume(translated);
    CHECK(exact>0);
    const double expected=static_cast<double>(exact)/6.0;
    CHECK(translatedProbe.volumeLower<=expected);
    CHECK(translatedProbe.volumeUpper>=expected);
    CHECK(translatedProbe.termCount<=algo::kT5SixVolumeExpansionTermsPerTet);

    const std::array<Vec3,4> hard{{
        {-0x1.f1b4bc6ee1736p-2,-0x1.50a02cc05e259p-1,
         -0x1.d62154f8e8101p-1},
        {-0x1.d71c68a96f7a7p-1, 0x1.c5a404d317909p-1,
          0x1.f38eb124853a8p-4},
        {-0x1.4f54dbd17be9cp-3,-0x1.02fca79e0b139p-1,
          0x1.d94eb1c8fd696p-2},
        {-0x1.5abfd59573ab1p-1,-0x1.2156129793eb1p-1,
          0x1.e2eaeb4a77124p-3}}};
    const auto hardProbe=algo::t5_detail::exactSixVolumeForTesting(hard);
    CHECK(hardProbe.complete);
    CHECK(hardProbe.termCount==6);
    CHECK(hardProbe.termCount==hardProbe.workUnits);
    CHECK(hardProbe.termCount<=algo::kT5SixVolumeExpansionTermsPerTet);
    CHECK(hardProbe.exactDigest==4163837141100653566ULL);
    CHECK(hardProbe.volumeLower==-0x1.41308de56a95ap-3);
    CHECK(hardProbe.volumeUpper==-0x1.41308de56a957p-3);
}

algo::T5Request singleRequest(const Ownership& provider,
                              const algo::T5RegionBinding& region) {
    algo::T5Request request;
    request.sourceRequest.expectedSurfaceSetStamp={15002,1};
    request.sourceRequest.expectedOwnershipProviderStamp=provider.current;
    request.regionBindings=std::span<const algo::T5RegionBinding>(&region,1);
    return request;
}

void requireCompleteXor(const algo::T5Outcome& outcome) {
    CHECK(outcome.report!=nullptr);
    if(!outcome.report)return;
    const auto& report=*outcome.report;
    CHECK(report.reportComplete);
    CHECK(report.processedUsed==report.processedStorage.size());
    CHECK(report.unattemptedUsed==report.unattemptedStorage.size());
    std::set<std::uint64_t> ordinals;
    for(const auto& row:report.processedStorage)
        CHECK(ordinals.insert(row.canonicalOrdinal).second);
    CHECK(ordinals.size()==report.processedStorage.size());
}

void singleRegionRoute() {
    TEST("t5 single-region authority publishes one audited mesh");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    Ownership provider;configure(*model,provider);
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(0.55);
    algo::T5RegionBinding region;
    region.regionTag=0;region.regionPersistentId=101;
    region.materialPersistentId=201;
    region.volumeClassification={3,0};
    auto request=singleRequest(provider,region);
    mesh::Mesh destination;
    auto outcome=algo::meshModelMultiRegionTet4(
        *model,spec,request,destination,provider);
    CHECK(outcome.published());
    CHECK(outcome.sourceReceipt!=nullptr);
    CHECK(outcome.retainedResult!=nullptr);
    CHECK(outcome.reportComplete);
    if(!outcome.retainedResult)return;
    const auto& result=*outcome.retainedResult;
    CHECK(result.audit.passed);
    requireCompleteXor(outcome);
    CHECK(result.regions.size()==1);
    CHECK(result.materials.size()==1);
    CHECK(result.interfaceFacets.empty());
    CHECK(result.privateMesh==nullptr);
    CHECK(result.metrics.perInterfaceGlobalTetScans==0);
    CHECK(result.metrics.perRegionGlobalFaceScans==0);
    CHECK(result.metrics.contextRebuilds==0);
    CHECK(result.metrics.publicationTailAllocations==0);
    CHECK(result.metrics.totalWorkUnits>0);
    CHECK(outcome.sourceReceipt->constructionMetrics.runComplete);
    const auto& construction=
        outcome.sourceReceipt->constructionMetrics;
    CHECK(construction.totalWorkUnits>0);
    CHECK(construction.importAttempts>=1);
    CHECK(construction.importEdgeValidationCellVisits==
          construction.importCellVisits);
    CHECK(construction.importEdgeBuildCellVisits==
          construction.importCellVisits);
    CHECK(construction.importFaceValidationCellVisits==
          construction.importCellVisits);
    CHECK(construction.importFaceBuildCellVisits==
          construction.importCellVisits);
    CHECK(construction.importEdgeUses==6*construction.importCellVisits);
    CHECK(construction.importFaceUses==4*construction.importCellVisits);
    CHECK(construction.importEdgeAudits==2*construction.importAttempts);
    CHECK(construction.importFaceAudits==construction.importAttempts);
    CHECK(construction.importPerCellGlobalAudits==0);
    CHECK(construction.partialImportRowsCopied==0);
    CHECK(construction.partialImportGlobalSlotsInitialized==0);
    CHECK(construction.partialImportCandidateMappingVisits==0);
    CHECK(construction.partialImportEdgeValidationRows==0);
    CHECK(construction.partialImportFaceValidationRows==0);
    CHECK(construction.partialImportProtectedEndpointVisits==0);
    CHECK(construction.partialImportAuthorityPointValidations==0);
    CHECK(construction.partialImportDuplicateKeyComparisons==0);
    CHECK(construction.partialImportResolutionNodeVisits==0);
    CHECK(construction.partialImportEdgeResolutionLookups==0);
    CHECK(construction.partialImportFaceResolutionLookups==0);
    CHECK(construction.partialImportEdgeDataReads==0);
    CHECK(construction.partialImportFaceDataCopies==0);
    CHECK(construction.partialImportFaceIncidenceRowsCopied==0);
    CHECK(construction.partialImportSortComparisons==0);
    const auto& extraction=outcome.sourceReceipt->extractionMetrics;
    CHECK(extraction.edgeBuild.validationCellVisits==
          extraction.rebuildCellInsertions);
    CHECK(extraction.edgeBuild.buildCellVisits==
          extraction.rebuildCellInsertions);
    CHECK(extraction.faceBuild.validationCellVisits==
          extraction.rebuildCellInsertions);
    CHECK(extraction.faceBuild.buildCellVisits==
          extraction.rebuildCellInsertions);
    CHECK(extraction.edgeBuild.edgeUses==
          6*extraction.rebuildCellInsertions);
    CHECK(extraction.faceBuild.faceUses==
          4*extraction.rebuildCellInsertions);
    CHECK(extraction.edgeBuild.finalEdgeAudits==1);
    CHECK(extraction.faceBuild.suppliedEdgeAudits==1);
    CHECK(extraction.faceBuild.finalFaceAudits==1);
    CHECK(extraction.edgeBuild.perCellGlobalAudits==0);
    CHECK(extraction.faceBuild.perCellGlobalAudits==0);
    CHECK(result.strictHash==mesh::meshHash(destination));
    const auto semantic=mesh::semanticHashEx(destination);
    CHECK(result.semanticHash.value==semantic.value);
    CHECK(result.semanticHash.canonical==semantic.canonical);
    CHECK(!result.authorityCoreBytes.empty());
    std::uint64_t regions=0,facets=0;
    CHECK(algo::t5_detail::decodeAuthorityCoreForTesting(
        result.authorityCoreBytes,regions,facets));
    CHECK(regions==1);CHECK(facets==0);
    CHECK(result.sourceNodeMappings.size()==
          static_cast<std::size_t>(destination.nodeCount()));
    CHECK(result.sourceSurfaceCellMappings.size()+
          result.sourceTetMappings.size()==
          static_cast<std::size_t>(destination.cellCount()));
}

void twoRegionRoute() {
    TEST("t5 two-region interface has one identity and exact owners");
    LayeredModel model;Ownership provider;
    for(const auto* face:model.faces())
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(10.0);
    std::array<algo::T5RegionBinding,2> regions{{
        {0,1001,2001,{3,0}},{1,1002,2002,{3,1}}}};
    algo::T5Request request;
    request.sourceRequest.expectedSurfaceSetStamp={15003,1};
    request.sourceRequest.expectedOwnershipProviderStamp=provider.current;
    request.regionBindings=regions;
    mesh::Mesh destination;
    const auto outcome=algo::meshModelMultiRegionTet4(
        model,spec,request,destination,provider);
    if(!outcome.published())
        std::printf("  t5 two-region status=%d reason=%d phase=%d\n",
            static_cast<int>(outcome.status),
            static_cast<int>(outcome.failure.reason),
            static_cast<int>(outcome.failure.phase));
    CHECK(outcome.published());
    CHECK(outcome.retainedResult!=nullptr);
    if(!outcome.retainedResult)return;
    const auto& result=*outcome.retainedResult;
    CHECK(result.audit.passed);
    CHECK(result.regions.size()==2);
    CHECK(result.materials.size()==2);
    CHECK(!result.interfaceFacets.empty());
    CHECK(result.audit.prescribedInterface==result.interfaceFacets.size());
    CHECK(result.audit.matchedInterface==result.interfaceFacets.size());
    CHECK(result.interfaceNodeIncidentFacets.size()==
          3*result.interfaceFacets.size());
    for(const auto& row:result.interfaceFacets) {
        CHECK(row.firstRegion==1001);
        CHECK(row.secondRegion==1002);
        CHECK(row.firstOutwardSense==-row.secondOutwardSense);
        CHECK(row.publishedTri3Cell>=0);
    }
    for(const auto& row:result.interfaceNodes) {
        CHECK(row.incidentFacetOffset+row.incidentFacetCount<=
              result.interfaceNodeIncidentFacets.size());
        CHECK(row.incidentFacetCount>0);
    }
    for(const auto& row:result.interfaceEdges) {
        CHECK(row.required);
        CHECK(row.incidentArtifactFacetOffset+
              row.incidentArtifactFacetCount<=
              result.interfaceEdgeIncidentArtifactFacets.size());
    }
    std::size_t regionFacetReferences=0;
    for(const auto& row:result.regions) {
        CHECK(row.interfaceReferenceOffset+row.interfaceReferenceCount<=
              result.regionInterfaceFacets.size());
        regionFacetReferences+=static_cast<std::size_t>(
            row.interfaceReferenceCount);
    }
    CHECK(regionFacetReferences==2*result.interfaceFacets.size());
    requireCompleteXor(outcome);
}

void fineInterfaceRoute() {
    TEST("t5 fine split box publishes exactly thirty two interface facets");
    LayeredModel model;Ownership provider;
    for(const auto* face:model.faces())
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(0.5);
    std::array<algo::T5RegionBinding,2> regions{{
        {0,5101,6101,{3,0}},{1,5102,6102,{3,1}}}};
    algo::T5Request request;
    request.sourceRequest.expectedSurfaceSetStamp={15100,1};
    request.sourceRequest.expectedOwnershipProviderStamp=provider.current;
    request.regionBindings=regions;
    mesh::Mesh destination;
    auto outcome=algo::meshModelMultiRegionTet4(
        model,spec,request,destination,provider);
    if(!outcome.published())
        std::printf("  t5 fine status=%d reason=%d phase=%d\n",
            static_cast<int>(outcome.status),
            static_cast<int>(outcome.failure.reason),
            static_cast<int>(outcome.failure.phase));
    CHECK(outcome.published());
    CHECK(outcome.retainedResult!=nullptr);
    if(!outcome.retainedResult)return;
    CHECK(outcome.retainedResult->interfaceFacets.size()==32);
    CHECK(outcome.retainedResult->audit.passed);
    CHECK(outcome.retainedResult->metrics.perInterfaceGlobalTetScans==0);
    CHECK(outcome.retainedResult->metrics.perRegionGlobalFaceScans==0);
    CHECK(outcome.retainedResult->metrics.contextRebuilds==0);
}

void realLayerNonCycleRefusal() {
    TEST("t5 real eight-layer route pins nested non-cycle edge refusal");
    constexpr std::size_t regionCount=8;
    constexpr double extent=4.0;
    LayeredModel model(static_cast<int>(regionCount),false,extent,1.0);
    Ownership provider;
    for(const auto* face:model.faces())
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,extent);
    spec.threads=1;
    std::vector<algo::T5RegionBinding> regions;
    regions.reserve(regionCount);
    for(std::size_t i=0;i<regionCount;++i)
        regions.push_back({static_cast<int>(i),17000+i,19000+i,
                           {3,static_cast<int>(i)}});
    algo::T5Request request;
    request.sourceRequest.expectedSurfaceSetStamp={16800,1};
    request.sourceRequest.expectedOwnershipProviderStamp=provider.current;
    request.regionBindings=regions;

    mesh::Mesh destination;
    destination.reserveAppend(4,1,4);
    (void)destination.addNode({-7.0,-11.0,-13.0},{0,91});
    const auto beforeHash=mesh::meshHash(destination);
    const auto beforeCapacity=destination.capacitySnapshot();
    const auto beforeNodes=destination.nodeCount();
    const auto beforeCells=destination.cellCount();

    bnd::resetTetEdgeCavityDiagnosticsForTesting(
        true,bnd::VolumeArtifactEdgeKey::of(40,127));
    auto outcome=algo::meshModelMultiRegionTet4(
        model,spec,request,destination,provider);
    CHECK(outcome.status==algo::T5Status::Refused);
    CHECK(outcome.failure.reason==algo::T5Reason::SourceRefused);
    CHECK(outcome.failure.phase==algo::T5Phase::SourceConstruction);
    CHECK(outcome.failure.sourceStatus==
          algo::BoundaryConformingV3Status::Refused);
    const auto& failure=outcome.failure.sourceFailure;
    CHECK(failure.reason==algo::BoundaryConformingV3Reason::NonStarEdgeCavity);
    CHECK(failure.phase==algo::BoundaryConformingV3Phase::EdgeDiscovery);
    CHECK(failure.scheduleOrdinal==0);
    const std::array<Index,2> expectedEdge{{40,127}};
    CHECK(failure.edgeKey.canonicalGlobalNodes==expectedEdge);
    CHECK(failure.exactWitnessDigest==7551663973499184525ULL);
    CHECK(failure.edgeOutcome.status==bnd::TetEdgeRecoveryStatus::Refused);
    CHECK(failure.edgeOutcome.reason==bnd::TetEdgeRecoveryReason::NonCycleLink);
    CHECK(failure.edgeOutcome.phase==bnd::TetEdgeRecoveryPhase::Discovery);
    CHECK(outcome.sourceOutcome!=nullptr);
    CHECK(outcome.sourceReceipt==nullptr);
    CHECK(outcome.retainedResult==nullptr);
    CHECK(outcome.destinationUnchanged);
    CHECK(outcome.sourceOutcome&&outcome.sourceOutcome->destinationUnchanged);
    CHECK(outcome.sourceOutcome&&outcome.sourceOutcome->reportComplete);
    CHECK(outcome.sourceOutcome&&outcome.sourceOutcome->report!=nullptr);
    const auto& cavity=bnd::tetEdgeCavityDiagnosticsForTesting();
    CHECK(cavity.enabled);
    CHECK(cavity.complete);
    CHECK(!cavity.allocationFailed);
    CHECK(cavity.targetGlobal==failure.edgeKey);
    CHECK(cavity.targetDense==(std::array<std::uint32_t,2>{{44,48}}));
    CHECK(!cavity.initiallySelected.empty());
    CHECK(cavity.selected.size()>=cavity.initiallySelected.size());
    CHECK(cavity.selected.size()==
          cavity.initiallySelected.size()+cavity.closureAdded.size());
    CHECK(cavity.initiallySelected.size()==4);
    CHECK(cavity.closureAdded.empty());
    CHECK(cavity.selected.size()==4);
    CHECK(cavity.boundaryFaces.size()==10);
    CHECK(cavity.boundaryFacesByTargetEndpoints==
          (std::array<std::size_t,3>{{4,6,0}}));
    CHECK(cavity.ringEdges.size()==6);
    CHECK(cavity.linkVertices.size()==5);
    CHECK(cavity.linkComponents==1);
    CHECK(!cavity.linkIsSingleCycle);
    CHECK(cavity.ballValid);
    CHECK(cavity.ballReason==bnd::TetFacetGeometryReason::None);
    CHECK(cavity.protectedContactEdges.empty());
    CHECK(cavity.protectedBoundaryFaces.empty());
    CHECK(cavity.regionBoundaryFaces.empty());
    const std::array<std::uint32_t,7> expectedDense{{24,44,45,48,59,60,63}};
    const std::array<std::uint64_t,7> expectedSequence{{24,44,45,48,59,60,63}};
    const std::array<Vec3,7> expectedPosition{{
        {1,-2.6000000000000001,-2},{2,-4,-4},
        {2,-4,-1.3333333333333335},{2,-1.9000000000000004,-1},
        {3,-4,-4},{3,-4,-1.3333333333333335},
        {3,-1.9000000000000004,-1}}};
    CHECK(cavity.initialVertices.size()==expectedDense.size());
    for(std::size_t i=0;i<expectedDense.size();++i) {
        CHECK(cavity.initialVertices[i].denseVertex==expectedDense[i]);
        CHECK(cavity.initialVertices[i].sequence==expectedSequence[i]);
        CHECK(cavity.initialVertices[i].position.x==expectedPosition[i].x);
        CHECK(cavity.initialVertices[i].position.y==expectedPosition[i].y);
        CHECK(cavity.initialVertices[i].position.z==expectedPosition[i].z);
    }
    const std::array<std::uint64_t,4> expectedCellSequence{{260,262,263,266}};
    const std::array<std::array<std::uint32_t,4>,4> expectedCells{{
        {{24,59,60,44}},{{24,45,63,48}},
        {{24,60,63,45}},{{59,60,63,24}}}};
    CHECK(cavity.initialCells.size()==expectedCells.size());
    for(std::size_t i=0;i<expectedCells.size();++i) {
        CHECK(cavity.initialCells[i].sequence==expectedCellSequence[i]);
        CHECK(cavity.initialCells[i].denseNodes==expectedCells[i]);
    }
    const std::array<std::array<std::uint32_t,3>,10> expectedBoundary{{
        {{59,60,44}},{{24,44,60}},{{24,59,44}},{{45,63,48}},
        {{24,48,63}},{{24,45,48}},{{60,63,45}},{{24,60,45}},
        {{59,24,63}},{{59,63,60}}}};
    CHECK(cavity.boundaryFaces.size()==expectedBoundary.size());
    for(std::size_t i=0;i<expectedBoundary.size();++i)
        CHECK(cavity.boundaryFaces[i].directedDenseNodes==expectedBoundary[i]);
    const std::array<std::array<std::uint32_t,2>,6> expectedRing{{
        {{24,45}},{{24,59}},{{24,60}},{{24,63}},{{45,63}},{{59,60}}}};
    for(std::size_t i=0;i<expectedRing.size();++i) {
        CHECK(cavity.ringEdges[i].denseNodes==expectedRing[i]);
        CHECK(cavity.ringEdges[i].multiplicity==1);
    }
    const std::array<std::pair<std::uint32_t,std::size_t>,5> expectedLink{{
        {24,4},{45,2},{59,2},{60,2},{63,2}}};
    for(std::size_t i=0;i<expectedLink.size();++i) {
        CHECK(cavity.linkVertices[i].denseVertex==expectedLink[i].first);
        CHECK(cavity.linkVertices[i].degree==expectedLink[i].second);
        CHECK(cavity.linkVertices[i].component==0);
    }
    std::set<topo::CellHandle> diagnosticCells(
        cavity.initiallySelected.begin(),cavity.initiallySelected.end());
    CHECK(diagnosticCells.size()==cavity.initiallySelected.size());
    for(const auto cell:cavity.closureAdded)
        CHECK(diagnosticCells.insert(cell).second);
    const std::set<topo::CellHandle> selectedDiagnosticCells(
        cavity.selected.begin(),cavity.selected.end());
    CHECK(diagnosticCells==selectedDiagnosticCells);

    CHECK(!cavity.completionAvailable);
    CHECK(cavity.completionAdded.size()==746);
    CHECK(cavity.completionSelected.size()==750);
    std::set<topo::CellHandle> completionCells(
        cavity.selected.begin(),cavity.selected.end());
    for(const auto cell:cavity.completionAdded)
        CHECK(completionCells.insert(cell).second);
    CHECK(completionCells==std::set<topo::CellHandle>(
          cavity.completionSelected.begin(),cavity.completionSelected.end()));
    CHECK(cavity.completionPasses==15);
    CHECK(cavity.completionWorkUnits==16856);
    CHECK(cavity.completionBoundaryFacesByTargetEndpoints==
          (std::array<std::size_t,3>{{12,0,0}}));
    CHECK(cavity.completionLinkComponents==0);
    CHECK(!cavity.completionLinkIsSingleCycle);
    CHECK(cavity.completionBallValid);
    CHECK(cavity.completionBallReason==bnd::TetFacetGeometryReason::None);
    CHECK(cavity.completionBlockedByBoundary);
    CHECK(!cavity.completionBlockedByProtectedFace);
    CHECK(!cavity.completionBlockedByRegion);
    if(outcome.sourceOutcome&&outcome.sourceOutcome->report) {
        const auto& report=*outcome.sourceOutcome->report;
        CHECK(report.reportComplete);
        CHECK(report.processed.size()==1);
        CHECK(!report.unattempted.empty());
        std::set<std::tuple<int,Index,Index,
                            algo::FacePersistentId,std::size_t>> keys;
        std::size_t committed=0;
        const auto insertKey=[&](const algo::BoundaryConformingV3WorkKey& key) {
            const auto value=std::make_tuple(
                static_cast<int>(key.kind),
                key.edgeKey.canonicalGlobalNodes[0],
                key.edgeKey.canonicalGlobalNodes[1],
                key.facetKey.facePersistentId,
                key.facetKey.normalizedLocalCellOrdinal);
            CHECK(keys.insert(value).second);
        };
        for(const auto& row:report.processed) {
            insertKey(row.key);
            committed+=row.status==
                algo::BoundaryConformingV3WorkStatus::Marked||
                row.status==algo::BoundaryConformingV3WorkStatus::Recovered;
        }
        for(const auto& key:report.unattempted)insertKey(key);
        CHECK(keys.size()==report.processed.size()+report.unattempted.size());
        CHECK(committed==0);
        CHECK(report.processed.front().failure.edgeOutcome.reason==
              bnd::TetEdgeRecoveryReason::NonCycleLink);
    }
    CHECK(mesh::meshHash(destination)==beforeHash);
    CHECK(destination.nodeCount()==beforeNodes);
    CHECK(destination.cellCount()==beforeCells);
    const auto afterCapacity=destination.capacitySnapshot();
    CHECK(afterCapacity.nodeCoordinateCapacity==
          beforeCapacity.nodeCoordinateCapacity);
    CHECK(afterCapacity.nodeClassificationCapacity==
          beforeCapacity.nodeClassificationCapacity);
    CHECK(afterCapacity.cellCapacity==beforeCapacity.cellCapacity);
    CHECK(afterCapacity.connectivityCapacity==
          beforeCapacity.connectivityCapacity);
    CHECK(afterCapacity.dedupBucketCount==beforeCapacity.dedupBucketCount);
    CHECK(afterCapacity.dedupEntryCount==beforeCapacity.dedupEntryCount);
    CHECK(afterCapacity.dedupIndexCapacity==beforeCapacity.dedupIndexCapacity);
    CHECK(afterCapacity.complete==beforeCapacity.complete);
    bnd::resetTetEdgeCavityDiagnosticsForTesting(false);
}

void realLayerV4ProtectedSubdivision() {
    TEST("t2c v4 passes dim2 edge then preserves exact facet conflict cycle");
    constexpr std::size_t regionCount=8;
    constexpr double extent=4.0;
    LayeredModel model(static_cast<int>(regionCount),false,extent,1.0);
    Ownership ownership;
    for(const auto* face:model.faces())
        ownership.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    CurvePlacement placement;
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,extent);
    spec.threads=1;
    algo::BoundaryConformingV4Request request;
    request.sourceRequest.expectedSurfaceSetStamp={16800,1};
    request.sourceRequest.expectedOwnershipProviderStamp=ownership.current;
    request.expectedPlacementProviderStamp=placement.current;
    request.limits={128,129,128,1000000,UINT64_C(1)<<34,
                    UINT64_C(1)<<50,128,4096,4096};
    mesh::Mesh baseEdges,refinedEdges;
    algo::meshEdges(model,spec,baseEdges);
    algo::t2c_detail::PrivateBuildTrace edgeTrace;
    algo::t2c_detail::meshEdgesWithProtectedSplits(
        model,spec,{},refinedEdges,edgeTrace);
    CHECK(mesh::meshHash(baseEdges)==mesh::meshHash(refinedEdges));
    mesh::Mesh destination;
    (void)destination.addNode({-7,-11,-13},{0,91});
    const auto before=mesh::meshHash(destination);
    bnd::resetTetEdgeTransactionDiagnosticsForTesting(
        true,bnd::VolumeArtifactEdgeKey::of(44,131));
    bnd::resetTetFacetConflictCycleDiagnosticsForTesting(true);
    auto outcome=algo::meshModelBoundaryConformingV4(
        model,spec,request,destination,ownership,placement);
    std::printf("  t2c v4 status=%d reason=%d phase=%d splits=%zu rounds=%llu "
                "frontier=%llu unique=%llu deferred=%llu surface=%llu "
                "tets=%llu\n",
        static_cast<int>(outcome.status),static_cast<int>(outcome.reason),
        static_cast<int>(outcome.phase),
        outcome.retainedResult?outcome.retainedResult->subdivisionLedger.size():0,
        static_cast<unsigned long long>(outcome.retainedResult?
            outcome.retainedResult->metrics.rounds:0),
        static_cast<unsigned long long>(outcome.retainedResult?
            outcome.retainedResult->metrics.frontierRows:0),
        static_cast<unsigned long long>(outcome.retainedResult?
            outcome.retainedResult->metrics.uniquePlacements:0),
        static_cast<unsigned long long>(outcome.retainedResult?
            outcome.retainedResult->metrics.deferredFrontierRows:0),
        static_cast<unsigned long long>(outcome.retainedResult?
            outcome.retainedResult->metrics.rebuiltSurfaceNodes:0),
        static_cast<unsigned long long>(outcome.retainedResult?
            outcome.retainedResult->metrics.candidateTetCells:0));
    const auto* terminalSource=outcome.lastSourceOutcome
        ? outcome.lastSourceOutcome.get()
        : outcome.retainedResult&&outcome.retainedResult->sourceOutcome
            ? outcome.retainedResult->sourceOutcome.get():nullptr;
    if(terminalSource)
        std::printf("    terminal source status=%d reason=%d phase=%d "
                    "kind=%d edge=(%lld,%lld) facet=(%llu,%zu) "
                    "schedule=%zu witness=%llu nested=%d nested-phase=%d\n",
            static_cast<int>(terminalSource->status),
            static_cast<int>(terminalSource->failure.reason),
            static_cast<int>(terminalSource->failure.phase),
            static_cast<int>(terminalSource->failure.workKind),
            static_cast<long long>(terminalSource->failure.edgeKey.
                canonicalGlobalNodes[0]),
            static_cast<long long>(terminalSource->failure.edgeKey.
                canonicalGlobalNodes[1]),
            static_cast<unsigned long long>(
                terminalSource->failure.facetKey.facePersistentId),
            terminalSource->failure.facetKey.normalizedLocalCellOrdinal,
            terminalSource->failure.scheduleOrdinal,
            static_cast<unsigned long long>(
                terminalSource->failure.exactWitnessDigest),
            static_cast<int>(terminalSource->failure.edgeOutcome.reason),
            static_cast<int>(terminalSource->failure.edgeOutcome.phase));
    const algo::BoundaryConformingV3ReportRow* repairedEdgeRow=nullptr;
    if(terminalSource&&terminalSource->report)
        for(const auto& row:terminalSource->report->processed)
            if(row.key.kind==algo::VolumeWorkKind::ArtifactEdge&&
               row.key.edgeKey.canonicalGlobalNodes==
                   (std::array<Index,2>{{44,131}})) {
                repairedEdgeRow=&row;
                std::printf("    repaired edge row status=%d reason=%d "
                            "phase=%d schedule=%zu\n",
                    static_cast<int>(row.status),
                    static_cast<int>(row.failure.edgeOutcome.reason),
                    static_cast<int>(row.failure.edgeOutcome.phase),
                    row.scheduleOrdinal);
            }
    if(outcome.lastSourceOutcome)
        std::printf("    source status=%d reason=%d phase=%d edge=(%lld,%lld) "
                    "nested=%d\n",
            static_cast<int>(outcome.lastSourceOutcome->status),
            static_cast<int>(outcome.lastSourceOutcome->failure.reason),
            static_cast<int>(outcome.lastSourceOutcome->failure.phase),
            static_cast<long long>(outcome.lastSourceOutcome->failure.
                edgeKey.canonicalGlobalNodes[0]),
            static_cast<long long>(outcome.lastSourceOutcome->failure.
                edgeKey.canonicalGlobalNodes[1]),
            static_cast<int>(outcome.lastSourceOutcome->failure.
                edgeOutcome.reason));
    if(outcome.lastSourceOutcome&&outcome.lastSourceOutcome->surfaceOutcome&&
       outcome.lastSourceOutcome->surfaceOutcome->retainedArtifact) {
        const auto edge=outcome.lastSourceOutcome->failure.edgeKey.
            canonicalGlobalNodes;
        for(const auto& face:outcome.lastSourceOutcome->surfaceOutcome->
                retainedArtifact->faces)
            if(face)for(const auto& persistent:face->persistentEdges)
                if(persistent.publishedGlobalEndpoints) {
                    auto endpoints=*persistent.publishedGlobalEndpoints;
                    if(endpoints[1]<endpoints[0])std::swap(endpoints[0],endpoints[1]);
                    if(endpoints==edge)
                        std::printf("    matching surface edge class=(%d,%d) "
                                    "face=%llu chain=%llu ordinal=%zu\n",
                            persistent.classification.dim,
                            persistent.classification.tag,
                            static_cast<unsigned long long>(
                                persistent.faceStamp.persistentId),
                            static_cast<unsigned long long>(
                                persistent.chainPersistentId),
                            persistent.originalSegmentOrdinal);
                }
    }
    CHECK(!outcome.published());
    CHECK(outcome.status==algo::BoundaryConformingV4Status::Refused);
    CHECK(outcome.reason==
          algo::BoundaryConformingV4Reason::ChildFacetRecoveryRefused);
    CHECK(outcome.phase==
          algo::BoundaryConformingV4Phase::ConstraintRecovery);
    CHECK(outcome.destinationUnchanged);
    CHECK(outcome.retainedResult!=nullptr);
    CHECK(terminalSource!=nullptr);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->subdivisionLedger.size()==63);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->metrics.rounds==10);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->metrics.frontierRows==67);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->metrics.uniquePlacements==63);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->metrics.deferredFrontierRows==3);
    if(outcome.retainedResult) {
        const auto& cacheReport=outcome.retainedResult->surfaceCacheReport;
        std::printf("    cache report complete=%d failed=%d face=%llu "
                    "processed=%zu unattempted=%zu\n",
            cacheReport.reportComplete?1:0,cacheReport.hasFailedFace?1:0,
            static_cast<unsigned long long>(
                cacheReport.failedFacePersistentId),
            cacheReport.processed.size(),cacheReport.unattempted.size());
        for(const auto& row:cacheReport.processed)
            if(row.surfaceStatus!=algo::SurfaceFaceStatus::Ready)
                std::printf("    cache stopped face=%llu disposition=%d "
                            "status=%d reason=%d allocation=%d\n",
                    static_cast<unsigned long long>(row.facePersistentId),
                    static_cast<int>(row.disposition),
                    static_cast<int>(row.surfaceStatus),
                    static_cast<int>(row.originalSurfaceFailure.reason),
                    static_cast<int>(row.originalSurfaceFailure.allocation));
        CHECK(cacheReport.reportComplete);
        CHECK(!cacheReport.hasFailedFace);
        CHECK(cacheReport.unattempted.empty());
        CHECK(cacheReport.processed.size()==model.faces().size());
        requireSurfaceCacheReportXor(model,cacheReport);
        std::set<algo::FacePersistentId> cacheFaces;
        std::size_t cleanRows=0,dirtyRows=0;
        for(const auto& row:cacheReport.processed) {
            CHECK(row.surfaceStatus==algo::SurfaceFaceStatus::Ready);
            CHECK(row.originalSurfaceFailure.reason==
                  algo::SurfaceReason::None);
            CHECK(cacheFaces.insert(row.facePersistentId).second);
            cleanRows+=row.disposition==
                algo::SurfaceCacheFaceDisposition::CleanReplay;
            dirtyRows+=row.disposition==
                algo::SurfaceCacheFaceDisposition::DirtyBuild;
        }
        CHECK(cacheFaces.size()==model.faces().size());
        CHECK(cleanRows>0);
        CHECK(dirtyRows>0);
        CHECK(outcome.retainedResult->metrics.surfaceCacheFacesBuilt>0);
        CHECK(outcome.retainedResult->metrics.surfaceCacheFacesReused>0);
        CHECK(outcome.retainedResult->metrics.surfaceCacheFacesReplayed==
              outcome.retainedResult->metrics.surfaceCacheFacesReused);
        CHECK(outcome.retainedResult->metrics.surfaceCacheFacesBuilt+
                  outcome.retainedResult->metrics.surfaceCacheFacesReused==
              outcome.retainedResult->metrics.rounds*model.faces().size());
        CHECK(outcome.retainedResult->metrics.
              surfaceCacheSizingReplayCalls>0);
        CHECK(outcome.retainedResult->metrics.
              surfaceCachePreparationWorkUnits>0);
        CHECK(outcome.retainedResult->metrics.
              surfaceCacheReplayWorkUnits>0);
        CHECK(outcome.retainedResult->metrics.
              surfaceCacheSetAuditWorkUnits>0);
        CHECK(outcome.retainedResult->metrics.surfaceCacheCommits==10);
        CHECK(outcome.retainedResult->metrics.surfaceCacheCommits==
              outcome.retainedResult->metrics.rounds);
        CHECK(outcome.retainedResult->metrics.surfaceCacheCleanNodeRows>0);
        CHECK(outcome.retainedResult->metrics.surfaceCacheDirtyNodeRows>0);
        CHECK(outcome.retainedResult->metrics.surfaceCacheCleanCellRows>0);
        CHECK(outcome.retainedResult->metrics.surfaceCacheDirtyCellRows>0);
        CHECK(outcome.retainedResult->metrics.surfaceCacheAccountedBytes>0);
        CHECK(outcome.retainedResult->metrics.
              surfaceCacheAllFaceScansPerDirtyFace==0);
    }
    CHECK(placement.placementCalls==63);
    CHECK(repairedEdgeRow!=nullptr);
    CHECK(repairedEdgeRow&&
          repairedEdgeRow->status==
              algo::BoundaryConformingV3WorkStatus::Recovered);
    CHECK(repairedEdgeRow&&repairedEdgeRow->scheduleOrdinal==199);
    CHECK(repairedEdgeRow&&
          repairedEdgeRow->failure.edgeOutcome.reason==
              bnd::TetEdgeRecoveryReason::None);
    CHECK(terminalSource&&
          terminalSource->status==
              algo::BoundaryConformingV3Status::Refused);
    CHECK(terminalSource&&
          terminalSource->failure.reason==
              algo::BoundaryConformingV3Reason::RecoveryConflictCycle);
    CHECK(terminalSource&&
          terminalSource->failure.phase==
              algo::BoundaryConformingV3Phase::FacetRouting);
    CHECK(terminalSource&&
          terminalSource->failure.workKind==algo::VolumeWorkKind::Facet);
    CHECK(terminalSource&&
          terminalSource->failure.facetKey==
              (bnd::VolumeFacetKey{40,4}));
    CHECK(terminalSource&&
          terminalSource->failure.scheduleOrdinal==1109);
    CHECK(terminalSource&&
          terminalSource->failure.exactWitnessDigest==
              UINT64_C(12033182254251637294));
    CHECK(terminalSource&&
          terminalSource->failure.facetOutcome.status==
              bnd::TetFacetRecoveryStatus::Refused);
    CHECK(terminalSource&&
          terminalSource->failure.facetOutcome.reason==
              bnd::TetFacetRecoveryReason::RecoveryConflictCycle);
    CHECK(terminalSource&&
          terminalSource->failure.facetOutcome.phase==
              bnd::TetFacetRecoveryPhase::Routing);
    const auto& facetCycleDiagnostic=
        bnd::tetFacetConflictCycleDiagnosticsForTesting();
    CHECK(facetCycleDiagnostic.enabled);
    CHECK(facetCycleDiagnostic.captured);
    CHECK(facetCycleDiagnostic.complete);
    CHECK(facetCycleDiagnostic.facetCount==462);
    CHECK(facetCycleDiagnostic.pendingCount==2);
    CHECK(facetCycleDiagnostic.memberCount==2);
    CHECK(facetCycleDiagnostic.dependencyCount==2);
    CHECK(facetCycleDiagnostic.members[0]==
          (bnd::VolumeFacetKey{40,4}));
    CHECK(facetCycleDiagnostic.members[1]==
          (bnd::VolumeFacetKey{41,4}));
    CHECK(facetCycleDiagnostic.dependencies[0][0]==
          (bnd::VolumeFacetKey{40,4}));
    CHECK(facetCycleDiagnostic.dependencies[0][1]==
          (bnd::VolumeFacetKey{41,4}));
    CHECK(facetCycleDiagnostic.dependencies[1][0]==
          (bnd::VolumeFacetKey{41,4}));
    CHECK(facetCycleDiagnostic.dependencies[1][1]==
          (bnd::VolumeFacetKey{40,4}));
    CHECK(facetCycleDiagnostic.schedulerMetrics.routeAssessments==2241);
    CHECK(facetCycleDiagnostic.schedulerMetrics.dependencyEdges==9937);
    CHECK(facetCycleDiagnostic.schedulerMetrics.dependencyInvalidations==21);
    CHECK(facetCycleDiagnostic.schedulerMetrics.selectedRebuilds==1109);
    CHECK(facetCycleDiagnostic.schedulerMetrics.conflictSccVisits==2);
    CHECK(facetCycleDiagnostic.schedulerMetrics.fullStageRescans==0);
    CHECK(facetCycleDiagnostic.pendingAssessmentCounters.cellsScanned==280);
    CHECK(facetCycleDiagnostic.pendingAssessmentCounters.
          crossedEdgeStars==5);
    CHECK(facetCycleDiagnostic.pendingAssessmentCounters.fixedPointPasses==2);
    CHECK(facetCycleDiagnostic.pendingAssessmentCounters.apexAttempts==7);
    CHECK(facetCycleDiagnostic.pendingAssessmentCounters.
          tetTetIntersectionTests==2508);
    CHECK(facetCycleDiagnostic.pendingRetiredCellRows==12);
    CHECK(facetCycleDiagnostic.allAssessmentCounters.cellsScanned==1164);
    CHECK(facetCycleDiagnostic.allAssessmentCounters.crossedEdgeStars==14);
    CHECK(facetCycleDiagnostic.allAssessmentCounters.fixedPointPasses==11);
    CHECK(facetCycleDiagnostic.allAssessmentCounters.apexAttempts==25);
    CHECK(facetCycleDiagnostic.allAssessmentCounters.
          tetTetIntersectionTests==7220);
    CHECK(facetCycleDiagnostic.allRetiredCellRows==39);
    CHECK(facetCycleDiagnostic.maximumRetiredCells==9);
    const auto transactionDiagnostic=
        bnd::tetEdgeTransactionDiagnosticsForTesting();
    CHECK(transactionDiagnostic.enabled);
    CHECK(transactionDiagnostic.targetGlobal==
          bnd::VolumeArtifactEdgeKey::of(44,131));
    CHECK(transactionDiagnostic.transactionRequestStamp==
          transactionDiagnostic.transactionCurrentStamp);
    CHECK(transactionDiagnostic.transactionRemovedCount!=0);
    CHECK(transactionDiagnostic.transactionCandidateCount!=0);
    CHECK(transactionDiagnostic.transactionRemovedDigest!=0);
    CHECK(transactionDiagnostic.transactionCandidateDigest!=0);
    CHECK(transactionDiagnostic.transactionCheckpoint==
          bnd::TetEdgeTransactionDiagnosticCheckpoint::StageEntered);
    CHECK(transactionDiagnostic.transactionClassification.dim==2);
    CHECK(transactionDiagnostic.transactionClassification.tag>=0);
    CHECK(transactionDiagnostic.transactionRequiredConstraint);
    CHECK(transactionDiagnostic.transactionRegion.dim==3);
    CHECK(transactionDiagnostic.transactionRegion.tag>=0);
    CHECK(mesh::meshHash(destination)==before);
    bnd::resetTetFacetConflictCycleDiagnosticsForTesting(false);
    bnd::resetTetEdgeTransactionDiagnosticsForTesting(false);
    if(!outcome.published()||!outcome.retainedResult) {
        return;
    }
    CHECK(!outcome.retainedResult->subdivisionLedger.empty());
    CHECK(outcome.retainedResult->sourceOutcome!=nullptr);
    CHECK(outcome.retainedResult->sourceOutcome->retainedResult!=nullptr);
    if(outcome.retainedResult->sourceOutcome->retainedResult) {
        const auto& audit=outcome.retainedResult->sourceOutcome->
            retainedResult->boundaryAudit;
        CHECK(audit.missingInterface==0);
        CHECK(audit.unexpectedInterface==0);
        CHECK(audit.orphanExterior==0);
        CHECK(audit.unexpectedExterior==0);
        CHECK(audit.passed);
    }
    CHECK(placement.placementCalls==
          outcome.retainedResult->subdivisionLedger.size());
    CHECK(placement.stampCalls==
          outcome.retainedResult->metrics.providerStampCalls);
    CHECK(placement.placementCalls==
          outcome.retainedResult->metrics.providerPlacementCalls);
    CHECK(outcome.retainedResult->metrics.coordinatorWorkUnits==
          outcome.retainedResult->metrics.rounds+
          outcome.retainedResult->metrics.providerStampCalls+
          outcome.retainedResult->metrics.providerPlacementCalls);
    for(const auto& row:outcome.retainedResult->subdivisionLedger) {
        CHECK(row.classification.dim==(row.edgeTag>=0?1:2));
        CHECK(row.classification.tag==(row.edgeTag>=0?row.edgeTag:row.faceTag));
        if(row.edgeTag>=0) {
            CHECK(row.midpointParameter>row.parentParameterLo);
            CHECK(row.midpointParameter<row.parentParameterHi);
            CHECK(row.coordinate.x==model.edge(row.edgeTag)->
                  pointAt(row.midpointParameter).x);
        } else {
            CHECK(row.coordinate.x==model.face(row.faceTag)->
                  pointAt(row.faceParameter).x);
        }
        CHECK(row.providerStamp==placement.current);
    }
    const auto referenceStrict=mesh::meshHash(destination);
    const auto referenceSemantic=mesh::semanticHashEx(destination);
    const auto replay=[&](unsigned threads,bool reversed,bool churn) {
        if(churn) {
            LayeredModel dummy(2);Ownership dummyOwnership;
            for(const auto* face:dummy.faces())
                dummyOwnership.rows[
                    static_cast<algo::FacePersistentId>(face->tag())+1]=
                    dummy.owners(face->tag());
            sizing::MeshSpec dummySpec;
            dummySpec.field=std::make_shared<sizing::ConstantField>(10.0);
            algo::BoundaryConformingV3Request dummyRequest;
            dummyRequest.expectedSurfaceSetStamp={17000,1};
            dummyRequest.expectedOwnershipProviderStamp=
                dummyOwnership.current;
            mesh::Mesh ignored;
            (void)algo::meshModelBoundaryConforming(
                dummy,dummySpec,dummyRequest,ignored,dummyOwnership);
        }
        LayeredModel replayModel(static_cast<int>(regionCount),false,extent,
                                 1.0,reversed);
        Ownership replayOwnership;
        for(const auto* face:replayModel.faces())
            replayOwnership.rows[
                static_cast<algo::FacePersistentId>(face->tag())+1]=
                replayModel.owners(face->tag());
        std::vector<algo::t2c_detail::DyadicEdgeSplit> replayEdges;
        std::vector<algo::t2c_detail::PrivateFacePoint> replayFaces;
        for(const auto& row:outcome.retainedResult->subdivisionLedger) {
            if(row.edgeTag>=0) {
                replayEdges.push_back({row.edgeTag,
                    row.childIdentity.parentOrdinal,row.midpointParameter,
                    replayModel.edge(row.edgeTag)->pointAt(row.midpointParameter),
                    row.childIdentity});
            } else {
                const auto* face=replayModel.face(row.faceTag);
                const Vec2 uv=reversed?face->paramAt(row.coordinate):
                    row.faceParameter;
                replayFaces.push_back({
                    static_cast<algo::FacePersistentId>(row.faceTag)+1,
                    uv,face->pointAt(uv),row.childIdentity});
            }
        }
        sizing::MeshSpec replaySpec=spec;replaySpec.threads=threads;
        algo::BoundaryConformingV3Request replayRequest=
            request.sourceRequest;
        replayRequest.expectedOwnershipProviderStamp=replayOwnership.current;
        mesh::Mesh replayMesh;
        algo::BoundaryConformingV3Outcome replayOutcome;
        std::unique_ptr<algo::BoundaryConformingV4Outcome> reversedOutcome;
        if(reversed) {
            CurvePlacement replayPlacement;
            algo::BoundaryConformingV4Request v4Request=request;
            v4Request.sourceRequest=replayRequest;
            v4Request.expectedPlacementProviderStamp=replayPlacement.current;
            reversedOutcome=std::make_unique<algo::BoundaryConformingV4Outcome>(
                algo::meshModelBoundaryConformingV4(
                    replayModel,replaySpec,v4Request,replayMesh,
                    replayOwnership,replayPlacement));
        } else {
            algo::t2c_detail::PrivateBuildTrace trace;
            replayOutcome=algo::t2c_detail::
                meshModelBoundaryConformingWithProtectedSplits(
                    replayModel,replaySpec,replayRequest,replayMesh,
                    replayOwnership,replayEdges,replayFaces,&trace);
        }
        const bool replayPublished=reversed
            ? reversedOutcome&&reversedOutcome->published()
            : replayOutcome.published();
        std::printf("    replay-result threads=%u reversed=%d churn=%d "
                    "published=%d strict=%llu semantic=%llu canonical=%d\n",
            threads,reversed?1:0,churn?1:0,replayPublished?1:0,
            static_cast<unsigned long long>(mesh::meshHash(replayMesh)),
            static_cast<unsigned long long>(
                mesh::semanticHashEx(replayMesh).value),
            mesh::semanticHashEx(replayMesh).canonical?1:0);
        std::fflush(stdout);
        if(!replayPublished)
            std::printf("    replay threads=%u reversed=%d churn=%d status=%d "
                        "reason=%d phase=%d\n",threads,reversed?1:0,churn?1:0,
                reversed?static_cast<int>(reversedOutcome->status):
                         static_cast<int>(replayOutcome.status),
                reversed?static_cast<int>(reversedOutcome->reason):
                         static_cast<int>(replayOutcome.failure.reason),
                reversed?static_cast<int>(reversedOutcome->phase):
                         static_cast<int>(replayOutcome.failure.phase));
        CHECK(replayPublished);
        if(reversed) {
            CHECK(reversedOutcome&&reversedOutcome->retainedResult);
            CHECK(reversedOutcome&&reversedOutcome->retainedResult&&
                  reversedOutcome->retainedResult->sourceOutcome&&
                  reversedOutcome->retainedResult->sourceOutcome->
                      retainedResult&&
                  reversedOutcome->retainedResult->sourceOutcome->
                      retainedResult->boundaryAudit.passed);
        } else {
            CHECK(replayOutcome.reportComplete);
            CHECK(replayOutcome.retainedResult&&
                  replayOutcome.retainedResult->boundaryAudit.passed);
        }
        CHECK(mesh::meshHash(replayMesh)==referenceStrict);
        const auto semantic=mesh::semanticHashEx(replayMesh);
        CHECK(semantic.value==referenceSemantic.value);
        CHECK(semantic.canonical==referenceSemantic.canonical);
    };
    replay(4,false,false);
    replay(0,false,true);
    replay(1,true,false);
}

void realLayerV4FacetConflictSubdivision() {
    TEST("t2c opt-in canonical facet SCC subdivision publishes real layer");
    constexpr std::size_t regionCount=8;
    constexpr double extent=4.0;
    LayeredModel model(static_cast<int>(regionCount),false,extent,1.0);
    Ownership ownership;
    for(const auto* face:model.faces())
        ownership.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    CurvePlacement placement;
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,extent);
    spec.threads=1;
    algo::BoundaryConformingV4Request request;
    request.sourceRequest.expectedSurfaceSetStamp={16800,1};
    request.sourceRequest.expectedOwnershipProviderStamp=ownership.current;
    request.expectedPlacementProviderStamp=placement.current;
    request.limits={128,129,128,1000000,UINT64_C(1)<<34,
                    UINT64_C(1)<<50,128,4096,4096};
    request.facetConflictMode=algo::FacetConflictRefinementMode::
        BoundedCanonicalSccSubdivision;
    request.facetConflictLimits={2,2,1,7,UINT64_C(1)<<30};
    request.recoverySchedulingMode=algo::TetRecoverySchedulingMode::
        MixedExactPresentBatching;
    {
        CurvePlacement cappedPlacement;
        auto cappedRequest=request;
        cappedRequest.expectedPlacementProviderStamp=cappedPlacement.current;
        cappedRequest.facetConflictLimits.maxSccMembers=1;
        mesh::Mesh cappedDestination;
        (void)cappedDestination.addNode({-7,-11,-13},{0,91});
        const auto cappedBefore=mesh::meshHash(cappedDestination);
        const auto capped=algo::meshModelBoundaryConformingV4(
            model,spec,cappedRequest,cappedDestination,ownership,
            cappedPlacement);
        CHECK(capped.status==algo::BoundaryConformingV4Status::Incomplete);
        CHECK(capped.reason==
              algo::BoundaryConformingV4Reason::FacetConflictBudgetReached);
        CHECK(capped.phase==
              algo::BoundaryConformingV4Phase::ConstraintRecovery);
        CHECK(capped.destinationUnchanged);
        CHECK(mesh::meshHash(cappedDestination)==cappedBefore);
        CHECK(cappedPlacement.placementCalls==63);
        CHECK(cappedPlacement.requests.size()==63);
        CHECK(std::none_of(cappedPlacement.requests.begin(),
                           cappedPlacement.requests.end(),
            [](const auto& row) {
                return row.kind==
                    algo::ProtectedFeatureEntityKind::FacetInterior;
            }));
        CHECK(capped.retainedResult!=nullptr);
        CHECK(capped.retainedResult&&
              capped.retainedResult->subdivisionLedger.size()==63);
        CHECK(capped.retainedResult&&
              capped.retainedResult->surfaceCacheReport.reportComplete);
        CHECK(capped.retainedResult&&
              capped.retainedResult->metrics.facetConflictSccMembers==2);
        CHECK(capped.retainedResult&&
              capped.retainedResult->metrics.facetConflictDependencyRows==2);
        CHECK(capped.lastSourceOutcome!=nullptr||
              (capped.retainedResult&&
               capped.retainedResult->sourceOutcome!=nullptr));
    }
    mesh::Mesh destination;
    (void)destination.addNode({-7,-11,-13},{0,91});
    const auto before=mesh::meshHash(destination);
    const auto outcome=algo::meshModelBoundaryConformingV4(
        model,spec,request,destination,ownership,placement);
    CHECK(outcome.published());
    CHECK(outcome.status==algo::BoundaryConformingV4Status::Published);
    CHECK(outcome.reason==algo::BoundaryConformingV4Reason::None);
    CHECK(outcome.phase==algo::BoundaryConformingV4Phase::Publication);
    CHECK(!outcome.destinationUnchanged);
    CHECK(mesh::meshHash(destination)!=before);
    CHECK(destination.nodeCount()==198);
    CHECK(destination.countCells(mesh::CellType::Tet4)==777);
    CHECK(mesh::meshHash(destination)==UINT64_C(2737597274546495434));
    const auto semantic=mesh::semanticHashEx(destination);
    CHECK(semantic.canonical);
    CHECK(semantic.value==UINT64_C(16753220727311122020));
    CHECK(outcome.retainedResult!=nullptr);
    if(!outcome.retainedResult)return;
    const auto& retained=*outcome.retainedResult;
    CHECK(retained.subdivisionLedger.size()==65);
    CHECK(retained.metrics.rounds==11);
    CHECK(retained.metrics.edgeSplits==63);
    CHECK(retained.metrics.facetSplits==2);
    CHECK(retained.metrics.uniquePlacements==65);
    CHECK(retained.metrics.facetConflictRounds==1);
    CHECK(retained.metrics.facetConflictSccMembers==2);
    CHECK(retained.metrics.facetConflictDependencyRows==2);
    CHECK(retained.metrics.facetConflictProviderCalls==7);
    CHECK(retained.metrics.facetConflictCoordinatorWorkUnits>7);
    const auto& mixed=retained.metrics.mixedPresentScheduling;
    CHECK(mixed.complete);
    CHECK(mixed.presentBindingBatches>0);
    CHECK(mixed.presentBindingEntities>0);
    CHECK(mixed.stateCopies==mixed.presentBindingBatches);
    CHECK(mixed.exactPresentSelectedRebuilds==0);
    CHECK(mixed.exactPresentDependencyTests==0);
    CHECK(mixed.schedulerBindingAuditRows>=mixed.presentBindingEntities);
    const auto& latest=
        retained.metrics.latestSourceMixedPresentScheduling;
    CHECK(latest.complete);
    CHECK(latest.nontrivialRecoveries==18);
    CHECK(mixed.nontrivialRecoveries>=latest.nontrivialRecoveries);
    CHECK(latest.exactPresentSelectedRebuilds==0);
    CHECK(latest.exactPresentDependencyTests==0);
    CHECK(retained.metrics.surfaceCacheCommits==11);
    CHECK(retained.surfaceCacheReport.reportComplete);
    CHECK(!retained.surfaceCacheReport.hasFailedFace);
    CHECK(retained.surfaceCacheReport.unattempted.empty());
    CHECK(retained.surfaceCacheReport.processed.size()==model.faces().size());
    requireSurfaceCacheReportXor(model,retained.surfaceCacheReport);
    CHECK(retained.sourceOutcome!=nullptr);
    if(!retained.sourceOutcome)return;
    CHECK(retained.sourceOutcome->published());
    CHECK(retained.sourceOutcome->reportComplete);
    CHECK(retained.sourceOutcome->retainedResult!=nullptr);
    if(!retained.sourceOutcome->retainedResult)return;
    const auto& audit=retained.sourceOutcome->retainedResult->boundaryAudit;
    CHECK(audit.passed);
    CHECK(audit.orphanExterior==0);
    CHECK(audit.unexpectedExterior==0);
    CHECK(audit.missingInterface==0);
    CHECK(audit.unexpectedInterface==0);
    CHECK(placement.placementCalls==65);
    CHECK(placement.requests.size()==65);
    CHECK(placement.stampCalls==retained.metrics.providerStampCalls);
    CHECK(placement.placementCalls==
          retained.metrics.providerPlacementCalls);
    const auto& firstFacetRequest=placement.requests[63];
    const auto& secondFacetRequest=placement.requests[64];
    CHECK(firstFacetRequest.kind==
          algo::ProtectedFeatureEntityKind::FacetInterior);
    CHECK(secondFacetRequest.kind==
          algo::ProtectedFeatureEntityKind::FacetInterior);
    CHECK(firstFacetRequest.entityStamp==(algo::ArtifactStamp{40,1}));
    CHECK(secondFacetRequest.entityStamp==(algo::ArtifactStamp{41,1}));
    CHECK(firstFacetRequest.requestedIdentity.parentFacetGlobalNodes==
          (std::array<Index,3>{{77,169,171}}));
    CHECK(secondFacetRequest.requestedIdentity.parentFacetGlobalNodes==
          (std::array<Index,3>{{89,180,183}}));
    for(const auto* requestRow:{&firstFacetRequest,&secondFacetRequest}) {
        CHECK(requestRow->requestedIdentity.domain==
              algo::SurfaceNodeIdentityDomain::RefinedFacetInterior);
        CHECK(requestRow->requestedIdentity.parentOrdinal==0);
        CHECK(requestRow->requestedIdentity.parentSubordinal==0);
        CHECK(requestRow->requestedIdentity.dyadicNumerator==0);
        CHECK(requestRow->requestedIdentity.dyadicLevel==0);
        CHECK(requestRow->requestedIdentity.parameterUBits==
              UINT64_C(0x3fc5eeeeeeeeeeee));
        CHECK(requestRow->requestedIdentity.parameterVBits==
              UINT64_C(0x3fdabb67ae8584ca));
    }
    const auto& firstFacetLedger=retained.subdivisionLedger[63];
    const auto& secondFacetLedger=retained.subdivisionLedger[64];
    CHECK(firstFacetLedger.childIdentity==
          firstFacetRequest.requestedIdentity);
    CHECK(secondFacetLedger.childIdentity==
          secondFacetRequest.requestedIdentity);
    CHECK(firstFacetLedger.faceTag==39);
    CHECK(secondFacetLedger.faceTag==40);
    CHECK(firstFacetLedger.classification.dim==2);
    CHECK(firstFacetLedger.classification.tag==39);
    CHECK(secondFacetLedger.classification.dim==2);
    CHECK(secondFacetLedger.classification.tag==40);
    const auto firstExpected=model.face(39)->pointAt(
        firstFacetLedger.faceParameter);
    const auto secondExpected=model.face(40)->pointAt(
        secondFacetLedger.faceParameter);
    CHECK(firstFacetLedger.coordinate.x==firstExpected.x&&
          firstFacetLedger.coordinate.y==firstExpected.y&&
          firstFacetLedger.coordinate.z==firstExpected.z);
    CHECK(secondFacetLedger.coordinate.x==secondExpected.x&&
          secondFacetLedger.coordinate.y==secondExpected.y&&
          secondFacetLedger.coordinate.z==secondExpected.z);

    const auto referenceStrict=mesh::meshHash(destination);
    const auto referenceSemantic=mesh::semanticHashEx(destination);
    const auto replay=[&](unsigned threads,bool reversed,bool churn) {
        if(churn) {
            // Consume allocator/pool identities through a complete disjoint
            // route. Those runtime tokens must not leak into retained surface
            // or protected-child identity.
            LayeredModel dummy(2);Ownership dummyOwnership;
            for(const auto* face:dummy.faces())
                dummyOwnership.rows[
                    static_cast<algo::FacePersistentId>(face->tag())+1]=
                    dummy.owners(face->tag());
            sizing::MeshSpec dummySpec;
            dummySpec.field=std::make_shared<sizing::ConstantField>(10.0);
            algo::BoundaryConformingV3Request dummyRequest;
            dummyRequest.expectedSurfaceSetStamp={17003,1};
            dummyRequest.expectedOwnershipProviderStamp=
                dummyOwnership.current;
            mesh::Mesh ignored;
            (void)algo::meshModelBoundaryConforming(
                dummy,dummySpec,dummyRequest,ignored,dummyOwnership);
        }
        LayeredModel replayModel(static_cast<int>(regionCount),false,extent,
                                 1.0,reversed);
        Ownership replayOwnership;
        for(const auto* face:replayModel.faces())
            replayOwnership.rows[
                static_cast<algo::FacePersistentId>(face->tag())+1]=
                replayModel.owners(face->tag());
        CurvePlacement replayPlacement;
        sizing::MeshSpec replaySpec;
        replaySpec.field=std::make_shared<LayerQualificationField>(1.0,extent);
        replaySpec.threads=threads;
        auto replayRequest=request;
        replayRequest.sourceRequest.expectedOwnershipProviderStamp=
            replayOwnership.current;
        replayRequest.expectedPlacementProviderStamp=replayPlacement.current;
        mesh::Mesh replayDestination;
        (void)replayDestination.addNode({-7,-11,-13},{0,91});
        const auto replayOutcome=algo::meshModelBoundaryConformingV4(
            replayModel,replaySpec,replayRequest,replayDestination,
            replayOwnership,replayPlacement);
        CHECK(replayOutcome.published());
        CHECK(replayOutcome.status==outcome.status);
        CHECK(replayOutcome.reason==outcome.reason);
        CHECK(replayOutcome.phase==outcome.phase);
        CHECK(replayOutcome.destinationUnchanged==outcome.destinationUnchanged);
        CHECK(mesh::meshHash(replayDestination)==referenceStrict);
        const auto replaySemantic=mesh::semanticHashEx(replayDestination);
        CHECK(replaySemantic.value==referenceSemantic.value);
        CHECK(replaySemantic.canonical==referenceSemantic.canonical);
        CHECK(replayOutcome.retainedResult!=nullptr);
        if(!replayOutcome.retainedResult)return;
        const auto& replayRetained=*replayOutcome.retainedResult;
        CHECK(sameSubdivisionLedger(retained.subdivisionLedger,
                                    replayRetained.subdivisionLedger));
        CHECK(sameSurfaceCacheReport(retained.surfaceCacheReport,
                                     replayRetained.surfaceCacheReport));
        requireSurfaceCacheReportXor(replayModel,
                                     replayRetained.surfaceCacheReport);
        CHECK(sameV4Metrics(retained.metrics,replayRetained.metrics));
        CHECK(replayPlacement.stampCalls==placement.stampCalls);
        CHECK(replayPlacement.placementCalls==placement.placementCalls);
        CHECK(sameProviderRequests(placement.requests,
                                   replayPlacement.requests));
        CHECK(replayRetained.sourceOutcome!=nullptr);
        if(!replayRetained.sourceOutcome)return;
        const auto& referenceSource=*retained.sourceOutcome;
        const auto& replaySource=*replayRetained.sourceOutcome;
        CHECK(replaySource.status==referenceSource.status);
        CHECK(sameV3Failure(replaySource.failure,referenceSource.failure));
        CHECK(replaySource.destinationUnchanged==
              referenceSource.destinationUnchanged);
        CHECK(sameCompleteV3Report(referenceSource,replaySource));
        CHECK(static_cast<bool>(referenceSource.surfaceOutcome)==
              static_cast<bool>(replaySource.surfaceOutcome));
        if(referenceSource.surfaceOutcome&&replaySource.surfaceOutcome)
            CHECK(sameRetainedSurfaceSet(*referenceSource.surfaceOutcome,
                                         *replaySource.surfaceOutcome));
        CHECK(static_cast<bool>(referenceSource.retainedResult)==
              static_cast<bool>(replaySource.retainedResult));
        if(referenceSource.retainedResult&&replaySource.retainedResult) {
            CHECK(sameRecoveryMetrics(
                referenceSource.retainedResult->recoveryMetrics,
                replaySource.retainedResult->recoveryMetrics));
            CHECK(referenceSource.retainedResult->globalByAuthoritativeNode==
                  replaySource.retainedResult->globalByAuthoritativeNode);
            CHECK(referenceSource.retainedResult->globalByAuthoritativeCell==
                  replaySource.retainedResult->globalByAuthoritativeCell);
        }
    };
    replay(4,false,false);
    replay(0,true,false);
    replay(1,false,true);
}

void facetConflictRequestPreflight() {
    TEST("t2c facet-conflict opt-in requires positive finite limits");
    LayeredModel model(2,false,4.0,1.0);
    Ownership ownership;
    for(const auto* face:model.faces())
        ownership.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,4.0);
    CurvePlacement placement;
    algo::BoundaryConformingV4Request request;
    request.sourceRequest.expectedSurfaceSetStamp={16800,1};
    request.sourceRequest.expectedOwnershipProviderStamp=ownership.current;
    request.expectedPlacementProviderStamp=placement.current;
    request.limits={128,129,128,1000000,UINT64_C(1)<<34,
                    UINT64_C(1)<<50,128,4096,4096};
    request.facetConflictMode=algo::FacetConflictRefinementMode::
        BoundedCanonicalSccSubdivision;
    request.facetConflictLimits={2,2,1,7,1024};
    for(int field=0;field<5;++field) {
        auto invalid=request;
        if(field==0)invalid.facetConflictLimits.maxSccMembers=0;
        if(field==1)invalid.facetConflictLimits.maxFacetInteriorPoints=0;
        if(field==2)invalid.facetConflictLimits.maxFacetConflictRounds=0;
        if(field==3)invalid.facetConflictLimits.maxProviderCalls=0;
        if(field==4)invalid.facetConflictLimits.maxCoordinatorWorkUnits=0;
        mesh::Mesh destination;
        const auto outcome=algo::meshModelBoundaryConformingV4(
            model,spec,invalid,destination,ownership,placement);
        CHECK(outcome.status==algo::BoundaryConformingV4Status::Rejected);
        CHECK(outcome.reason==algo::BoundaryConformingV4Reason::InvalidRequest);
        CHECK(outcome.phase==algo::BoundaryConformingV4Phase::Preflight);
        CHECK(outcome.destinationUnchanged);
        CHECK(placement.stampCalls==0);
        CHECK(placement.placementCalls==0);
    }
    auto defaultRequest=request;
    defaultRequest.facetConflictMode=algo::FacetConflictRefinementMode::Refuse;
    defaultRequest.facetConflictLimits={};
    defaultRequest.limits.maxRounds=1;
    mesh::Mesh defaultDestination;
    const auto defaultOutcome=algo::meshModelBoundaryConformingV4(
        model,spec,defaultRequest,defaultDestination,ownership,placement);
    CHECK(defaultOutcome.status!=algo::BoundaryConformingV4Status::Rejected||
          defaultOutcome.reason!=algo::BoundaryConformingV4Reason::
              InvalidRequest);
}

void facetConflictProviderFailures() {
    TEST("t2c facet-conflict provider protocol is atomic at every SCC position");
    constexpr std::size_t regionCount=8;
    constexpr double extent=4.0;
    LayeredModel model(static_cast<int>(regionCount),false,extent,1.0);
    Ownership ownership;
    for(const auto* face:model.faces())
        ownership.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,extent);
    spec.threads=1;
    const auto makeRequest=[&](const CurvePlacement& placement) {
        algo::BoundaryConformingV4Request request;
        request.sourceRequest.expectedSurfaceSetStamp={16800,1};
        request.sourceRequest.expectedOwnershipProviderStamp=ownership.current;
        request.expectedPlacementProviderStamp=placement.current;
        request.limits={128,129,128,1000000,UINT64_C(1)<<34,
                        UINT64_C(1)<<50,128,4096,4096};
        request.facetConflictMode=algo::FacetConflictRefinementMode::
            BoundedCanonicalSccSubdivision;
        request.facetConflictLimits={2,2,1,7,UINT64_C(1)<<30};
        return request;
    };
    struct Row {
        CurvePlacement::FacetFault fault;
        std::size_t ordinal;
        algo::BoundaryConformingV4Status status;
        algo::BoundaryConformingV4Reason reason;
    };
    const std::array rows{
        Row{CurvePlacement::FacetFault::WrongEntity,1,
            algo::BoundaryConformingV4Status::Failed,
            algo::BoundaryConformingV4Reason::PlacementContractViolation},
        Row{CurvePlacement::FacetFault::WrongParameter,2,
            algo::BoundaryConformingV4Status::Failed,
            algo::BoundaryConformingV4Reason::PlacementContractViolation},
        Row{CurvePlacement::FacetFault::OffEntity,2,
            algo::BoundaryConformingV4Status::Failed,
            algo::BoundaryConformingV4Reason::PlacementContractViolation},
        Row{CurvePlacement::FacetFault::Rejected,2,
            algo::BoundaryConformingV4Status::Rejected,
            algo::BoundaryConformingV4Reason::ProviderRejected},
        Row{CurvePlacement::FacetFault::Throw,2,
            algo::BoundaryConformingV4Status::Failed,
            algo::BoundaryConformingV4Reason::ProviderFailed},
        Row{CurvePlacement::FacetFault::BadAlloc,2,
            algo::BoundaryConformingV4Status::Failed,
            algo::BoundaryConformingV4Reason::AllocationFailed}};
    for(const auto& row:rows) {
        CurvePlacement placement;
        placement.facetFault=row.fault;
        placement.facetFaultAt=row.ordinal;
        mesh::Mesh destination;
        (void)destination.addNode({-7,-11,-13},{0,91});
        const auto before=mesh::meshHash(destination);
        const auto outcome=algo::meshModelBoundaryConformingV4(
            model,spec,makeRequest(placement),destination,ownership,placement);
        CHECK(outcome.status==row.status);
        CHECK(outcome.reason==row.reason);
        CHECK(outcome.phase==algo::BoundaryConformingV4Phase::EdgePlacement);
        CHECK(outcome.destinationUnchanged);
        CHECK(mesh::meshHash(destination)==before);
        CHECK(placement.facetPlacementCalls==row.ordinal);
        CHECK(placement.placementCalls==63+row.ordinal);
        CHECK(outcome.retainedResult!=nullptr);
        if(!outcome.retainedResult)continue;
        const auto& retained=*outcome.retainedResult;
        CHECK(retained.subdivisionLedger.size()==63);
        CHECK(retained.metrics.edgeSplits==63);
        CHECK(retained.metrics.facetSplits==0);
        CHECK(retained.metrics.facetConflictRounds==1);
        CHECK(retained.surfaceCacheReport.reportComplete);
        requireSurfaceCacheReportXor(model,retained.surfaceCacheReport);
        CHECK(retained.sourceOutcome!=nullptr||outcome.lastSourceOutcome!=nullptr);
        CHECK(std::none_of(retained.subdivisionLedger.begin(),
                           retained.subdivisionLedger.end(),[](const auto& item) {
            return item.childIdentity.domain==
                algo::SurfaceNodeIdentityDomain::RefinedFacetInterior;
        }));
    }

    // Both the outer attempt needed to audit the proposed surface and the
    // complete two-row provider protocol (including its final checkpoint)
    // are secured before the first facet placement.
    for(int budgetKind=0;budgetKind<2;++budgetKind) {
        CurvePlacement placement;
        auto limited=makeRequest(placement);
        if(budgetKind==0)limited.limits.maxRounds=10;
        else limited.facetConflictLimits.maxProviderCalls=6;
        mesh::Mesh destination;
        (void)destination.addNode({-7,-11,-13},{0,91});
        const auto before=mesh::meshHash(destination);
        const auto outcome=algo::meshModelBoundaryConformingV4(
            model,spec,limited,destination,ownership,placement);
        CHECK(outcome.status==algo::BoundaryConformingV4Status::Incomplete);
        CHECK(outcome.reason==
              algo::BoundaryConformingV4Reason::FacetConflictBudgetReached);
        CHECK(outcome.phase==
              algo::BoundaryConformingV4Phase::ConstraintRecovery);
        CHECK(outcome.destinationUnchanged);
        CHECK(mesh::meshHash(destination)==before);
        CHECK(placement.facetPlacementCalls==0);
        CHECK(placement.placementCalls==63);
        CHECK(outcome.retainedResult!=nullptr);
        CHECK(outcome.retainedResult&&
              outcome.retainedResult->subdivisionLedger.size()==63);
        CHECK(outcome.retainedResult&&
              outcome.retainedResult->metrics.facetSplits==0);
        CHECK(outcome.retainedResult&&
              outcome.retainedResult->surfaceCacheReport.reportComplete);
        if(outcome.retainedResult)
            requireSurfaceCacheReportXor(
                model,outcome.retainedResult->surfaceCacheReport);
    }

    // An identical authoritative facet child is structural no-progress, not
    // another placement attempt and not a budget result.  Seed only the
    // fixed-value testing discriminator; the production ledger remains
    // untouched and the real SCC must resolve to this exact identity.
    algo::SurfaceNodeIdentity repeatedIdentity;
    repeatedIdentity.domain=
        algo::SurfaceNodeIdentityDomain::RefinedFacetInterior;
    repeatedIdentity.parentArtifact={40,1};
    repeatedIdentity.parentFacetGlobalNodes={{77,169,171}};
    repeatedIdentity.parameterUBits=UINT64_C(0x3fc5eeeeeeeeeeee);
    repeatedIdentity.parameterVBits=UINT64_C(0x3fdabb67ae8584ca);
    algo::t2c_detail::setFacetConflictRepeatedIdentityForTesting(
        repeatedIdentity);
    CurvePlacement repeatedPlacement;
    mesh::Mesh repeatedDestination;
    (void)repeatedDestination.addNode({-7,-11,-13},{0,91});
    const auto repeatedBefore=mesh::meshHash(repeatedDestination);
    const auto repeatedOutcome=algo::meshModelBoundaryConformingV4(
        model,spec,makeRequest(repeatedPlacement),repeatedDestination,
        ownership,repeatedPlacement);
    algo::t2c_detail::resetFacetConflictRepeatedIdentityForTesting();
    CHECK(repeatedOutcome.status==algo::BoundaryConformingV4Status::Incomplete);
    CHECK(repeatedOutcome.reason==
          algo::BoundaryConformingV4Reason::NoProgressCycle);
    CHECK(repeatedOutcome.phase==
          algo::BoundaryConformingV4Phase::ConstraintRecovery);
    CHECK(repeatedOutcome.destinationUnchanged);
    CHECK(mesh::meshHash(repeatedDestination)==repeatedBefore);
    CHECK(repeatedPlacement.facetPlacementCalls==0);
    CHECK(repeatedPlacement.placementCalls==63);
    CHECK(repeatedOutcome.retainedResult!=nullptr);
    CHECK(repeatedOutcome.retainedResult&&
          repeatedOutcome.retainedResult->subdivisionLedger.size()==63);
    CHECK(repeatedOutcome.retainedResult&&
          repeatedOutcome.retainedResult->metrics.facetSplits==0);
    CHECK(repeatedOutcome.retainedResult&&
          repeatedOutcome.retainedResult->surfaceCacheReport.reportComplete);
    if(repeatedOutcome.retainedResult)
        requireSurfaceCacheReportXor(
            model,repeatedOutcome.retainedResult->surfaceCacheReport);

    // The post-audit checkpoint is the final provider event.  It retains the
    // complete private result and both facet rows but cannot publish it under
    // a stale authority revision.
    CurvePlacement finalPlacement;
    finalPlacement.mutateAtStampCall=131;
    mesh::Mesh finalDestination;
    (void)finalDestination.addNode({-7,-11,-13},{0,91});
    const auto finalBefore=mesh::meshHash(finalDestination);
    const auto finalOutcome=algo::meshModelBoundaryConformingV4(
        model,spec,makeRequest(finalPlacement),finalDestination,ownership,
        finalPlacement);
    CHECK(finalOutcome.status==algo::BoundaryConformingV4Status::Failed);
    CHECK(finalOutcome.reason==
          algo::BoundaryConformingV4Reason::ProviderRevisionChanged);
    CHECK(finalOutcome.phase==algo::BoundaryConformingV4Phase::Publication);
    CHECK(finalOutcome.destinationUnchanged);
    CHECK(mesh::meshHash(finalDestination)==finalBefore);
    CHECK(finalPlacement.facetPlacementCalls==2);
    CHECK(finalPlacement.placementCalls==65);
    CHECK(finalPlacement.stampCalls==131);
    CHECK(finalOutcome.retainedResult!=nullptr);
    CHECK(finalOutcome.retainedResult&&
          finalOutcome.retainedResult->subdivisionLedger.size()==65);
    CHECK(finalOutcome.retainedResult&&
          finalOutcome.retainedResult->metrics.facetSplits==2);
    CHECK(finalOutcome.retainedResult&&
          finalOutcome.retainedResult->sourceOutcome&&
          finalOutcome.retainedResult->sourceOutcome->published());
}

void allDirtyCacheMatchesUncachedSurface() {
    TEST("t2c all-dirty cache build exactly matches uncached surface");
    LayeredModel model(2,false,4.0,1.0);
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,4.0);
    spec.threads=1;
    constexpr algo::ArtifactStamp stamp{16790,1};
    mesh::Mesh uncachedMesh,cachedMesh;
    algo::meshEdges(model,spec,uncachedMesh);
    algo::meshEdges(model,spec,cachedMesh);
    const auto uncached=algo::detail::buildProductionSurfaceSetForV3(
        model,spec,stamp,uncachedMesh);

    algo::t2c_detail::SurfaceFaceCacheSession cache;
    std::vector<algo::SurfaceCacheFaceReportRow> processed(
        model.faces().size());
    std::vector<algo::FacePersistentId> unattempted(model.faces().size());
    algo::t2c_detail::SurfaceCacheDiagnosticChannel channel{
        processed,unattempted};
    cache.attachDiagnosticChannel(&channel);
    const auto cached=algo::detail::buildProductionSurfaceSetForV3(
        model,spec,stamp,cachedMesh,{}, {},&cache);

    CHECK(uncached.surface.setStatus==algo::SurfaceSetStatus::Published);
    CHECK(cached.surface.setStatus==algo::SurfaceSetStatus::Published);
    CHECK(uncached.surface.retainedArtifact!=nullptr);
    CHECK(cached.surface.retainedArtifact!=nullptr);
    CHECK(uncached.surface.reportComplete==cached.surface.reportComplete);
    CHECK(uncached.surface.setStatus==cached.surface.setStatus);
    CHECK(sameSurfaceFailure(uncached.surface.setFailure,
                             cached.surface.setFailure));
    CHECK(uncached.surface.processed.size()==cached.surface.processed.size());
    CHECK(uncached.surface.unattempted==cached.surface.unattempted);
    for(std::size_t i=0;i<uncached.surface.processed.size()&&
                            i<cached.surface.processed.size();++i) {
        const auto& x=uncached.surface.processed[i];
        const auto& y=cached.surface.processed[i];
        CHECK(x.facePersistentId==y.facePersistentId);
        CHECK(x.status==y.status);
        CHECK(sameSurfaceFailure(x.failure,y.failure));
        CHECK(static_cast<bool>(x.readyResult)==
              static_cast<bool>(y.readyResult));
    }
    CHECK(mesh::meshHash(uncachedMesh)==mesh::meshHash(cachedMesh));
    CHECK(mesh::semanticHashEx(uncachedMesh).value==
          mesh::semanticHashEx(cachedMesh).value);
    CHECK(channel.reportComplete);
    CHECK(channel.processedUsed==model.faces().size());
    CHECK(channel.unattemptedUsed==0);
    for(std::size_t i=0;i<channel.processedUsed;++i) {
        CHECK(channel.processedStorage[i].disposition==
              algo::SurfaceCacheFaceDisposition::DirtyBuild);
        CHECK(channel.processedStorage[i].surfaceStatus==
              algo::SurfaceFaceStatus::Ready);
        CHECK(channel.processedStorage[i].originalSurfaceFailure.reason==
              algo::SurfaceReason::None);
    }
    CHECK(cache.metrics().facesBuilt==model.faces().size());
    CHECK(cache.metrics().facesReused==0);
    CHECK(cache.metrics().cacheCommits==1);
    if(!uncached.surface.retainedArtifact||
       !cached.surface.retainedArtifact)return;
    const auto& a=*uncached.surface.retainedArtifact;
    const auto& b=*cached.surface.retainedArtifact;
    CHECK(a.setStamp==b.setStamp);
    CHECK(a.faces.size()==b.faces.size());
    CHECK(sameSharedLedger(a.sharedLedger,b.sharedLedger));
    CHECK(a.publishedGlobalByRefinedIdentity==
          b.publishedGlobalByRefinedIdentity);
    const auto slots=cache.slotsForTesting();
    const auto slotFaces=cache.canonicalFacesForTesting();
    CHECK(slots.size()==a.faces.size());
    CHECK(slotFaces.size()==slots.size());
    for(std::size_t i=0;i<a.faces.size()&&i<b.faces.size();++i) {
        CHECK(a.faces[i]!=nullptr); CHECK(b.faces[i]!=nullptr);
        if(!a.faces[i]||!b.faces[i])continue;
        CHECK(a.faces[i]->faceStamp==b.faces[i]->faceStamp);
        CHECK(a.faces[i]->localNodes.size()==b.faces[i]->localNodes.size());
        CHECK(a.faces[i]->localCells.size()==b.faces[i]->localCells.size());
        CHECK(a.faces[i]->boundaryLoops==b.faces[i]->boundaryLoops);
        CHECK(a.faces[i]->publishedGlobalByLocalNode==
              b.faces[i]->publishedGlobalByLocalNode);
        CHECK(a.faces[i]->publishedGlobalCellByLocalCell==
              b.faces[i]->publishedGlobalCellByLocalCell);
        CHECK(a.faces[i]->refinedIdentityByLocalNode==
              b.faces[i]->refinedIdentityByLocalNode);
        CHECK(a.faces[i]->existingGlobalIndex==b.faces[i]->existingGlobalIndex);
        CHECK(a.faces[i]->persistentEdges.size()==
              b.faces[i]->persistentEdges.size());
        for(std::size_t edge=0;edge<a.faces[i]->persistentEdges.size()&&
                                    edge<b.faces[i]->persistentEdges.size();
                                    ++edge)
            CHECK(a.faces[i]->persistentEdges[edge].publishedGlobalEndpoints==
                  b.faces[i]->persistentEdges[edge].publishedGlobalEndpoints);

        const auto faceId=a.faces[i]->faceStamp.persistentId;
        const auto slotPos=std::lower_bound(
            slotFaces.begin(),slotFaces.end(),faceId);
        CHECK(slotPos!=slotFaces.end()&&*slotPos==faceId);
        if(slotPos==slotFaces.end()||*slotPos!=faceId)continue;
        const auto slot=static_cast<std::size_t>(slotPos-slotFaces.begin());
        CHECK(slot<slots.size());CHECK(slots[slot]!=nullptr);
        if(slot>=slots.size()||!slots[slot])continue;
        algo::ProductionSurfaceFaceMetrics metric{};
        const auto metricRow=std::find_if(uncached.metrics.begin(),
            uncached.metrics.end(),[&](const auto& row) {
                return row.facePersistentId==faceId;
            });
        CHECK(metricRow!=uncached.metrics.end());
        if(metricRow!=uncached.metrics.end())metric=*metricRow;
        const auto uncachedSnapshot=
            algo::t2c_detail::snapshotSurfaceFaceForCache(
                *a.faces[i],slots[slot]->authorityBytes,
                slots[slot]->sizingSamples,metric,a.sharedLedger);
        CHECK(uncachedSnapshot!=nullptr);
        CHECK(uncachedSnapshot&&sameCachedArtifact(
              *uncachedSnapshot,*slots[slot]));
    }
}

void cachedArtifactCorruptionIsRejected() {
    TEST("t2c cache rejects every retained-artifact corruption class");
    LayeredModel model(2,false,4.0,1.0);
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,4.0);
    spec.threads=1;
    constexpr algo::ArtifactStamp stamp{16792,1};
    mesh::Mesh destination;
    algo::meshEdges(model,spec,destination);
    algo::t2c_detail::SurfaceFaceCacheSession cache;
    std::vector<algo::SurfaceCacheFaceReportRow> processed(
        model.faces().size());
    std::vector<algo::FacePersistentId> unattempted(model.faces().size());
    algo::t2c_detail::SurfaceCacheDiagnosticChannel channel{
        processed,unattempted};
    cache.attachDiagnosticChannel(&channel);
    const auto built=algo::detail::buildProductionSurfaceSetForV3(
        model,spec,stamp,destination,{}, {},&cache);
    CHECK(built.surface.setStatus==algo::SurfaceSetStatus::Published);
    CHECK(built.surface.retainedArtifact!=nullptr);
    const auto slots=cache.slotsForTesting();
    const auto selected=std::find_if(slots.begin(),slots.end(),
        [](const auto& entry) {
            return entry&&!entry->localNodes.empty()&&
                !entry->localCells.empty()&&!entry->persistentEdges.empty()&&
                !entry->sharedKeysByLocalNode.empty()&&
                !entry->boundaryLoops.empty()&&
                !entry->localMappingParams.empty();
        });
    CHECK(selected!=slots.end());
    if(selected==slots.end()||!(*selected)||
       !built.surface.retainedArtifact)return;
    const auto& ledger=built.surface.retainedArtifact->sharedLedger;
    const auto rejects=[&](auto mutate) {
        auto corrupted=std::make_shared<
            algo::t2c_detail::CachedSurfaceFaceArtifact>(**selected);
        mutate(*corrupted);
        try {
            auto restored=algo::t2c_detail::rehydrateSurfaceFaceFromCache(
                *corrupted,ledger,destination);
            return !restored||!algo::t2c_detail::
                auditRehydratedSurfaceFace(*restored,destination);
        } catch(...) { return true; }
    };

    CHECK(rejects([](auto& row) {
        const auto local=row.sharedKeysByLocalNode.front().first;
        row.localNodes[local].coord.x+=1.0;
    }));
    CHECK(rejects([](auto& row) {
        std::swap(row.localCells.front().nodes[1],
                  row.localCells.front().nodes[2]);
    }));
    CHECK(rejects([](auto& row) {
        row.persistentEdges.front().required=false;
    }));
    CHECK(rejects([](auto& row) {
        ++row.localNodes.front().cls.tag;
    }));
    CHECK(rejects([](auto& row) {
        row.sharedKeysByLocalNode.front().second.nodeOrdinal+=1000000U;
    }));
    CHECK(rejects([](auto& row) {
        std::reverse(row.boundaryLoops.front().begin(),
                     row.boundaryLoops.front().end());
    }));
    CHECK(rejects([](auto& row) {
        row.localMappingParams.front().u=std::nextafter(
            row.localMappingParams.front().u,
            (std::numeric_limits<double>::infinity)());
    }));
}

void sharedCurveDirtyPropagationIsExact() {
    TEST("t2c shared-curve split dirties exactly its incident faces");
    LayeredModel model(2,false,4.0,1.0);
    Ownership ownership;
    for(const auto* face:model.faces())
        ownership.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    sizing::MeshSpec spec;
    // Keep each source curve as one segment so the dyadic midpoint is a
    // genuinely new child rather than an already-present natural node.
    spec.field=std::make_shared<sizing::ConstantField>(10.0);
    spec.threads=1;
    algo::BoundaryConformingV3Request request;
    request.expectedSurfaceSetStamp={16791,1};
    request.expectedOwnershipProviderStamp=ownership.current;

    const geo::GEdge* selected=nullptr;
    std::size_t incidentFaces=0;
    for(const auto* edge:model.edges()) {
        if(!edge)continue;
        std::size_t count=0;
        for(const auto* face:model.faces()) {
            if(!face)continue;
            const auto tags=face->boundaryEdgeTags();
            count+=static_cast<std::size_t>(
                std::find(tags.begin(),tags.end(),edge->tag())!=tags.end());
        }
        if(count>incidentFaces) { selected=edge; incidentFaces=count; }
    }
    CHECK(selected!=nullptr); CHECK(incidentFaces>1);
    if(!selected||incidentFaces<=1)return;
    std::printf("  selected shared curve %d on faces",selected->tag());
    for(const auto* face:model.faces()) {
        if(!face)continue;
        const auto tags=face->boundaryEdgeTags();
        if(std::find(tags.begin(),tags.end(),selected->tag())!=tags.end())
            std::printf(" %d",face->tag()+1);
    }
    std::printf("\n");

    algo::t2c_detail::SurfaceFaceCacheSession cache;
    std::vector<algo::SurfaceCacheFaceReportRow> processed(
        model.faces().size());
    std::vector<algo::FacePersistentId> unattempted(model.faces().size());
    algo::t2c_detail::SurfaceCacheDiagnosticChannel channel{
        processed,unattempted};
    cache.attachDiagnosticChannel(&channel);
    algo::t2c_detail::PrivateBuildTrace firstTrace;
    firstTrace.surfaceCacheSession=&cache;
    firstTrace.surfaceCacheChannel=&channel;
    mesh::Mesh firstDestination;
    const auto first=algo::t2c_detail::
        meshModelBoundaryConformingWithProtectedSplits(
            model,spec,request,firstDestination,ownership,{}, {},&firstTrace);
    CHECK(first.published());
    const auto firstBuilt=cache.metrics().facesBuilt;
    const auto firstReused=cache.metrics().facesReused;

    algo::SurfaceNodeIdentity identity;
    identity.domain=algo::SurfaceNodeIdentityDomain::RefinedCurve;
    identity.parentArtifact={
        static_cast<std::uint64_t>(selected->tag())+1U,1};
    identity.parentOrdinal=0;
    identity.parameterUBits=std::bit_cast<std::uint64_t>(0.5);
    identity.dyadicNumerator=1;
    identity.dyadicLevel=1;
    const algo::t2c_detail::DyadicEdgeSplit split{
        selected->tag(),0,0.5,selected->pointAt(0.5),identity};
    const std::span<const algo::t2c_detail::DyadicEdgeSplit> splitSpan(
        &split,1);
    algo::t2c_detail::PrivateBuildTrace secondTrace;
    secondTrace.surfaceCacheSession=&cache;
    secondTrace.surfaceCacheChannel=&channel;
    secondTrace.proposedEdgeSplits=splitSpan;
    mesh::Mesh secondDestination;
    const auto second=algo::t2c_detail::
        meshModelBoundaryConformingWithProtectedSplits(
            model,spec,request,secondDestination,ownership,{}, {},
            &secondTrace);
    std::printf("  shared-curve edge=%d incident=%zu second status=%d reason=%d phase=%d "
                "cache=(%d,%d,%d face=%llu) built_delta=%llu reused_delta=%llu "
                "report=%zu+%zu\n",
        selected->tag(),incidentFaces,
        static_cast<int>(second.status),static_cast<int>(second.failure.reason),
        static_cast<int>(second.failure.phase),
        static_cast<int>(secondTrace.surfaceCacheDiagnostic.status),
        static_cast<int>(secondTrace.surfaceCacheDiagnostic.reason),
        static_cast<int>(secondTrace.surfaceCacheDiagnostic.phase),
        static_cast<unsigned long long>(
            secondTrace.surfaceCacheDiagnostic.failedFace),
        static_cast<unsigned long long>(
            cache.metrics().facesBuilt-firstBuilt),
        static_cast<unsigned long long>(
            cache.metrics().facesReused-firstReused),
        channel.processedUsed,channel.unattemptedUsed);
    CHECK(second.published());
    CHECK(channel.reportComplete);
    CHECK(channel.processedUsed==model.faces().size());
    CHECK(channel.unattemptedUsed==0);
    CHECK(cache.metrics().facesBuilt-firstBuilt==incidentFaces);
    CHECK(cache.metrics().facesReused-firstReused==
          model.faces().size()-incidentFaces);
    std::size_t dirtyRows=0,cleanRows=0;
    for(std::size_t i=0;i<channel.processedUsed;++i) {
        dirtyRows+=static_cast<std::size_t>(
            channel.processedStorage[i].disposition==
                algo::SurfaceCacheFaceDisposition::DirtyBuild);
        cleanRows+=static_cast<std::size_t>(
            channel.processedStorage[i].disposition==
                algo::SurfaceCacheFaceDisposition::CleanReplay);
    }
    CHECK(dirtyRows==incidentFaces);
    CHECK(cleanRows==model.faces().size()-incidentFaces);
    CHECK(cache.metrics().allFaceScansPerDirtyFace==0);
}

void checkSinglePassMeasuredWork(
        const algo::BoundaryConformingV4Result& result,
        std::uint64_t referenceTotal,std::uint64_t msvcTotal,
        std::uint64_t nonComparatorTotal) {
    const auto& metrics=result.metrics;
    CHECK(result.sourceOutcome&&result.sourceOutcome->retainedResult);
    if(!result.sourceOutcome||!result.sourceOutcome->retainedResult)return;
    const auto& construction=
        result.sourceOutcome->retainedResult->constructionMetrics;
    const auto& surface=construction.surface.scalable;
    const auto& candidate=construction.tetCandidate;
    // The original-coordinate positivity audit visits all 959 returned
    // tetrahedra. This fixture needs no discarded construction or retry.
    const auto orientationAuditWork=
        candidate.construction.originalOrientationAuditTests;
    CHECK(orientationAuditWork==UINT64_C(959));
    // The twelve enclosing-box triangles support all eight box corners.
    // The hull audit can certify the input after those 96 orientation tests.
    CHECK(candidate.uniquePoints==UINT64_C(189));
    CHECK(candidate.boundaryFaces==UINT64_C(12));
    CHECK(candidate.validationHullTests==UINT64_C(12)*UINT64_C(8));
    CHECK(candidate.discardedCandidateWork==0);
    CHECK(candidate.construction.discardedConstructionWork==0);
    CHECK(!candidate.originalCoordinateRetryUsed);
    CHECK(!candidate.construction.originalCoordinateRetryUsed);
    // These five counters charge comparator invocations made by std::sort
    // and std::set. Their invocation counts depend on the standard library;
    // the remaining work in this fixed dyadic fixture is invariant.
    const auto surfaceComparisons=surface.coordinateSortComparisons+
                                 surface.insertionOrderComparisons;
    const auto importComparisons=construction.importSelectionSortComparisons+
        construction.importGlobalMapComparisons+
        construction.partialImportSortComparisons;
    const auto comparisons=surfaceComparisons+importComparisons;
    CHECK(construction.surfaceWorkUnits>=surfaceComparisons);
    CHECK(construction.importWorkUnits>=importComparisons);
    CHECK(construction.surfaceWorkUnits-surfaceComparisons==UINT64_C(16051));
    CHECK(construction.tetCandidateWorkUnits==
          UINT64_C(62502)+orientationAuditWork);
    CHECK(construction.importWorkUnits-importComparisons==UINT64_C(56545));
    CHECK(construction.extractionWorkUnits==UINT64_C(18634));
    CHECK(construction.totalWorkUnits==construction.surfaceWorkUnits+
        construction.tetCandidateWorkUnits+construction.importWorkUnits+
        construction.extractionWorkUnits);
    CHECK(construction.runComplete);
    CHECK(metrics.nestedRecoveryWorkUnits==UINT64_C(11749));
    CHECK(metrics.coordinatorWorkUnits==UINT64_C(208));
    const auto cacheWork=metrics.surfaceCachePreparationWorkUnits+
        metrics.surfaceCacheReplayWorkUnits+metrics.surfaceCacheSetAuditWorkUnits;
    CHECK(metrics.nestedConstructionWorkUnits==
          construction.totalWorkUnits+cacheWork);
    CHECK(metrics.totalMeasuredWorkUnits==metrics.nestedConstructionWorkUnits+
        metrics.nestedRecoveryWorkUnits+metrics.coordinatorWorkUnits+cacheWork);
    CHECK(metrics.totalMeasuredWorkUnits>=comparisons);
    CHECK(metrics.totalMeasuredWorkUnits-comparisons==
          nonComparatorTotal+orientationAuditWork);
    std::printf("  measured comparator work: %llu; invariant other work: %llu\n",
        static_cast<unsigned long long>(comparisons),
        static_cast<unsigned long long>(metrics.totalMeasuredWorkUnits-comparisons));

    // Keep exact totals for qualified release environments only. Other
    // compiler/STL versions still enforce every exact invariant above.
#if defined(__linux__) && defined(__x86_64__) && defined(__GNUC__) && \
    !defined(__clang__) && __GNUC__==16 && __GNUC_MINOR__==2 && \
    __GNUC_PATCHLEVEL__==0 && defined(_GLIBCXX_RELEASE) && \
    _GLIBCXX_RELEASE==16 && !defined(_GLIBCXX_DEBUG) && \
    !defined(_GLIBCXX_ASSERTIONS) && !defined(_GLIBCXX_PARALLEL) && \
    !defined(_DEBUG)
    CHECK(metrics.totalMeasuredWorkUnits==referenceTotal+orientationAuditWork);
    (void)msvcTotal;
#elif defined(_WIN64) && defined(_M_X64) && !defined(__clang__) && \
    defined(_MSC_FULL_VER) && _MSC_FULL_VER==194435207 && \
    defined(_MSC_BUILD) && _MSC_BUILD==1 && defined(_MSVC_STL_VERSION) && \
    _MSVC_STL_VERSION==143 && defined(_MSVC_STL_UPDATE) && \
    _MSVC_STL_UPDATE==202503L && defined(_ITERATOR_DEBUG_LEVEL) && \
    _ITERATOR_DEBUG_LEVEL==0 && !defined(_DEBUG) && \
    defined(_M_FP_STRICT) && _M_FP_STRICT==1
    CHECK(metrics.totalMeasuredWorkUnits==msvcTotal+orientationAuditWork);
    (void)referenceTotal;
#else
    (void)referenceTotal;
    (void)msvcTotal;
    std::puts("  no qualified compiler/STL total golden; exact non-comparator checks active");
#endif
}

algo::BoundaryConformingV4Request makeSinglePassPrefillRequest(
        std::size_t regions,const Ownership& ownership,
        const CurvePlacement& placement,std::uint32_t grid=3,
        algo::ProtectedFeatureExecutionStrategy strategy=algo::
            ProtectedFeatureExecutionStrategy::SinglePassInterfacePrefill) {
    algo::BoundaryConformingV4Request request;
    request.sourceRequest.expectedSurfaceSetStamp={17000,1};
    request.sourceRequest.expectedOwnershipProviderStamp=ownership.current;
    request.expectedPlacementProviderStamp=placement.current;
    const auto splitLimit=static_cast<std::uint64_t>(regions*16U+16U);
    request.limits={splitLimit,1,splitLimit,
                    static_cast<std::uint64_t>(regions*5000U+5000U),
                    8ULL*1024ULL*1024ULL*1024ULL,UINT64_C(1)<<50,
                    128,4096,4096};
    request.executionStrategy=strategy;
    request.interfacePrefill={grid,grid,
        static_cast<std::uint64_t>(regions),splitLimit,3.5};
    request.facetConflictMode=algo::FacetConflictRefinementMode::Refuse;
    request.recoverySchedulingMode=algo::TetRecoverySchedulingMode::
        MixedExactPresentBatching;
    return request;
}

void singlePassInterfacePrefillRoute() {
    TEST("t2c single-pass interface prefill publishes one audited explicit mesh");
    constexpr std::size_t regions=8;
    LayeredModel model(static_cast<int>(regions),false,4.0,1.0);
    Ownership ownership;
    for(const auto* face:model.faces())
        ownership.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    CurvePlacement placement;
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,4.0);
    spec.threads=1;
    mesh::Mesh destination;
    const auto outcome=algo::meshModelBoundaryConformingV4(
        model,spec,makeSinglePassPrefillRequest(
            regions,ownership,placement),destination,ownership,placement);
    CHECK(outcome.published());
    CHECK(outcome.status==algo::BoundaryConformingV4Status::Published);
    CHECK(outcome.reason==algo::BoundaryConformingV4Reason::None);
    CHECK(outcome.phase==algo::BoundaryConformingV4Phase::Publication);
    CHECK(!outcome.destinationUnchanged);
    CHECK(outcome.retainedResult!=nullptr);
    if(!outcome.retainedResult)return;
    const auto& result=*outcome.retainedResult;
    CHECK(result.executionStrategy==algo::
          ProtectedFeatureExecutionStrategy::SinglePassInterfacePrefill);
    CHECK(result.metrics.rounds==1);
    CHECK(result.metrics.interfacePrefillCandidateFaces==7);
    CHECK(result.metrics.interfacePrefillFaces==6);
    CHECK(result.metrics.interfacePrefillDenseFacesSkipped==1);
    CHECK(result.metrics.interfacePrefillPoints==54);
    CHECK(result.metrics.interfacePrefillSkippedParameters==0);
    CHECK(result.metrics.interfacePrefillProviderCalls==54);
    CHECK(result.metrics.facetSplits==54);
    CHECK(result.metrics.edgeSplits==0);
    CHECK(result.metrics.uniquePlacements==54);
    CHECK(result.metrics.frontierRows==0);
    CHECK(result.metrics.latestSourceMixedPresentScheduling.
          nontrivialRecoveries==0);
    std::printf("  single-pass measured work: %llu\n",
        static_cast<unsigned long long>(result.metrics.totalMeasuredWorkUnits));
    checkSinglePassMeasuredWork(result,UINT64_C(1420812),
        UINT64_C(1415944),UINT64_C(1389519));
    CHECK(result.metrics.preparedSurfaceRecoveryFaces==model.faces().size());
    CHECK(result.metrics.directSurfaceBuildAttempts==0);
    CHECK(result.metrics.directSurfaceBuildFaces==0);
    CHECK(result.metrics.runComplete);
    CHECK(result.subdivisionLedger.size()==54);
    CHECK(placement.placementCalls==54);
    CHECK(placement.facetPlacementCalls==54);
    CHECK(placement.requests.size()==54);
    for(const auto& row:placement.requests) {
        CHECK(row.kind==algo::ProtectedFeatureEntityKind::FacetInterior);
        CHECK(row.requestedIdentity.domain==
              algo::SurfaceNodeIdentityDomain::InterfacePrefill);
    }
    CHECK(destination.nodeCount()==181);
    CHECK(destination.countCells(mesh::CellType::Tet4)==657);
    CHECK(mesh::meshHash(destination)==UINT64_C(1950448103972936315));
    const auto semantic=mesh::semanticHashEx(destination);
    CHECK(semantic.canonical);
    CHECK(semantic.value==UINT64_C(18033145305553960814));
    CHECK(result.sourceOutcome!=nullptr);
    if(result.sourceOutcome&&result.sourceOutcome->retainedResult) {
        const auto& source=*result.sourceOutcome->retainedResult;
        CHECK(source.boundaryAudit.passed);
        CHECK(source.boundaryAudit.orphanExterior==0);
        CHECK(source.boundaryAudit.unexpectedExterior==0);
        CHECK(source.boundaryAudit.missingInterface==0);
        CHECK(source.boundaryAudit.unexpectedInterface==0);
        CHECK(source.recoveryMetrics.edgeCavityCommits==0);
        CHECK(source.recoveryMetrics.facetCommits==0);
        CHECK(source.recoveryMetrics.productionGlobalScanFallbacks==0);
        CHECK(source.recoveryMetrics.fullStageRescans==0);
    }
    const auto quality=mesh::analyzeTets(destination);
    CHECK(quality.nTets==657);
    CHECK(quality.nInverted==0);
    CHECK(quality.nNonFinite==0);
    CHECK(quality.radiusRatio.min>0.19);
    CHECK(quality.scaledJacobian.min>0.17);
    CHECK(quality.minDihedralDeg>11.7);
    CHECK(quality.maxDihedralDeg<159.0);
    CHECK(std::abs(quality.totalVolume-512.0)<1.0e-9);
}

void singlePassInterfacePrefillDirectRoute() {
    TEST("t2c direct single-pass bypasses cache with exact v2 output");
    constexpr std::size_t regions=8;
    LayeredModel model(static_cast<int>(regions),false,4.0,1.0);
    Ownership ownership;
    for(const auto* face:model.faces())
        ownership.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    CurvePlacement placement;
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,4.0);
    spec.threads=1;
    auto request=makeSinglePassPrefillRequest(
        regions,ownership,placement,3,
        algo::ProtectedFeatureExecutionStrategy::
            SinglePassInterfacePrefillDirect);
    mesh::Mesh destination;
    const auto outcome=algo::meshModelBoundaryConformingV4(
        model,spec,request,destination,ownership,placement);
    CHECK(outcome.published());
    CHECK(outcome.retainedResult!=nullptr);
    if(!outcome.retainedResult)return;
    const auto& result=*outcome.retainedResult;
    CHECK(result.executionStrategy==algo::
          ProtectedFeatureExecutionStrategy::
              SinglePassInterfacePrefillDirect);
    CHECK(result.metrics.rounds==1);
    CHECK(result.metrics.directSurfaceBuildAttempts==1);
    CHECK(result.metrics.directSurfaceBuildFaces==model.faces().size());
    CHECK(result.metrics.preparedSurfaceRecoveryFaces==model.faces().size());
    CHECK(result.metrics.surfaceCacheFacesBuilt==0);
    CHECK(result.metrics.surfaceCacheFacesReused==0);
    CHECK(result.metrics.surfaceCacheFacesReplayed==0);
    CHECK(result.metrics.surfaceCachePreparationWorkUnits==0);
    CHECK(result.metrics.surfaceCacheReplayWorkUnits==0);
    CHECK(result.metrics.surfaceCacheSetAuditWorkUnits==0);
    CHECK(result.metrics.surfaceCacheAccountedBytes==0);
    CHECK(result.surfaceCacheReport.reportComplete);
    CHECK(!result.surfaceCacheReport.hasFailedFace);
    CHECK(result.surfaceCacheReport.processed.size()==model.faces().size());
    CHECK(result.surfaceCacheReport.unattempted.empty());
    for(const auto& row:result.surfaceCacheReport.processed) {
        CHECK(row.disposition==algo::SurfaceCacheFaceDisposition::DirtyBuild);
        CHECK(row.surfaceStatus==algo::SurfaceFaceStatus::Ready);
    }
    CHECK(destination.nodeCount()==181);
    CHECK(destination.countCells(mesh::CellType::Tet4)==657);
    CHECK(mesh::meshHash(destination)==UINT64_C(1950448103972936315));
    const auto semantic=mesh::semanticHashEx(destination);
    CHECK(semantic.canonical);
    CHECK(semantic.value==UINT64_C(18033145305553960814));
    const auto quality=mesh::analyzeTets(destination);
    CHECK(quality.nTets==657);
    CHECK(quality.nInverted==0);
    CHECK(quality.radiusRatio.min>0.19);
    CHECK(quality.scaledJacobian.min>0.17);
}

void singlePassPreparedConstraintPrefillRoute() {
    TEST("t2c v4 prepared constraints preserve semantic mesh at lower work");
    constexpr std::size_t regions=8;
    LayeredModel model(static_cast<int>(regions),false,4.0,1.0);
    Ownership ownership;
    for(const auto* face:model.faces())
        ownership.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    CurvePlacement placement;
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,4.0);
    spec.threads=1;
    const auto request=makeSinglePassPrefillRequest(
        regions,ownership,placement,3,
        algo::ProtectedFeatureExecutionStrategy::
            SinglePassPreparedConstraintPrefill);
    mesh::Mesh destination;
    const auto outcome=algo::meshModelBoundaryConformingV4(
        model,spec,request,destination,ownership,placement);
    CHECK(outcome.published());
    CHECK(outcome.retainedResult!=nullptr);
    if(!outcome.retainedResult)return;
    const auto& result=*outcome.retainedResult;
    CHECK(result.executionStrategy==algo::
          ProtectedFeatureExecutionStrategy::
              SinglePassPreparedConstraintPrefill);
    CHECK(result.metrics.rounds==1);
    CHECK(result.metrics.interfacePrefillFaces==6);
    CHECK(result.metrics.interfacePrefillPoints==54);
    CHECK(result.metrics.directSurfaceBuildAttempts==1);
    CHECK(result.metrics.directSurfaceBuildFaces==model.faces().size());
    CHECK(result.metrics.preparedSurfaceRecoveryFaces==model.faces().size());
    CHECK(result.metrics.surfaceCacheFacesBuilt==0);
    CHECK(result.metrics.surfaceCacheFacesReused==0);
    CHECK(result.metrics.surfaceCachePreparationWorkUnits==0);
    CHECK(result.metrics.surfaceCacheAccountedBytes==0);
    std::printf("  prepared-constraint measured work: %llu\n",
        static_cast<unsigned long long>(result.metrics.totalMeasuredWorkUnits));
    checkSinglePassMeasuredWork(result,UINT64_C(196982),
        UINT64_C(192114),UINT64_C(165689));
    CHECK(result.metrics.runComplete);
    CHECK(destination.nodeCount()==181);
    CHECK(destination.countCells(mesh::CellType::Tet4)==657);
    CHECK(mesh::meshHash(destination)==UINT64_C(2499501827151892571));
    const auto semantic=mesh::semanticHashEx(destination);
    CHECK(semantic.canonical);
    CHECK(semantic.value==UINT64_C(18033145305553960814));
    CHECK(result.sourceOutcome!=nullptr);
    if(result.sourceOutcome&&result.sourceOutcome->retainedResult) {
        const auto& source=*result.sourceOutcome->retainedResult;
        CHECK(source.boundaryAudit.passed);
        CHECK(source.boundaryAudit.orphanExterior==0);
        CHECK(source.boundaryAudit.unexpectedExterior==0);
        CHECK(source.boundaryAudit.missingInterface==0);
        CHECK(source.boundaryAudit.unexpectedInterface==0);
        CHECK(source.recoveryMetrics.edgeCavityCommits==0);
        CHECK(source.recoveryMetrics.facetCommits==0);
        CHECK(source.recoveryMetrics.productionGlobalScanFallbacks==0);
        CHECK(source.recoveryMetrics.fullStageRescans==0);
    }
    const auto quality=mesh::analyzeTets(destination);
    CHECK(quality.nTets==657);
    CHECK(quality.nInverted==0);
    CHECK(quality.nNonFinite==0);
    CHECK(quality.radiusRatio.min>0.19);
    CHECK(quality.scaledJacobian.min>0.17);
    CHECK(quality.minDihedralDeg>11.7);
    CHECK(quality.maxDihedralDeg<159.0);
    CHECK(std::abs(quality.totalVolume-512.0)<1.0e-9);
}

void sizingDrivenSinglePassPrefillRoute() {
    TEST("t2c v6 compiled controls drive bounded one-pass interface density");
    constexpr std::size_t regions=8;
    LayeredModel model(static_cast<int>(regions),false,4.0,1.0);
    Ownership ownership;
    for(const auto* face:model.faces())
        ownership.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    CurvePlacement placement;
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,4.0);
    spec.threads=1;
    sizing::MeshControlSpecification controlSpec;
    controlSpec.minimumSize=0.05;
    controlSpec.maximumSize=4.0;
    controlSpec.maximumGrowthRate=1.3;
    controlSpec.entitySizes={{3,0,4.0}};
    const auto compiled=sizing::compileMeshControls(controlSpec);
    CHECK(compiled.compiled);
    auto request=makeSinglePassPrefillRequest(
        regions,ownership,placement,3,
        algo::ProtectedFeatureExecutionStrategy::
            SinglePassPreparedConstraintPrefill);
    request.interfacePrefill.densityMode=
        algo::InterfacePrefillDensityMode::SizingDriven;
    request.interfacePrefill.minimumInteriorU=3;
    request.interfacePrefill.minimumInteriorV=3;
    request.interfacePrefill.maximumInteriorU=3;
    request.interfacePrefill.maximumInteriorV=3;
    request.interfacePrefill.sizingScale=1.0;
    request.meshControls=std::make_shared<const sizing::CompiledMeshControls>(
        compiled.controls);
    mesh::Mesh destination;
    const auto outcome=algo::meshModelBoundaryConformingV4(
        model,spec,request,destination,ownership,placement);
    CHECK(outcome.published());
    CHECK(outcome.retainedResult!=nullptr);
    if(!outcome.retainedResult)return;
    const auto& result=*outcome.retainedResult;
    CHECK(result.metrics.interfacePrefillSizingDrivenFaces==6);
    CHECK(result.metrics.interfacePrefillControlForcedFaces==0);
    CHECK(result.metrics.interfacePrefillMaximumPlannedU==3);
    CHECK(result.metrics.interfacePrefillMaximumPlannedV==3);
    CHECK(result.metrics.interfacePrefillPoints==54);
    CHECK(result.metrics.rounds==1);
    CHECK(result.metrics.runComplete);
    CHECK(destination.nodeCount()==181);
    CHECK(destination.countCells(mesh::CellType::Tet4)==657);
    const auto semantic=mesh::semanticHashEx(destination);
    CHECK(semantic.canonical);
    CHECK(semantic.value==UINT64_C(18033145305553960814));
}


void sizingDrivenFinerBodyControlRoute() {
    TEST("t2c v6 finer body control changes one-pass interface density");
    constexpr std::size_t regions=8;
    LayeredModel model(static_cast<int>(regions),false,4.0,1.0);
    Ownership ownership;
    for(const auto* face:model.faces())
        ownership.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    CurvePlacement placement;
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,4.0);
    spec.threads=1;
    sizing::MeshControlSpecification controlSpec;
    controlSpec.minimumSize=0.05;
    controlSpec.maximumSize=4.0;
    controlSpec.maximumGrowthRate=1.3;
    controlSpec.entitySizes={{3,0,1.5}};
    const auto compiled=sizing::compileMeshControls(controlSpec);
    CHECK(compiled.compiled);
    auto request=makeSinglePassPrefillRequest(
        regions,ownership,placement,3,
        algo::ProtectedFeatureExecutionStrategy::
            SinglePassPreparedConstraintPrefill);
    request.interfacePrefill.densityMode=
        algo::InterfacePrefillDensityMode::SizingDriven;
    request.interfacePrefill.minimumInteriorU=3;
    request.interfacePrefill.minimumInteriorV=3;
    request.interfacePrefill.maximumInteriorU=6;
    request.interfacePrefill.maximumInteriorV=6;
    request.interfacePrefill.sizingScale=1.0;
    request.meshControls=std::make_shared<const sizing::CompiledMeshControls>(
        compiled.controls);
    mesh::Mesh destination;
    const auto outcome=algo::meshModelBoundaryConformingV4(
        model,spec,request,destination,ownership,placement);
    CHECK(outcome.published());
    CHECK(outcome.retainedResult!=nullptr);
    if(!outcome.retainedResult)return;
    const auto& result=*outcome.retainedResult;
    CHECK(result.metrics.rounds==1);
    CHECK(result.metrics.interfacePrefillCandidateFaces==7);
    CHECK(result.metrics.interfacePrefillFaces==7);
    CHECK(result.metrics.interfacePrefillDenseFacesSkipped==0);
    CHECK(result.metrics.interfacePrefillSizingDrivenFaces==7);
    CHECK(result.metrics.interfacePrefillControlForcedFaces==1);
    CHECK(result.metrics.interfacePrefillMaximumPlannedU==5);
    CHECK(result.metrics.interfacePrefillMaximumPlannedV==5);
    CHECK(result.metrics.interfacePrefillPoints==79);
    CHECK(result.metrics.facetSplits==79);
    CHECK(result.metrics.runComplete);
    CHECK(destination.countCells(mesh::CellType::Tet4)>0);
    const auto quality=mesh::analyzeTets(destination);
    CHECK(quality.nTets==destination.countCells(mesh::CellType::Tet4));
    CHECK(quality.nInverted==0);
    CHECK(quality.nNonFinite==0);
    CHECK(quality.totalVolume>0.0);
    const auto semantic=mesh::semanticHashEx(destination);
    CHECK(semantic.canonical);
}

void singlePassPreparedConstraintDeterminism() {
    TEST("t2c v4 prepared constraints are invariant to threads and presentation");
    struct Receipt {
        std::uint64_t strict=0,semantic=0,work=0;
        Index nodes=0,tets=0;
    };
    const auto run=[](std::size_t threads,bool reverse) {
        constexpr std::size_t regions=8;
        LayeredModel model(static_cast<int>(regions),false,4.0,1.0,reverse);
        Ownership ownership;
        for(const auto* face:model.faces())
            ownership.rows[
                static_cast<algo::FacePersistentId>(face->tag())+1]=
                    model.owners(face->tag());
        CurvePlacement placement;
        sizing::MeshSpec spec;
        spec.field=std::make_shared<LayerQualificationField>(1.0,4.0);
        spec.threads=threads;
        mesh::Mesh destination;
        const auto outcome=algo::meshModelBoundaryConformingV4(
            model,spec,makeSinglePassPrefillRequest(
                regions,ownership,placement,3,
                algo::ProtectedFeatureExecutionStrategy::
                    SinglePassPreparedConstraintPrefill),
            destination,ownership,placement);
        Receipt receipt;
        if(!outcome.published()||!outcome.retainedResult)return receipt;
        receipt.strict=mesh::meshHash(destination);
        const auto semantic=mesh::semanticHashEx(destination);
        receipt.semantic=semantic.canonical?semantic.value:0;
        receipt.work=outcome.retainedResult->metrics.totalMeasuredWorkUnits;
        receipt.nodes=destination.nodeCount();
        receipt.tets=destination.countCells(mesh::CellType::Tet4);
        return receipt;
    };
    const auto one=run(1,false);
    const auto four=run(4,false);
    const auto automatic=run(0,true);
    CHECK(one.strict==UINT64_C(2499501827151892571));
    CHECK(one.strict==four.strict);
    CHECK(one.strict==automatic.strict);
    CHECK(one.semantic==four.semantic);
    CHECK(one.semantic==automatic.semantic);
    CHECK(one.work==four.work);
    CHECK(one.work==automatic.work);
    CHECK(one.nodes==four.nodes&&one.nodes==automatic.nodes);
    CHECK(one.tets==four.tets&&one.tets==automatic.tets);
}

void singlePassInterfacePrefillDeterminism() {
    TEST("t2c single-pass prefill is invariant to threads and presentation");
    struct Receipt {
        std::uint64_t strict=0,semantic=0,work=0;
        Index nodes=0,tets=0;
        std::vector<algo::SurfaceNodeIdentity> identities;
    };
    const auto run=[](std::size_t threads,bool reverse) {
        constexpr std::size_t regions=8;
        LayeredModel model(static_cast<int>(regions),false,4.0,1.0,reverse);
        Ownership ownership;
        for(const auto* face:model.faces())
            ownership.rows[
                static_cast<algo::FacePersistentId>(face->tag())+1]=
                    model.owners(face->tag());
        CurvePlacement placement;
        sizing::MeshSpec spec;
        spec.field=std::make_shared<LayerQualificationField>(1.0,4.0);
        spec.threads=threads;
        mesh::Mesh destination;
        const auto outcome=algo::meshModelBoundaryConformingV4(
            model,spec,makeSinglePassPrefillRequest(
                regions,ownership,placement),destination,
            ownership,placement);
        Receipt receipt;
        if(!outcome.published()||!outcome.retainedResult)return receipt;
        receipt.strict=mesh::meshHash(destination);
        const auto semantic=mesh::semanticHashEx(destination);
        receipt.semantic=semantic.canonical?semantic.value:0;
        receipt.work=outcome.retainedResult->metrics.totalMeasuredWorkUnits;
        receipt.nodes=destination.nodeCount();
        receipt.tets=destination.countCells(mesh::CellType::Tet4);
        for(const auto& row:outcome.retainedResult->subdivisionLedger)
            receipt.identities.push_back(row.childIdentity);
        return receipt;
    };
    const auto one=run(1,false);
    const auto four=run(4,false);
    const auto automatic=run(0,true);
    CHECK(one.strict!=0);
    CHECK(one.strict==four.strict);
    CHECK(one.strict==automatic.strict);
    CHECK(one.semantic==four.semantic);
    CHECK(one.semantic==automatic.semantic);
    CHECK(one.work==four.work);
    CHECK(one.work==automatic.work);
    CHECK(one.nodes==four.nodes);
    CHECK(one.nodes==automatic.nodes);
    CHECK(one.tets==four.tets);
    CHECK(one.tets==automatic.tets);
    CHECK(one.identities==four.identities);
    CHECK(one.identities==automatic.identities);
}

void singlePassInterfacePrefillRefusesIncompleteMesh() {
    TEST("t2c single-pass underfill returns incomplete without publication");
    constexpr std::size_t regions=8;
    LayeredModel model(static_cast<int>(regions),false,4.0,1.0);
    Ownership ownership;
    for(const auto* face:model.faces())
        ownership.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    CurvePlacement placement;
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,4.0);
    spec.threads=1;
    mesh::Mesh destination;
    destination.reserveAppend(3,1,4);
    (void)destination.addNode({-17,-19,-23},{0,73});
    const auto beforeHash=mesh::meshHash(destination);
    const auto beforeCapacity=destination.capacitySnapshot();
    const auto outcome=algo::meshModelBoundaryConformingV4(
        model,spec,makeSinglePassPrefillRequest(
            regions,ownership,placement,2),destination,ownership,placement);
    CHECK(!outcome.published());
    CHECK(outcome.status==algo::BoundaryConformingV4Status::Incomplete);
    CHECK(outcome.reason==algo::BoundaryConformingV4Reason::
          SinglePassRecoveryIncomplete);
    CHECK(outcome.phase==algo::BoundaryConformingV4Phase::ConstraintRecovery);
    CHECK(outcome.destinationUnchanged);
    CHECK(mesh::meshHash(destination)==beforeHash);
    CHECK(destination.nodeCount()==1);
    CHECK(destination.cellCount()==0);
    const auto afterCapacity=destination.capacitySnapshot();
    CHECK(afterCapacity.nodeCoordinateCapacity==
          beforeCapacity.nodeCoordinateCapacity);
    CHECK(afterCapacity.nodeClassificationCapacity==
          beforeCapacity.nodeClassificationCapacity);
    CHECK(afterCapacity.cellCapacity==beforeCapacity.cellCapacity);
    CHECK(afterCapacity.connectivityCapacity==
          beforeCapacity.connectivityCapacity);
    CHECK(outcome.retainedResult!=nullptr);
    if(!outcome.retainedResult)return;
    const auto& result=*outcome.retainedResult;
    CHECK(result.metrics.rounds==1);
    CHECK(result.metrics.interfacePrefillCandidateFaces==7);
    CHECK(result.metrics.interfacePrefillFaces==6);
    CHECK(result.metrics.interfacePrefillDenseFacesSkipped==1);
    CHECK(result.metrics.interfacePrefillPoints==24);
    CHECK(result.metrics.uniquePlacements==24);
    CHECK(result.metrics.runComplete==false);
    CHECK(result.sourceOutcome!=nullptr);
    CHECK(result.sourceOutcome&&
          result.sourceOutcome->status==
              algo::BoundaryConformingV3Status::Refused);

    // The direct no-cache route has the same safe one-pass refusal contract.
    // It must never publish a partial mesh merely because its reusable cache
    // is intentionally disabled.
    mesh::Mesh directDestination;
    directDestination.reserveAppend(3,1,4);
    (void)directDestination.addNode({-41,-43,-47},{0,75});
    const auto directBeforeHash=mesh::meshHash(directDestination);
    const auto directBeforeCapacity=directDestination.capacitySnapshot();
    const auto directOutcome=algo::meshModelBoundaryConformingV4(
        model,spec,makeSinglePassPrefillRequest(
            regions,ownership,placement,2,
            algo::ProtectedFeatureExecutionStrategy::
                SinglePassInterfacePrefillDirect),
        directDestination,ownership,placement);
    CHECK(!directOutcome.published());
    CHECK(directOutcome.status==algo::BoundaryConformingV4Status::Incomplete);
    CHECK(directOutcome.reason==algo::BoundaryConformingV4Reason::
          SinglePassRecoveryIncomplete);
    CHECK(directOutcome.phase==
          algo::BoundaryConformingV4Phase::ConstraintRecovery);
    CHECK(directOutcome.destinationUnchanged);
    CHECK(mesh::meshHash(directDestination)==directBeforeHash);
    CHECK(directDestination.nodeCount()==1);
    CHECK(directDestination.cellCount()==0);
    const auto directAfterCapacity=directDestination.capacitySnapshot();
    CHECK(directAfterCapacity.nodeCoordinateCapacity==
          directBeforeCapacity.nodeCoordinateCapacity);
    CHECK(directAfterCapacity.nodeClassificationCapacity==
          directBeforeCapacity.nodeClassificationCapacity);
    CHECK(directAfterCapacity.cellCapacity==directBeforeCapacity.cellCapacity);
    CHECK(directAfterCapacity.connectivityCapacity==
          directBeforeCapacity.connectivityCapacity);
    CHECK(directOutcome.retainedResult!=nullptr);
    if(!directOutcome.retainedResult)return;
    const auto& directResult=*directOutcome.retainedResult;
    CHECK(directResult.executionStrategy==algo::
          ProtectedFeatureExecutionStrategy::
              SinglePassInterfacePrefillDirect);
    CHECK(directResult.metrics.rounds==1);
    CHECK(directResult.metrics.interfacePrefillPoints==24);
    CHECK(directResult.metrics.uniquePlacements==24);
    CHECK(directResult.metrics.directSurfaceBuildAttempts==1);
    CHECK(directResult.metrics.directSurfaceBuildFaces==model.faces().size());
    CHECK(directResult.metrics.runComplete==false);
    CHECK(directResult.sourceOutcome!=nullptr);
    CHECK(directResult.sourceOutcome&&
          directResult.sourceOutcome->status==
              algo::BoundaryConformingV3Status::Refused);

    mesh::Mesh v4Destination;
    v4Destination.reserveAppend(3,1,4);
    (void)v4Destination.addNode({-53,-59,-61},{0,76});
    const auto v4BeforeHash=mesh::meshHash(v4Destination);
    const auto v4BeforeCapacity=v4Destination.capacitySnapshot();
    const auto v4Outcome=algo::meshModelBoundaryConformingV4(
        model,spec,makeSinglePassPrefillRequest(
            regions,ownership,placement,2,
            algo::ProtectedFeatureExecutionStrategy::
                SinglePassPreparedConstraintPrefill),
        v4Destination,ownership,placement);
    CHECK(!v4Outcome.published());
    CHECK(v4Outcome.status==algo::BoundaryConformingV4Status::Incomplete);
    CHECK(v4Outcome.reason==algo::BoundaryConformingV4Reason::
          SinglePassRecoveryIncomplete);
    CHECK(v4Outcome.phase==algo::BoundaryConformingV4Phase::
          ConstraintRecovery);
    CHECK(v4Outcome.destinationUnchanged);
    CHECK(mesh::meshHash(v4Destination)==v4BeforeHash);
    CHECK(v4Destination.nodeCount()==1);
    CHECK(v4Destination.cellCount()==0);
    const auto v4AfterCapacity=v4Destination.capacitySnapshot();
    CHECK(v4AfterCapacity.nodeCoordinateCapacity==
          v4BeforeCapacity.nodeCoordinateCapacity);
    CHECK(v4AfterCapacity.nodeClassificationCapacity==
          v4BeforeCapacity.nodeClassificationCapacity);
    CHECK(v4AfterCapacity.cellCapacity==v4BeforeCapacity.cellCapacity);
    CHECK(v4AfterCapacity.connectivityCapacity==
          v4BeforeCapacity.connectivityCapacity);
    CHECK(v4Outcome.retainedResult!=nullptr);
    if(!v4Outcome.retainedResult)return;
    const auto& v4Result=*v4Outcome.retainedResult;
    CHECK(v4Result.executionStrategy==algo::
          ProtectedFeatureExecutionStrategy::
              SinglePassPreparedConstraintPrefill);
    CHECK(v4Result.metrics.rounds==1);
    CHECK(v4Result.metrics.interfacePrefillPoints==24);
    CHECK(v4Result.metrics.directSurfaceBuildAttempts==1);
    CHECK(v4Result.metrics.directSurfaceBuildFaces==model.faces().size());
    CHECK(!v4Result.metrics.runComplete);
}

void singlePassInterfacePrefillFailureContracts() {
    TEST("t2c single-pass prefill pins authority contract and budget stops");
    constexpr std::size_t regions=8;
    LayeredModel model(static_cast<int>(regions),false,4.0,1.0);
    const auto loadOwnership=[&](Ownership& ownership) {
        for(const auto* face:model.faces())
            ownership.rows[
                static_cast<algo::FacePersistentId>(face->tag())+1]=
                    model.owners(face->tag());
    };
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,4.0);
    spec.threads=1;
    mesh::Mesh destination;
    (void)destination.addNode({-29,-31,-37},{0,74});
    const auto digest=mesh::meshHash(destination);
    {
        Ownership ownership;loadOwnership(ownership);
        ownership.mutateAtStampCall=2;
        CurvePlacement placement;
        const auto out=algo::meshModelBoundaryConformingV4(
            model,spec,makeSinglePassPrefillRequest(
                regions,ownership,placement),destination,
            ownership,placement);
        CHECK(out.status==algo::BoundaryConformingV4Status::Failed);
        CHECK(out.reason==algo::BoundaryConformingV4Reason::
              InterfacePrefillAuthorityChanged);
        CHECK(out.phase==algo::BoundaryConformingV4Phase::InterfacePrefill);
        CHECK(out.destinationUnchanged);
        CHECK(out.retainedResult!=nullptr);
        CHECK(out.retainedResult&&out.retainedResult->metrics.rounds==0);
        CHECK(placement.placementCalls==0);
        CHECK(mesh::meshHash(destination)==digest);
    }
    {
        Ownership ownership;loadOwnership(ownership);
        CurvePlacement placement;placement.wrongParameter=true;
        const auto out=algo::meshModelBoundaryConformingV4(
            model,spec,makeSinglePassPrefillRequest(
                regions,ownership,placement),destination,
            ownership,placement);
        CHECK(out.status==algo::BoundaryConformingV4Status::Failed);
        CHECK(out.reason==algo::BoundaryConformingV4Reason::
              PlacementContractViolation);
        CHECK(out.phase==algo::BoundaryConformingV4Phase::InterfacePrefill);
        CHECK(out.destinationUnchanged);
        CHECK(placement.placementCalls==1);
        CHECK(mesh::meshHash(destination)==digest);
    }
    {
        Ownership ownership;loadOwnership(ownership);
        CurvePlacement placement;
        auto request=makeSinglePassPrefillRequest(
            regions,ownership,placement);
        request.interfacePrefill.maxPoints=53;
        request.limits.maxAddedNodes=53;
        const auto out=algo::meshModelBoundaryConformingV4(
            model,spec,request,destination,ownership,placement);
        CHECK(out.status==algo::BoundaryConformingV4Status::Incomplete);
        CHECK(out.reason==algo::BoundaryConformingV4Reason::
              InterfacePrefillBudgetReached);
        CHECK(out.phase==algo::BoundaryConformingV4Phase::InterfacePrefill);
        CHECK(out.destinationUnchanged);
        CHECK(out.retainedResult!=nullptr);
        CHECK(out.retainedResult&&
              out.retainedResult->metrics.interfacePrefillPoints==53);
        CHECK(out.retainedResult&&out.retainedResult->metrics.rounds==0);
        CHECK(mesh::meshHash(destination)==digest);
    }
    {
        Ownership ownership;loadOwnership(ownership);
        CurvePlacement placement;
        auto request=makeSinglePassPrefillRequest(
            regions,ownership,placement);
        request.facetConflictMode=algo::FacetConflictRefinementMode::
            BoundedCanonicalSccSubdivision;
        request.facetConflictLimits={64,128,128,385,UINT64_C(1)<<34};
        const auto out=algo::meshModelBoundaryConformingV4(
            model,spec,request,destination,ownership,placement);
        CHECK(out.status==algo::BoundaryConformingV4Status::Rejected);
        CHECK(out.reason==algo::BoundaryConformingV4Reason::InvalidRequest);
        CHECK(out.phase==algo::BoundaryConformingV4Phase::Preflight);
        CHECK(out.destinationUnchanged);
        CHECK(mesh::meshHash(destination)==digest);
    }
}

void surfaceCacheBudgetTransportAndDigest() {
    TEST("t2c cache budget retains complete xor report and destination");
    constexpr std::size_t regionCount=2;
    constexpr double extent=4.0;
    LayeredModel model(static_cast<int>(regionCount),false,extent,1.0);
    Ownership ownership;
    for(const auto* face:model.faces())
        ownership.rows[
            static_cast<algo::FacePersistentId>(face->tag())+1U]=
                model.owners(face->tag());
    CurvePlacement placement;
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,extent);
    spec.threads=1;
    algo::BoundaryConformingV4Request request;
    request.sourceRequest.expectedSurfaceSetStamp={16800,1};
    request.sourceRequest.expectedOwnershipProviderStamp=ownership.current;
    request.expectedPlacementProviderStamp=placement.current;
    request.limits={128,129,128,1000000,UINT64_C(1)<<34,
                    UINT64_C(1)<<50,128,4096,4096};
    request.surfaceCacheLimits.maxDirtyFacesPerBatch=1;

    mesh::Mesh destination;
    destination.reserveAppend(5,2,7);
    (void)destination.addNode({-7,-11,-13},{0,91});
    const auto beforeHash=mesh::meshHash(destination);
    const auto beforeCapacity=destination.capacitySnapshot();
    const auto beforeNodes=destination.nodeCount();
    const auto beforeCells=destination.cellCount();

    const auto outcome=algo::meshModelBoundaryConformingV4(
        model,spec,request,destination,ownership,placement);
    CHECK(outcome.status==algo::BoundaryConformingV4Status::Incomplete);
    CHECK(outcome.reason==
          algo::BoundaryConformingV4Reason::SurfaceCacheBudgetReached);
    CHECK(outcome.phase==
          algo::BoundaryConformingV4Phase::SurfaceCachePreparation);
    CHECK(!outcome.published());
    CHECK(outcome.destinationUnchanged);
    CHECK(outcome.retainedResult!=nullptr);
    CHECK(outcome.lastSourceOutcome!=nullptr);
    CHECK(placement.placementCalls==0);
    CHECK(placement.stampCalls==1);
    CHECK(mesh::meshHash(destination)==beforeHash);
    CHECK(destination.nodeCount()==beforeNodes);
    CHECK(destination.cellCount()==beforeCells);
    const auto afterCapacity=destination.capacitySnapshot();
    CHECK(afterCapacity.nodeCoordinateCapacity==
          beforeCapacity.nodeCoordinateCapacity);
    CHECK(afterCapacity.nodeClassificationCapacity==
          beforeCapacity.nodeClassificationCapacity);
    CHECK(afterCapacity.cellCapacity==beforeCapacity.cellCapacity);
    CHECK(afterCapacity.connectivityCapacity==
          beforeCapacity.connectivityCapacity);
    CHECK(afterCapacity.dedupBucketCount==beforeCapacity.dedupBucketCount);
    CHECK(afterCapacity.dedupEntryCount==beforeCapacity.dedupEntryCount);
    CHECK(afterCapacity.dedupIndexCapacity==beforeCapacity.dedupIndexCapacity);
    CHECK(afterCapacity.complete==beforeCapacity.complete);
    if(!outcome.retainedResult)return;

    const auto& report=outcome.retainedResult->surfaceCacheReport;
    requireSurfaceCacheReportXor(model,report);
    CHECK(!report.hasFailedFace);
    CHECK(report.processed.size()==model.faces().size());
    CHECK(report.unattempted.empty());
    CHECK(report.unattempted.size()+report.processed.size()==
          model.faces().size());
    for(const auto& row:report.processed) {
        CHECK(row.disposition==algo::SurfaceCacheFaceDisposition::Stopped);
        CHECK(row.surfaceStatus==algo::SurfaceFaceStatus::Failed);
        CHECK(row.originalSurfaceFailure.reason==
              algo::SurfaceReason::InternalError);
    }
    const auto& metrics=outcome.retainedResult->metrics;
    CHECK(metrics.rounds==1);
    CHECK(metrics.surfaceCacheFacesBuilt==0);
    CHECK(metrics.surfaceCacheFacesReused==0);
    CHECK(metrics.surfaceCacheFacesReplayed==0);
    CHECK(metrics.surfaceCacheCommits==0);
    CHECK(metrics.surfaceCacheAllFaceScansPerDirtyFace==0);
    CHECK(metrics.providerPlacementCalls==0);
    CHECK(metrics.providerStampCalls==1);
}

void v4ProviderBudgetAndDigestContracts() {
    TEST("t2c v4 pins provider budget allocation and digest outcomes");
    constexpr std::size_t regionCount=8;
    constexpr double extent=4.0;
    LayeredModel model(static_cast<int>(regionCount),false,extent,1.0);
    auto makeOwnership=[&]() {
        Ownership provider;
        for(const auto* face:model.faces())
            provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
                model.owners(face->tag());
        return provider;
    };
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,extent);
    spec.threads=1;
    const auto makeRequest=[](const Ownership& ownership,
                              const CurvePlacement& placement) {
        algo::BoundaryConformingV4Request request;
        request.sourceRequest.expectedSurfaceSetStamp={16800,1};
        request.sourceRequest.expectedOwnershipProviderStamp=ownership.current;
        request.expectedPlacementProviderStamp=placement.current;
        request.limits={2,3,2,1000000,UINT64_C(1)<<34,
                        UINT64_C(1)<<50,1,4096,4096};
        return request;
    };
    mesh::Mesh destination;
    (void)destination.addNode({-7,-11,-13},{0,91});
    const auto digest=mesh::meshHash(destination);

    {
        auto ownership=makeOwnership();CurvePlacement placement;
        auto request=makeRequest(ownership,placement);request.limits={};
        const auto out=algo::meshModelBoundaryConformingV4(
            model,spec,request,destination,ownership,placement);
        CHECK(out.status==algo::BoundaryConformingV4Status::Rejected);
        CHECK(out.reason==algo::BoundaryConformingV4Reason::InvalidRequest);
        CHECK(placement.stampCalls==0);CHECK(placement.placementCalls==0);
        CHECK(mesh::meshHash(destination)==digest);
    }
    {
        auto ownership=makeOwnership();CurvePlacement placement;
        placement.status=algo::ProtectedFeaturePlacementStatus::Rejected;
        const auto out=algo::meshModelBoundaryConformingV4(
            model,spec,makeRequest(ownership,placement),destination,
            ownership,placement);
        CHECK(out.status==algo::BoundaryConformingV4Status::Rejected);
        CHECK(out.reason==algo::BoundaryConformingV4Reason::ProviderRejected);
        CHECK(placement.placementCalls==1);
        CHECK(mesh::meshHash(destination)==digest);
    }
    {
        auto ownership=makeOwnership();CurvePlacement placement;
        placement.allocatePlacement=true;
        const auto out=algo::meshModelBoundaryConformingV4(
            model,spec,makeRequest(ownership,placement),destination,
            ownership,placement);
        CHECK(out.status==algo::BoundaryConformingV4Status::Failed);
        CHECK(out.reason==algo::BoundaryConformingV4Reason::AllocationFailed);
        CHECK(mesh::meshHash(destination)==digest);
    }
    for(int corruption=0;corruption<3;++corruption) {
        auto ownership=makeOwnership();CurvePlacement placement;
        placement.wrongEntity=corruption==0;
        placement.wrongParameter=corruption==1;
        placement.offEntity=corruption==2;
        const auto out=algo::meshModelBoundaryConformingV4(
            model,spec,makeRequest(ownership,placement),destination,
            ownership,placement);
        CHECK(out.status==algo::BoundaryConformingV4Status::Failed);
        CHECK(out.reason==
              algo::BoundaryConformingV4Reason::PlacementContractViolation);
        CHECK(out.retainedResult!=nullptr);
        CHECK(placement.placementCalls==1);
        CHECK(mesh::meshHash(destination)==digest);
    }
    {
        auto ownership=makeOwnership();CurvePlacement placement;
        placement.throwPlacement=true;
        const auto out=algo::meshModelBoundaryConformingV4(
            model,spec,makeRequest(ownership,placement),destination,
            ownership,placement);
        CHECK(out.status==algo::BoundaryConformingV4Status::Failed);
        CHECK(out.reason==algo::BoundaryConformingV4Reason::ProviderFailed);
        CHECK(out.retainedResult!=nullptr);
        CHECK(placement.placementCalls==1);
        CHECK(mesh::meshHash(destination)==digest);
    }
    {
        auto ownership=makeOwnership();CurvePlacement placement;
        placement.mutateAtStampCall=1;
        const auto out=algo::meshModelBoundaryConformingV4(
            model,spec,makeRequest(ownership,placement),destination,
            ownership,placement);
        CHECK(out.status==algo::BoundaryConformingV4Status::Failed);
        CHECK(out.reason==
              algo::BoundaryConformingV4Reason::ProviderRevisionChanged);
        CHECK(placement.placementCalls==0);
        CHECK(mesh::meshHash(destination)==digest);
    }
    {
        auto ownership=makeOwnership();CurvePlacement placement;
        auto request=makeRequest(ownership,placement);
        request.limits.maxEdgeSplits=1;
        request.limits.maxAddedNodes=1;
        const auto out=algo::meshModelBoundaryConformingV4(
            model,spec,request,destination,ownership,placement);
        CHECK(out.status==algo::BoundaryConformingV4Status::Incomplete);
        CHECK(out.reason==
              algo::BoundaryConformingV4Reason::SubdivisionBudgetReached);
        CHECK(out.retainedResult!=nullptr);
        CHECK(out.retainedResult&&
              out.retainedResult->subdivisionLedger.size()==1);
        CHECK(mesh::meshHash(destination)==digest);
    }
}

void privateCadEdgeSubdivision() {
    TEST("t2c private CAD edge subdivision preserves shared identity");
    LayeredModel model(2);
    Ownership ownership;
    for(const auto* face:model.faces())
        ownership.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(10.0);
    spec.threads=1;
    const auto* edge=model.edges().front();
    CHECK(edge!=nullptr);
    if(!edge)return;
    algo::SurfaceNodeIdentity edgeIdentity;
    edgeIdentity.domain=algo::SurfaceNodeIdentityDomain::RefinedCurve;
    edgeIdentity.parentArtifact={
        static_cast<std::uint64_t>(edge->tag())+1U,1};
    edgeIdentity.parentOrdinal=0;
    edgeIdentity.parameterUBits=std::bit_cast<std::uint64_t>(0.5);
    edgeIdentity.dyadicNumerator=1;
    edgeIdentity.dyadicLevel=1;
    const algo::t2c_detail::DyadicEdgeSplit split{
        edge->tag(),0,0.5,edge->pointAt(0.5),edgeIdentity};
    mesh::Mesh refined;
    algo::t2c_detail::PrivateBuildTrace trace;
    algo::t2c_detail::meshEdgesWithProtectedSplits(
        model,spec,std::span<const algo::t2c_detail::DyadicEdgeSplit>(&split,1),
        refined,trace);
    std::vector<Index> midsides;
    for(Index node=0;node<refined.nodeCount();++node)
        if(refined.nodeClass(node).dim==1&&
           refined.nodeClass(node).tag==edge->tag()&&
           refined.node(node).x==split.coordinate.x&&
           refined.node(node).y==split.coordinate.y&&
           refined.node(node).z==split.coordinate.z)
            midsides.push_back(node);
    CHECK(midsides.size()==1);
    std::size_t childSegments=0;
    for(Index cell=0;cell<refined.cellCount();++cell)
        if(refined.cell(cell).type==mesh::CellType::Line2&&
           refined.cell(cell).onEntity.tag==edge->tag())++childSegments;
    CHECK(childSegments==2);
    CHECK(trace.segments.size()==
          static_cast<std::size_t>(refined.cellCount()));

    algo::BoundaryConformingV3Request request;
    request.expectedSurfaceSetStamp={16900,1};
    request.expectedOwnershipProviderStamp=ownership.current;
    mesh::Mesh destination;
    auto outcome=algo::t2c_detail::
        meshModelBoundaryConformingWithProtectedSplits(
            model,spec,request,destination,ownership,
            std::span<const algo::t2c_detail::DyadicEdgeSplit>(&split,1),{},
            &trace);
    CHECK(outcome.published());
    CHECK(outcome.surfaceOutcome!=nullptr);
    CHECK(outcome.surfaceOutcome&&outcome.surfaceOutcome->retainedArtifact);
    if(!outcome.surfaceOutcome||!outcome.surfaceOutcome->retainedArtifact)return;
    std::set<Index> sharedIndices;
    std::size_t faceUses=0;
    for(const auto& face:outcome.surfaceOutcome->retainedArtifact->faces)
        if(face)for(const auto& [local,global]:face->existingGlobalIndex) {
            (void)local;
            if(destination.nodeClass(global).dim==1&&
               destination.nodeClass(global).tag==edge->tag()&&
               destination.node(global).x==split.coordinate.x&&
               destination.node(global).y==split.coordinate.y&&
               destination.node(global).z==split.coordinate.z) {
                sharedIndices.insert(global);
                ++faceUses;
            }
        }
    CHECK(sharedIndices.size()==1);
    CHECK(faceUses>=2);
    CHECK(outcome.retainedResult&&outcome.retainedResult->boundaryAudit.passed);
}

void privateInterfaceFacetSubdivision() {
    TEST("t2c private interface facet point preserves both owner sides");
    LayeredModel model(2);
    Ownership ownership;
    const geo::GFace* interfaceFace=nullptr;
    for(const auto* face:model.faces()) {
        const auto& owners=model.owners(face->tag());
        ownership.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            owners;
        if(owners.size()==2)interfaceFace=face;
    }
    CHECK(interfaceFace!=nullptr);
    if(!interfaceFace)return;
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(10.0);
    spec.threads=1;
    const Vec2 uv{0.375,0.625};
    algo::SurfaceNodeIdentity faceIdentity;
    faceIdentity.domain=
        algo::SurfaceNodeIdentityDomain::RefinedFacetInterior;
    faceIdentity.parentArtifact={
        static_cast<algo::FacePersistentId>(interfaceFace->tag())+1U,1};
    // This private helper fixture supplies an authored canonical parent
    // token; the public coordinator derives this triple from its artifact.
    faceIdentity.parentFacetGlobalNodes={0,1,2};
    faceIdentity.parameterUBits=std::bit_cast<std::uint64_t>(uv.u);
    faceIdentity.parameterVBits=std::bit_cast<std::uint64_t>(uv.v);
    const algo::t2c_detail::PrivateFacePoint point{
        static_cast<algo::FacePersistentId>(interfaceFace->tag())+1,
        uv,interfaceFace->pointAt(uv),faceIdentity};
    algo::BoundaryConformingV3Request request;
    request.expectedSurfaceSetStamp={16901,1};
    request.expectedOwnershipProviderStamp=ownership.current;
    mesh::Mesh destination;
    algo::t2c_detail::PrivateBuildTrace trace;
    auto outcome=algo::t2c_detail::
        meshModelBoundaryConformingWithProtectedSplits(
            model,spec,request,destination,ownership,{},
            std::span<const algo::t2c_detail::PrivateFacePoint>(&point,1),
            &trace);
    CHECK(outcome.published());
    CHECK(outcome.retainedResult!=nullptr);
    if(!outcome.retainedResult)return;
    CHECK(outcome.retainedResult->boundaryAudit.passed);
    CHECK(outcome.retainedResult->boundaryAudit.missingInterface==0);
    CHECK(outcome.retainedResult->boundaryAudit.unexpectedInterface==0);
    std::size_t identityCount=0;
    for(Index node=0;node<destination.nodeCount();++node)
        if(destination.node(node).x==point.coordinate.x&&
           destination.node(node).y==point.coordinate.y&&
           destination.node(node).z==point.coordinate.z&&
           destination.nodeClass(node).dim==2&&
           destination.nodeClass(node).tag==interfaceFace->tag())
            ++identityCount;
    CHECK(identityCount==1);
    const auto& state=outcome.retainedResult->constraintState;
    std::size_t ownerFacets=0;
    for(const auto& [key,binding]:state.bindings) {
        (void)binding;
        if(key.facePersistentId==point.facePersistentId)++ownerFacets;
    }
    CHECK(ownerFacets>=3);
}

#if defined(PAACS_T5_SCALING_MAIN)
struct ScalingRow {
    std::size_t regions=0,interfaceFacets=0;
    std::uint64_t work=0,strictHash=0,semanticHash=0;
    std::array<std::uint8_t,32> authority{};
    std::array<std::uint64_t,static_cast<std::size_t>(
        algo::boundary_conforming_detail::BoundaryConformingWallPhase::Count)>
        phaseNanoseconds{};
    algo::TetCandidateMetrics tetCandidateMetrics{};
};

ScalingRow runScalingRow(std::size_t regionCount,unsigned threads,
                         bool independentPairs=true,
                         bool qualificationLayers=false) {
    const double layerSpacing=1.0;
    const double effectiveExtent=qualificationLayers
        ? static_cast<double>(regionCount)*0.5 : 1.0;
    LayeredModel model(static_cast<int>(regionCount),independentPairs,
                       effectiveExtent,layerSpacing);
    Ownership provider;
    for(const auto* face:model.faces())
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    sizing::MeshSpec spec;
    spec.field=qualificationLayers
        ? std::shared_ptr<sizing::SizeField>(
              std::make_shared<LayerQualificationField>(
                  layerSpacing,effectiveExtent))
        : std::shared_ptr<sizing::SizeField>(
              std::make_shared<sizing::ConstantField>(0.5));
    spec.threads=threads;
    std::vector<algo::T5RegionBinding> regions;
    regions.reserve(regionCount);
    for(std::size_t i=0;i<regionCount;++i)
        regions.push_back({static_cast<int>(i),7000+i,9000+i,
                           {3,static_cast<int>(i)}});
    algo::T5Request request;
    request.sourceRequest.expectedSurfaceSetStamp={16000+regionCount,1};
    request.sourceRequest.expectedOwnershipProviderStamp=provider.current;
    request.regionBindings=regions;
    mesh::Mesh destination;
    auto& wallTrace=algo::boundary_conforming_detail::
        boundaryConformingWallTraceForTesting();
    wallTrace.reset(true);
    auto outcome=algo::meshModelMultiRegionTet4(
        model,spec,request,destination,provider);
    const auto phaseNanoseconds=wallTrace.nanoseconds;
    const auto tetCandidateMetrics=wallTrace.tetCandidateMetrics;
    wallTrace.reset(false);
    if(!outcome.published())
        std::printf("  t5 scaling regions=%zu status=%d reason=%d phase=%d source_status=%d source_reason=%d source_phase=%d source_edge=(%lld,%lld) schedule=%zu witness=%llu\n",
            regionCount,static_cast<int>(outcome.status),
            static_cast<int>(outcome.failure.reason),
            static_cast<int>(outcome.failure.phase),
            static_cast<int>(outcome.failure.sourceStatus),
            static_cast<int>(outcome.failure.sourceFailure.reason),
            static_cast<int>(outcome.failure.sourceFailure.phase),
            static_cast<long long>(outcome.failure.sourceFailure.edgeKey.
                canonicalGlobalNodes[0]),
            static_cast<long long>(outcome.failure.sourceFailure.edgeKey.
                canonicalGlobalNodes[1]),
            outcome.failure.sourceFailure.scheduleOrdinal,
            static_cast<unsigned long long>(
                outcome.failure.sourceFailure.exactWitnessDigest));
#ifdef PAACS_T5_SCALING_MAIN
    if(!outcome.published()&&
       outcome.failure.sourceFailure.reason==
           algo::BoundaryConformingV3Reason::TetConstructionFailed) {
        mesh::Mesh surface;
        algo::meshEdges(model,spec,surface);
        const auto surfaceOutcome=algo::meshFacesConstrained(model,spec,surface);
        const auto candidate=algo::buildTetCandidate(surface.nodes());
        std::printf("  t5 diagnostic surface_status=%d nodes=%lld cells=%lld candidate_status=%d candidate_reason=%d\n",
            static_cast<int>(surfaceOutcome.status),
            static_cast<long long>(surface.nodeCount()),
            static_cast<long long>(surface.cellCount()),
            static_cast<int>(candidate.status),
            static_cast<int>(candidate.reason));
    }
#endif
    CHECK(outcome.published());
    CHECK(outcome.retainedResult!=nullptr);
    ScalingRow row;
    if(!outcome.retainedResult)return row;
    const auto& result=*outcome.retainedResult;
    row.regions=result.regions.size();
    row.interfaceFacets=result.interfaceFacets.size();
    row.work=result.metrics.totalWorkUnits;
    row.strictHash=result.strictHash;
    row.semanticHash=result.semanticHash.value;
    row.authority=result.semanticAuthorityCoreDigest;
    row.phaseNanoseconds=phaseNanoseconds;
    row.tetCandidateMetrics=tetCandidateMetrics;
    CHECK(result.audit.passed);
    CHECK(result.metrics.perInterfaceGlobalTetScans==0);
    CHECK(result.metrics.perRegionGlobalFaceScans==0);
    CHECK(result.metrics.contextRebuilds==0);
    CHECK(result.metrics.publicationTailAllocations==0);
    return row;
}

void publicScalingRoute() {
    TEST("t5 independent-pair phase diagnostic has exact rows");
    constexpr std::array<std::size_t,3> regionCounts{{2,8,32}};
    std::array<ScalingRow,3> rows{};
    for(std::size_t i=0;i<regionCounts.size();++i) {
        const auto started=std::chrono::steady_clock::now();
        rows[i]=runScalingRow(regionCounts[i],1);
        const auto elapsed=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-started).count();
        std::printf("  t5 scaling regions=%zu facets=%zu work=%llu seconds=%.3f\n",
            rows[i].regions,rows[i].interfaceFacets,
            static_cast<unsigned long long>(rows[i].work),elapsed);
        std::printf("  t5 phases_ms surface=%.3f artifact=%.3f envelope=%.3f tet=%.3f import=%.3f context=%.3f recovery=%.3f extraction=%.3f audit_publish=%.3f\n",
            rows[i].phaseNanoseconds[0]/1.0e6,
            rows[i].phaseNanoseconds[1]/1.0e6,
            rows[i].phaseNanoseconds[2]/1.0e6,
            rows[i].phaseNanoseconds[3]/1.0e6,
            rows[i].phaseNanoseconds[4]/1.0e6,
            rows[i].phaseNanoseconds[5]/1.0e6,
            rows[i].phaseNanoseconds[6]/1.0e6,
            rows[i].phaseNanoseconds[7]/1.0e6,
            rows[i].phaseNanoseconds[8]/1.0e6);
        const auto& construction=rows[i].tetCandidateMetrics.construction;
        std::printf("  t5 tet_stats points=%llu cells=%llu walk_cells=%llu conflict_queries=%llu cavity_cells=%llu fallback=%llu\n",
            static_cast<unsigned long long>(construction.inputPoints),
            static_cast<unsigned long long>(construction.outputCells),
            static_cast<unsigned long long>(construction.walkCellVisits),
            static_cast<unsigned long long>(construction.conflictQueries),
            static_cast<unsigned long long>(construction.cavityCellVisits),
            static_cast<unsigned long long>(construction.fallbackCellTests));
        std::fflush(stdout);
        CHECK(rows[i].regions==regionCounts[i]);
        CHECK(rows[i].interfaceFacets==16*regionCounts[i]);
        CHECK(rows[i].work>0);
        if(i>0) {
            CHECK(rows[i].work>rows[i-1].work);
            CHECK(rows[i].work<12*rows[i-1].work);
        }
    }
    const auto one=runScalingRow(2,1);
    const auto four=runScalingRow(2,4);
    const auto automatic=runScalingRow(2,0);
    CHECK(one.strictHash==four.strictHash);
    CHECK(one.strictHash==automatic.strictHash);
    CHECK(one.semanticHash==four.semanticHash);
    CHECK(one.semanticHash==automatic.semanticHash);
    CHECK(one.authority==four.authority);
    CHECK(one.authority==automatic.authority);
}
#endif

void fourRegionSharedMaterial() {
    TEST("t5 four pairwise regions permit a deliberately shared material");
    LayeredModel model(4);Ownership provider;
    for(const auto* face:model.faces())
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(10.0);
    std::array<algo::T5RegionBinding,4> regions{{
        {0,3101,4101,{3,0}},{1,3102,4101,{3,1}},
        {2,3103,4102,{3,2}},{3,3104,4103,{3,3}}}};
    algo::T5Request request;
    request.sourceRequest.expectedSurfaceSetStamp={15005,1};
    request.sourceRequest.expectedOwnershipProviderStamp=provider.current;
    request.regionBindings=regions;
    mesh::Mesh destination;
    auto outcome=algo::meshModelMultiRegionTet4(
        model,spec,request,destination,provider);
    CHECK(outcome.published());CHECK(outcome.retainedResult!=nullptr);
    if(!outcome.retainedResult)return;
    const auto& result=*outcome.retainedResult;
    CHECK(result.regions.size()==4);
    CHECK(result.materials.size()==3);
    CHECK(result.interfaceFacets.size()==6);
    CHECK(result.regionInterfaceFacets.size()==12);
    CHECK(result.audit.passed);
    std::set<std::pair<algo::RegionPersistentId,
                       algo::RegionPersistentId>> pairs;
    for(const auto& facet:result.interfaceFacets)
        pairs.emplace(facet.firstRegion,facet.secondRegion);
    const std::set<std::pair<algo::RegionPersistentId,
                             algo::RegionPersistentId>> expectedPairs{
        {3101,3102},{3102,3103},{3103,3104}};
    CHECK(pairs==expectedPairs);
    const auto shared=std::find_if(result.materials.begin(),
        result.materials.end(),[](const auto& row){
            return row.materialPersistentId==4101;});
    CHECK(shared!=result.materials.end());
    CHECK(shared!=result.materials.end()&&shared->regionReferenceCount==2);
    requireCompleteXor(outcome);
}

void preflightAndPostSourceBudgets() {
    TEST("t5 preflight rejects before providers and post-source budgets retain handoff");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(0.55);
    Ownership preflightProvider;configure(*model,preflightProvider);
    std::array<algo::T5RegionBinding,2> duplicate{{
        {0,101,201,{3,0}},{0,102,202,{3,0}}}};
    algo::T5Request invalid;
    invalid.sourceRequest.expectedSurfaceSetStamp={15004,1};
    invalid.sourceRequest.expectedOwnershipProviderStamp=
        preflightProvider.current;
    invalid.regionBindings=duplicate;
    mesh::Mesh destination;
    destination.addNode({-1,-2,-3},{0,99});
    const auto before=mesh::meshHash(destination);
    auto rejected=algo::meshModelMultiRegionTet4(
        *model,spec,invalid,destination,preflightProvider);
    CHECK(rejected.status==algo::T5Status::Rejected);
    CHECK(rejected.failure.reason==algo::T5Reason::DuplicateRegionTag);
    CHECK(preflightProvider.stampCalls==0);
    CHECK(preflightProvider.ownershipCalls==0);
    CHECK(preflightProvider.envelopeCalls==0);
    CHECK(mesh::meshHash(destination)==before);

    const int rounding=std::fegetround();
    if(std::fesetround(FE_DOWNWARD)==0) {
        Ownership fpProvider;configure(*model,fpProvider);
        const algo::T5RegionBinding fpRegion{0,301,401,{3,0}};
        auto fpRequest=singleRequest(fpProvider,fpRegion);
        auto unsupported=algo::meshModelMultiRegionTet4(
            *model,spec,fpRequest,destination,fpProvider);
        CHECK(unsupported.status==algo::T5Status::Unsupported);
        CHECK(unsupported.failure.reason==
              algo::T5Reason::FloatingPointEnvironmentUnsupported);
        CHECK(fpProvider.stampCalls==0);
        CHECK(mesh::meshHash(destination)==before);
        CHECK(std::fesetround(rounding)==0);
    }

    Ownership budgetProvider;configure(*model,budgetProvider);
    const algo::T5RegionBinding region{0,101,201,{3,0}};
    auto limited=singleRequest(budgetProvider,region);
    limited.maxPublicationNodes=0;
    auto incomplete=algo::meshModelMultiRegionTet4(
        *model,spec,limited,destination,budgetProvider);
    CHECK(incomplete.status==algo::T5Status::Incomplete);
    CHECK(incomplete.failure.reason==algo::T5Reason::ResultBudgetReached);
    CHECK(incomplete.failure.phase==algo::T5Phase::MeshBuild);
    CHECK(incomplete.sourceOutcome==nullptr);
    CHECK(incomplete.sourceReceipt!=nullptr);
    CHECK(incomplete.retainedResult!=nullptr);
    CHECK(incomplete.retainedResult&&
          incomplete.retainedResult->privateMesh==nullptr);
    CHECK(mesh::meshHash(destination)==before);

    Ownership workProvider;configure(*model,workProvider);
    auto workLimited=singleRequest(workProvider,region);
    workLimited.maxAuditWork=0;
    auto workStop=algo::meshModelMultiRegionTet4(
        *model,spec,workLimited,destination,workProvider);
    CHECK(workStop.status==algo::T5Status::Incomplete);
    CHECK(workStop.failure.reason==algo::T5Reason::WorkBudgetReached);
    CHECK(workStop.failure.phase==algo::T5Phase::InterfaceAudit);
    CHECK(mesh::meshHash(destination)==before);

    Ownership retryProvider;configure(*model,retryProvider);
    auto cleanRequest=singleRequest(retryProvider,region);
    auto retry=algo::meshModelMultiRegionTet4(
        *model,spec,cleanRequest,destination,retryProvider);
    CHECK(retry.published());
}

void providerFinalCheckpoint() {
    TEST("t5 provider change at the final checkpoint cannot publish stale authority");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(0.55);
    const algo::T5RegionBinding region{0,501,601,{3,0}};
    Ownership reference;configure(*model,reference);
    auto referenceRequest=singleRequest(reference,region);
    mesh::Mesh referenceMesh;
    auto clean=algo::meshModelMultiRegionTet4(
        *model,spec,referenceRequest,referenceMesh,reference);
    CHECK(clean.published());
    CHECK(reference.stampCalls>1);

    Ownership changing;configure(*model,changing);
    changing.mutateAtStampCall=reference.stampCalls;
    auto changingRequest=singleRequest(changing,region);
    mesh::Mesh destination;destination.addNode({-7,-11,-13},{0,88});
    const auto before=mesh::meshHash(destination);
    auto stale=algo::meshModelMultiRegionTet4(
        *model,spec,changingRequest,destination,changing);
    CHECK(stale.status==algo::T5Status::Failed);
    CHECK(stale.failure.reason==algo::T5Reason::StaleOwnershipProvider);
    CHECK(stale.failure.phase==algo::T5Phase::InterfaceAudit);
    CHECK(stale.sourceOutcome!=nullptr);
    CHECK(stale.sourceReceipt==nullptr);
    CHECK(mesh::meshHash(destination)==before);
}

void exportBindingRoundTrip() {
    TEST("t5 VTK and MSH exports bind independently to one authority core");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    Ownership provider;configure(*model,provider);
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(0.55);
    const algo::T5RegionBinding region{0,101,201,{3,0}};
    auto request=singleRequest(provider,region);
    mesh::Mesh destination;
    auto outcome=algo::meshModelMultiRegionTet4(
        *model,spec,request,destination,provider);
    CHECK(outcome.published());
    CHECK(outcome.retainedResult!=nullptr);
    if(!outcome.retainedResult)return;
    const auto nonce=std::chrono::steady_clock::now().time_since_epoch().count();
    const auto root=std::filesystem::temp_directory_path()/
        ("paacs_t5_export_"+std::to_string(nonce));
    std::filesystem::create_directories(root);
    const auto vtk=(root/"mesh.vtk").string();
    const auto vtkCore=(root/"mesh.vtk.t5core").string();
    const auto vtkBinding=(root/"mesh.vtk.t5binding").string();
    const auto msh=(root/"mesh.msh").string();
    const auto mshCore=(root/"mesh.msh.t5core").string();
    const auto mshBinding=(root/"mesh.msh.t5binding").string();
    const auto vtkOut=algo::writeT5VtkExport(
        destination,*outcome.retainedResult,vtk,vtkCore,vtkBinding);
    const auto mshOut=algo::writeT5MshExport(
        destination,*outcome.retainedResult,msh,mshCore,mshBinding);
    CHECK(vtkOut.complete);CHECK(mshOut.complete);
    CHECK(vtkOut.binding.exportBindingDigest!=
          mshOut.binding.exportBindingDigest);
    const auto vtkReload=algo::reloadT5ExportBinding(
        vtk,vtkCore,vtkBinding);
    const auto mshReload=algo::reloadT5ExportBinding(
        msh,mshCore,mshBinding);
    CHECK(vtkReload.passed);CHECK(mshReload.passed);
    CHECK(vtkReload.meshHashesMatch);
    CHECK(mshReload.meshHashesMatch);
    CHECK(vtkReload.authorityProjectionValid);
    CHECK(mshReload.authorityProjectionValid);
    CHECK(vtkReload.authorityMappingsMatch);
    CHECK(mshReload.authorityMappingsMatch);
    CHECK(vtkReload.strictHash==mesh::meshHash(destination));
    CHECK(mshReload.strictHash==vtkReload.strictHash);
    CHECK(mshReload.semanticHash.value==vtkReload.semanticHash.value);
    CHECK(vtkReload.nodeCount==static_cast<std::uint64_t>(destination.nodeCount()));
    CHECK(vtkReload.cellCount==static_cast<std::uint64_t>(destination.cellCount()));
    CHECK(mshReload.nodeCount==vtkReload.nodeCount);
    CHECK(mshReload.cellCount==vtkReload.cellCount);
    {
        std::fstream corrupt(vtkBinding,
            std::ios::binary|std::ios::in|std::ios::out);
        corrupt.seekp(-1,std::ios::end);char value=0;
        corrupt.write(&value,1);
    }
    const auto corrupted=algo::reloadT5ExportBinding(
        vtk,vtkCore,vtkBinding);
    CHECK(!corrupted.passed);
    std::filesystem::remove_all(root);
}

struct DeterministicSignature {
    std::uint64_t meshStrict=0,meshSemantic=0;
    std::array<std::uint8_t,32> semanticAuthority{};
    std::vector<algo::T5ProcessedRow> processed;
    std::vector<algo::T5InterfaceFacetRow> facets;
};

DeterministicSignature runTwoRegion(unsigned threads) {
    LayeredModel model;Ownership provider;
    for(const auto* face:model.faces())
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(10.0);
    spec.threads=threads;
    std::array<algo::T5RegionBinding,2> regions{{
        {0,1001,2001,{3,0}},{1,1002,2002,{3,1}}}};
    algo::T5Request request;
    request.sourceRequest.expectedSurfaceSetStamp={15003,1};
    request.sourceRequest.expectedOwnershipProviderStamp=provider.current;
    request.regionBindings=regions;
    mesh::Mesh mesh;
    auto outcome=algo::meshModelMultiRegionTet4(
        model,spec,request,mesh,provider);
    CHECK(outcome.published());
    CHECK(outcome.retainedResult!=nullptr);CHECK(outcome.report!=nullptr);
    DeterministicSignature result;
    if(!outcome.retainedResult||!outcome.report)return result;
    result.meshStrict=outcome.retainedResult->strictHash;
    result.meshSemantic=outcome.retainedResult->semanticHash.value;
    result.semanticAuthority=
        outcome.retainedResult->semanticAuthorityCoreDigest;
    result.processed=outcome.report->processedStorage;
    result.facets=outcome.retainedResult->interfaceFacets;
    return result;
}

bool sameWorkKey(const algo::T5WorkKey& a,const algo::T5WorkKey& b) noexcept {
    return a.kind==b.kind&&a.region==b.region&&a.material==b.material&&
        a.globalNode==b.globalNode&&a.globalEdge==b.globalEdge&&
        a.facet==b.facet&&a.surfaceFace==b.surfaceFace&&
        a.normalizedLocalCellOrdinal==b.normalizedLocalCellOrdinal&&
        a.sourceSlot==b.sourceSlot;
}

void threadDeterminism() {
    TEST("t5 1 4 auto routes preserve mesh report and semantic authority");
    const auto one=runTwoRegion(1),four=runTwoRegion(4),automatic=runTwoRegion(0);
    CHECK(one.meshStrict==four.meshStrict);
    CHECK(one.meshStrict==automatic.meshStrict);
    CHECK(one.meshSemantic==four.meshSemantic);
    CHECK(one.meshSemantic==automatic.meshSemantic);
    CHECK(one.semanticAuthority==four.semanticAuthority);
    CHECK(one.semanticAuthority==automatic.semanticAuthority);
    CHECK(one.processed.size()==four.processed.size());
    CHECK(one.processed.size()==automatic.processed.size());
    for(std::size_t i=0;i<one.processed.size()&&i<four.processed.size()&&
        i<automatic.processed.size();++i) {
        CHECK(sameWorkKey(one.processed[i].key,four.processed[i].key));
        CHECK(sameWorkKey(one.processed[i].key,automatic.processed[i].key));
        CHECK(one.processed[i].status==four.processed[i].status);
        CHECK(one.processed[i].status==automatic.processed[i].status);
        CHECK(one.processed[i].canonicalOrdinal==
              four.processed[i].canonicalOrdinal);
        CHECK(one.processed[i].canonicalOrdinal==
              automatic.processed[i].canonicalOrdinal);
    }
    CHECK(one.facets.size()==four.facets.size());
    CHECK(one.facets.size()==automatic.facets.size());
    for(std::size_t i=0;i<one.facets.size()&&i<four.facets.size()&&
        i<automatic.facets.size();++i) {
        CHECK(one.facets[i].facet==four.facets[i].facet);
        CHECK(one.facets[i].facet==automatic.facets[i].facet);
        CHECK(one.facets[i].canonicalPublishedNodes==
              four.facets[i].canonicalPublishedNodes);
        CHECK(one.facets[i].canonicalPublishedNodes==
              automatic.facets[i].canonicalPublishedNodes);
        CHECK(one.facets[i].firstRegion==four.facets[i].firstRegion);
        CHECK(one.facets[i].secondRegion==four.facets[i].secondRegion);
    }
}

} // namespace

namespace paacs_tests {
struct Entry { const char* name;void(*fn)(); };
const Entry* t5Tests(std::size_t& count) {
    static const Entry tests[]{{"t5.capacity_snapshot",capacitySnapshotContract},
        {"t5.capacity_dedup",capacitySnapshotDedupAndOverflow},
        {"t5.sha256",sha256KnownAnswer},
        {"t5.exact_volume",exactVolumeContract},
        {"t5.single_region",singleRegionRoute},
        {"t5.two_region",twoRegionRoute},
        {"t5.fine_interface",fineInterfaceRoute},
        {"t5.real_layer_noncycle",realLayerNonCycleRefusal},
        {"t2c.real_layer_v4",realLayerV4ProtectedSubdivision},
        {"t2c.facet_conflict_v4",realLayerV4FacetConflictSubdivision},
        {"t2c.facet_conflict_preflight",facetConflictRequestPreflight},
        {"t2c.facet_conflict_provider",facetConflictProviderFailures},
        {"t2c.cache_all_dirty_equivalence",allDirtyCacheMatchesUncachedSurface},
        {"t2c.cache_corruption",cachedArtifactCorruptionIsRejected},
        {"t2c.cache_shared_curve_dirty",sharedCurveDirtyPropagationIsExact},
        {"t2c.cache_budget_transport",surfaceCacheBudgetTransportAndDigest},
        {"t2c.fast_single_pass",singlePassInterfacePrefillRoute},
        {"t2c.fast_direct",singlePassInterfacePrefillDirectRoute},
        {"t2c.v4_prepared",singlePassPreparedConstraintPrefillRoute},
        {"t2c.v4_determinism",singlePassPreparedConstraintDeterminism},
        {"t2c.v6_sizing_controls",sizingDrivenSinglePassPrefillRoute},
        {"t2c.v6_finer_body_control",sizingDrivenFinerBodyControlRoute},
        {"t2c.fast_determinism",singlePassInterfacePrefillDeterminism},
        {"t2c.fast_incomplete",singlePassInterfacePrefillRefusesIncompleteMesh},
        {"t2c.fast_contracts",singlePassInterfacePrefillFailureContracts},
        {"t2c.contracts",v4ProviderBudgetAndDigestContracts},
        {"t2c.private_edge",privateCadEdgeSubdivision},
        {"t2c.private_facet",privateInterfaceFacetSubdivision},
        {"t5.four_region_shared_material",fourRegionSharedMaterial},
        {"t5.budgets",preflightAndPostSourceBudgets},
        {"t5.provider_final_checkpoint",providerFinalCheckpoint},
        {"t5.export_binding",exportBindingRoundTrip},
        {"t5.thread_determinism",threadDeterminism}};
    count=std::size(tests);return tests;
}
} // namespace paacs_tests

#ifdef PAACS_T5_SCALING_MAIN
int gFailures=0;
int main(int argc,char** argv) {
    const bool layered=argc>=3&&std::string_view(argv[1])=="--layers";
    if(argc==2||layered) {
        const char* countText=layered?argv[2]:argv[1];
        const auto regions=static_cast<std::size_t>(
            std::strtoull(countText,nullptr,10));
        TEST(layered?"t5 selected contiguous-layer qualification row":
                     "t5 selected independent-pair performance diagnostic row");
        const auto started=std::chrono::steady_clock::now();
        const auto row=runScalingRow(
            regions,1,!layered,layered);
        const auto elapsed=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-started).count();
        std::printf("  t5 scaling regions=%zu facets=%zu work=%llu seconds=%.3f\n",
            row.regions,row.interfaceFacets,
            static_cast<unsigned long long>(row.work),elapsed);
        std::printf("  t5 phases_ms surface=%.3f artifact=%.3f envelope=%.3f tet=%.3f import=%.3f context=%.3f recovery=%.3f extraction=%.3f audit_publish=%.3f\n",
            row.phaseNanoseconds[0]/1.0e6,
            row.phaseNanoseconds[1]/1.0e6,
            row.phaseNanoseconds[2]/1.0e6,
            row.phaseNanoseconds[3]/1.0e6,
            row.phaseNanoseconds[4]/1.0e6,
            row.phaseNanoseconds[5]/1.0e6,
            row.phaseNanoseconds[6]/1.0e6,
            row.phaseNanoseconds[7]/1.0e6,
            row.phaseNanoseconds[8]/1.0e6);
        const auto& construction=row.tetCandidateMetrics.construction;
        std::printf("  t5 tet_stats points=%llu cells=%llu walk_cells=%llu conflict_queries=%llu cavity_cells=%llu fallback=%llu\n",
            static_cast<unsigned long long>(construction.inputPoints),
            static_cast<unsigned long long>(construction.outputCells),
            static_cast<unsigned long long>(construction.walkCellVisits),
            static_cast<unsigned long long>(construction.conflictQueries),
            static_cast<unsigned long long>(construction.cavityCellVisits),
            static_cast<unsigned long long>(construction.fallbackCellTests));
        CHECK(row.regions==regions);
        CHECK(row.interfaceFacets==16*regions);
    } else {
        publicScalingRoute();
    }
    if(gFailures) {
        std::printf("\n%d FAILURE(S)\n",gFailures);
        return 1;
    }
    std::printf("\nT5 PUBLIC SCALING PASSED\n");
    return 0;
}
#endif

#ifdef PAACS_T2C_SCALING_MAIN
int gFailures=0;
int main(int argc,char** argv) {
    if(argc<2||argc>4) {
        std::printf("usage: paacs_t2c_scaling <regions> [batch-cap] "
                    "[iterative|fast|direct|v4|v6]\n");
        return 2;
    }
    const auto regions=static_cast<std::size_t>(
        std::strtoull(argv[1],nullptr,10));
    if(regions!=2&&regions!=8&&regions!=32) {
        std::printf("regions must be one of 2, 8, 32\n");
        return 2;
    }
    const bool diagnostic=argc==3;
    const auto batchCap=argc>=3
        ? static_cast<std::uint64_t>(std::strtoull(argv[2],nullptr,10))
        : UINT64_C(128);
    const std::string_view strategy=argc==4
        ? std::string_view(argv[3]):std::string_view("iterative");
    const bool direct=strategy=="direct";
    const bool v4=strategy=="v4";
    const bool v6=strategy=="v6";
    const bool fast=direct||v4||v6||strategy=="fast";
    if(!fast&&strategy!="iterative") {
        std::printf(
            "strategy must be iterative, fast, direct, v4, or v6\n");
        return 2;
    }
    if(batchCap==0) {
        std::printf("batch-cap must be positive\n");
        return 2;
    }
    LayeredModel model(static_cast<int>(regions),false,4.0,1.0);
    Ownership ownership;
    for(const auto* face:model.faces())
        ownership.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
    CurvePlacement placement;
    sizing::MeshSpec spec;
    spec.field=std::make_shared<LayerQualificationField>(1.0,4.0);
    spec.threads=1;
    if(const char* text=std::getenv("PAACS_T2C_THREADS")) {
        const auto parsed=std::strtoull(text,nullptr,10);
        if(parsed>1024) {
            std::printf("PAACS_T2C_THREADS must be in [0,1024]\n");
            return 2;
        }
        spec.threads=static_cast<std::size_t>(parsed);
    }
    algo::BoundaryConformingV4Request request;
    request.sourceRequest.expectedSurfaceSetStamp={17000,1};
    request.sourceRequest.expectedOwnershipProviderStamp=ownership.current;
    request.expectedPlacementProviderStamp=placement.current;
    const auto splitLimit=static_cast<std::uint64_t>(regions*16U+16U);
    request.limits={splitLimit,splitLimit+1U,splitLimit,
                    static_cast<std::uint64_t>(regions*5000U+5000U),
                    8ULL*1024ULL*1024ULL*1024ULL,UINT64_C(1)<<50,
                    batchCap,4096,4096};
    if(fast) {
        std::uint32_t prefillGrid=3;
        if(const char* text=std::getenv("PAACS_T2C_PREFILL_GRID")) {
            const auto parsed=std::strtoull(text,nullptr,10);
            if(parsed==0||parsed>1000) {
                std::printf("PAACS_T2C_PREFILL_GRID must be in [1,1000]\n");
                return 2;
            }
            prefillGrid=static_cast<std::uint32_t>(parsed);
        }
        request.executionStrategy=(v4||v6)
            ? algo::ProtectedFeatureExecutionStrategy::
                SinglePassPreparedConstraintPrefill
            : direct
                ? algo::ProtectedFeatureExecutionStrategy::
                    SinglePassInterfacePrefillDirect
                : algo::ProtectedFeatureExecutionStrategy::
                    SinglePassInterfacePrefill;
        request.interfacePrefill={prefillGrid,prefillGrid,
            static_cast<std::uint64_t>(regions),splitLimit,3.5};
        if(v6) {
            sizing::MeshControlSpecification controlSpec;
            controlSpec.minimumSize=0.05;
            controlSpec.maximumSize=4.0;
            controlSpec.maximumGrowthRate=1.3;
            controlSpec.entitySizes.reserve(regions);
            for(std::size_t solid=0;solid<regions;++solid)
                controlSpec.entitySizes.push_back(
                    {3,static_cast<int>(solid),4.0});
            const auto compiled=sizing::compileMeshControls(controlSpec);
            if(!compiled.compiled) {
                std::printf("failed to compile v6 mesh controls\n");
                return 2;
            }
            request.interfacePrefill.densityMode=
                algo::InterfacePrefillDensityMode::SizingDriven;
            request.interfacePrefill.minimumInteriorU=prefillGrid;
            request.interfacePrefill.minimumInteriorV=prefillGrid;
            request.interfacePrefill.maximumInteriorU=prefillGrid;
            request.interfacePrefill.maximumInteriorV=prefillGrid;
            request.interfacePrefill.sizingScale=1.0;
            request.meshControls=
                std::make_shared<const sizing::CompiledMeshControls>(
                    compiled.controls);
        }
        request.limits.maxRounds=1;
        request.facetConflictMode=
            algo::FacetConflictRefinementMode::Refuse;
    } else {
        request.facetConflictMode=algo::FacetConflictRefinementMode::
            BoundedCanonicalSccSubdivision;
        request.facetConflictLimits={64,splitLimit,splitLimit,
            splitLimit*3U+1U,UINT64_C(1)<<34};
    }
    request.recoverySchedulingMode=algo::TetRecoverySchedulingMode::
        MixedExactPresentBatching;
    mesh::Mesh destination;
    auto& wallTrace=algo::boundary_conforming_detail::
        boundaryConformingWallTraceForTesting();
    wallTrace.reset(true);
    const auto started=std::chrono::steady_clock::now();
    const auto outcome=algo::meshModelBoundaryConformingV4(
        model,spec,request,destination,ownership,placement);
    const auto seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-started).count();
    wallTrace.reset(false);
    const auto* result=outcome.retainedResult.get();
    const auto* source=result&&result->sourceOutcome
        ? result->sourceOutcome->retainedResult.get():nullptr;
    std::printf("t2c_scaling regions=%zu strategy=%s published=%d "
                "status=%d reason=%d batch=%llu edge_splits=%llu "
                "facet_splits=%llu prefill_candidates=%llu "
                "prefill_faces=%llu prefill_dense_skipped=%llu "
                "prefill_points=%llu prefill_work=%llu attempts=%llu "
                "frontier=%llu "
                "unique=%llu deferred=%llu work=%llu coordinator=%llu "
                "nodes=%lld tets=%lld strict=%llu semantic=%llu seconds=%.3f\n",
        regions,v6?"v6":v4?"v4":direct?"direct":fast?"fast":"iterative",
        outcome.published()?1:0,
        static_cast<int>(outcome.status),static_cast<int>(outcome.reason),
        static_cast<unsigned long long>(batchCap),
        static_cast<unsigned long long>(result?result->metrics.edgeSplits:0),
        static_cast<unsigned long long>(result?result->metrics.facetSplits:0),
        static_cast<unsigned long long>(result?
            result->metrics.interfacePrefillCandidateFaces:0),
        static_cast<unsigned long long>(result?
            result->metrics.interfacePrefillFaces:0),
        static_cast<unsigned long long>(result?
            result->metrics.interfacePrefillDenseFacesSkipped:0),
        static_cast<unsigned long long>(result?
            result->metrics.interfacePrefillPoints:0),
        static_cast<unsigned long long>(result?
            result->metrics.interfacePrefillWorkUnits:0),
        static_cast<unsigned long long>(result?result->metrics.rounds:0),
        static_cast<unsigned long long>(result?result->metrics.frontierRows:0),
        static_cast<unsigned long long>(result?result->metrics.uniquePlacements:0),
        static_cast<unsigned long long>(result?result->metrics.deferredFrontierRows:0),
        static_cast<unsigned long long>(result?
            result->metrics.totalMeasuredWorkUnits:0),
        static_cast<unsigned long long>(result?
            result->metrics.coordinatorWorkUnits:0),
        static_cast<long long>(destination.nodeCount()),
        static_cast<long long>(destination.countCells(mesh::CellType::Tet4)),
        static_cast<unsigned long long>(mesh::meshHash(destination)),
        static_cast<unsigned long long>(mesh::semanticHashEx(destination).value),
        seconds);
    const auto quality=mesh::analyzeTets(destination);
    std::printf("t2c_quality rr_min=%.17g rr_mean=%.17g sj_min=%.17g "
                "skew_max=%.17g dih_min=%.17g dih_max=%.17g "
                "inverted=%lld volume=%.17g\n",
        quality.radiusRatio.min,quality.radiusRatio.mean,
        quality.scaledJacobian.min,quality.skewness.max,
        quality.minDihedralDeg,quality.maxDihedralDeg,
        static_cast<long long>(quality.nInverted),quality.totalVolume);
    const auto* terminal=outcome.lastSourceOutcome
        ? outcome.lastSourceOutcome.get()
        : result&&result->sourceOutcome?result->sourceOutcome.get():nullptr;
    if(terminal)
        std::printf("t2c_terminal source_status=%d source_reason=%d "
                    "source_phase=%d edge=(%lld,%lld) nested_reason=%d "
                    "nested_phase=%d\n",
            static_cast<int>(terminal->status),
            static_cast<int>(terminal->failure.reason),
            static_cast<int>(terminal->failure.phase),
            static_cast<long long>(terminal->failure.edgeKey.
                canonicalGlobalNodes[0]),
            static_cast<long long>(terminal->failure.edgeKey.
                canonicalGlobalNodes[1]),
            static_cast<int>(terminal->failure.edgeOutcome.reason),
            static_cast<int>(terminal->failure.edgeOutcome.phase));
    std::printf("t2c_phases_ms surface=%.3f artifact=%.3f envelope=%.3f "
                "tet=%.3f import=%.3f context=%.3f recovery=%.3f "
                "extraction=%.3f audit_publish=%.3f\n",
        (result?result->metrics.nestedPhaseNanoseconds[0]:0)/1.0e6,
        (result?result->metrics.nestedPhaseNanoseconds[1]:0)/1.0e6,
        (result?result->metrics.nestedPhaseNanoseconds[2]:0)/1.0e6,
        (result?result->metrics.nestedPhaseNanoseconds[3]:0)/1.0e6,
        (result?result->metrics.nestedPhaseNanoseconds[4]:0)/1.0e6,
        (result?result->metrics.nestedPhaseNanoseconds[5]:0)/1.0e6,
        (result?result->metrics.nestedPhaseNanoseconds[6]:0)/1.0e6,
        (result?result->metrics.nestedPhaseNanoseconds[7]:0)/1.0e6,
        (result?result->metrics.nestedPhaseNanoseconds[8]:0)/1.0e6);
    std::printf("t2c_cache built=%llu reused=%llu rejected=%llu "
                "key_rows=%llu dirty_visits=%llu prep_work=%llu "
                "replay_work=%llu audit_work=%llu copied_bytes=%llu "
                "high_slots=%llu high_rows=%llu accounted=%zu "
                "prepared_recovery_faces=%llu direct_attempts=%llu "
                "direct_faces=%llu\n",
        static_cast<unsigned long long>(result?
            result->metrics.surfaceCacheFacesBuilt:0),
        static_cast<unsigned long long>(result?
            result->metrics.surfaceCacheFacesReused:0),
        static_cast<unsigned long long>(result?
            result->metrics.surfaceCacheFacesRejected:0),
        static_cast<unsigned long long>(result?
            result->metrics.surfaceCacheKeyRowsVisited:0),
        static_cast<unsigned long long>(result?
            result->metrics.surfaceCacheDirtyPropagationVisits:0),
        static_cast<unsigned long long>(result?
            result->metrics.surfaceCachePreparationWorkUnits:0),
        static_cast<unsigned long long>(result?
            result->metrics.surfaceCacheReplayWorkUnits:0),
        static_cast<unsigned long long>(result?
            result->metrics.surfaceCacheSetAuditWorkUnits:0),
        static_cast<unsigned long long>(result?
            result->metrics.surfaceCacheCopiedBytes:0),
        static_cast<unsigned long long>(result?
            result->metrics.surfaceCacheHighWaterSlots:0),
        static_cast<unsigned long long>(result?
            result->metrics.surfaceCacheHighWaterRows:0),
        result?result->metrics.surfaceCacheAccountedBytes:0,
        static_cast<unsigned long long>(result?
            result->metrics.preparedSurfaceRecoveryFaces:0),
        static_cast<unsigned long long>(result?
            result->metrics.directSurfaceBuildAttempts:0),
        static_cast<unsigned long long>(result?
            result->metrics.directSurfaceBuildFaces:0));
    if(source) {
        const auto& surface=source->constructionMetrics.surface;
        std::printf("t2c_surface_metrics faces=%llu boundary_nodes=%llu "
                    "interior=%llu spacing_queries=%llu spacing_buckets=%llu "
                    "spacing_tests=%llu constraint_tests=%llu input=%llu "
                    "coord_sort=%llu order_sort=%llu location_tests=%llu "
                    "walk=%llu cavity_tests=%llu cavity_members=%llu "
                    "boundary_edges=%llu triangles=%llu slots=%llu "
                    "fallbacks=%llu complete=%d\n",
            static_cast<unsigned long long>(surface.faceRows),
            static_cast<unsigned long long>(surface.boundaryNodes),
            static_cast<unsigned long long>(surface.interiorCandidates),
            static_cast<unsigned long long>(surface.spacingQueries),
            static_cast<unsigned long long>(surface.spacingBucketsVisited),
            static_cast<unsigned long long>(surface.spacingDistanceTests),
            static_cast<unsigned long long>(surface.constraintPointTests),
            static_cast<unsigned long long>(surface.scalable.inputPoints),
            static_cast<unsigned long long>(
                surface.scalable.coordinateSortComparisons),
            static_cast<unsigned long long>(
                surface.scalable.insertionOrderComparisons),
            static_cast<unsigned long long>(
                surface.scalable.locationTriangleTests),
            static_cast<unsigned long long>(surface.scalable.walkSteps),
            static_cast<unsigned long long>(
                surface.scalable.cavityPredicateTests),
            static_cast<unsigned long long>(surface.scalable.cavityMembers),
            static_cast<unsigned long long>(surface.scalable.boundaryEdges),
            static_cast<unsigned long long>(surface.scalable.trianglesCommitted),
            static_cast<unsigned long long>(
                surface.scalable.materializedTriangleSlots),
            static_cast<unsigned long long>(surface.scalable.fullScanFallbacks),
            surface.runComplete?1:0);
        const auto& metrics=source->recoveryMetrics;
        std::printf("t2c_recovery assessments=%llu edge_queries=%llu "
                    "facet_queries=%llu index_nodes=%llu returned=%llu "
                    "segment_tests=%llu target_tests=%llu survivor_tests=%llu "
                    "star_queries=%llu apex_pairs=%llu edge_commits=%llu "
                    "facet_commits=%llu dependency_edges=%llu invalidations=%llu "
                    "selected_rebuilds=%llu batch_preflight=%llu "
                    "scc_visits=%llu stage_rescans=%llu global_fallbacks=%llu\n",
            static_cast<unsigned long long>(metrics.routeAssessments),
            static_cast<unsigned long long>(metrics.edgeAabbQueries),
            static_cast<unsigned long long>(metrics.facetAabbQueries),
            static_cast<unsigned long long>(metrics.indexNodesVisited),
            static_cast<unsigned long long>(metrics.returnedTetCandidates),
            static_cast<unsigned long long>(metrics.segmentTetTests),
            static_cast<unsigned long long>(metrics.targetTetTests),
            static_cast<unsigned long long>(metrics.candidateSurvivorTests),
            static_cast<unsigned long long>(metrics.completeStarQueries),
            static_cast<unsigned long long>(metrics.apexPairs),
            static_cast<unsigned long long>(metrics.edgeCavityCommits),
            static_cast<unsigned long long>(metrics.facetCommits),
            static_cast<unsigned long long>(metrics.dependencyEdges),
            static_cast<unsigned long long>(metrics.dependencyInvalidations),
            static_cast<unsigned long long>(metrics.selectedRebuilds),
            static_cast<unsigned long long>(metrics.batchPreflightAssessments),
            static_cast<unsigned long long>(metrics.conflictSccVisits),
            static_cast<unsigned long long>(metrics.fullStageRescans),
            static_cast<unsigned long long>(
                metrics.productionGlobalScanFallbacks));
        const auto& construction=source->constructionMetrics;
        std::printf("t2c_partial_import attempts=%llu rows=%llu slots=%llu "
                    "mapping=%llu edge_rows=%llu face_rows=%llu "
                    "endpoints=%llu authority=%llu duplicates=%llu "
                    "resolution_nodes=%llu edge_lookups=%llu "
                    "face_lookups=%llu edge_data=%llu face_data=%llu "
                    "face_incidence=%llu sort=%llu import_work=%llu\n",
            static_cast<unsigned long long>(construction.importAttempts),
            static_cast<unsigned long long>(
                construction.partialImportRowsCopied),
            static_cast<unsigned long long>(
                construction.partialImportGlobalSlotsInitialized),
            static_cast<unsigned long long>(
                construction.partialImportCandidateMappingVisits),
            static_cast<unsigned long long>(
                construction.partialImportEdgeValidationRows),
            static_cast<unsigned long long>(
                construction.partialImportFaceValidationRows),
            static_cast<unsigned long long>(
                construction.partialImportProtectedEndpointVisits),
            static_cast<unsigned long long>(
                construction.partialImportAuthorityPointValidations),
            static_cast<unsigned long long>(
                construction.partialImportDuplicateKeyComparisons),
            static_cast<unsigned long long>(
                construction.partialImportResolutionNodeVisits),
            static_cast<unsigned long long>(
                construction.partialImportEdgeResolutionLookups),
            static_cast<unsigned long long>(
                construction.partialImportFaceResolutionLookups),
            static_cast<unsigned long long>(
                construction.partialImportEdgeDataReads),
            static_cast<unsigned long long>(
                construction.partialImportFaceDataCopies),
            static_cast<unsigned long long>(
                construction.partialImportFaceIncidenceRowsCopied),
            static_cast<unsigned long long>(
                construction.partialImportSortComparisons),
            static_cast<unsigned long long>(construction.importWorkUnits));
    }
    if(result) {
        const auto& mixed=result->metrics.mixedPresentScheduling;
        std::printf("t2c_mixed batches=%llu present=%llu copies=%llu "
                    "copied_rows=%llu audits=%llu pending=%llu pairs=%llu "
                    "nontrivial=%llu latest_nontrivial=%llu "
                    "present_rebuilds=%llu "
                    "present_dependencies=%llu complete=%d\n",
            static_cast<unsigned long long>(mixed.presentBindingBatches),
            static_cast<unsigned long long>(mixed.presentBindingEntities),
            static_cast<unsigned long long>(mixed.stateCopies),
            static_cast<unsigned long long>(mixed.copiedRows),
            static_cast<unsigned long long>(
                mixed.schedulerBindingAuditRows),
            static_cast<unsigned long long>(mixed.pendingRowsVisited),
            static_cast<unsigned long long>(mixed.dependencyPairsTested),
            static_cast<unsigned long long>(mixed.nontrivialRecoveries),
            static_cast<unsigned long long>(result->metrics.
                latestSourceMixedPresentScheduling.nontrivialRecoveries),
            static_cast<unsigned long long>(
                mixed.exactPresentSelectedRebuilds),
            static_cast<unsigned long long>(
                mixed.exactPresentDependencyTests),mixed.complete?1:0);
        std::printf("t2c_mixed_status edge=%llu/%llu/%llu/%llu/%llu/%llu "
                    "facet=%llu/%llu/%llu/%llu/%llu/%llu\n",
            static_cast<unsigned long long>(
                mixed.edgeAssessmentEventStatusCounts[0]),
            static_cast<unsigned long long>(
                mixed.edgeAssessmentEventStatusCounts[1]),
            static_cast<unsigned long long>(
                mixed.edgeAssessmentEventStatusCounts[2]),
            static_cast<unsigned long long>(
                mixed.edgeAssessmentEventStatusCounts[3]),
            static_cast<unsigned long long>(
                mixed.edgeAssessmentEventStatusCounts[4]),
            static_cast<unsigned long long>(
                mixed.edgeAssessmentEventStatusCounts[5]),
            static_cast<unsigned long long>(
                mixed.facetAssessmentEventStatusCounts[0]),
            static_cast<unsigned long long>(
                mixed.facetAssessmentEventStatusCounts[1]),
            static_cast<unsigned long long>(
                mixed.facetAssessmentEventStatusCounts[2]),
            static_cast<unsigned long long>(
                mixed.facetAssessmentEventStatusCounts[3]),
            static_cast<unsigned long long>(
                mixed.facetAssessmentEventStatusCounts[4]),
            static_cast<unsigned long long>(
                mixed.facetAssessmentEventStatusCounts[5]));
        const auto& assessed=mixed.assessmentBindingWork;
        const auto& executed=mixed.recoveryBindingWork;
        std::printf("t2c_mixed_binding prep_allocs=%llu "
                    "prep_logical_bytes=%llu "
                    "assess_edge=%llu/%llu/%llu "
                    "assess_facet=%llu/%llu data=%llu/%llu/%llu "
                    "assess_state=%llu/%llu "
                    "recover_edge=%llu/%llu/%llu "
                    "recover_facet=%llu/%llu data=%llu/%llu/%llu "
                    "recover_state=%llu/%llu\n",
            static_cast<unsigned long long>(
                mixed.modeledBatchPreparationAllocationRequests),
            static_cast<unsigned long long>(
                mixed.modeledBatchPreparationLogicalBytes),
            static_cast<unsigned long long>(assessed.edgeBindingAuditCalls),
            static_cast<unsigned long long>(assessed.edgeBindingRowsVisited),
            static_cast<unsigned long long>(assessed.edgeAliasRowsVisited),
            static_cast<unsigned long long>(assessed.facetBindingAuditCalls),
            static_cast<unsigned long long>(assessed.facetBindingRowsVisited),
            static_cast<unsigned long long>(assessed.edgeDataCopyRows),
            static_cast<unsigned long long>(assessed.faceDataCopyRows),
            static_cast<unsigned long long>(
                assessed.faceIncidenceRowsCopied),
            static_cast<unsigned long long>(assessed.stateCopies),
            static_cast<unsigned long long>(assessed.stateCopyRows),
            static_cast<unsigned long long>(executed.edgeBindingAuditCalls),
            static_cast<unsigned long long>(executed.edgeBindingRowsVisited),
            static_cast<unsigned long long>(executed.edgeAliasRowsVisited),
            static_cast<unsigned long long>(executed.facetBindingAuditCalls),
            static_cast<unsigned long long>(executed.facetBindingRowsVisited),
            static_cast<unsigned long long>(executed.edgeDataCopyRows),
            static_cast<unsigned long long>(executed.faceDataCopyRows),
            static_cast<unsigned long long>(
                executed.faceIncidenceRowsCopied),
            static_cast<unsigned long long>(executed.stateCopies),
            static_cast<unsigned long long>(executed.stateCopyRows));
    }
    if(!diagnostic)CHECK(outcome.published());
    CHECK(result!=nullptr);
    if(!diagnostic)CHECK(source!=nullptr);
    if(source) {
        CHECK(source->boundaryAudit.passed);
        CHECK(source->boundaryAudit.missingInterface==0);
        CHECK(source->boundaryAudit.unexpectedInterface==0);
        CHECK(source->recoveryMetrics.productionGlobalScanFallbacks==0);
        CHECK(source->recoveryMetrics.fullStageRescans==0);
        const auto& construction=source->constructionMetrics;
        CHECK(construction.partialImportRowsCopied>0);
        CHECK(construction.partialImportGlobalSlotsInitialized>0);
        CHECK(construction.partialImportCandidateMappingVisits>0);
        CHECK(construction.partialImportEdgeValidationRows>0);
        CHECK(construction.partialImportFaceValidationRows>0);
        CHECK(construction.partialImportProtectedEndpointVisits>0);
        CHECK(construction.partialImportAuthorityPointValidations>0);
        CHECK(construction.partialImportDuplicateKeyComparisons>0);
        CHECK(construction.partialImportResolutionNodeVisits>0);
        CHECK(construction.partialImportEdgeResolutionLookups>0);
        CHECK(construction.partialImportFaceResolutionLookups>0);
        CHECK(construction.partialImportEdgeDataReads>0);
        CHECK(construction.partialImportFaceDataCopies>0);
        CHECK(construction.partialImportFaceIncidenceRowsCopied>0);
        CHECK(construction.partialImportSortComparisons>0);
    }
    if(outcome.published())CHECK(mesh::semanticHashEx(destination).canonical);
    if(gFailures)return 1;
    return 0;
}
#endif
