// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Compact, allocation-stable copy of the real contiguous-layer geometry used
// by the T2-C facet-conflict production-route allocation probe.  The provider
// keeps its complete call trace in caller-owned fixed storage so an injected
// allocation can never be confused with probe bookkeeping.

#include "paacs/algo/ProtectedFeatureRefinement.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/sizing/SizeField.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <span>
#include <utility>
#include <vector>

namespace paacs::tests::t2c_facet_probe {

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
    int startVertexTag() const override { return a_->tag(); }
    int endVertexTag() const override { return b_->tag(); }
private:
    const geo::GVertex* a_=nullptr;
    const geo::GVertex* b_=nullptr;
};

class InterfaceFace final : public geo::GFace {
public:
    InterfaceFace(int tag,Vec3 origin,Vec3 u,Vec3 v,std::vector<int> loop)
        : GFace(tag),origin_(origin),u_(u),v_(v),loop_(std::move(loop)) {}
    Vec3 pointAt(const Vec2& p) const override {
        return origin_+u_*p.u+v_*p.v;
    }
    Vec2 paramAt(const Vec3& p) const override {
        const auto d=p-origin_;
        const double uu=u_.norm2(),uv=u_.dot(v_),vv=v_.norm2();
        const double determinant=uu*vv-uv*uv;
        return {(d.dot(u_)*vv-d.dot(v_)*uv)/determinant,
                (d.dot(v_)*uu-d.dot(u_)*uv)/determinant};
    }
    Vec3 normalAt(const Vec2&) const override {
        return u_.cross(v_).normalized();
    }
    void paramBounds(Vec2& lo,Vec2& hi) const override {
        lo={0,0}; hi={1,1};
    }
    std::vector<int> boundaryEdgeTags() const override { return loop_; }
    double curvatureAt(const Vec2&) const override { return 0.0; }
    Vec2 metricScale() const override { return {u_.norm(),v_.norm()}; }
private:
    Vec3 origin_{},u_{},v_{};
    std::vector<int> loop_;
};

class LayeredModel final : public geo::Model {
public:
    explicit LayeredModel(int regionCount=8,double halfExtent=4.0,
                          double layerSpacing=1.0)
        : regionCount_(regionCount) {
        for(int x=0;x<=regionCount_;++x)
            for(int z=0;z<2;++z)
                for(int y=0;y<2;++y) {
                    const Vec3 point{
                        layerSpacing*static_cast<double>(x),
                        halfExtent*static_cast<double>(2*y-1),
                        halfExtent*static_cast<double>(2*z-1)};
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

    const std::vector<bnd::VolumeRegionUse>& owners(int faceTag) const {
        return owners_.at(static_cast<std::size_t>(faceTag-1));
    }
    std::vector<const geo::GVertex*> vertices() const override {
        std::vector<const geo::GVertex*> out;
        out.reserve(vertices_.size());
        for(const auto& row:vertices_)out.push_back(row.get());
        return out;
    }
    std::vector<const geo::GEdge*> edges() const override {
        std::vector<const geo::GEdge*> out;
        out.reserve(edges_.size());
        for(const auto& row:edges_)out.push_back(row.get());
        return out;
    }
    std::vector<const geo::GFace*> faces() const override {
        std::vector<const geo::GFace*> out;
        out.reserve(faces_.size());
        for(const auto& row:faces_)out.push_back(row.get());
        return out;
    }
    std::vector<const geo::GSolid*> solids() const override { return {}; }
    const geo::GVertex* vertex(int tag) const override {
        for(const auto& row:vertices_)if(row->tag()==tag)return row.get();
        return nullptr;
    }
    const geo::GEdge* edge(int tag) const override {
        for(const auto& row:edges_)if(row->tag()==tag)return row.get();
        return nullptr;
    }
    const geo::GFace* face(int tag) const override {
        for(const auto& row:faces_)if(row->tag()==tag)return row.get();
        return nullptr;
    }
    const geo::GSolid* solid(int) const override { return nullptr; }

private:
    std::size_t v(int x,int y,int z) const noexcept {
        return static_cast<std::size_t>(x*4+z*2+y);
    }
    void addFace(std::array<std::size_t,4> loop,
                 std::vector<bnd::VolumeRegionUse> owners) {
        std::vector<int> faceEdges;
        faceEdges.reserve(4);
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

    std::vector<std::unique_ptr<geo::GVertex>> vertices_;
    std::vector<std::unique_ptr<InterfaceEdge>> edges_;
    std::vector<std::unique_ptr<InterfaceFace>> faces_;
    std::vector<std::vector<bnd::VolumeRegionUse>> owners_;
    std::map<std::array<std::size_t,2>,int> edgeTags_;
    int regionCount_=8;
};

class LayerField final : public sizing::SizeField {
public:
    double eval(const Vec3& point) const override {
        return 4.0*(std::abs(point.x-1.0)<=0.125?0.5:0.75);
    }
};

class Ownership final : public algo::VolumeBoundaryOwnershipProvider {
public:
    algo::VolumeBoundaryProviderStamp current{15001,1};
    std::map<algo::FacePersistentId,std::vector<bnd::VolumeRegionUse>> rows;
    algo::VolumeBoundaryProviderStamp stamp() const noexcept override {
        return current;
    }
    algo::VolumeFaceOwnership ownership(
            algo::FacePersistentId face) const override {
        algo::VolumeFaceOwnership out;
        out.providerStamp=current;
        out.facePersistentId=face;
        const auto found=rows.find(face);
        if(found!=rows.end())out.owners=found->second;
        return out;
    }
    algo::ProviderRecoveryEnvelope recoveryEnvelope(
            std::span<const Vec3>) const override {
        algo::ProviderRecoveryEnvelope out;
        out.providerStamp=current;
        out.availability=
            algo::VolumeProviderAvailability::RecoveryEnvelopeUnavailable;
        return out;
    }
};

enum class ProviderCallKind : std::uint8_t { Stamp, Placement };
struct ProviderCall {
    ProviderCallKind kind=ProviderCallKind::Stamp;
    algo::ProtectedFeatureEntityKind entityKind=
        algo::ProtectedFeatureEntityKind::Curve;
    algo::ArtifactStamp entityStamp{};
    algo::SurfaceNodeIdentity identity{};
    friend bool operator==(const ProviderCall&,const ProviderCall&)=default;
};

class Placement final : public algo::ProtectedFeaturePlacementProvider {
public:
    algo::ProtectedFeatureProviderStamp current{15002,1};
    mutable std::array<ProviderCall,256> calls{};
    mutable std::size_t callsUsed=0;
    mutable bool overflow=false;
    bool* providerStarted=nullptr;

    algo::ProtectedFeatureProviderStamp stamp() const noexcept override {
        note({ProviderCallKind::Stamp});
        return current;
    }
    algo::ProtectedFeaturePlacement placeOnEdge(
            const geo::GEdge& edge,
            const algo::ProtectedFeaturePlacementRequest& request) const
            override {
        note({ProviderCallKind::Placement,request.kind,
              request.entityStamp,request.requestedIdentity});
        if(providerStarted)*providerStarted=true;
        return {current,request.entityStamp,request.requestedIdentity,
            algo::ProtectedFeaturePlacementStatus::Placed,
            edge.pointAt(request.curveParameter)};
    }
    algo::ProtectedFeaturePlacement placeOnFace(
            const geo::GFace& face,
            const algo::ProtectedFeaturePlacementRequest& request) const
            override {
        note({ProviderCallKind::Placement,request.kind,
              request.entityStamp,request.requestedIdentity});
        if(providerStarted)*providerStarted=true;
        return {current,request.entityStamp,request.requestedIdentity,
            algo::ProtectedFeaturePlacementStatus::Placed,
            face.pointAt(request.faceParameter)};
    }
    std::span<const ProviderCall> trace() const noexcept {
        return {calls.data(),callsUsed};
    }
private:
    void note(const ProviderCall& row) const noexcept {
        if(callsUsed<calls.size())calls[callsUsed++]=row;
        else overflow=true;
    }
};

inline algo::BoundaryConformingV4Request requestFor(
        const Ownership& ownership,const Placement& placement) {
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
}

inline void bindOwnership(const LayeredModel& model,Ownership& ownership) {
    for(const auto* face:model.faces())
        ownership.rows[
            static_cast<algo::FacePersistentId>(face->tag())+1]=
            model.owners(face->tag());
}

} // namespace paacs::tests::t2c_facet_probe
