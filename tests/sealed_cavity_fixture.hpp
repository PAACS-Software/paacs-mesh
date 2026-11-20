// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Shared-shape style fixture used by the T3-A retained-boundary proof.  This
// intentionally supplies no GSolid: the authoritative ownership provider is
// the region authority for the outer shell plus sealed void component.

#include "paacs/geometry/Kernel.hpp"

#include <algorithm>
#include <array>
#include <map>
#include <memory>
#include <vector>

namespace paacs_tests::fixtures {

class ShellEdge final : public paacs::geo::GEdge {
public:
    ShellEdge(int tag,const paacs::geo::GVertex* a,
              const paacs::geo::GVertex* b)
        : GEdge(tag),a_(a),b_(b) {}
    paacs::Vec3 pointAt(double t) const override {
        return a_->position()+(b_->position()-a_->position())*t;
    }
    paacs::Vec3 tangentAt(double) const override {
        return (b_->position()-a_->position()).normalized();
    }
    double length() const override {
        return (b_->position()-a_->position()).norm();
    }
    int startVertexTag() const override { return a_->tag(); }
    int endVertexTag() const override { return b_->tag(); }
private:
    const paacs::geo::GVertex* a_=nullptr;
    const paacs::geo::GVertex* b_=nullptr;
};

class ShellFace final : public paacs::geo::GFace {
public:
    ShellFace(int tag,paacs::Vec3 origin,paacs::Vec3 uAxis,
              paacs::Vec3 vAxis,std::vector<int> loop)
        : GFace(tag),origin_(origin),uAxis_(uAxis),vAxis_(vAxis),
          loop_(std::move(loop)) {}
    paacs::Vec3 pointAt(const paacs::Vec2& uv) const override {
        return origin_+uAxis_*uv.u+vAxis_*uv.v;
    }
    paacs::Vec2 paramAt(const paacs::Vec3& point) const override {
        const paacs::Vec3 d=point-origin_;
        return {d.dot(uAxis_)/uAxis_.norm2(),
                d.dot(vAxis_)/vAxis_.norm2()};
    }
    paacs::Vec3 normalAt(const paacs::Vec2&) const override {
        return uAxis_.cross(vAxis_).normalized();
    }
    void paramBounds(paacs::Vec2& lo,paacs::Vec2& hi) const override {
        lo={0,0}; hi={1,1};
    }
    std::vector<int> boundaryEdgeTags() const override { return loop_; }
    double curvatureAt(const paacs::Vec2&) const override { return 0.0; }
    paacs::Vec2 metricScale() const override {
        return {uAxis_.norm(),vAxis_.norm()};
    }
private:
    paacs::Vec3 origin_{},uAxis_{},vAxis_{};
    std::vector<int> loop_;
};

class SealedCavityModel final : public paacs::geo::Model {
public:
    SealedCavityModel() {
        addShell({-2,-2,-2},{2,2,2},false);
        addShell({-0.5,-0.5,-0.5},{0.5,0.5,0.5},true);
    }
    std::vector<const paacs::geo::GVertex*> vertices() const override {
        std::vector<const paacs::geo::GVertex*> result;
        for (const auto& value:vertices_) result.push_back(value.get());
        return result;
    }
    std::vector<const paacs::geo::GEdge*> edges() const override {
        std::vector<const paacs::geo::GEdge*> result;
        for (const auto& value:edges_) result.push_back(value.get());
        return result;
    }
    std::vector<const paacs::geo::GFace*> faces() const override {
        std::vector<const paacs::geo::GFace*> result;
        for (const auto& value:faces_) result.push_back(value.get());
        return result;
    }
    std::vector<const paacs::geo::GSolid*> solids() const override { return {}; }
    const paacs::geo::GVertex* vertex(int tag) const override {
        for (const auto& value:vertices_) if (value->tag()==tag) return value.get();
        return nullptr;
    }
    const paacs::geo::GEdge* edge(int tag) const override {
        for (const auto& value:edges_) if (value->tag()==tag) return value.get();
        return nullptr;
    }
    const paacs::geo::GFace* face(int tag) const override {
        for (const auto& value:faces_) if (value->tag()==tag) return value.get();
        return nullptr;
    }
    const paacs::geo::GSolid* solid(int) const override { return nullptr; }
private:
    void addShell(const paacs::Vec3& lo,const paacs::Vec3& hi,bool inward) {
        const std::size_t first=vertices_.size();
        for (std::size_t i=0;i<8;++i)
            vertices_.push_back(std::make_unique<paacs::geo::GVertex>(
                static_cast<int>(vertices_.size()+1),
                paacs::Vec3{(i&1U)?hi.x:lo.x,(i&2U)?hi.y:lo.y,
                            (i&4U)?hi.z:lo.z}));
        const std::array<std::array<std::size_t,4>,6> outward{{
            {{0,4,6,2}},{{1,3,7,5}},{{0,1,5,4}},
            {{2,6,7,3}},{{0,2,3,1}},{{4,5,7,6}}}};
        std::map<std::array<std::size_t,2>,int> edgeTags;
        for (auto loop:outward) {
            if (inward) std::reverse(loop.begin()+1,loop.end());
            std::vector<int> faceEdges;
            for (std::size_t i=0;i<4;++i) {
                const std::size_t a=first+loop[i];
                const std::size_t b=first+loop[(i+1)%4];
                auto key=std::array<std::size_t,2>{a,b};
                if (key[1]<key[0]) std::swap(key[0],key[1]);
                auto found=edgeTags.find(key);
                if (found==edgeTags.end()) {
                    const int tag=static_cast<int>(edges_.size()+1);
                    edges_.push_back(std::make_unique<ShellEdge>(
                        tag,vertices_[key[0]].get(),vertices_[key[1]].get()));
                    found=edgeTags.emplace(key,tag).first;
                }
                faceEdges.push_back(found->second);
            }
            const paacs::Vec3 p0=vertices_[first+loop[0]]->position();
            const paacs::Vec3 p1=vertices_[first+loop[1]]->position();
            const paacs::Vec3 p3=vertices_[first+loop[3]]->position();
            faces_.push_back(std::make_unique<ShellFace>(
                static_cast<int>(faces_.size()+1),p0,p1-p0,p3-p0,
                std::move(faceEdges)));
        }
    }
    std::vector<std::unique_ptr<paacs::geo::GVertex>> vertices_;
    std::vector<std::unique_ptr<ShellEdge>> edges_;
    std::vector<std::unique_ptr<ShellFace>> faces_;
};

} // namespace paacs_tests::fixtures
