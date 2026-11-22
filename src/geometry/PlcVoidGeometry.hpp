// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Private source-PLC queries for initial point placement. Numerical clearance
// ranks proposals; exact ray parity separately proves strict source inclusion.
#pragma once

#include "../bnd/TetPlcGiftWrapping.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <utility>

namespace paacs::geo::detail {

// Numerical queue metric only. Strict source inclusion below uses exact signs;
// it never infers inclusion from a positive distance or a native CAD sample.
inline double pointTriangleDistanceSquared(Vec3 p, Vec3 a, Vec3 b, Vec3 c) {
    using L=long double;
    using V=std::array<L,3>;
    const auto sub=[](Vec3 x,Vec3 y){return V{L(x.x)-L(y.x),L(x.y)-L(y.y),L(x.z)-L(y.z)};};
    const auto dot=[](V x,V y){return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];};
    auto u=sub(b,a),v=sub(c,a),w=sub(p,a);
    L scale=0;for(auto x:{u,v,w})for(auto y:x)scale=std::max(scale,std::abs(y));
    if(!(scale>0)||!std::isfinite(scale))return 0.;
    for(unsigned i=0;i<3;++i){u[i]/=scale;v[i]/=scale;w[i]/=scale;}
    const auto finish=[&](V q){
        V d{};for(unsigned i=0;i<3;++i)d[i]=w[i]-q[i];
        const L squared=std::max(L(0),dot(d,d))*scale*scale;
        return double(squared);
    };
    const V zero{};
    const L d1=dot(u,w),d2=dot(v,w);
    if(d1<=0&&d2<=0)return finish(zero);
    V wb{};for(unsigned i=0;i<3;++i)wb[i]=w[i]-u[i];
    const L d3=dot(u,wb),d4=dot(v,wb);
    if(d3>=0&&d4<=d3)return finish(u);
    const L vc=d1*d4-d3*d2;
    if(vc<=0&&d1>=0&&d3<=0){
        const L t=d1/(d1-d3);V q{};for(unsigned i=0;i<3;++i)q[i]=t*u[i];return finish(q);
    }
    V wc{};for(unsigned i=0;i<3;++i)wc[i]=w[i]-v[i];
    const L d5=dot(u,wc),d6=dot(v,wc);
    if(d6>=0&&d5<=d6)return finish(v);
    const L vb=d5*d2-d1*d6;
    if(vb<=0&&d2>=0&&d6<=0){
        const L t=d2/(d2-d6);V q{};for(unsigned i=0;i<3;++i)q[i]=t*v[i];return finish(q);
    }
    const L va=d3*d6-d5*d4;
    if(va<=0&&(d4-d3)>=0&&(d5-d6)>=0){
        const L t=(d4-d3)/((d4-d3)+(d5-d6));V q{};
        for(unsigned i=0;i<3;++i)q[i]=u[i]+t*(v[i]-u[i]);return finish(q);
    }
    const L denominator=va+vb+vc;
    if(!(denominator>0)||!std::isfinite(denominator))return 0.;
    const L s=vb/denominator,t=vc/denominator;V q{};
    for(unsigned i=0;i<3;++i)q[i]=s*u[i]+t*v[i];return finish(q);
}

struct PlcVoidGeometry {
    using Front=bnd::detail::PlcGiftWrapping;
    using Box=Front::Box;
    enum class Inclusion {Inside,Outside,Ambiguous};
    Front& front;
    std::function<void(std::uint64_t)> charge;
    Box bounds{};

    explicit PlcVoidGeometry(Front& source,std::function<void(std::uint64_t)> account={})
        :front(source),charge(std::move(account)) {
        tick();
        if(!front.initialized||front.tree.empty()||front.boundary.empty())return;
        bounds=front.tree.front().bounds;
        valid_=finite(bounds.lo)&&finite(bounds.hi)&&bounds.lo.x<=bounds.hi.x&&
            bounds.lo.y<=bounds.hi.y&&bounds.lo.z<=bounds.hi.z;
    }

    // Clamped numerical nearest-triangle distance, used only to rank voids.
    double clearance(Vec3 p,double upperBound)const {
        tick();
        if(!valid_||!finite(p)||!(upperBound>0.))return 0.;
        double best=upperBound;
        std::array<unsigned,128> stack{};std::size_t count=1,visits=0;stack[0]=0;
        while(count){
            const auto id=stack[--count];tick();
            if(id>=front.tree.size()||++visits>front.tree.size())return 0.;
            const auto& node=front.tree[id];
            if(boxDistanceLower(p,node.bounds)>=best)continue;
            if(node.left||node.right){
                if(!node.left||!node.right||node.left>=front.tree.size()||node.right>=front.tree.size()||count+2>stack.size())return 0.;
                auto near=node.left,far=node.right;
                if(boxDistanceLower(p,front.tree[far].bounds)<boxDistanceLower(p,front.tree[near].bounds))std::swap(near,far);
                stack[count++]=far;stack[count++]=near;
            }else{
                if(node.first>node.last||node.last>front.order.size())return 0.;
                for(unsigned j=node.first;j<node.last;++j){
                    tick(16);const auto row=front.order[j];
                    if(row>=front.boundary.size())return 0.;
                    const auto f=front.boundary[row];
                    const auto d2=pointTriangleDistanceSquared(p,front.vertices[f[0]].position,
                        front.vertices[f[1]].position,front.vertices[f[2]].position);
                    if(!(d2>=0.))return 0.;
                    best=std::min(best,std::sqrt(d2));
                    if(best==0.)return 0.;
                }
            }
        }
        return best;
    }

    Inclusion inside(Vec3 p)const {
        tick();
        if(!valid_||!finite(p))return Inclusion::Ambiguous;
        if(outside(p))return Inclusion::Outside;
        const auto span=std::max({bounds.hi.x-bounds.lo.x,bounds.hi.y-bounds.lo.y,bounds.hi.z-bounds.lo.z});
        if(!(span>0.)||!std::isfinite(span))return Inclusion::Ambiguous;
        // Four unrelated directions; endpoints are verified strictly outside
        // the complete source AABB. An ambiguous edge/vertex/coplanar ray is
        // retried, never assigned a parity by a floating tolerance.
        const std::array<Vec3,4> ends{{
            {bounds.hi.x+span,std::lerp(bounds.lo.y,bounds.hi.y,3./11.),std::lerp(bounds.lo.z,bounds.hi.z,7./13.)},
            {bounds.lo.x-span,std::lerp(bounds.lo.y,bounds.hi.y,8./17.),std::lerp(bounds.lo.z,bounds.hi.z,5./19.)},
            {std::lerp(bounds.lo.x,bounds.hi.x,7./23.),bounds.hi.y+span,std::lerp(bounds.lo.z,bounds.hi.z,11./29.)},
            {std::lerp(bounds.lo.x,bounds.hi.x,13./31.),std::lerp(bounds.lo.y,bounds.hi.y,17./37.),bounds.hi.z+span}
        }};
        for(const auto q:ends){
            tick();if(!finite(q)||!outside(q)){continue;}
            const auto result=ray(p,q);
            if(result!=Inclusion::Ambiguous)return result;
            
        }
        return Inclusion::Ambiguous;
    }

private:
    bool valid_=false;
    void tick(std::uint64_t n=1)const{if(charge)charge(n);}
    static bool finite(Vec3 p){return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);}
    bool outside(Vec3 p)const{return p.x<bounds.lo.x||p.x>bounds.hi.x||p.y<bounds.lo.y||p.y>bounds.hi.y||p.z<bounds.lo.z||p.z>bounds.hi.z;}
    static double coordinate(Vec3 p,unsigned i){return i==0?p.x:i==1?p.y:p.z;}
    // L-infinity distance is a conservative lower bound on Euclidean distance.
    static double boxDistanceLower(Vec3 p,Box b){
        double result=0.;for(unsigned i=0;i<3;++i){const auto x=coordinate(p,i),lo=coordinate(b.lo,i),hi=coordinate(b.hi,i);
            if(x<lo)result=std::max(result,std::nextafter(lo-x,0.));
            if(x>hi)result=std::max(result,std::nextafter(x-hi,0.));}
        return result;
    }
    struct Interval{double lo,hi;};
    static Interval quotient(double numerator,double denominator){
        const auto infinity=std::numeric_limits<double>::infinity();
        if(!std::isfinite(numerator)||!std::isfinite(denominator))return {-infinity,infinity};
        const Interval n{std::nextafter(numerator,-infinity),std::nextafter(numerator,infinity)};
        const Interval d{std::nextafter(denominator,-infinity),std::nextafter(denominator,infinity)};
        // Finite extrema can expand to infinity. Retain the BVH branch rather
        // than allow an inf/inf NaN to make a false-negative slab decision.
        if(!std::isfinite(n.lo)||!std::isfinite(n.hi)||!std::isfinite(d.lo)||!std::isfinite(d.hi)||
            (d.lo<=0&&d.hi>=0))return {-infinity,infinity};
        const std::array values{n.lo/d.lo,n.lo/d.hi,n.hi/d.lo,n.hi/d.hi};
        const auto [lo,hi]=std::minmax_element(values.begin(),values.end());
        return {std::nextafter(*lo,-infinity),std::nextafter(*hi,infinity)};
    }
    static bool rayBox(Vec3 p,Vec3 q,Box b){
        double enter=0.,leave=1.;
        for(unsigned i=0;i<3;++i){
            const auto x=coordinate(p,i),y=coordinate(q,i),lo=coordinate(b.lo,i),hi=coordinate(b.hi,i);
            if(x==y){if(x<lo||x>hi)return false;continue;}
            const auto a=quotient(lo-x,y-x),c=quotient(hi-x,y-x);
            enter=std::max(enter,std::min(a.lo,c.lo));leave=std::min(leave,std::max(a.hi,c.hi));
            if(enter>leave)return false;
        }return true;
    }
    Inclusion ray(Vec3 p,Vec3 q)const {
        using predicates::Sign;
        bool odd=false;
        std::array<unsigned,128> stack{};std::size_t count=1,visits=0;stack[0]=0;
        while(count){
            const auto id=stack[--count];tick();
            if(id>=front.tree.size()||++visits>front.tree.size())return Inclusion::Ambiguous;
            const auto& node=front.tree[id];if(!rayBox(p,q,node.bounds))continue;
            if(node.left||node.right){
                if(!node.left||!node.right||node.left>=front.tree.size()||node.right>=front.tree.size()||count+2>stack.size())return Inclusion::Ambiguous;
                stack[count++]=node.left;stack[count++]=node.right;
            }else{
                if(node.first>node.last||node.last>front.order.size())return Inclusion::Ambiguous;
                for(unsigned j=node.first;j<node.last;++j){
                    tick(6);const auto row=front.order[j];
                    if(row>=front.boundary.size())return Inclusion::Ambiguous;
                    const auto f=front.boundary[row];
                    const auto a=front.vertices[f[0]].position,b=front.vertices[f[1]].position,c=front.vertices[f[2]].position;
                    const auto sp=predicates::orient3d(a,b,c,p),sq=predicates::orient3d(a,b,c,q);
                    if(sp==Sign::Zero){
                        if(sq==Sign::Zero)return Inclusion::Ambiguous;
                        const auto projection=bnd::tet_facet_detail::projectionForTriangle(a,b,c);
                        if(projection==bnd::tet_facet_detail::Projection::Invalid)return Inclusion::Ambiguous;
                        if(bnd::tet_facet_detail::pointInTriangle2(bnd::tet_facet_detail::project(p,projection),
                            bnd::tet_facet_detail::project(a,projection),bnd::tet_facet_detail::project(b,projection),
                            bnd::tet_facet_detail::project(c,projection))!=bnd::tet_facet_detail::PointTri2::Outside)return Inclusion::Ambiguous;
                        continue; // Start is only on the supporting-plane extension.
                    }
                    if(sq==Sign::Zero||sp==sq)continue; // Outside endpoint cannot be on a source triangle.
                    const std::array edges{predicates::orient3d(p,q,a,b),predicates::orient3d(p,q,b,c),predicates::orient3d(p,q,c,a)};
                    bool positive=false,negative=false,zero=false;
                    for(const auto sign:edges){positive|=sign==Sign::Positive;negative|=sign==Sign::Negative;zero|=sign==Sign::Zero;}
                    if(positive&&negative)continue;
                    if(zero)return Inclusion::Ambiguous;
                    odd=!odd;
                }
            }
        }
        return odd?Inclusion::Inside:Inclusion::Outside;
    }
};

} // namespace paacs::geo::detail
