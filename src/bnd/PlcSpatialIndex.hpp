// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/core/Types.hpp"
#include <algorithm>
#include <numeric>
#include <vector>
namespace paacs::bnd::detail {
class PlcGrowingAabbIndex {
public:
    struct Box {Vec3 lo,hi;};
private:
    struct Node{Box bounds;std::size_t begin=0,end=0,left=0,right=0;};
    std::vector<Node> nodes_;std::vector<std::size_t> order_;std::vector<std::size_t> roots_;std::size_t indexed_=0;
    static Box united(Box a,Box b){return {{std::min(a.lo.x,b.lo.x),std::min(a.lo.y,b.lo.y),std::min(a.lo.z,b.lo.z)},{std::max(a.hi.x,b.hi.x),std::max(a.hi.y,b.hi.y),std::max(a.hi.z,b.hi.z)}};}
    static bool overlaps(Box a,Box b){return a.lo.x<=b.hi.x&&b.lo.x<=a.hi.x&&a.lo.y<=b.hi.y&&b.lo.y<=a.hi.y&&a.lo.z<=b.hi.z&&b.lo.z<=a.hi.z;}
    template<class Bounds,class Charge>std::size_t build(std::size_t a,std::size_t b,Bounds bounds,Charge charge){charge(1);const auto i=nodes_.size();nodes_.push_back({});Box box=bounds(order_[a]);for(auto j=a+1;j<b;++j){charge(1);box=united(box,bounds(order_[j]));}nodes_[i]={box,a,b,0,0};if(b-a<=8)return i;const auto extent=box.hi-box.lo;const unsigned axis=extent.x>=extent.y&&extent.x>=extent.z?0:extent.y>=extent.z?1:2;const auto coordinate=[&](std::size_t id){auto p=bounds(id);return axis==0?p.lo.x+p.hi.x:axis==1?p.lo.y+p.hi.y:p.lo.z+p.hi.z;};const auto mid=(a+b)/2;std::nth_element(order_.begin()+a,order_.begin()+mid,order_.begin()+b,[&](auto x,auto y){charge(1);const auto p=coordinate(x),q=coordinate(y);return p!=q?p<q:x<y;});const auto left=build(a,mid,bounds,charge),right=build(mid,b,bounds,charge);nodes_[i].left=left;nodes_[i].right=right;return i;}
    template<class Bounds,class Charge,class Visitor>bool visit(std::size_t id,Box box,Bounds bounds,Charge charge,Visitor visitor)const{charge(1);const auto& n=nodes_[id];if(!overlaps(box,n.bounds))return true;if(n.left||n.right)return visit(n.left,box,bounds,charge,visitor)&&visit(n.right,box,bounds,charge,visitor);for(auto j=n.begin;j<n.end;++j){charge(1);const auto k=order_[j];if(overlaps(box,bounds(k))&&!visitor(k))return false;}return true;}
    template<class Bounds,class Charge,class Visitor>
    bool visitAll(std::size_t count,Box box,Bounds bounds,Charge charge,Visitor visitor)const{
        for(auto root:roots_)if(!visit(root,box,bounds,charge,visitor))return false;
        for(auto i=indexed_;i<count;++i){
            charge(1);
            if(overlaps(box,bounds(i))&&!visitor(i))return false;
        }
        return true;
    }
public:
    template<class Bounds,class Charge>
    void rebuild(std::size_t count,Bounds bounds,Charge charge){
        indexed_=0;
        roots_.clear();
        nodes_.clear();
        order_.resize(count);
        std::iota(order_.begin(),order_.end(),0);
        nodes_.reserve(count*2);
        if(count)roots_.push_back(build(0,count,bounds,charge));
        indexed_=count;
    }
    template<class Bounds,class Charge,class Visitor>
    bool all(std::size_t count,Box box,Bounds bounds,Charge charge,Visitor visitor,std::size_t rebuildInterval=256){
        if(count<indexed_||count-indexed_>=rebuildInterval)rebuild(count,bounds,charge);
        return visitAll(count,box,bounds,charge,visitor);
    }
    // The live tetrahedron IDs are append-only. Index short batches without
    // rebuilding the old tree, then consolidate after enough batches accrue.
    template<class Bounds,class Charge,class Visitor>
    bool allAppending(std::size_t count,Box box,Bounds bounds,Charge charge,Visitor visitor){
        constexpr std::size_t batch=128, maximumRoots=32;
        if(count<indexed_||(count-indexed_>=batch&&(roots_.empty()||roots_.size()>=maximumRoots))){
            rebuild(count,bounds,charge);
        }else if(count-indexed_>=batch){
            const auto previousNodes=nodes_.size();
            order_.resize(count);
            std::iota(order_.begin()+indexed_,order_.end(),indexed_);
            try{
                roots_.push_back(build(indexed_,count,bounds,charge));
                indexed_=count;
            }catch(...){
                nodes_.resize(previousNodes);
                order_.resize(indexed_);
                throw;
            }
        }
        return visitAll(count,box,bounds,charge,visitor);
    }
};
}
