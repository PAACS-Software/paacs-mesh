// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <map>
#include <vector>
namespace paacs::geo::detail {
// Cancel only an exactly opposite pair on the same original CAD face.
// Each removal preserves the oriented boundary and signed volume exactly.
template<class Charge,class Eligible>
std::size_t cancelOcctOppositeFacePairs(std::vector<std::array<int,3>>& triangles,
    std::vector<int>& faceTags,Charge charge,Eligible eligible){
    using Key=std::array<int,4>;constexpr auto none=(std::numeric_limits<std::size_t>::max)();
    struct Heads {std::array<std::size_t,2> value{(std::numeric_limits<std::size_t>::max)(),(std::numeric_limits<std::size_t>::max)()};};
    const auto less=[&](const Key& a,const Key& b){charge(1);return a<b;};
    std::map<Key,Heads,decltype(less)> pending(less);
    std::vector<std::size_t> next(triangles.size(),none);std::vector<bool> removed(triangles.size());
    std::size_t pairs=0;
    for(std::size_t i=0;i<triangles.size();++i){
        charge(1);if(!eligible(faceTags[i]))continue;const auto t=triangles[i];auto sorted=t;std::sort(sorted.begin(),sorted.end());
        const unsigned parity=unsigned((t[0]>t[1])+(t[0]>t[2])+(t[1]>t[2]))&1;
        auto& heads=pending[{faceTags[i],sorted[0],sorted[1],sorted[2]}].value;
        auto& opposite=heads[1-parity];
        if(opposite!=none){removed[i]=removed[opposite]=true;opposite=next[opposite];++pairs;}
        else{next[i]=heads[parity];heads[parity]=i;}
    }
    if(!pairs)return 0;
    std::size_t kept=0;for(std::size_t i=0;i<triangles.size();++i){charge(1);if(!removed[i]){triangles[kept]=triangles[i];faceTags[kept++]=faceTags[i];}}
    triangles.resize(kept);faceTags.resize(kept);return pairs;
}
template<class Charge>
std::size_t cancelOcctOppositeFacePairs(std::vector<std::array<int,3>>& triangles,
    std::vector<int>& faceTags,Charge charge){
    return cancelOcctOppositeFacePairs(triangles,faceTags,charge,[](int){return true;});
}
}
