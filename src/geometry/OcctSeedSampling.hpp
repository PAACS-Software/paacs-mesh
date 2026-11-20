// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// This sample proposes optional interior seeds only. The mesher retains every
// boundary vertex and facet in its exact geometric admission checks.
#include "paacs/core/Types.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <span>
#include <utility>
#include <vector>
namespace paacs::geo::detail {
template<class Charge>
std::vector<Index> occtSeedBoundarySample(std::span<const Vec3> points,Vec3 lower,Vec3 upper,Charge charge){
    constexpr std::size_t fullLimit=8192,sampleLimit=2048;
    std::vector<Index> selected;
    if(points.size()<=fullLimit){selected.resize(points.size());std::iota(selected.begin(),selected.end(),Index(0));return selected;}
    std::vector<std::pair<std::uint64_t,Index>> spatial;spatial.reserve(points.size());
    const std::array lo{lower.x,lower.y,lower.z},hi{upper.x,upper.y,upper.z};
    for(std::size_t i=0;i<points.size();++i){
        charge(1);const std::array p{points[i].x,points[i].y,points[i].z};std::array<std::uint32_t,3> q{};
        for(unsigned axis=0;axis<3;++axis)if(hi[axis]>lo[axis])q[axis]=std::uint32_t(std::clamp((p[axis]-lo[axis])/(hi[axis]-lo[axis]),0.,1.)*2097151.);
        std::uint64_t code=0;for(unsigned bit=0;bit<21;++bit)for(unsigned axis=0;axis<3;++axis)code|=std::uint64_t((q[axis]>>bit)&1U)<<(3*bit+axis);
        spatial.emplace_back(code,Index(i));
    }
    std::sort(spatial.begin(),spatial.end(),[&](const auto& a,const auto& b){charge(1);return a<b;});
    selected.reserve(sampleLimit);for(std::size_t i=0;i<sampleLimit;++i){charge(1);selected.push_back(spatial[i*(spatial.size()-1)/(sampleLimit-1)].second);}
    return selected;
}
}
