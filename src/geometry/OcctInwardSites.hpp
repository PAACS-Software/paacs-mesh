// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/mesh/Mesh.hpp"
#include <algorithm>
#include <array>
#include <map>
#include <span>
#include <vector>

namespace paacs::geo::detail {

// Collect interior support for the CAD faces next to rejected cells.
// Keep the complete original shell and append only interior sites.
template<class AdmitInterior, class AdmitCentroid, class Charge>
std::vector<Vec3> occtInwardSites(const mesh::Mesh& boundary, const mesh::Mesh& filled,
    std::span<const unsigned char> admitted, double chord, std::size_t maximumSites,
    AdmitInterior admitInterior, AdmitCentroid admitCentroid, Charge charge) {
    using Face=std::array<Index,3>;
    const auto key=[](Face face){std::ranges::sort(face);return face;};
    std::map<Face,Face> source;
    for(Index i=0;i<boundary.cellCount();++i){
        charge(1);
        const auto* ids=boundary.cellNodes(i);
        const Face face{ids[0],ids[1],ids[2]};
        source.emplace(key(face),face);
    }
    constexpr unsigned maximumProposals=384;
    unsigned proposals=0;
    std::vector<Vec3> sites;
    sites.reserve(maximumSites);
    std::size_t tet=0;
    for(Index i=0;i<filled.cellCount()&&sites.size()<maximumSites&&proposals<maximumProposals;++i){
        if(filled.cell(i).type!=mesh::CellType::Tet4)continue;
        charge(1);
        if(admitted[tet++])continue;
        const auto* ids=filled.cellNodes(i);
        for(unsigned side=0;side<4&&sites.size()<maximumSites&&proposals<maximumProposals;++side){
            const auto found=source.find(key({ids[(side+1)%4],ids[(side+2)%4],ids[(side+3)%4]}));
            if(found==source.end())continue;
            const auto [a,b,c]=found->second;
            source.erase(found);
            const auto p=boundary.node(a),q=boundary.node(b),r=boundary.node(c);
            const auto center=(p+q+r)/3.;
            const auto normal=(q-p).cross(r-p);
            const auto length=normal.norm();
            if(!(length>0.)||!std::isfinite(length))continue;
            const auto inward=normal*(-1./length);
            const auto edge=std::max({(q-p).norm(),(r-p).norm(),(r-q).norm()});
            for(const auto factor:{.4,.24,.65,.12,.08}){
                charge(sites.size()+1);
                const auto site=center+inward*std::clamp(factor*edge,.18*chord,10.*chord);
                const auto separation=std::max(.18*chord,.18*edge);
                if(std::ranges::any_of(sites,[&](Vec3 old){return (old-site).norm()<separation;}))continue;
                if(proposals==maximumProposals)break;
                ++proposals;
                if(admitInterior(site)&&admitCentroid((p+q+r+site)*.25)){
                    sites.push_back(site);
                    break;
                }
            }
        }
    }
    return sites;
}

}
