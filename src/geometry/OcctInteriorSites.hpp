// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/mesh/Mesh.hpp"
#include "paacs/sizing/RegionSizeField.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace paacs::geo::detail {

struct OcctInteriorSiteReport {
    std::size_t boxCandidates=0, faceCandidates=0;
    std::size_t boxSites=0, faceSites=0;
    std::size_t spacingRejects=0, nativeRejects=0, invalidCandidates=0;
};

// Supply the first constrained Delaunay construction with interior spacing
// sites. The original nodes, classifications and oriented Tri3 cells stay
// intact. These are proposals only: strict native interior admission is
// required for every appended node, and all normal recovery/audits still run.
// SizeField is OcctSpatialSizing (or an equivalent bounded size-field view).
template<class SizeField, class AdmitNode, class Charge>
mesh::Mesh appendOcctInteriorSites(const mesh::Mesh& boundary, const SizeField& field,
    mesh::Classification region, std::size_t maximumSites, AdmitNode admitNode,
    Charge charge, OcctInteriorSiteReport* report=nullptr, bool sampleControls=true) {
    mesh::Mesh result=boundary;
    if (report) *report={};
    if (!maximumSites || !boundary.nodeCount() || region.dim!=3) return result;
    maximumSites=std::min<std::size_t>(maximumSites,1024);
    const auto baseSize=field.boxes().baseSize();
    const auto bucketWidth=.65*baseSize;
    if (!(bucketWidth>0.) || !std::isfinite(bucketWidth)) return result;
    OcctInteriorSiteReport counts;
    const auto finite=[](Vec3 p) {
        return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);
    };
    const auto sizeAt=[&](Vec3 p, int face=-1) {
        charge(field.queryWork());
        return field.eval(p,face);
    };
    using Key=std::array<std::int64_t,3>;
    struct Hash {
        std::size_t operator()(const Key& key) const noexcept {
            std::uint64_t value=1469598103934665603ULL;
            for (auto component:key) {
                value^=std::uint64_t(component);
                value*=1099511628211ULL;
            }
            return std::size_t(value);
        }
    };
    struct Site { Vec3 point; double size; };
    std::vector<Site> known;
    known.reserve(std::size_t(boundary.nodeCount())+maximumSites);
    std::unordered_map<Key,std::vector<std::size_t>,Hash> buckets;
    const auto origin=boundary.node(0);
    const auto keyAt=[&](Vec3 point,Key& key) {
        const auto relative=(point-origin)/bucketWidth;
        const std::array values{relative.x,relative.y,relative.z};
        // Leave room for the fixed neighbouring-bucket offsets below.
        const auto bound=double(std::numeric_limits<std::int64_t>::max()/4);
        for (unsigned axis=0;axis<3;++axis) {
            if (!std::isfinite(values[axis]) || std::abs(values[axis])>bound) return false;
            key[axis]=std::int64_t(std::floor(values[axis]));
        }
        return true;
    };
    const auto remember=[&](Vec3 p,double h,const Key& key) {
        buckets[key].push_back(known.size());
        known.push_back({p,h});
    };
    for (Index i=0;i<boundary.nodeCount();++i) {
        charge(1);
        const auto p=boundary.node(i);
        const auto classification=boundary.nodeClass(i);
        const auto h=sizeAt(p,classification.dim==2?classification.tag:-1);
        Key key{};
        if (!finite(p) || !(h>0.) || !std::isfinite(h) || !keyAt(p,key)) return result;
        remember(p,h,key);
    }
    const auto propose=[&](Vec3 p,bool box) {
        if (counts.boxSites+counts.faceSites>=maximumSites) return false;
        if (box) ++counts.boxCandidates; else ++counts.faceCandidates;
        charge(1);
        if (!finite(p)) {++counts.invalidCandidates;return false;}
        const auto h=sizeAt(p);
        Key key{};
        if (!(h>0.) || !std::isfinite(h) || h>baseSize || !keyAt(p,key)) {
            ++counts.invalidCandidates;
            return false;
        }
        // The field never exceeds baseSize, so every violating neighbour is
        // in these 27 buckets. Check its cached size rather than projecting
        // onto the CAD model again for each candidate/neighbour pair.
        for (int x=-1;x<=1;++x) for (int y=-1;y<=1;++y) for (int z=-1;z<=1;++z) {
            const auto found=buckets.find({key[0]+x,key[1]+y,key[2]+z});
            if (found==buckets.end()) continue;
            for (const auto id:found->second) {
                charge(1);
                const auto separation=.65*std::min(h,known[id].size);
                const auto delta=known[id].point-p;
                if (std::hypot(delta.x,delta.y,delta.z)<separation) {
                    ++counts.spacingRejects;
                    return false;
                }
            }
        }
        if (!admitNode(p)) {++counts.nativeRejects;return false;}
        result.addNode(p,region);
        remember(p,h,key);
        if (box) ++counts.boxSites; else ++counts.faceSites;
        return true;
    };

    // Enclosed controls must influence the first triangulation, even when no
    // boundary node intersects the box. Reserve most sites for surface support.
    const auto controls=field.boxes().controls();
    const auto boxCount=std::count_if(controls.begin(),controls.end(),[&](const auto& box) {
        return box.kind==sizing::SizeRegionKind::Box && box.targetSize<baseSize;
    });
    const auto boxBudget=std::max<std::size_t>(1,maximumSites/4);
    const auto share=std::max<std::size_t>(1,
        boxBudget/std::max<std::size_t>(1,std::size_t(boxCount)));
    for (const auto& box:controls) {
        charge(1);
        if (!sampleControls) break;
        if (box.kind!=sizing::SizeRegionKind::Box || box.targetSize>=baseSize) continue;
        const auto h=std::max(field.boxes().minimumSize(),box.targetSize);
        const auto extent=box.upper-box.lower;
        const std::array lengths{extent.x,extent.y,extent.z};
        std::array<unsigned,3> count{};
        for (unsigned axis=0;axis<3;++axis) count[axis]=unsigned(std::clamp(
            std::ceil(lengths[axis]/h),1.,double(share)));
        while (std::size_t(count[0])*count[1]*count[2]>share) {
            auto largest=std::max_element(count.begin(),count.end());
            *largest=(*largest+1)/2;
        }
        for (unsigned i=0;i<count[0]&&counts.boxSites<boxBudget;++i)
            for (unsigned j=0;j<count[1]&&counts.boxSites<boxBudget;++j)
                for (unsigned k=0;k<count[2]&&counts.boxSites<boxBudget;++k)
                    propose({std::lerp(box.lower.x,box.upper.x,(i+.5)/count[0]),
                        std::lerp(box.lower.y,box.upper.y,(j+.5)/count[1]),
                        std::lerp(box.lower.z,box.upper.z,(k+.5)/count[2])},true);
    }

    using Edge=std::tuple<Index,Index,int>;
    std::map<Edge,double> edgeSizes;
    const auto edgeSize=[&](Index a,Index b,int face) {
        const auto key=Edge{std::min(a,b),std::max(a,b),face};
        if (const auto found=edgeSizes.find(key);found!=edgeSizes.end()) return found->second;
        const auto h=sizeAt((boundary.node(a)+boundary.node(b))*.5,face);
        edgeSizes.emplace(key,h);
        return h;
    };
    for (Index i=0;i<boundary.cellCount()&&counts.boxSites+counts.faceSites<maximumSites;++i) {
        charge(1);
        if (boundary.cell(i).type!=mesh::CellType::Tri3) continue;
        const auto* v=boundary.cellNodes(i);
        const auto a=boundary.node(v[0]),b=boundary.node(v[1]),c=boundary.node(v[2]);
        const auto u=b-a,w=c-a,normal=u.cross(w);
        const auto n2=normal.norm2();
        if (!(n2>0.) || !std::isfinite(n2)) {++counts.invalidCandidates;continue;}
        const auto face=boundary.cell(i).onEntity.tag;
        const auto h=(edgeSize(v[0],v[1],face)+edgeSize(v[1],v[2],face)+
            edgeSize(v[2],v[0],face))/3.;
        if (!(h>0.) || !std::isfinite(h)) {++counts.invalidCandidates;continue;}
        auto center=(a+b+c)/3.;
        // An acute face's circumcentre lies strictly inside the face. For a
        // right/obtuse face, the centroid is only a spacing proposal and has
        // no offcentre quality guarantee.
        if (u.dot(w)>0. && (a-b).dot(c-b)>0. && (a-c).dot(b-c)>0.)
            center=a+(w.cross(normal)*u.norm2()+normal.cross(u)*w.norm2())/(2*n2);
        const auto r2=std::max({(center-a).norm2(),(center-b).norm2(),(center-c).norm2()});
        if (!(r2<h*h)) {++counts.invalidCandidates;continue;}
        // Height inspired by 2015 IMR paper 26, eq. 2 (DOI
        // 10.1016/j.proeng.2015.10.143), with fixed face-edge field samples:
        // https://www.sandia.gov/imr/2015%20IMR%20Papers/26.pdf
        // a=min(sqrt(hbar^2-rface^2),sqrt(2/3)*hbar). This is a proposal
        // heuristic; no coupled sizing iteration or complete refinement proof.
        const auto height=std::min(std::sqrt(h*h-r2),std::sqrt(2./3.)*h);
        propose(center-normal*(height/std::sqrt(n2)),false);
    }
    if (report) *report=counts;
    return result;
}

} // namespace paacs::geo::detail
