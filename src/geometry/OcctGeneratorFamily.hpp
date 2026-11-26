// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/geometry/OcctMesher.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string_view>
#include <vector>

namespace paacs::geo::detail {

// This is a route eligibility hint, not a sweep or material certificate.
// It preserves the existing audited legacy/partition constructor for simple
// analytic generator families. Explicit controls retain their usual dispatch.
inline bool independentGeneratorHintAllowed(const NativeCadMeshOptions& o,std::size_t originalBodies) noexcept {
    return originalBodies>1&&o.assemblyPolicy==NativeCadAssemblyPolicy::IndependentParts&&
        o.elements==algo::PrimitiveElements::Tet4&&!o.useAdvancingInteriorFront.has_value()&&
        o.construction.resourcePolicy==CadResourcePolicy::Normal&&
        o.geometry.strategy==algo::GeometryStrategy::Auto&&o.geometry.bodyControls.empty()&&
        !o.geometry.spatialField&&o.geometry.maximumElementSizeRatio==0&&o.geometry.warningElementSizeRatio==0&&
        o.faceSizes.empty()&&o.sizeBoxes.empty()&&!o.useSpatialSizing&&o.minimumSize==0&&
        !o.improveNativeSurfaceSpacing&&!o.useRegularSurfaceFront&&!o.regularMeshQualityEnabled()&&
        !o.curvaturePointsPerTurn&&!o.improveNativeTetQuality&&!o.useLocalDelaunayRefinement&&!o.useLiveMaterialRecovery&&
        !o.initialTopologyQuality.has_value()&&!o.optimizeProfileInterior&&
        o.construction.effort==CadConstructionEffort::Legacy&&
        !o.construction.surfacePointsPerFace&&!o.construction.surfaceChordPointsPerFace&&
        !o.construction.interiorCandidates&&!o.construction.interiorSites;
}

inline bool generatorLegacyRetryReason(NativeCadMeshReason reason) noexcept {
    return reason==NativeCadMeshReason::GeometryError||reason==NativeCadMeshReason::UnsupportedBody||
        reason==NativeCadMeshReason::RecognitionFailed;
}
inline bool independentAutoMayRetryLegacy(const NativeCadMeshOptions& options,std::size_t originalBodies,
                                        const NativeCadMeshOutcome& first) noexcept {
    return independentGeneratorHintAllowed(options,originalBodies)&&!first.published&&first.destinationUnchanged&&
        generatorLegacyRetryReason(first.reason)&&!first.routeResolution.wholeBodyVerified&&
        first.routeResolution.count==1&&first.routeResolution.candidate&&first.routeResolution.attempts[0].method&&
        std::string_view(first.routeResolution.candidate)=="native_front"&&
        std::string_view(first.routeResolution.attempts[0].method)=="native_front";
}
inline std::uint64_t generatorRetryWork(std::uint64_t original,std::uint64_t spent,std::uint64_t sharedRemaining) noexcept {
    return spent<original?std::min(original-spent,sharedRemaining):0;
}
inline bool accumulateIndependentAttemptMetrics(NativeCadMeshOutcome& parent,const NativeCadMeshOutcome& child) noexcept {
    auto p=std::array{parent.nativeProjectionQueries,parent.nativeProjectionInitializations,parent.nativeProjectionReuses,
        parent.nativeTrimmedProjectionQueries,parent.nativeTrimmedProjectionMemoHits};
    const auto c=std::array{child.nativeProjectionQueries,child.nativeProjectionInitializations,child.nativeProjectionReuses,
        child.nativeTrimmedProjectionQueries,child.nativeTrimmedProjectionMemoHits};
    auto seconds=std::array{parent.preparationSeconds,parent.generationSeconds,parent.validationSeconds,parent.auditSeconds,parent.qualitySeconds};
    const auto added=std::array{child.preparationSeconds,child.generationSeconds,child.validationSeconds,child.auditSeconds,child.qualitySeconds};
    for(std::size_t i=0;i<p.size();++i){if(c[i]>UINT64_MAX-p[i])return false;p[i]+=c[i];
        if(!std::isfinite(added[i])||added[i]<0||!std::isfinite(seconds[i]+added[i]))return false;seconds[i]+=added[i];}
    parent.nativeProjectionQueries=p[0];parent.nativeProjectionInitializations=p[1];parent.nativeProjectionReuses=p[2];
    parent.nativeTrimmedProjectionQueries=p[3];parent.nativeTrimmedProjectionMemoHits=p[4];
    parent.preparationSeconds=seconds[0];parent.generationSeconds=seconds[1];parent.validationSeconds=seconds[2];
    parent.auditSeconds=seconds[3];parent.qualitySeconds=seconds[4];return true;
}

template<class Tick>
bool remapIndependentSourceShell(NativeCadSourceShellReport& shell,const std::array<std::vector<int>,4>& tags,Tick&& tick){
    if(shell.count>shell.attempts.size())return false;
    const auto map=[&](int& tag,unsigned dimension){
        if(tag<=0){tag=-1;return true;}
        if(static_cast<std::size_t>(tag)>=tags[dimension].size()||tags[dimension][tag]<=0)return false;
        tag=tags[dimension][tag];return true;
    };
    for(unsigned i=0;i<shell.count;++i){
        tick();auto witness=shell.attempts[i].witness;
        if(!map(witness.edgeTag,1)||!map(witness.faceTag,2))return false;
        // Edge-use IDs enumerate face/wire occurrences in one import. A part
        // copy has a different enumeration; entity history cannot map them.
        witness.useId.reset();shell.attempts[i].witness=witness;
    }
    return true;
}

template<class Tick>
bool remapIndependentChildDiagnostics(NativeCadMeshOutcome& child,const std::array<std::vector<int>,4>& tags,Tick&& tick){
    if(!remapIndependentSourceShell(child.initialFrontFallback.sourceShell,tags,tick))return false;
    for(auto& body:child.bodies){
        tick();
        if(!remapIndependentSourceShell(body.sourceShell,tags,tick)||
           !remapIndependentSourceShell(body.initialFrontFallback.sourceShell,tags,tick))return false;
        for(auto& spacing:body.surfaceSpacing){
            tick();const auto tag=spacing.faceTag;
            if(tag<=0||static_cast<std::size_t>(tag)>=tags[2].size()||tags[2][tag]<=0)return false;
            spacing.faceTag=tags[2][tag];
        }
    }
    return true;
}

enum class GeneratorFaceKind { Other,Plane,Cylinder };
struct GeneratorFace {
    GeneratorFaceKind kind=GeneratorFaceKind::Other;
    Vec3 origin,axis;
    double tolerance=0;
};

// Two scans, constant auxiliary storage. faceAt and tick are supplied by the
// Native wrapper so every geometry query uses its shared work/cancellation
// admission. Cap levels use the reported face tolerances, without changing
// any tolerances or the final original-CAD acceptance checks.
template<class FaceAt,class Tick>
bool analyticGeneratorFamily(std::size_t count,FaceAt&& faceAt,Tick&& tick) {
    if(count<4)return false;
    const auto finite=[](Vec3 p){return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);};
    const auto valid=[&](const GeneratorFace& f){return finite(f.origin)&&finite(f.axis)&&
        std::isfinite(f.tolerance)&&f.tolerance>=0&&std::abs(f.axis.norm2()-1)<=1e-10;};
    // Analytic axes returned by OCCT are unit directions. This tolerance only
    // chooses a constructor; it does not certify approximate extrusion.
    constexpr double angular=1e-10;
    Vec3 axis,origin;bool cylinder=false;
    for(std::size_t i=0;i<count;++i){
        tick();const auto f=faceAt(i);
        if(f.kind==GeneratorFaceKind::Other||!valid(f))return false;
        if(f.kind==GeneratorFaceKind::Cylinder){
            if(!cylinder){axis=f.axis;origin=f.origin;cylinder=true;}
            else if(axis.cross(f.axis).norm()>angular)return false;
        }
    }
    if(!cylinder)return false;
    unsigned levels=0;bool side=false;
    double lower[2]{},upper[2]{},extentLower[2]{},extentUpper[2]{};
    for(std::size_t i=0;i<count;++i){
        tick();const auto f=faceAt(i);
        if(f.kind!=GeneratorFaceKind::Plane)continue;
        if(axis.cross(f.axis).norm()<=angular){
            const double offset=(f.origin-origin).dot(axis);
            if(!std::isfinite(offset))return false;
            const double lo=offset-f.tolerance,hi=offset+f.tolerance;
            if(!std::isfinite(lo)||!std::isfinite(hi))return false;
            unsigned match=2;
            for(unsigned level=0;level<levels;++level)if(lo<=upper[level]&&lower[level]<=hi){
                // Do not choose an arbitrary cap when a broad source
                // tolerance overlaps both. A conservative refusal is safe.
                if(match!=2)return false;match=level;
            }
            if(match==2){if(levels==2)return false;lower[levels]=extentLower[levels]=lo;upper[levels]=extentUpper[levels]=hi;++levels;}
            else{lower[match]=std::max(lower[match],lo);upper[match]=std::min(upper[match],hi);
                extentLower[match]=std::min(extentLower[match],lo);extentUpper[match]=std::max(extentUpper[match],hi);}
        }else if(std::abs(axis.dot(f.axis))<=angular)side=true;
        else return false;
    }
    if(!side||levels!=2||!(upper[0]<lower[1]||upper[1]<lower[0]))return false;
    const unsigned left=upper[0]<lower[1]?0:1,right=1-left;
    // A broad interval seen before the other cap also must not overlap that
    // final cap. Retained extrema detect this without a third Native scan.
    return extentUpper[left]<lower[right]&&extentLower[right]>upper[left];
}

} // namespace paacs::geo::detail
