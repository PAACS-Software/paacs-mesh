// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctSourceFallbackPolicy.hpp"
#include "paacs/sizing/BoundedSizeField.hpp"
#include <array>
#include <iostream>
#include <stdexcept>
#include <string_view>

using namespace paacs;
using namespace paacs::geo;
using namespace paacs::geo::detail;
namespace {
void check(bool condition,const char* detail){if(!condition)throw std::runtime_error(detail);}
struct Field final:sizing::BoundedSizeField {
    double eval(const Vec3&)const override{return .5;}
    bool valid()const noexcept override{return true;}
    double minimumOnBox(const Vec3&,const Vec3&)const noexcept override{return .5;}
};
void admission(){
    auto options=normalNativeCadMeshOptions();
    OcctBoundaryConstructionWitness witness{"wire_chart","pcurve_position_mismatch",4,7,3,true};
    const auto decide=[&](OcctBoundaryReason reason=OcctBoundaryReason::NonconformingEdge,
                          bool budget=false,bool available=true){
        return nativeSourceFallbackDecision(options,1,true,true,reason,witness,budget,available);
    };
    check(!decide().permitted,"provider fallback changed its default");
    options.allowNativeSourceFallback=true;
    for(const auto reason:std::array{OcctBoundaryReason::NonconformingEdge,OcctBoundaryReason::OpenBoundary,
                                    OcctBoundaryReason::MissingTriangulation})
        check(decide(reason).permitted,"named witnessed geometry cannot switch provider");
    for(const auto reason:std::array{OcctBoundaryReason::None,OcctBoundaryReason::InvalidRequest,
        OcctBoundaryReason::UnsupportedTopology,OcctBoundaryReason::GeometryError,
        OcctBoundaryReason::InvalidSurface,OcctBoundaryReason::Nonmanifold,
        OcctBoundaryReason::BudgetExceeded,OcctBoundaryReason::Cancelled,
        OcctBoundaryReason::AllocationFailed,OcctBoundaryReason::KernelFailure,
        OcctBoundaryReason::SourceChanged,OcctBoundaryReason::InternalError})
        check(!decide(reason).permitted,"terminal/invalid/unknown source reason bypassed its gate");
    check(!decide(OcctBoundaryReason::NonconformingEdge,true).permitted,
          "construction-budget refusal was hidden by a topology enum");
    check(!decide(OcctBoundaryReason::NonconformingEdge,false,false).permitted,
          "source-provider fallback exceeded the shared attempt bound");
    witness.geometricRefusal=false;witness.detail="chart_node_limit";
    check(!decide().permitted,"resource witness was reinterpreted as geometry");
    witness={};check(!decide().permitted,"an absent witness permitted a retry");
    witness={"shared_edge_sampling","curve_parameter_resolution",4,-1,{},true};
    check(decide().permitted,"curve-resolution witness cannot change provider");
    for(auto solids:std::array<std::size_t,3>{0,2,9})
        check(!nativeSourceFallbackDecision(options,solids,true,true,
            OcctBoundaryReason::NonconformingEdge,witness,false,true).permitted,
            "provider switch crossed the single-body scope");
    for(auto state:std::array{std::array{false,true},std::array{true,false}})
        check(!nativeSourceFallbackDecision(options,1,state[0],state[1],
            OcctBoundaryReason::NonconformingEdge,witness,false,true).permitted,
            "provider toggle silently selected a different constructor");
    options.useNativeSurfaceFallback=false;check(!decide().permitted,"disabled native provider was used");
    options.useNativeSurfaceFallback=true;
    std::atomic_bool cancellation=true;options.geometry.cancellation=&cancellation;
    check(!decide().permitted&&std::string_view(decide().diagnostic)=="cancelled",
          "provider substitution ignored cancellation");
}
void controls(){
    auto options=normalNativeCadMeshOptions();options.allowNativeSourceFallback=true;
    const OcctBoundaryConstructionWitness witness{"wire_chart","pcurve_fit_failed",4,7,3,true};
    const auto allowed=[&]{return nativeSourceFallbackDecision(options,1,true,true,
        OcctBoundaryReason::NonconformingEdge,witness,false,true).permitted;};
    options.faceSizes.push_back({7,.1,1.3});check(!allowed(),"native provider bypassed face sizing");
    options.faceSizes.clear();options.sizeBoxes.push_back({{0,0,0},{1,1,1},.1,1.3});
    check(!allowed(),"native provider bypassed a volume control");options.sizeBoxes.clear();
    algo::GeometryBodyControl body;body.solidTag=1;options.geometry.bodyControls.push_back(body);
    check(!allowed(),"native provider bypassed body controls");options.geometry.bodyControls.clear();
    const Field field;options.geometry.spatialField=&field;
    check(!allowed(),"native provider bypassed an arbitrary spatial field");options.geometry.spatialField=nullptr;
    options.geometry.maximumElementSizeRatio=1.4;check(!allowed(),"size-ratio requirement was bypassed");
    options.geometry.maximumElementSizeRatio=0;options.geometry.warningElementSizeRatio=1.4;
    check(!allowed(),"size-ratio warning scope was changed");options.geometry.warningElementSizeRatio=0;
    options.construction.surfacePointsPerFace=10;check(!allowed(),"explicit surface point cap was bypassed");
    options.construction.surfacePointsPerFace=0;options.construction.surfaceChordPointsPerFace=10;
    check(!allowed(),"explicit chord point cap was bypassed");options.construction.surfaceChordPointsPerFace=0;
    options.minimumSize=.1;
    options.construction.interiorSites=100;options.construction.interiorCandidates=1000;
    check(allowed(),"unchanged live-constructor controls incorrectly disable the native source");
}
}
int main(){try{admission();controls();
    std::cout<<"native source-provider geometry/resource admission and control preservation passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
